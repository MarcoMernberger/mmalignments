"""Wrapper for mmfqcount — a fast Rust-based FASTQ read-sequence counter.

mmfqcount has two subcommands:

``count``
    Iterates one (single-end) or two (paired-end) FASTQ files and writes a TSV
    table of every unique sequence (or sequence pair) with its count and one
    representative read name.  Supports optional adapter trimming.

``match``
    Reads a counts TSV produced by ``count`` and a predefined-sequences TSV and
    produces two output files: matched sequences (with Read Count + Frequency
    columns appended) and unmatched sequences (counted but not predefined).

The class mirrors the ``FastGrab`` pattern:
- each step exposes a high-level ``@element`` method that returns an
  :class:`~mmalignments.models.elements.Element` (with full dependency
  tracking, signature-based caching, etc.)
- and a low-level ``@subroutine`` method that builds the subprocess command.

Parameters are described inline (no separate JSON file needed — we own the
tool and its CLI is stable).
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Callable, Mapping

import numpy as np
import pandas as pd
from pandas import DataFrame, Series

from mmalignments.models.artifacts import (
    ArtifactSet,
    FileSpec,
    OutputSpec,
)
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FileSource,
    NextGenSample,
    TableElement,
    element,
)
from mmalignments.models.parameters import (
    ParamRegistry,
    Params,
    ParamSet,
    ParamSpec,
    render_value,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.dependencies import function_hash
from mmalignments.services.io import parents, read_frame, write_frames

from ..externals import (
    External,
    ExternalRunConfig,
    Runnable,
    SubroutineIn,
    subroutine,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameter definitions (inlined — no JSON file needed for own Rust tools)
# ---------------------------------------------------------------------------


def _build_param_registry() -> ParamRegistry:
    """Parameter registry for mmfqcount.

    Defined here in Python because we own and maintain the Rust source —
    the CLI is stable and there is no need to parse --help output or maintain
    a separate JSON file.  If the CLI ever grows new flags, add them here.

    count subcommand
    ----------------
    --trim-start   str   Trim read from the first occurrence of this k-mer (inclusive).
    --trim-stop    str   Trim read up to (exclusive) the last occurrence of this k-mer.
    --trim-length  int   Keep at most this many bases after adapter trimming.
    --split-by     str   Split counts by this tag in read names (e.g. "sgRNAid").

    match subcommand
    ----------------
    --seq-col      str   Column holding the R1 sequence.   default: "Sequence"
    --r2-col       str   Column holding the R2 sequence (paired mode).
    --id-col       str   Column holding the sequence identifier.  default: "Name"
    """
    count_specs: dict[str, ParamSpec] = {
        "trim_start": ParamSpec(
            "trim_start",
            "--trim-start",
            str,
            render=render_value,
            description="Trim from first occurrence of k-mer (inclusive).",
        ),
        "trim_stop": ParamSpec(
            "trim_stop",
            "--trim-stop",
            str,
            render=render_value,
            description="Trim up to (exclusive) last occurrence of k-mer.",
        ),
        "trim_length": ParamSpec(
            "trim_length",
            "--trim-length",
            int,
            render=render_value,
            description="Keep at most this many bases after trimming.",
        ),
        "split_by": ParamSpec(
            "split_by",
            "--split-by",
            str,
            render=render_value,
            description="Split counts by this read-name tag (e.g. 'sgRNAid').",
        ),
    }
    match_specs: dict[str, ParamSpec] = {
        "seq_col": ParamSpec(
            "seq_col",
            "--seq-col",
            str,
            default="Sequence",
            render=render_value,
            description="Column holding the R1 sequence.",
        ),
        "r2_col": ParamSpec(
            "r2_col",
            "--r2-col",
            str,
            render=render_value,
            description="Column holding the R2 sequence (paired mode).",
        ),
        "id_col": ParamSpec(
            "id_col",
            "--id-col",
            str,
            default="Name",
            render=render_value,
            description="Column holding the sequence identifier.",
        ),
    }
    return ParamRegistry(
        default=ParamSet({}, "mmfqcount", "default"),
        by_subcommand={
            "count": ParamSet(count_specs, "mmfqcount", "count"),
            "match": ParamSet(match_specs, "mmfqcount", "match"),
        },
    )


# ---------------------------------------------------------------------------
# MmFqCount class
# ---------------------------------------------------------------------------


class MmFqCount(External):
    """Wrapper for the ``mmfqcount`` Rust CLI tool.

    ``mmfqcount`` counts FASTQ read sequences by identity and can match
    counted sequences against a predefined set.  It is fast (AHash, parallel
    I/O), supports plain and gzip-compressed FASTQs, and handles both
    single-end and paired-end reads.

    Typical usage::

        counter = MmFqCount()

        # Step 1 — count all unique sequences in a sample
        count_el = counter.count(sample)

        # Step 2 — match against known sequences
        match_el = counter.match(
            count_el,
            predefined=predefined_tsv_element,
            seq_col="Sequence",
            id_col="Name",
        )

    With trimming::

        count_el = counter.count(
            sample,
            params=Params(trim_start="ACGT", trim_length=20),
        )
    """

    def __init__(
        self,
        name: str = "mmfqcount",
        primary_binary: str = "mmfqcount",
        version: str | None = None,
        source: str = "https://github.com/MarcoMernberger/mmfqcount",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialise the MmFqCount wrapper.

        Parameters
        ----------
        name : str
            Logical tool name, default ``"mmfqcount"``.
        primary_binary : str
            Executable name on ``$PATH``, default ``"mmfqcount"``.
        version : str | None
            Optional version string override; auto-detected when *None*.
        source : str
            Documentation / source URL.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Parameter sets for invocations.  When *None*, the built-in
            parameter registry (derived from the stable CLI) is used.
        """
        # Build parameters from the known CLI instead of a JSON file
        resolved_parameters: Mapping[str, ParamSet] | ParamSet = (
            parameters if parameters is not None else _build_param_registry_as_mapping()
        )
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=resolved_parameters,
        )
        # Overwrite param_registry with our richer version (includes
        # per-subcommand ParamSets) regardless of how External.__init__
        # processes the mapping.
        self.param_registry = _build_param_registry()

    # -----------------------------------------------------------------------
    # Version detection
    # -----------------------------------------------------------------------

    def get_version(self, fallback: str | None = None) -> str | None:
        """Return the mmfqcount version string.

        Runs ``mmfqcount --version`` and extracts the semver token.  The
        output format is stable because we own the Rust source
        (clap auto-generates ``mmfqcount X.Y.Z``).

        Parameters
        ----------
        fallback : str | None
            Value returned when the version cannot be determined.

        Returns
        -------
        str | None
            Version string (e.g. ``"0.1.0"``) or *fallback*.
        """
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if cp.returncode == 0 and cp.stdout:
                # clap output: "mmfqcount 0.1.0\n"
                m = re.search(r"(\d+\.\d+\.\d+)", cp.stdout)
                if m:
                    return m.group(1)
        except Exception:
            pass
        return fallback

    # -----------------------------------------------------------------------
    # Default paths
    # -----------------------------------------------------------------------

    def default_output_dir(self, sample_name: str) -> Path:
        """Return the default output directory for a given sample."""
        return Path("results") / "counts" / self.version_name / sample_name

    def default_output(
        self, sample_name: str = "", outdir: Path | None = None, ext: str = "tsv"
    ) -> OutputSpec:
        """Return the default output spec for a given sample."""
        return OutputSpec(
            outdir=outdir or self.default_output_dir(sample_name),
            ext=ext,
        )

    # -----------------------------------------------------------------------
    # count — high-level @element
    # -----------------------------------------------------------------------

    @element
    def count(
        self,
        sample: Element | NextGenSample,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Count all unique read sequences (or read pairs) in a sample.

        Runs ``mmfqcount count`` on the FASTQ file(s) of *sample* and writes
        a TSV table with columns:

        * single-end: ``R1``, ``Count``, ``R1 Name``
        * paired-end: ``R1``, ``R2``, ``Count``, ``R1 Name``, ``R2 Name``

        The output is sorted by ``Count`` descending.

        Parameters
        ----------
        sample : Element | NextGenSample
            Input sample. Both single-end and paired-end samples are handled
            automatically.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None
            Output specification; defaults to
            ``results/counts/<version>/<sample.name>`` with extension ``tsv``.
        params : Params | None
            Trimming parameters.  Recognised keys:

            ``trim_start`` (str)
                Trim from the first occurrence of this k-mer (inclusive).
                Reads without the k-mer are discarded.
            ``trim_stop`` (str)
                Trim up to (exclusive) the last occurrence of this k-mer.
            ``trim_length`` (int)
                Keep at most this many bases after trimming.

        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        Element
            Element whose artifact ``"tsv"`` is the path to the counts TSV.
        """
        fastq_r1 = sample.primary.r1
        fastq_r2 = sample.primary.r2 if hasattr(sample.primary, "r2") else None

        tag = from_prior(
            sample.tag,
            tag,
            stage=Stage.QUANT,
            method=Method.MMFQCOUNT,
            state=State.COUNT,
            ext="tsv",
        )

        spec = self.default_output(sample.root).merge(outspec)
        output_file = spec.path(tag.default_name)

        runner = self.run_count(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            output_file=output_file,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(tag, "count")
        determinants = self.signature_determinants(params, subroutine="count")
        inputs = (fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,)
        pres = (sample,) if isinstance(sample, Element) else sample.pres
        artifacts = ArtifactSet(output_file, primary_name=spec.ext)
        return Element(
            key,
            runner,
            tag=tag,
            determinants=determinants,
            artifacts=artifacts,
            inputs=inputs,
            pres=pres,
            name=name,
        )

    # -----------------------------------------------------------------------
    # count — low-level @subroutine
    # -----------------------------------------------------------------------

    @subroutine
    def run_count(
        self,
        fastq_r1: Path,
        fastq_r2: Path | None,
        output_file: Path,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mmfqcount count``.

        Parameters
        ----------
        fastq_r1 : Path | str
            R1 FASTQ file (plain or gzip).
        fastq_r2 : Path | str | None
            R2 FASTQ file (plain or gzip), or *None* for single-end.
        output_file : Path | str
            Destination path for the counts TSV.
        params : Params | None
            Trimming parameters (``trim_start``, ``trim_stop``,
            ``trim_length``).
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        SubroutineIn
            Tuple consumed by the ``@subroutine`` decorator.
        """
        arguments = [
            "count",
            "--r1",
            str(fastq_r1),
            "--output",
            str(output_file),
        ]
        if fastq_r2 is not None:
            arguments += ["--r2", str(Path(fastq_r2).resolve())]
        in_paths = [fastq_r1] + ([Path(fastq_r2).resolve()] if fastq_r2 else [])
        out_paths = [output_file]

        return (
            arguments,
            "count",
            in_paths,
            out_paths,
            None,  # no piped output
            None,  # no pre-hook
            None,  # no post-hook
        )

    # -----------------------------------------------------------------------
    # match — high-level @element
    # -----------------------------------------------------------------------

    @element
    def match(
        self,
        counts: Element,
        predefined: Element | FileSource,
        *,
        seq_col: str = "Sequence",
        r2_col: str | None = None,
        id_col: str = "Name",
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Match counted sequences against a predefined-sequences table.

        Runs ``mmfqcount match`` on the counts TSV produced by :meth:`count`
        and writes two output files:

        * **matched TSV** — all rows from the predefined table, with
          ``Read Count`` and ``Frequency`` columns appended.
        * **unmatched TSV** — all counted sequences that had no predefined
          match, sorted by count descending.

        Parameters
        ----------
        counts : Element
            Element produced by :meth:`count`; its ``"tsv"`` artifact is used
            as the counts input.
        predefined : Element | FileSource
            Either a pipeline Element whose ``"tsv"`` (or ``"path"``) artifact
            points to the predefined-sequences TSV, or a plain FileSource.
        seq_col : str
            Column in the predefined TSV holding the R1 sequence.
            Default: ``"Sequence"``.
        r2_col : str | None
            Column in the predefined TSV holding the R2 sequence
            (paired mode).  *None* for single-end matching.
        id_col : str
            Column in the predefined TSV holding the sequence identifier.
            Default: ``"Name"``.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None
            Output specification.
        params : Params | None
            Additional ``match`` parameters (``seq_col``, ``r2_col``,
            ``id_col`` can also be supplied here as ``Params`` overrides
            instead of keyword arguments).
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        Element
            Element with artifacts:

            ``"matched"``   → matched sequences TSV
            ``"unmatched"`` → unmatched sequences TSV
        """
        # count file
        counts_file = counts.primary.resolve()
        # Resolve predefined path
        pred_path = predefined.primary.resolve()
        if isinstance(predefined, Element):
            pred_pres: tuple = (predefined,)
        elif isinstance(predefined, FileSource):
            pred_pres = ()
        else:
            raise TypeError(
                f"predefined must be an Element or FileSource, got {type(predefined)}"
            )
        # Merge column settings: explicit kwargs take priority over params
        merged_params = Params(
            seq_col=seq_col,
            id_col=id_col,
            **({"r2_col": r2_col} if r2_col else {}),
        )
        if params:
            merged_params = merged_params.override(**params.to_dict())
        tag = from_prior(
            counts.tag,
            tag,
            stage=Stage.QUANT,
            method=Method.MMFQCOUNT,
            state=State.ANNOTATED,
        )
        spec = self.default_output(outdir=counts_file.parent).merge(outspec)
        spec = spec.add_output(
            "unmatched", FileSpec(tag.default_name + ".unmatched", ext=spec.ext)
        )
        artifacts = ArtifactSet.generate(tag, spec=spec)
        runner = self.run_match(
            counts_tsv=counts_file,
            predefined_tsv=pred_path,
            matched_tsv=artifacts.primary.resolve(),
            unmatched_tsv=artifacts["unmatched"].resolve(),
            params=merged_params,
            cfg=cfg,
        )

        key, name = self.build_element_name(
            tag, "match", seq_col=seq_col, id_col=id_col
        )
        determinants = self.signature_determinants(merged_params, subroutine="match")
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(counts_file, pred_path),
            pres=(counts,) + pred_pres,
            name=name,
        )

    # -----------------------------------------------------------------------
    # match — low-level @subroutine
    # -----------------------------------------------------------------------

    @subroutine
    def run_match(
        self,
        counts_tsv: Path | str,
        predefined_tsv: Path | str,
        matched_tsv: Path | str,
        unmatched_tsv: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``mmfqcount match``.

        Parameters
        ----------
        counts_tsv : Path | str
            Counts TSV produced by ``mmfqcount count``.
        predefined_tsv : Path | str
            Predefined-sequences TSV.
        matched_tsv : Path | str
            Output path for matched sequences.
        unmatched_tsv : Path | str
            Output path for unmatched sequences.
        params : Params | None
            Match parameters (``seq_col``, ``r2_col``, ``id_col``).
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        SubroutineIn
            Tuple consumed by the ``@subroutine`` decorator.
        """
        counts_tsv = Path(counts_tsv).absolute()
        predefined_tsv = Path(predefined_tsv).absolute()
        matched_tsv = Path(matched_tsv).absolute()
        unmatched_tsv = Path(unmatched_tsv).absolute()

        arguments = [
            "match",
            "--counts",
            str(counts_tsv),
            "--predefined",
            str(predefined_tsv),
            "--output",
            str(matched_tsv),
            "--unmatched",
            str(unmatched_tsv),
        ]

        # Append column flags from params
        params = params or Params()
        # cli_extras = self.to_cli(params, subroutine="match")
        # arguments += cli_extras
        return (
            arguments,
            "match",
            [counts_tsv, predefined_tsv],
            [matched_tsv, unmatched_tsv],
            None,
            None,
            None,
        )

    # -----------------------------------------------------------------------
    # match — high-level @element
    # -----------------------------------------------------------------------

    @element
    def compare(
        self,
        compare: TableElement,
        against: TableElement,
        *,
        score: Callable[[Series, Series], Series] | None = None,
        keys: list[str] | None = ["R2", "Annotation"],
        filter: Callable[[pd.DataFrame], pd.DataFrame] | bool = True,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
        mode: str = "both",
    ) -> Element:
        """
        compare the counts of two count results and output a tsv with the counts of
        both samples and a predefined score between them.
        Both input elements should be TableElements with a "tsv" artifact. The
        TSVs must have the same counts columns (e.g. "Count") and the same
        sequence identifier columns (e.g. "Sequence").

        The result will produce 2 tsv files: a score TSV with all sequences
        from the "compare" table, with the counts from both inputs and the
        score; and a filtered score TSV with only sequences that meet the score
        threshold (if exclude_score is set).

        Parameters
        ----------
        compare : NextGenSampleElement
            Element containing the counts TSV to compare (condition).
        against : NextGenSampleElement
            Element containing the counts TSV to compare against (control).
        score : Callable[[int, int], float]] | None, optional
            Function that takes two counts (compare, against) and returns a score.
            Default is log2 fold change with inf if denominator is zero.
        key_columns : list[str] | None, optional
            Column names in the input TSVs that contain the sequence
            identifiers. Default is ["R2", "Annotation"].
        count_column : str, optional
            Column name in the input TSVs that contains the counts. Default is "Count".
        id_column : str, optional
            Column name in the input TSVs that contains the identifiers. Default is
            "Annotation".
        exclude_score : float | None | Callable, optional
            If set, sequences with a score equal to this will be excluded from
            the filtered score TSV. If a callable is provided, it will be called
            with the score and should return True if the sequence should be
            excluded.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None, optional
            Optional output specification override.
        params : Params | None, optional
            Additional parameters, by default None
        cfg : ExternalRunConfig | None, optional
            Runtime configuration, by default None

        Returns
        -------
        TableElement
            Element with artifacts:

            "score" → TSV with all sequences and their scores
            "filtered_score" → TSV with sequences that meet the score threshold
        """

        def default_score(count1: Series, count2: Series) -> Series:
            offset = 1
            c1 = np.asarray(count1, dtype=np.float64)
            c2 = np.asarray(count2, dtype=np.float64)
            c1 = c1 + offset
            c2 = c2 + offset
            np.divide(c1, c1.sum(), out=c1)  # normalize to lib size
            np.divide(c2, c2.sum(), out=c2)
            np.divide(c1, c2, out=c1)  # ratio
            np.log2(c1, out=c1)  # log fold change
            return pd.Series(c1, index=count1.index)

        def get_filter(
            filter_column: str, filter: Callable[[pd.DataFrame], pd.DataFrame] | bool
        ) -> Callable[[DataFrame], DataFrame] | None:
            def filter_func(df: DataFrame) -> DataFrame:
                df = df[df[filter_column] > 0]
                return df

            if callable(filter):
                return filter
            elif filter is True:
                return filter_func
            else:
                return None

        tag = from_prior(
            compare.tag,
            tag,
            stage=Stage.QUANT,
            method=Method.MMFQCOUNT,
            state=State.SCORE,
            ext="tsv",
        )
        params = Params(
            count_column="Count",
            freq_column="Frequency",
            annotation_column="Annotation",
            seq_column="R2",
            score_name="Score" if score else "Log2 Relative Enrichment Score",
        ).update(params)

        keys = keys or ["R2", "Annotation"]
        filter_func = get_filter(f"{params.count_column} ({compare.tag.root})", filter)
        spec = self.default_output(outdir=compare.file.parent).merge(outspec)
        score_tsv = spec.path(tag.default_name)
        pres = (compare, against)
        inputs = (compare.file, against.file)
        determinants = [function_hash(score), str(params)] + keys
        determinants += (
            [function_hash(filter_func)]
            if callable(filter_func)
            else [str(f"filter={filter}")]
        )
        runner = self.compare_counts(
            compare=compare.tag.root,
            against=against.tag.root,
            compare_file=compare.file,
            against_file=against.file,
            out_tsv=score_tsv,
            score=score or default_score,
            keys=keys,
            filter=filter_func,
            mode=mode,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(
            tag,
            "compare",
            param_str="_against=" + against.tag.root,
        )
        artifacts = {"tsv": score_tsv}
        if filter:
            artifacts["filtered"] = score_tsv.with_suffix(".filtered.tsv")

        return Element(
            key,
            runner,
            tag=tag,
            determinants=tuple(determinants),
            inputs=inputs,
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    # -----------------------------------------------------------------------
    # compare_counts — low-level helper
    # -----------------------------------------------------------------------

    def compare_counts(
        self,
        compare: str,
        against: str,
        compare_file: Path | str,
        against_file: Path | str,
        out_tsv: Path | str,
        score: Callable[[Series, Series], Series],
        *,
        keys: list[str] | None = ["R2", "Annotation"],
        filter: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
        mode: str = "both",
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:
        """Compare two count TSVs and assign a score to each sequence."""
        params = Params(
            count_column="Count",
            freq_column="Frequency",
            annotation_column="Annotation",
            seq_column="R2",
            score_name="Score" if score else "Log2 Relative Enrichment Score",
            how="left",
        ).update(params)
        count_column = params.count_column
        freq_column = params.freq_column

        def reduce_to_keys(
            key_columms: list[str],
        ) -> Callable[[pd.DataFrame], pd.DataFrame]:
            def __reduce(df: pd.DataFrame) -> DataFrame:
                drop_cols = [
                    c
                    for c in df.columns
                    if c
                    not in key_columms
                    + [count_column, freq_column, "R1 Name", "R2 Name"]
                ]
                df = df.drop(columns=drop_cols)
                df = df.groupby(
                    key_columms,
                    as_index=False,
                ).agg(
                    {
                        count_column: "sum",
                        freq_column: "sum",
                        "R1 Name": "first",
                        "R2 Name": "first",
                    }
                )
                return df

            return __reduce

        def merge_frames(
            compare_df: pd.DataFrame, against_df: pd.DataFrame
        ) -> pd.DataFrame:
            merged = compare_df.merge(
                against_df,
                on=keys,
                how=params.how,
                suffixes=(f" ({compare})", f" ({against})"),
            )
            merged = merged.fillna(
                dict.fromkeys(
                    [f"{count_column} ({against})", f"{freq_column} ({against})"], 0
                )
            )
            return merged

        compare_path = Path(compare_file)
        against_path = Path(against_file)
        score_path = Path(out_tsv)

        def _runner() -> None:

            parents(score_path)
            compare_df = read_frame(compare_path)
            against_df = read_frame(against_path)
            if keys:
                compare_df = reduce_to_keys(keys)(compare_df)
                against_df = reduce_to_keys(keys)(against_df)
            score_name = params.get("score_name", "score")
            compare_df = merge_frames(compare_df, against_df)
            comp_count_column = f"{count_column} ({compare})"
            against_count_column = f"{count_column} ({against})"
            for col in [f"{freq_column} ({compare})", f"{freq_column} ({against})"]:
                if col not in compare_df.columns:
                    compare_df[col] = compare_df.fillna(0)[col].astype("float64")
            for col in [comp_count_column, against_count_column]:
                compare_df[col] = compare_df[col].fillna(0).astype("int64")
            compare_df[score_name] = score(
                compare_df[comp_count_column], compare_df[against_count_column]
            )
            compare_df = compare_df.sort_values(by=score_name, ascending=False)
            write_frames(compare_df, [score_path])
            # filter out zero counts
            if filter:
                filtered_path = Path(score_path.with_suffix(".filtered.tsv"))
                filtered = filter(compare_df)
                write_frames(filtered, [filtered_path])

        callspec = CallSpec(
            path=("compare_counts",),
            kwargs={
                "compare": compare_file,
                "against": against_file,
                "compare_file": compare_file,
                "against_file": against_file,
                "out_tsv": score_path,
                "score": score,
                "keys": keys,
                "filter": filter,
                "params": params,
            },
        ).render()
        return Runnable(
            _runner,
            display=callspec,
        )

    # def a_better_compare_counts(
    #     self,
    #     compare: str,
    #     against: str,
    #     compare_file: Path | str,
    #     against_file: Path | str,
    #     out_tsv: Path | str,
    #     score: Callable[[Series, Series], Series],
    #     *,
    #     params: Params | None = None,
    #     cfg: ExternalRunConfig | None = None,
    # ) -> Runnable:
    #     """Compare two count TSVs and assign a score to each sequence."""
    #     keys: list[str] | None = ["R2", "Annotation"]
    #     filter: Callable[[pd.DataFrame], pd.DataFrame] | None = None
    #     params = Params(
    #         count_column="Count",
    #         freq_column="Frequency",
    #         annotation_column="Annotation",
    #         seq_column="R2",
    #         score_name="Score" if score else "Log2 Relative Enrichment Score",
    #         how="left",
    #     ).update(params)
    #     count_column = params.count_column
    #     freq_column = params.freq_column

    #     def reduce_to_keys(
    #         key_columms: list[str],
    #     ) -> Callable[[pd.DataFrame], pd.DataFrame]:
    #         def __reduce(df: pd.DataFrame) -> DataFrame:
    #             drop_cols = [
    #                 c
    #                 for c in df.columns
    #                 if c
    #                 not in key_columms
    #                 + [count_column, freq_column, "R1 Name", "R2 Name"]
    #             ]
    #             df = df.drop(columns=drop_cols)
    #             df = df.groupby(
    #                 key_columms,
    #                 as_index=False,
    #             ).agg(
    #                 {
    #                     count_column: "sum",
    #                     freq_column: "sum",
    #                     "R1 Name": "first",
    #                     "R2 Name": "first",
    #                 }
    #             )
    #             return df

    #         return __reduce

    #     def merge_frames(
    #         compare_df: pd.DataFrame, against_df: pd.DataFrame
    #     ) -> pd.DataFrame:
    #         print(params.how)
    #         merged = compare_df.merge(
    #             against_df,
    #             on=keys,
    #             how=params.how,
    #             suffixes=(f" ({compare})", f" ({against})"),
    #         )
    #         merged = merged.fillna(
    #             dict.fromkeys(
    #                 [f"{count_column} ({against})", f"{freq_column} ({against})"], 0
    #             )
    #         )
    #         return merged

    #     compare_path = Path(compare_file)
    #     against_path = Path(against_file)
    #     score_path = Path(out_tsv)

    #     def _runner() -> None:

    #         parents(score_path)
    #         compare_df = read_frame(compare_path)
    #         against_df = read_frame(against_path)
    #         if keys:
    #             compare_df = reduce_to_keys(keys)(compare_df)
    #             against_df = reduce_to_keys(keys)(against_df)

    #         df = (
    #             merge_frames(compare_df, against_df)
    #             .pipe(add(), score)  #  score_name = params.get("score_name", "score")
    #             .pipe(fill(), columns)
    #             .pipe(types())
    #             .pipe(sort())
    #         )
    #         compare_df = merge_frames(compare_df, against_df)
    #         comp_count_column = f"{count_column} ({compare})"
    #         against_count_column = f"{count_column} ({against})"
    #         for col in [f"{freq_column} ({compare})", f"{freq_column} ({against})"]:
    #             if col not in compare_df.columns:
    #                 compare_df[col] = compare_df.fillna(0)[col].astype("float64")
    #         for col in [comp_count_column, against_count_column]:
    #             compare_df[col] = compare_df[col].fillna(0).astype("int64")
    #         compare_df[score_name] = score(
    #             compare_df[comp_count_column], compare_df[against_count_column]
    #         )
    #         compare_df = compare_df.sort_values(by=score_name, ascending=False)
    #         write_frames(compare_df, score_path, mode)
    #         # filter out zero counts
    #         if filter:
    #             filtered_path = Path(score_path.with_suffix(".filtered.tsv"))
    #             filtered = filter(compare_df)
    #             write_frames(filtered, filtered_path, mode)

    #     callspec = CallSpec(
    #         path=("compare_counts",),
    #         kwargs={
    #             "compare": compare_file,
    #             "against": against_file,
    #             "compare_file": compare_file,
    #             "against_file": against_file,
    #             "out_tsv": score_path,
    #             "score": score,
    #             "keys": keys,
    #             "filter": filter,
    #             "params": params,
    #         },
    #     ).render()
    #     return Runnable(
    #         _runner,
    #         display=callspec,
    #     )

    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------

    def count_and_match(
        self,
        sample: NextGenSample,
        predefined: Element | FileSource,
        *,
        seq_col: str = "Sequence",
        r2_col: str | None = None,
        id_col: str = "Name",
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        count_params: Params | None = None,
        match_params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[Element, Element]:
        """Count reads and immediately match against predefined sequences.

        Convenience wrapper that chains :meth:`count` → :meth:`match`.

        Parameters
        ----------
        sample : NextGenSample
            Input sample.
        predefined : Element | Path | str
            Predefined-sequences TSV (Element or plain path).
        seq_col : str
            Sequence column in the predefined TSV. Default: ``"Sequence"``.
        r2_col : str | None
            R2-sequence column for paired matching. *None* for single-end.
        id_col : str
            ID column in the predefined TSV. Default: ``"Name"``.
        tag : PartialElementTag | ElementTag | None
            Tag override (propagated to both elements).
        outspec : OutputSpec | None
            Output specification overrides.
        count_params : Params | None
            Trimming parameters for the ``count`` step.
        match_params : Params | None
            Column parameters for the ``match`` step.
        cfg : ExternalRunConfig | None
            Subprocess configuration (shared by both steps).

        Returns
        -------
        tuple[Element, Element]
            ``(match_element, count_element)``
        """
        spec = self.default_output(sample.root).merge(outspec)
        count_el = self.count(
            sample,
            tag=tag,
            outspec=spec,
            params=count_params,
            cfg=cfg,
        )
        match_el = self.match(
            counts=count_el,
            predefined=predefined,
            seq_col=seq_col,
            r2_col=r2_col,
            id_col=id_col,
            tag=tag,
            outspec=spec,
            params=match_params,
            cfg=cfg,
        )
        return match_el, count_el

    def samplecount(
        self,
        samples: Mapping[str, NextGenSample],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> dict[str, Element]:
        """
        samplecount is a convenience method that takes a mapping of sample names
        to NextGenSampleElements and runs the count method on each sample,
        returning a mapping of sample names to their corresponding count
        Elements. This allows for easy batch processing of multiple samples
        without having to manually call count for each one.

        Parameters
        ----------
        samples : Mapping[str, NextGenSample]
            A mapping of sample names to NextGenSample instances to be counted.
        tag : PartialElementTag | ElementTag | None, optional
            Optional tag override applied to all count elements, by default None.
        outspec : OutputSpec | None, optional
            Output specification for all count elements, by default None. If None,
            defaults to results/counts/<version>/<sample_name> for each sample.
        params : Params | None, optional
            Parameters for the count method, by default None.
        cfg : ExternalRunConfig | None, optional
            Configuration for external run, by default None.

        Returns
        -------
        dict[str, Element]
            A mapping of sample names to their corresponding count Elements.
        """
        count_elements = {}
        for sample_name, sample in samples.items():
            spec = self.default_output(sample.root).merge(outspec)
            counted = self.count(
                sample,
                tag=tag,
                outspec=spec,
                params=params,
                cfg=cfg,
            )
            count_elements[sample_name] = counted
        return count_elements


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _build_param_registry_as_mapping() -> dict:
    """Return a plain dict so External.__init__ can accept it.

    The actual ParamRegistry is set directly in MmFqCount.__init__ after the
    super().__init__() call.
    """
    return {}
