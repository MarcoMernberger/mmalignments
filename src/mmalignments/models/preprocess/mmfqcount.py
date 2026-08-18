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

import matplotlib.pyplot as plt  # type: ignore[import]
import numpy as np  # type: ignore[import]
import pandas as pd  # type: ignore[import]
from matplotlib.patches import Rectangle
from pandas import DataFrame, Series  # type: ignore[import]

from mmalignments.core.biopython import (
    analyze_read_with_predefined,
    # build_mutation_lookup,
    build_sequence_profile_lookup,
    create_aligner,
    dataframe_to_igv_bam,
)
from mmalignments.models.artifacts import (
    ArtifactSet,
    FileSpec,
    OutputSpec,
    TableArtifact,
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
from mmalignments.services.io import (
    parents,
    read_frame,
    save_figure,
    write_frame,
    write_frames,
)

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
        ----------st du
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

    @element
    def profile(
        self,
        predefined: Element | FileSource,
        wt: Element | FileSource,
        *,
        id_column: str = "mut_ID",
        sequence_column: str = "full_sequence",
        exon_column: str = "exon",
        profile_column: str = "profile",
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        Creates mutation profile for predefined signatures later used for classification
        of the reads.

        Parameters
        ----------
        predefined : Element | FileSource
            Element or file containing the predefined sequences.
        wt : Element | FileSource
            Element or file containing the wild-type sequence.
        id_column : str, optional
            Column name in the predefined TSV that contains the mutation identifiers.
            Default is "mut_ID".
        sequence_column : str, optional
            Column name in the predefined TSV that contains the full sequences.
            Default is "full_sequence".
        exon_column : str, optional
            Column name in the predefined TSV that contains the exon identifiers.
            Default is "exon".
        profile_column : str
            Column name in the output TSV that will contain the mutation profiles.
            Default is "profile".
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
            Element with the alignment results as a TSV file.
        """

        tag = from_prior(
            predefined.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MMFQCOUNT,
            state=State.MAP,
        )
        params = Params().update(params)

        pres = ()
        predefined_path = predefined.primary.resolve()
        wt_path = wt.primary.resolve()
        outfile = (
            outspec.path(tag.default_name)
            if outspec
            else predefined_path.parent / f"{tag.default_name}.tsv"
        )
        runner = self.create_predefined_profile(
            wt_path,
            predefined_path,
            outfile,
            id_column=id_column,
            exon_column=exon_column,
            sequence_column=sequence_column,
            profile_column=profile_column,
        )
        determinants = (id_column, sequence_column) + params.determinants()
        key, name = self.build_element_name(
            tag,
            "profile",
        )
        artifacts = ArtifactSet(
            TableArtifact(outfile),
            primary_name=outfile.suffix,
        )

        return Element(
            key,
            runner,
            tag=tag,
            determinants=tuple(determinants),
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    def create_predefined_profile(
        self,
        wt_path: Path,
        predefined_path: Path,
        output_path: Path,
        id_column: str = "mut_ID",
        exon_column: str = "exon",
        sequence_column: str = "full_sequence",
        profile_column: str = "profile",
    ) -> Runnable:
        def call():
            predefined = read_frame(predefined_path)
            wt = read_frame(wt_path)
            sequence_profile_lookup = {}
            for exon, group in wt.groupby(exon_column):
                wt_sequence = group[sequence_column].values[0]
                if len(group) > 1:
                    raise ValueError(
                        f"Multiple wild-type sequences found for exon {exon}. "
                        "Please ensure that the wild-type sequence file contains "
                        "only one sequence per exon."
                    )
                predefined_group = predefined[predefined[exon_column] == exon]
                predefined_sequences = dict(
                    zip(
                        predefined_group[id_column].astype(str),
                        predefined_group[sequence_column].astype(str),
                    )
                )
                sequence_profile_lookup.update(
                    build_sequence_profile_lookup(predefined_sequences, wt_sequence)
                )
            predefined[profile_column] = predefined["full_sequence"].map(
                sequence_profile_lookup
            )
            write_frame(predefined, output_path)

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "MmFqCount",
                    "create_predefined_profile",
                ),
            ).render(),
        )

    @element
    def aligncount(
        self,
        count: TableElement,
        wt: Element | FileSource,
        predefined: Element | FileSource,
        exon: str,
        *,
        top: int = 50,
        sequence_column: str = "full_sequence",
        id_column: str = "mut_ID",
        min_coverage: float = 0.8,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        Selects the top N sequences from a count table and aligns them to a
        wild-type sequence, in order to display what sort of sequences we
        actually have in the top N sequences. This is useful for visualizing the
        results of our counting. If a read corresponds to a predefined sequence,
        it will be annotated with the predefined sequence name, no alignment
        needed.

        Parameters
        ----------
        count : TableElement
            Element containing the counts TSV to process.
        wt : Element | FileSource
            Element or file containing the wild-type sequence.
        predefined : Element | FileSource
            Element or file containing the predefined sequences.
        top : int, optional
            Number of top sequences to select. Default is 50.
        exon : str
            Exon to select wt sequence from. Default is "Ex5".
        sequence_column : str, optional
            Column name in the predefined TSV that contains the full sequences.
            Default is "full_sequence".
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None, optional
            Optional output specification override.
        params : Params | None, optional
            Additional parameters, by default None
        cfg : ExternalRunConfig | None, optional
            Runtime configuration, by default None
        mode : str, optional
            Mode of alignment, by default "both"
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
            Element with the alignment results as a TSV file.
        """

        tag = from_prior(
            count.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MMFQCOUNT,
            state=State.MAP,
        )
        params = Params().update(params)

        pres = (count, predefined)
        determinants = params.determinants()
        outfile = (
            outspec.path(tag.default_name)
            if outspec
            else count.file.parent / f"{tag.default_name}.tsv"
        )
        runner = self.align_counts(
            outpath=outfile,
            frame_path=count.primary.resolve(),
            wt_path=wt.primary.resolve(),
            exon=exon,
            predefined_path=predefined.primary.resolve(),
            sequence_column=sequence_column,
            id_column=id_column,
            min_coverage=min_coverage,
            top=top,
        )

        key, name = self.build_element_name(
            tag,
            "aligncount",
        )
        artifacts = ArtifactSet(
            TableArtifact(outfile),
            primary_name=outfile.suffix,
        )

        return Element(
            key,
            runner,
            tag=tag,
            determinants=tuple(determinants),
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    def align_counts(
        self,
        outpath: Path,
        frame_path: Path,
        wt_path: Path,
        exon: str,
        predefined_path: Path,
        sequence_column: str = "full_sequence",
        id_column: str = "mut_ID",
        min_coverage: float = 0.8,
        top: int = 50,
    ) -> Runnable:

        def call():
            df = read_frame(frame_path)
            if top:
                df = df.nlargest(top, "Count")
            wt_df = read_frame(wt_path)
            wt = wt_df[wt_df["exon"] == exon][sequence_column].values[0]
            predefined = read_frame(predefined_path)
            predefined_sequences = dict(
                zip(
                    predefined[id_column],
                    predefined[sequence_column],
                )
            )
            aligner = create_aligner()
            # mutation_lookup = build_mutation_lookup(predefined_sequences, wt)

            results = []

            for _, row in df.iterrows():

                read = row["R1"]

                result = analyze_read_with_predefined(
                    read=read,
                    wt=wt,
                    predefined_sequences=predefined_sequences,
                    aligner=aligner,
                    min_coverage=min_coverage,
                )

                result["count"] = row["Count"]

                # ssresult = classify_against_predefined(result, mutation_lookup)

                results.append(result)

            analysis_df = pd.DataFrame(results)

            write_frame(analysis_df, outpath)

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "MmFqCount",
                    "align_counts",
                ),
                kwargs={
                    "outpath": str(outpath),
                    "frame_path": str(frame_path),
                    "wt_path": str(wt_path),
                    "exon": exon,
                    "predefined_path": str(predefined_path),
                    "sequence_column": sequence_column,
                    "min_coverage": min_coverage,
                },
            ).render(),
        )

    @element
    def bamaligncount(
        self,
        count: TableElement,
        *,
        reference_name: str = "WT",
        min_mapq: int = 60,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        tag = from_prior(
            count.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MMFQCOUNT,
            state=State.PLOT,
        )
        outspec = OutputSpec(outdir=count.file.parent).merge(outspec)
        params = Params().update(params)

        pres = (count,)
        determinants = params.determinants()
        outfile = (
            outspec.path(tag.default_name)
            if outspec
            else count.file.parent / f"{tag.default_name}.svg"
        )

        def run():
            df = read_frame(count.primary.resolve())
            return dataframe_to_igv_bam(
                df=df,
                output_dir=outspec.outdir or count.file.parent,
                sample_name=tag.default_name,
                reference_name=reference_name,
                min_mapq=min_mapq,
            )

        key, name = self.build_element_name(
            tag,
            "bamaligncount",
        )
        bam_unsorted = (
            outspec.outdir or count.file.parent
        ) / f"{tag.default_name}.unsorted.bam"
        bam = (outspec.outdir or count.file.parent) / f"{tag.default_name}.bam"
        bai = Path(str(bam) + ".bai")

        artifacts = (
            ArtifactSet(bam, primary_name="bam")
            .with_extra("bai", bai)
            .with_extra("unsorted_bam", bam_unsorted)
        )

        return Element(
            key,
            run,
            tag=tag,
            determinants=tuple(determinants),
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    @element
    def plotaligncount(
        self,
        count: TableElement,
        *,
        sort_by: str = "count",
        figsize=None,
        show_mutation_labels: bool = True,
        show_count: bool = True,
        show_coverage: bool = True,
        show_grid: bool = True,
        show_classification: bool = True,
        fontsize: int = 9,
        start: int | None = None,
        end: int | None = None,
        show_position_labels: bool = True,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        tag = from_prior(
            count.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MMFQCOUNT,
            state=State.PLOT,
        )
        params = Params().update(params)

        pres = (count,)
        determinants = params.determinants()
        outfile = (
            outspec.path(tag.default_name)
            if outspec
            else count.file.parent / f"{tag.default_name}.svg"
        )
        runner = self.plot_read_alignments(
            outfile=outfile,
            analysis_df_path=count.primary.resolve(),
            sort_by=sort_by,
            figsize=figsize,
            fontsize=fontsize,
            show_mutation_labels=show_mutation_labels,
            show_count=show_count,
            show_classification=show_classification,
            show_coverage=show_coverage,
            show_position_labels=show_position_labels,
            show_grid=show_grid,
            start=start,
            end=end,
        )

        key, name = self.build_element_name(
            tag,
            "plotaligncount",
        )
        artifacts = ArtifactSet(outfile, primary_name="svg")

        return Element(
            key,
            runner,
            tag=tag,
            determinants=tuple(determinants),
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    def plot_read_alignments(
        self,
        outfile: Path,
        analysis_df_path: Path,
        sort_by: str = "count",
        figsize=None,
        start: int | None = None,
        end: int | None = None,
        fontsize: int = 9,
        show_mutation_labels: bool = True,
        show_count: bool = True,
        show_classification: bool = True,
        show_coverage: bool = False,
        show_position_labels: bool = True,
        show_grid: bool = True,
    ) -> Runnable:
        """
        Plot aligned reads base-by-base against the WT sequence.

        Each nucleotide is displayed as a colored cell:

            A = adenine
            C = cytosine
            G = guanine
            T = thymine
            - = gap/deletion

        Parameters
        ----------
        analysis_df_path : Path
            DataFrame containing at least:

                wt
                wt_aligned
                read_aligned
                count

            Optional:

                classification
                coverage
                mutation_profile

        sort_by : str
            Column used for sorting, usually "count".

        start : int or None
            First WT position to display (1-based).

        end : int or None
            Last WT position to display (1-based).

        figsize : tuple or None
            Figure size.

        fontsize : int
            Font size for nucleotide letters.

        show_mutation_labels : bool
            Annotate mutations above the corresponding base.

        show_count : bool
            Show read count.

        show_classification : bool
            Show read classification.

        show_coverage : bool
            Show coverage percentage.

        show_position_labels : bool
            Show WT coordinate labels.

        Returns
        -------
        fig, ax
        """

        def call(
            figsize=figsize, sort_by=sort_by, start=start, end=end, fontsize=fontsize
        ):

            df = read_frame(analysis_df_path)

            if len(df) == 0:
                raise ValueError("analysis_df is empty.")

            # =========================================================
            # Validate
            # =========================================================

            required = {
                "wt",
                "wt_aligned",
                "read_aligned",
            }

            missing = required - set(df.columns)

            if missing:
                raise ValueError(f"Missing columns: {sorted(missing)}")

            # =========================================================
            # Sort
            # =========================================================

            if sort_by is not None:

                if sort_by not in df.columns:
                    raise ValueError(f"Column {sort_by!r} not found.")

                df = df.sort_values(
                    sort_by,
                    ascending=False,
                    kind="stable",
                )

            # =========================================================
            # WT
            # =========================================================

            wt = str(df.iloc[0]["wt"]).upper()

            wt_length = len(wt)

            # User coordinates are 1-based
            if start is None:
                start = 1

            if end is None:
                end = wt_length

            if start < 1 or end > wt_length or start > end:
                raise ValueError(
                    f"Invalid range: {start}-{end}. " f"WT length = {wt_length}"
                )

            # Convert to zero-based slicing
            plot_start = start - 1
            plot_end = end

            n_positions = end - start + 1

            # =========================================================
            # Figure
            # =========================================================

            if figsize is None:

                width = max(
                    12,
                    n_positions * 0.32,
                )

                height = max(
                    3,
                    (len(df) + 2) * 0.42,
                )

                figsize = (
                    width,
                    height,
                )

            fig, ax = plt.subplots(figsize=figsize)

            # =========================================================
            # Base colors
            # =========================================================

            base_colors = {
                "A": "#4C78A8",
                "C": "#59A14F",
                "G": "#F28E2B",
                "T": "#E15759",
                "-": "#BDBDBD",
                "N": "#9C9C9C",
            }

            # =========================================================
            # Draw cells
            # =========================================================

            cell_width = 1.0
            cell_height = 0.85

            # WT is at top
            wt_y = len(df)

            # ---------------------------------------------------------
            # WT sequence
            # ---------------------------------------------------------

            for x, base in enumerate(wt[plot_start:plot_end]):

                color = base_colors.get(
                    base,
                    "#CCCCCC",
                )

                rect = Rectangle(
                    (
                        x,
                        wt_y - cell_height / 2,
                    ),
                    cell_width,
                    cell_height,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.5,
                )

                ax.add_patch(rect)

                ax.text(
                    x + 0.5,
                    wt_y,
                    base,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    fontfamily="monospace",
                    fontweight="bold",
                    color="white",
                )

            # WT label
            ax.text(
                -0.5,
                wt_y,
                "WT",
                ha="right",
                va="center",
                fontweight="bold",
            )

            # =========================================================
            # Helper: aligned strings -> WT-position mapping
            # =========================================================

            def alignment_to_positions(
                wt_aligned,
                read_aligned,
            ):
                """
                Convert alignment columns into WT coordinates.

                Returns
                -------
                dict
                    WT coordinate (0-based) ->
                    read base / '-' / insertion information
                """

                if len(wt_aligned) != len(read_aligned):
                    raise ValueError(
                        "wt_aligned and read_aligned " "must have identical length."
                    )

                mapping = {}

                wt_pos = 0

                for wt_base, read_base in zip(
                    wt_aligned.upper(),
                    read_aligned.upper(),
                ):

                    # Insertion relative to WT
                    if wt_base == "-":

                        # Attach insertion to the preceding WT position
                        if wt_pos > 0:
                            mapping.setdefault(wt_pos - 1, {}).setdefault(
                                "insertions", []
                            ).append(read_base)

                        continue

                    # Normal WT position
                    mapping[wt_pos] = {
                        "wt": wt_base,
                        "read": read_base,
                        "insertions": [],
                    }

                    wt_pos += 1

                return mapping

            # =========================================================
            # Mutation helper
            # =========================================================

            def get_mutations(row):

                profile = row.get(
                    "mutation_profile",
                    None,
                )

                if profile is None:
                    return []

                # Current MutationProfile object
                if hasattr(profile, "mutations"):
                    return list(profile.mutations)

                # Already a list
                if isinstance(
                    profile,
                    (list, tuple),
                ):
                    return list(profile)

                return []

            # =========================================================
            # Reads
            # =========================================================

            for i, (_, row) in enumerate(df.iterrows()):

                y = len(df) - i - 1

                wt_aligned = str(row["wt_aligned"]).upper()

                read_aligned = str(row["read_aligned"]).upper()

                mapping = alignment_to_positions(
                    wt_aligned,
                    read_aligned,
                )

                # -----------------------------------------------------
                # Read label
                # -----------------------------------------------------

                label_parts = []

                if show_classification:
                    label_parts.append(
                        str(
                            row.get(
                                "classification",
                                "read",
                            )
                        )
                    )

                if show_coverage and "coverage" in row:
                    if pd.notna(row["coverage"]):
                        label_parts.append(f"{row['coverage'] * 100:.1f}%")

                if show_count:
                    label_parts.append(f"n={int(row.get('count', 1))}")

                label = "  ".join(label_parts)

                ax.text(
                    -0.5,
                    y,
                    label,
                    ha="right",
                    va="center",
                    fontsize=9,
                )

                # -----------------------------------------------------
                # Draw bases
                # -----------------------------------------------------

                for wt_pos in range(
                    plot_start,
                    plot_end,
                ):

                    x = wt_pos - plot_start

                    info = mapping.get(wt_pos)

                    if info is None:
                        continue

                    read_base = info["read"]

                    # ---------------------------------------------
                    # Deletion
                    # ---------------------------------------------

                    if read_base == "-":

                        base = "-"

                    else:

                        base = read_base

                    color = base_colors.get(
                        base,
                        "#CCCCCC",
                    )

                    rect = Rectangle(
                        (
                            x,
                            y - cell_height / 2,
                        ),
                        cell_width,
                        cell_height,
                        facecolor=color,
                        edgecolor="white",
                        linewidth=0.5,
                    )

                    ax.add_patch(rect)

                    ax.text(
                        x + 0.5,
                        y,
                        base,
                        ha="center",
                        va="center",
                        fontsize=fontsize,
                        fontfamily="monospace",
                        fontweight="bold",
                        color="white",
                    )

                    # ---------------------------------------------
                    # Insertions
                    # ---------------------------------------------

                    insertions = info.get(
                        "insertions",
                        [],
                    )

                    if insertions:

                        insertion_text = "".join(insertions)

                        ax.text(
                            x + 0.5,
                            y + 0.42,
                            f"+{insertion_text}",
                            ha="center",
                            va="bottom",
                            fontsize=max(
                                6,
                                fontsize - 2,
                            ),
                            fontweight="bold",
                        )

                # =====================================================
                # Mutation labels
                # =====================================================

                if show_mutation_labels:

                    mutations = get_mutations(row)

                    for mutation in mutations:

                        mutation_type = mutation.type

                        if mutation_type == "SNV":

                            pos = mutation.position

                            # Assuming mutation positions are 1-based
                            x = pos - start

                            if 0 <= x < n_positions:

                                ax.text(
                                    x + 0.5,
                                    y + 0.48,
                                    mutation.description,
                                    ha="center",
                                    va="bottom",
                                    fontsize=7,
                                    rotation=45,
                                )

                        elif mutation_type == "INS":

                            pos = mutation.position

                            x = pos - start

                            if 0 <= x < n_positions:

                                ax.text(
                                    x + 0.5,
                                    y + 0.48,
                                    mutation.description,
                                    ha="center",
                                    va="bottom",
                                    fontsize=7,
                                    rotation=45,
                                )

                        elif mutation_type == "DEL":

                            pos = mutation.position

                            x = pos - start

                            if 0 <= x < n_positions:

                                end_pos = getattr(
                                    mutation,
                                    "end",
                                    pos,
                                )

                                x2 = end_pos - start

                                ax.plot(
                                    [
                                        x + 0.1,
                                        x2 + 0.9,
                                    ],
                                    [
                                        y,
                                        y,
                                    ],
                                    linewidth=3,
                                    color="black",
                                )

            # =========================================================
            # Position axis
            # =========================================================

            if show_position_labels:

                tick_step = max(
                    1,
                    n_positions // 20,
                )

                ticks = list(
                    range(
                        0,
                        n_positions,
                        tick_step,
                    )
                )

                ax.set_xticks([x + 0.5 for x in ticks])

                ax.set_xticklabels(
                    [str(start + x) for x in ticks],
                    fontsize=8,
                )

                ax.set_xlabel("WT position")

            else:

                ax.set_xticks([])

            # =========================================================
            # Limits
            # =========================================================

            ax.set_xlim(
                0,
                n_positions,
            )

            ax.set_ylim(
                -1,
                len(df) + 1.3,
            )

            ax.set_yticks([])

            # =========================================================
            # Title
            # =========================================================

            ax.set_title(
                f"Read alignments — " f"WT positions {start}-{end}",
                fontsize=13,
                fontweight="bold",
            )

            # =========================================================
            # Clean axes
            # =========================================================

            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)

            ax.tick_params(
                axis="x",
                length=0,
            )

            # =========================================================
            # Legend
            # =========================================================

            legend_elements = []

            for base in ["A", "C", "G", "T", "-"]:

                legend_elements.append(
                    Rectangle(
                        (0, 0),
                        1,
                        1,
                        facecolor=base_colors[base],
                        edgecolor="none",
                        label=base,
                    )
                )

            ax.legend(
                handles=legend_elements,
                title="Base",
                loc="upper right",
                ncol=5,
                frameon=False,
            )

            plt.tight_layout()
            save_figure(fig, outfile)
            return fig, ax

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "MmFqCount",
                    "plot_read_alignments",
                ),
                kwargs={
                    "outfile": str(outfile),
                    "analysis_df_path": str(analysis_df_path),
                    "sort_by": sort_by,
                    "figsize": figsize,
                    "show_mutation_labels": show_mutation_labels,
                    "show_coverage": show_coverage,
                    "show_grid": show_grid,
                },
            ).render(),
        )

    def plot_read_alignments_first(
        self,
        outfile: Path,
        analysis_df_path: Path,
        sort_by: str = "count",
        figsize=None,
        show_mutation_labels: bool = True,
        show_coverage: bool = True,
        show_grid: bool = True,
    ) -> Runnable:
        """
        Plot read alignments relative to WT coordinates.

        Parameters
        ----------
        analysis_df_path : Path
            Path to the DataFrame returned by the read-analysis workflow.

            Expected columns:
                - wt
                - wt_start / wt_end
                - read_start / read_end
                - count
                - coverage
                - classification
                - mutation_profile

        sort_by : str
            Column used for sorting. Usually "count".

        figsize : tuple or None
            Matplotlib figure size. If None, calculated automatically.

        show_mutation_labels : bool
            Show mutation descriptions next to mutations.

        show_coverage : bool
            Show coverage percentage in the read label.

        show_grid : bool
            Show vertical WT coordinate grid.
        """

        def call(figsize=figsize, sort_by=sort_by):
            df = read_frame(analysis_df_path)
            if len(df) == 0:
                raise ValueError("analysis_df is empty.")

            # ---------------------------------------------------------
            # Sort
            # ---------------------------------------------------------

            if sort_by is not None:
                if sort_by not in df.columns:
                    raise ValueError(f"Column {sort_by!r} not found in analysis_df.")

                df = df.sort_values(
                    sort_by,
                    ascending=False,
                    kind="stable",
                )

            # ---------------------------------------------------------
            # WT length
            # ---------------------------------------------------------

            wt_length = len(df.iloc[0]["wt"])

            # ---------------------------------------------------------
            # Figure size
            # ---------------------------------------------------------

            if figsize is None:
                figsize = (
                    16,
                    max(4, 0.55 * (len(df) + 2)),
                )

            fig, ax = plt.subplots(figsize=figsize)

            # ---------------------------------------------------------
            # Y positions
            # ---------------------------------------------------------

            y_wt = len(df)

            # ---------------------------------------------------------
            # WT reference
            # ---------------------------------------------------------

            ax.plot(
                [1, wt_length],
                [y_wt, y_wt],
                linewidth=7,
                solid_capstyle="butt",
            )

            ax.text(
                -0.02 * wt_length,
                y_wt,
                "WT",
                ha="right",
                va="center",
                fontweight="bold",
            )

            # ---------------------------------------------------------
            # Helper: extract mutations
            # ---------------------------------------------------------

            def get_mutations(row):

                profile = row.get("mutation_profile", None)

                if profile is None:
                    return []

                # Your current MutationProfile dataclass
                if hasattr(profile, "mutations"):
                    return list(profile.mutations)

                # In case the profile is already a list/tuple
                if isinstance(profile, (list, tuple)):
                    return list(profile)

                return []

            # ---------------------------------------------------------
            # Reads
            # ---------------------------------------------------------

            for i, (_, row) in enumerate(df.iterrows()):

                y = len(df) - i - 1

                wt_start = row.get("wt_start")
                wt_end = row.get("wt_end")

                if wt_start is None or wt_end is None:
                    continue

                # -----------------------------------------------------
                # Read coverage
                # -----------------------------------------------------

                ax.plot(
                    [wt_start, wt_end],
                    [y, y],
                    linewidth=7,
                    solid_capstyle="butt",
                    alpha=0.8,
                )

                # -----------------------------------------------------
                # Read label
                # -----------------------------------------------------

                classification = row.get(
                    "classification",
                    "",
                )

                count = row.get(
                    "count",
                    1,
                )

                if show_coverage and "coverage" in row:
                    coverage = row["coverage"] * 100

                    label = (
                        f"{classification}   "
                        f"coverage={coverage:.1f}%   "
                        f"n={count}"
                    )
                else:
                    label = f"{classification}   " f"n={count}"

                ax.text(
                    -0.02 * wt_length,
                    y,
                    label,
                    ha="right",
                    va="center",
                    fontsize=9,
                )

                # -----------------------------------------------------
                # Mutations
                # -----------------------------------------------------

                mutations = get_mutations(row)

                for mutation in mutations:

                    mutation_type = mutation.type

                    # -------------------------------------------------
                    # SNV
                    # -------------------------------------------------

                    if mutation_type == "SNV":

                        pos = mutation.position

                        ax.scatter(
                            pos,
                            y,
                            s=55,
                            zorder=5,
                        )

                        if show_mutation_labels:

                            ax.text(
                                pos,
                                y + 0.18,
                                mutation.description,
                                ha="center",
                                va="bottom",
                                fontsize=8,
                                rotation=45,
                            )

                    # -------------------------------------------------
                    # INSERTION
                    # -------------------------------------------------

                    elif mutation_type == "INS":

                        pos = mutation.position

                        ax.scatter(
                            pos,
                            y,
                            marker="^",
                            s=70,
                            zorder=5,
                        )

                        if show_mutation_labels:

                            ax.text(
                                pos,
                                y + 0.18,
                                mutation.description,
                                ha="center",
                                va="bottom",
                                fontsize=8,
                                rotation=45,
                            )

                    # -------------------------------------------------
                    # DELETION
                    # -------------------------------------------------

                    elif mutation_type == "DEL":

                        start = mutation.position
                        end = mutation.end

                        ax.plot(
                            [start, end],
                            [y, y],
                            linewidth=10,
                            alpha=0.5,
                            solid_capstyle="butt",
                            zorder=4,
                        )

                        if show_mutation_labels:

                            ax.text(
                                (start + end) / 2,
                                y + 0.18,
                                mutation.description,
                                ha="center",
                                va="bottom",
                                fontsize=8,
                            )

            # ---------------------------------------------------------
            # Axes
            # ---------------------------------------------------------

            ax.set_xlim(
                1,
                wt_length,
            )

            ax.set_ylim(
                -1,
                len(df) + 1,
            )

            ax.set_xlabel(
                "WT position",
                fontsize=11,
            )

            ax.set_yticks([])

            # ---------------------------------------------------------
            # Grid
            # ---------------------------------------------------------

            if show_grid:

                ax.grid(
                    axis="x",
                    linestyle=":",
                    alpha=0.4,
                )

            # ---------------------------------------------------------
            # Title
            # ---------------------------------------------------------

            ax.set_title(
                f"Read alignments ({len(df)} reads)",
                fontsize=13,
                fontweight="bold",
            )

            # ---------------------------------------------------------
            # Clean frame
            # ---------------------------------------------------------

            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)

            plt.tight_layout()
            save_figure(fig, outfile)
            return fig, ax

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "MmFqCount",
                    "plot_read_alignments",
                ),
                kwargs={
                    "outfile": str(outfile),
                    "analysis_df_path": str(analysis_df_path),
                    "sort_by": sort_by,
                    "figsize": figsize,
                    "show_mutation_labels": show_mutation_labels,
                    "show_coverage": show_coverage,
                    "show_grid": show_grid,
                },
            ).render(),
        )


# def align_reads_to_wt(
#     self,
#     df: DataFrame,
#     predefined: DataFrame,
#     wt_df: DataFrame,
# ) -> DataFrame:
#     """
#     Align all unique sequences in a dataframe.

#     Required columns:
#         R1
#         Count

#     Optional columns:
#         Frequency
#         Annotation
#     """
#     sequence_col = "R1"
#     results = []

#     for _, row in df.iterrows():
#         if row[sequence_col] in predefined["sequence"].values:
#             # If the sequence is in the predefined sequences, we can skip alignment
#             annotation = predefined.loc[
#                 predefined["sequence"] == row[sequence_col], "mut_id"
#             ].values[0]
#             results.append(
#                 {
#                     "sequence": row[sequence_col],
#                     "count": row["Count"],
#                     "frequency": row.get("Frequency", np.nan),
#                     "annotation": annotation,
#                     "aligned_read": row[sequence_col],
#                     "wt_start": 0,
#                     "wt_end": len(row[sequence_col]),
#                     "aligned_length": len(row[sequence_col]),
#                     "alignment_score": np.nan,
#                     "insertions": [],
#                 }
#             )
#             continue
#         if row[sequence_col] in wt_df["sequence"].values:
#             # If the sequence is in the wild-type sequences, we can skip alignment
#             wt_row = wt_df.loc[wt_df["sequence"] == row[sequence_col]].iloc[0]
#             results.append(
#                 {
#                     "sequence": row[sequence_col],
#                     "count": row["Count"],
#                     "frequency": row.get("Frequency", np.nan),
#                     "annotation": "WT",
#                     "aligned_read": row[sequence_col],
#                     "wt_start": 0,
#                     "wt_end": len(row[sequence_col]),
#                     "aligned_length": len(row[sequence_col]),
#                     "alignment_score": np.nan,
#                     "insertions": [],
#                 }
#             )
#             continue
#         best_score = -np.inf
#         result = None
#         for _, wtrow in wt_df.iterrows():
#             wt_sequence = wtrow["sequence"]
#             flank_5 = wtrow.get("5_flank", "")
#             flank_3 = wtrow.get("3_flank", "")
#             tmp = align_read_to_wt(
#                 read=row[sequence_col],
#                 wt_sequence=wt_sequence,
#                 flank_5=flank_5,
#                 flank_3=flank_3,
#                 name=row.get("name", ""),
#             )
#             if tmp["alignment_score"] > best_score:
#                 best_score = tmp["alignment_score"]
#                 result = tmp

#         results.append(
#             {
#                 "sequence": row[sequence_col],
#                 "count": row["Count"],
#                 "frequency": row.get("Frequency", np.nan),
#                 "annotation": result["annotation"],
#                 "aligned_read": result["aligned_read"],
#                 "wt_start": result["wt_start"],
#                 "wt_end": result["wt_end"],
#                 "aligned_length": result["aligned_length"],
#                 "alignment_score": result["alignment_score"],
#                 "insertions": result["insertions"],
#                 "has_flank_5": result["has_flank_5"],
#                 "has_flank_3": result["has_flank_3"],
#                 "flank_5_start": result["flank_5_start"],
#                 "flank_5_end": result["flank_5_end"],
#                 "flank_3_start": result["flank_3_start"],
#                 "flank_3_end": result["flank_3_end"],
#                 "flank_annotation": result["flank_annotation"],
#             }
#         )

#     wt_row = pd.DataFrame(
#         [
#             {
#                 "sequence": wt_sequence,
#                 "count": np.nan,
#                 "frequency": np.nan,
#                 "annotation": "WT",
#                 "aligned_read": wt_sequence,
#                 "wt_start": 0,
#                 "wt_end": len(wt_sequence),
#                 "aligned_length": len(wt_sequence),
#                 "alignment_score": np.nan,
#                 "insertions": [],
#             }
#         ]
#     )
#     alignment_df = pd.DataFrame(results)
#     alignment_df = pd.concat([wt_row, alignment_df], ignore_index=True)
#     return alignment_df

# def plot_sequence_alignment(alignment_df, wt_sequence, top_n=20, figsize=(14, 8)):
#     """
#     Plot sequencing variants aligned to a WT sequence.

#     Parameters
#     ----------
#     alignment_df : pandas.DataFrame
#         Output from align_reads_to_wt().
#     wt_sequence : str
#         WT/reference sequence.
#     top_n : int
#         Number of variants to display.
#     figsize : tuple
#         Figure size.
#     """

#     # Separate WT and reads
#     wt_length = len(wt_sequence)

#     reads = alignment_df[alignment_df["annotation"] != "WT"].copy()

#     # Sort by abundance
#     reads = reads.sort_values("count", ascending=False).head(top_n)

#     # WT first
#     plot_rows = []

#     plot_rows.append(
#         {
#             "label": "WT",
#             "sequence": wt_sequence,
#             "aligned_read": wt_sequence,
#             "count": np.nan,
#             "annotation": "WT",
#         }
#     )

#     for i, (_, row) in enumerate(reads.iterrows()):

#         plot_rows.append(
#             {
#                 "label": f"Variant {i+1}",
#                 "sequence": row["sequence"],
#                 "aligned_read": row["aligned_read"],
#                 "count": row["count"],
#                 "annotation": row["annotation"],
#             }
#         )

#     fig, ax = plt.subplots(figsize=figsize)

#     row_height = 0.6

#     for y, row in enumerate(plot_rows):

#         aligned = row["aligned_read"]

#         # Draw sequence positions
#         for x, base in enumerate(aligned):

#             if base == ".":
#                 continue

#             if base == "-":
#                 # deletion
#                 ax.plot([x, x + 1], [y, y], linewidth=3)

#             else:
#                 # sequence present
#                 ax.add_patch(Rectangle((x, y - row_height / 2), 1, row_height))

#     # Labels
#     labels = []

#     for row in plot_rows:

#         if pd.isna(row["count"]):
#             labels.append(row["label"])
#         else:
#             labels.append(f'{row["label"]} ({int(row["count"]):,})')

#     ax.set_yticks(range(len(plot_rows)))
#     ax.set_yticklabels(labels)

#     ax.set_xlim(0, wt_length)
#     ax.invert_yaxis()

#     ax.set_xlabel("Position in WT sequence")
#     ax.set_ylabel("Sequence")

#     ax.set_title("Sequencing variants aligned to WT")

#     plt.tight_layout()

#     return fig, ax


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _build_param_registry_as_mapping() -> dict:
    """Return a plain dict so External.__init__ can accept it.

    The actual ParamRegistry is set directly in MmFqCount.__init__ after the
    super().__init__() call.
    """
    return {}


class AmpliconQC(External):
    """Wrapper for the ``amplicon-wc`` Rust CLI tool.

    ``amplicon-wc`` is a Rust command-line tool for counting the lengths of
    inserts between flanking sequences in FASTQ files.
    """

    def __init__(
        self,
        name: str = "amplicon-qc",
        primary_binary: str = "amplicon-qc",
        version: str | None = None,
        source: str = "https://github.com/MarcoMernberger/ampolicon-qc",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialise the AmpliconQC wrapper.

        Parameters
        ----------
        name : str
            Logical tool name, default ``"amplicon-qc"``.
        primary_binary : str
            Executable name on ``$PATH``, default ``"amplicon-qc"``.
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
        self.param_registry = _build_param_registry_amplicon_qc()

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
        return Path("results") / "qc" / self.version_name / sample_name

    def default_output(
        self, sample_name: str = "", outdir: Path | None = None, ext: str = "tsv"
    ) -> OutputSpec:
        """Return the default output spec for a given sample."""
        return OutputSpec(
            outdir=outdir or self.default_output_dir(sample_name),
            ext=ext,
        )

    # -----------------------------------------------------------------------
    # amplicon-wc
    # -----------------------------------------------------------------------

    @element
    def flankinserts(
        self,
        sample: Element | NextGenSample,
        flank_start: FileSource,
        flank_end: FileSource,
        *,
        max_hamming: int = 1,
        limit: int = 10000,
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Counts the read lengths between flanks, in other words, what remains
        after trimming for trouvleshooting.

        Runs ``amplicon-qc`` on the FASTQ file(s) of *sample* and writes
        a TSV table with lengths and categories.

        Parameters
        ----------
        sample : Element | NextGenSample
            Input sample. Both single-end and paired-end samples are handled
            automatically.
        flank_start : FileSource
            5' flank sequence file (FASTA).
        flank_end : FileSource
            3' flank sequence file (FASTA).
        max_hamming : int
            Maximum Hamming distance for flank matching. Default: 1.
        limit : int
            Maximum number of reads to process. Default: 10000.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None
            Output specification; defaults to
            ``results/counts/<version>/<sample.name>`` with extension ``tsv``.
        params : Params | None
            Additional parameters for the count step. Default: None.

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
            stage=Stage.QC,
            method=Method.MMFQCOUNT,
            state=State.COUNT,
            ext="tsv",
        )

        spec = self.default_output(sample.root).merge(outspec)
        output_file = spec.path(tag.default_name)
        len_file = output_file.with_name(output_file.stem + "_lengths.tsv")
        categories_file = output_file.with_name(output_file.stem + "_categories.tsv")

        runner = self.check_flank_insert(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            flank_start=flank_start.artifacts.primary.resolve(),
            flank_end=flank_end.artifacts.primary.resolve(),
            max_hamming=max_hamming,
            limit=limit,
            output_file=output_file,
            len_file=len_file,
            categories_file=categories_file,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(tag, "count")
        determinants = self.signature_determinants(params, subroutine="count")
        inputs = (fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,)
        pres = (sample,) if isinstance(sample, Element) else sample.pres
        artifacts = ArtifactSet(output_file, primary_name=spec.ext)
        artifacts = artifacts.with_extra("lengths", value=len_file).with_extra(
            "categories", value=categories_file
        )
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
    def check_flank_insert(
        self,
        output_file: Path,
        fastq_r1: Path,
        fastq_r2: Path,
        flank_start: Path,
        flank_end: Path,
        *,
        len_file: Path | None = None,
        categories_file: Path | None = None,
        max_hamming: int = 1,
        limit: int = 10000,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``amplicon-qc``.

        Parameters
        ----------
        fastq_r1 : Path | str
            R1 FASTQ file (plain or gzip).
        fastq_r2 : Path | str | None
            R2 FASTQ file (plain or gzip), or *None* for single-end.
        flank_start : Path | str
            Path to the 5' flank sequence file (FASTA).
        flank_end : Path | str
            Path to the 3' flank sequence file (FASTA).
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
        len_file = len_file or output_file.with_name(output_file.stem + "_lengths.tsv")
        categories_file = categories_file or output_file.with_name(
            output_file.stem + "_categories.tsv"
        )
        arguments = [
            "--r1",
            str(fastq_r1),
            "--start-flanks",
            str(flank_start),
            "--end-flanks",
            str(flank_end),
            "--max-hamming",
            str(max_hamming),
            "--output",
            str(output_file),
            "--histogram",
            len_file,
            "--categories",
            categories_file,
            "--limit",
            str(limit),
        ]
        if fastq_r2 is not None:
            arguments += ["--r2", str(fastq_r2)]
        in_paths = [fastq_r1] + ([fastq_r2] if fastq_r2 else [])
        out_paths = [output_file, len_file, categories_file]

        return (
            arguments,
            "count",
            in_paths,
            out_paths,
            None,  # no piped output
            None,  # no pre-hook
            None,  # no post-hook
        )

    @element
    def plotflankdist(
        self,
        histogram: TableElement,
        *,
        figsize: dict[str, tuple[int, int]] = {
            "hist": (14, 7),
            "cat": (10, 6),
            "pos": (10, 6),
        },
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
    ) -> Element:
        """
        Plots the histogram calculated by amplicon-qc.

        Parameters
        ----------
        histogram : TableElement
            Element containing the histogram TSV to process.
        figsize : tuple[int, int], optional
            Size of the figure. Default is (14, 7).
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None, optional
            Optional output specification override.
        params : Params | None, optional
            Additional parameters, by default None

        Returns
        -------
        Element
            Element with the histogram plot as an SVG file.
        """
        outspec = OutputSpec(
            outdir=histogram.primary.resolve().parent, ext="svg"
        ).merge(outspec)
        tag = from_prior(
            histogram.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.MMFQCOUNT,
            state=State.PLOT,
        )
        params = Params().update(params)

        pres = (histogram,)
        determinants = params.determinants()
        outfile = (
            outspec.path(tag.default_name)
            if outspec
            else histogram.file.parent / f"{tag.default_name}.svg"
        )
        cat_file = outfile.with_name(outfile.stem + "_categories.svg")
        hist_file = outfile.with_name(outfile.stem + "_lengths.svg")
        runner = self.plot_flank_lengths(
            output_hist_file=hist_file,
            output_category_file=cat_file,
            output_position_file=outfile,
            position_file=histogram.primary.resolve(),
            histogram_path=histogram.artifacts["lengths"].resolve(),
            category_path=histogram.artifacts["categories"].resolve(),
            figsize=figsize,
        )

        key, name = self.build_element_name(
            tag,
            "plotflankdist",
        )
        artifacts = ArtifactSet(
            outfile,
            primary_name=outfile.suffix,
        ).with_extra("categories", value=cat_file)

        return Element(
            key,
            runner,
            tag=tag,
            determinants=tuple(determinants),
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    def plot_flank_lengths(
        self,
        output_hist_file: Path,
        output_category_file: Path,
        output_position_file: Path,
        position_file: Path,
        histogram_path: Path,
        category_path: Path,
        figsize: dict[str, tuple[int, int]] = {
            "hist": (14, 7),
            "cat": (10, 6),
            "pos": (10, 6),
        },
    ) -> Runnable:

        def call():

            plot_histogram(
                histogram_file=str(histogram_path),
                output_svg=str(output_hist_file),
                figsize=figsize["hist"],
            )
            plot_pair_categories(
                categories_file=str(category_path),
                output_svg=str(output_category_file),
                figsize=figsize["cat"],
            )
            plot_flank_positions(
                tsv_file=str(position_file),
                output_svg=str(output_position_file),
                figsize=figsize["pos"],
            )

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "MmFqCount",
                    "align_counts",
                ),
                kwargs={
                    "output_hist_file": str(output_hist_file),
                    "output_category_file": str(output_category_file),
                    "output_position_file": str(output_position_file),
                    "position_file": str(position_file),
                    "histogram_path": str(histogram_path),
                    "category_path": str(category_path),
                    "figsize": {
                        "hist": figsize["hist"],
                        "cat": figsize["cat"],
                        "pos": figsize["pos"],
                    },
                },
            ).render(),
        )


def _build_param_registry_amplicon_qc() -> ParamRegistry:
    """Parameter registry for amplicon-qc.

    Defined here in Python because we own and maintain the Rust source —
    the CLI is stable and there is no need to parse --help output or maintain
    a separate JSON file.  If the CLI ever grows new flags, add them here.
    """
    specs: dict[str, ParamSpec] = {
        "r1": ParamSpec(
            "r1",
            "--r1",
            str,
            render=render_value,
            description="R1 read file.",
        ),
        "r2": ParamSpec(
            "r2",
            "--r2",
            str,
            render=render_value,
            description="R2 read file.",
        ),
        "start_flanks": ParamSpec(
            "start-flanks",
            "--start-flanks",
            str,
            render=render_value,
            description="File with start flanking sequences.",
        ),
        "end_flanks": ParamSpec(
            "end-flanks",
            "--end-flanks",
            str,
            render=render_value,
            description="File with end flanking sequences.",
        ),
        "max_hamming": ParamSpec(
            "max-hamming",
            "--max-hamming",
            int,
            render=render_value,
            description="Maximum allowed Hamming distance.",
        ),
        "output": ParamSpec(
            "output",
            "--output",
            str,
            render=render_value,
            description="Per-read TSV output.",
        ),
        "histogram": ParamSpec(
            "histogram",
            "--histogram",
            str,
            render=render_value,
            description="Aggregated histogram output.",
        ),
        "categories": ParamSpec(
            "categories",
            "--categories",
            str,
            render=render_value,
            description="Pair-category count output.",
        ),
        "limit": ParamSpec(
            "limit",
            "--limit",
            int,
            render=render_value,
            description="Optional maximum number of read pairs to process.",
        ),
    }
    return ParamRegistry(default=ParamSet(specs, "amplicon-qc", "default"))


################################################################################


def plot_histogram(
    histogram_file: str | Path,
    output_svg: str | Path,
    figsize=(14, 7),
    include_none_none=False,
):
    """
    Plot paired R1/R2 length distributions from the pairwise
    histogram TSV.

    Input format:

        pair_category    r1_length    r2_length    count

    The input contains counts for individual R1/R2 length pairs.

    Example:

        end_only__end_only    23    22    10
        end_only__end_only    23    23    250
        end_only__end_only    24    23    40
        end_only__end_only    24    24    220

    The figure contains two panels:

        LEFT  = R1 length distribution
        RIGHT = R2 length distribution

    Each pair_category gets one color, which is used consistently
    in both panels.

    The R1 and R2 lengths are NOT added together.

    Since the input contains pairwise R1/R2 combinations, the
    marginal distributions are calculated by summing count over
    the opposite read:

        R1 distribution:
            sum over r2_length

        R2 distribution:
            sum over r1_length

    By default, the "none__none" category is excluded from the
    plotted distributions because it usually represents complete
    reads and can have a much larger count than the flank-positive
    categories. This prevents the other distributions from being
    visually compressed.

    Set include_none_none=True to include it.
    """

    histogram_file = Path(histogram_file)
    output_svg = Path(output_svg)

    # ------------------------------------------------------------
    # Read histogram
    # ------------------------------------------------------------

    df = pd.read_csv(
        histogram_file,
        sep="\t",
    )

    required_columns = {
        "pair_category",
        "r1_length",
        "r2_length",
        "count",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns in histogram file: {sorted(missing)}")

    # ------------------------------------------------------------
    # Ensure numeric columns have the correct type
    # ------------------------------------------------------------

    df["r1_length"] = pd.to_numeric(
        df["r1_length"],
        errors="raise",
    )

    df["r2_length"] = pd.to_numeric(
        df["r2_length"],
        errors="raise",
    )

    df["count"] = pd.to_numeric(
        df["count"],
        errors="raise",
    )

    # ------------------------------------------------------------
    # Categories
    # ------------------------------------------------------------

    all_categories = sorted(df["pair_category"].unique())

    # Keep track of none__none separately.
    none_none_category = "none__none"

    if include_none_none:
        categories = all_categories
    else:
        categories = [
            category for category in all_categories if category != none_none_category
        ]

    # ------------------------------------------------------------
    # One color per category
    # ------------------------------------------------------------

    cmap = plt.get_cmap("tab20")

    category_colors = {
        category: cmap(i % cmap.N) for i, category in enumerate(all_categories)
    }

    # ------------------------------------------------------------
    # Calculate marginal R1/R2 distributions
    # ------------------------------------------------------------

    r1_distribution = df.groupby(
        [
            "pair_category",
            "r1_length",
        ],
        as_index=False,
    )["count"].sum()

    r2_distribution = df.groupby(
        [
            "pair_category",
            "r2_length",
        ],
        as_index=False,
    )["count"].sum()

    # ------------------------------------------------------------
    # Calculate none__none information
    # ------------------------------------------------------------

    none_none_count = df.loc[
        df["pair_category"] == none_none_category,
        "count",
    ].sum()

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    fig, (ax_r1, ax_r2) = plt.subplots(
        1,
        2,
        figsize=figsize,
        sharey=True,
    )

    # ------------------------------------------------------------
    # R1
    # ------------------------------------------------------------

    for category in categories:

        subset = r1_distribution[
            r1_distribution["pair_category"] == category
        ].sort_values("r1_length")

        ax_r1.plot(
            subset["r1_length"],
            subset["count"],
            label=category,
            color=category_colors[category],
            linewidth=2,
        )

    ax_r1.set_title("R1")
    ax_r1.set_xlabel("R1 observed length")
    ax_r1.set_ylabel("Read-pair count")

    ax_r1.grid(
        True,
        alpha=0.25,
    )

    # ------------------------------------------------------------
    # R2
    # ------------------------------------------------------------

    for category in categories:

        subset = r2_distribution[
            r2_distribution["pair_category"] == category
        ].sort_values("r2_length")

        ax_r2.plot(
            subset["r2_length"],
            subset["count"],
            label=category,
            color=category_colors[category],
            linewidth=2,
        )

    ax_r2.set_title("R2")
    ax_r2.set_xlabel("R2 observed length")

    ax_r2.grid(
        True,
        alpha=0.25,
    )

    # ------------------------------------------------------------
    # Shared title
    # ------------------------------------------------------------

    if include_none_none:
        title = "Observed length distributions by paired-read category"
    else:
        title = (
            "Observed length distributions by paired-read category"
            "\n"
            f"none__none excluded from plot ({none_none_count:,} read pairs)"
        )

    fig.suptitle(
        title,
        fontsize=16,
    )

    # ------------------------------------------------------------
    # One shared legend
    # ------------------------------------------------------------

    handles, labels = ax_r1.get_legend_handles_labels()

    if handles:
        fig.legend(
            handles,
            labels,
            title="Pair category",
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
        )

    # ------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------

    fig.tight_layout(
        rect=(0, 0, 0.84, 0.92),
    )

    # ------------------------------------------------------------
    # Save SVG
    # ------------------------------------------------------------

    output_svg.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_svg,
        format="svg",
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Wrote plot: {output_svg}")


def plot_pair_categories(
    categories_file: str | Path,
    output_svg: str | Path,
    figsize=(10, 6),
):
    """
    Plot read-pair category counts as a bar chart.

    Input format:

        pair_category    count

    Example:

        end_only__end_only    5885
        end_only__none        2847

    The resulting figure is saved as SVG.
    """

    categories_file = Path(categories_file)
    output_svg = Path(output_svg)

    # ------------------------------------------------------------
    # Read category counts
    # ------------------------------------------------------------

    df = pd.read_csv(
        categories_file,
        sep="\t",
    )

    required_columns = {
        "pair_category",
        "count",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns in categories file: {sorted(missing)}")

    df["count"] = pd.to_numeric(
        df["count"],
        errors="raise",
    )

    # Sort by count, largest first.
    df = df.sort_values(
        "count",
        ascending=False,
    )

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=figsize,
    )

    ax.bar(
        df["pair_category"],
        df["count"],
    )

    ax.set_title(
        "Read-pair categories",
        fontsize=16,
    )

    ax.set_xlabel(
        "Pair category",
    )

    ax.set_ylabel(
        "Read-pair count",
    )

    ax.grid(
        axis="y",
        alpha=0.25,
    )

    # Rotate category labels because they can be long.
    plt.setp(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
    )

    # Add counts above bars.
    for i, count in enumerate(df["count"]):
        ax.text(
            i,
            count,
            f"{count:,}",
            ha="center",
            va="bottom",
        )

    fig.tight_layout()

    # ------------------------------------------------------------
    # Save SVG
    # ------------------------------------------------------------

    output_svg.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_svg,
        format="svg",
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Wrote plot: {output_svg}")


def plot_flank_positions(
    tsv_file: str | Path,
    output_svg: str | Path,
    figsize=(14, 10),
):
    """
    Plot flank positions from the per-read TSV.

    The figure contains four panels:

        top-left     = R1 start-flank positions
        top-right    = R2 start-flank positions
        bottom-left  = R1 end-flank positions
        bottom-right = R2 end-flank positions

    Each pair_category gets its own color, consistently across
    all four panels.

    Only reads where the corresponding flank was actually found
    are included in the respective panel.
    """

    tsv_file = Path(tsv_file)
    output_svg = Path(output_svg)

    # ------------------------------------------------------------
    # Read TSV
    # ------------------------------------------------------------

    df = pd.read_csv(
        tsv_file,
        sep="\t",
    )

    required_columns = {
        "pair_category",
        "r1_start_position",
        "r1_end_position",
        "r2_rev_end_position",
        "r2_rev_start_position",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns in TSV: {sorted(missing)}")

    # Positions are empty when the flank was not found.
    position_columns = [
        "r1_start_position",
        "r1_end_position",
        "r2_rev_end_position",
        "r2_rev_start_position",
    ]

    for column in position_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    # ------------------------------------------------------------
    # Categories
    # ------------------------------------------------------------

    categories = sorted(df["pair_category"].dropna().unique())

    cmap = plt.get_cmap("tab20")

    category_colors = {
        category: cmap(i % cmap.N) for i, category in enumerate(categories)
    }

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    fig, axes = plt.subplots(
        2,
        2,
        figsize=figsize,
        sharex=False,
        sharey=False,
    )

    ax_r1_start = axes[0, 0]
    ax_r2_start = axes[0, 1]
    ax_r1_end = axes[1, 0]
    ax_r2_end = axes[1, 1]

    # ------------------------------------------------------------
    # Helper
    # ------------------------------------------------------------

    def plot_position_distribution(
        ax,
        position_column,
        title,
    ):
        for category in categories:

            subset = df[df["pair_category"] == category][position_column].dropna()

            if subset.empty:
                continue

            counts = subset.value_counts().sort_index()

            ax.plot(
                counts.index,
                counts.values,
                marker="o",
                markersize=3,
                linewidth=1.8,
                color=category_colors[category],
            )

        ax.set_title(title)
        ax.set_xlabel("Flank position")
        ax.set_ylabel("Read count")

        ax.grid(
            True,
            alpha=0.25,
        )

    # ------------------------------------------------------------
    # R1 start
    # ------------------------------------------------------------

    plot_position_distribution(
        ax_r1_start,
        "r1_start_position",
        "R1 start-flank position",
    )

    # ------------------------------------------------------------
    # R2 start
    # ------------------------------------------------------------

    plot_position_distribution(
        ax_r2_start,
        "r2_rev_end_position",
        "R2 end-flank position",
    )

    # ------------------------------------------------------------
    # R1 end
    # ------------------------------------------------------------

    plot_position_distribution(
        ax_r1_end,
        "r1_end_position",
        "R1 end-flank position",
    )

    # ------------------------------------------------------------
    # R2 end
    # ------------------------------------------------------------

    plot_position_distribution(
        ax_r2_end,
        "r2_rev_start_position",
        "R2 start-flank position",
    )

    # ------------------------------------------------------------
    # Shared title
    # ------------------------------------------------------------

    fig.suptitle(
        "Flank position distributions by paired-read category",
        fontsize=16,
    )

    # ------------------------------------------------------------
    # Shared legend
    # ------------------------------------------------------------

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            color=category_colors[category],
            linewidth=1.8,
            marker="o",
            markersize=3,
            label=category,
        )
        for category in categories
    ]

    if legend_handles:
        fig.legend(
            handles=legend_handles,
            title="Pair category",
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
        )

    # ------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------

    fig.tight_layout(
        rect=(0, 0, 0.84, 0.95),
    )

    # ------------------------------------------------------------
    # Save SVG
    # ------------------------------------------------------------

    output_svg.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_svg,
        format="svg",
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Wrote flank-position plot: {output_svg}")
