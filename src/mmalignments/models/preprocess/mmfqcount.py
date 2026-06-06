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
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Mapping, Callable

from mmalignments.models.elements import (
    Element,
    NextGenSampleElement,
    TableElement,
    element,
    CallSpec
)
from mmalignments.services.io import parents
from mmalignments.models.parameters import (
    Params,
    ParamSet,
    ParamSpec,
    ParamRegistry,
    render_value,
    render_flag,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    Omics,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.dependencies import function_hash
from ..externals import External, ExternalRunConfig, Runnable, SubroutineIn, subroutine

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
        "trim_start":  ParamSpec("trim_start",  "--trim-start",  str, render=render_value,
                                 description="Trim from first occurrence of k-mer (inclusive)."),
        "trim_stop":   ParamSpec("trim_stop",   "--trim-stop",   str, render=render_value,
                                 description="Trim up to (exclusive) last occurrence of k-mer."),
        "trim_length": ParamSpec("trim_length", "--trim-length", int, render=render_value,
                                 description="Keep at most this many bases after trimming."),
        "split_by":    ParamSpec("split_by",    "--split-by",    str, render=render_value,
                                 description="Split counts by this read-name tag (e.g. 'sgRNAid')."),
    }
    match_specs: dict[str, ParamSpec] = {
        "seq_col": ParamSpec("seq_col", "--seq-col", str, default="Sequence", render=render_value,
                             description="Column holding the R1 sequence."),
        "r2_col":  ParamSpec("r2_col",  "--r2-col",  str, render=render_value,
                             description="Column holding the R2 sequence (paired mode)."),
        "id_col":  ParamSpec("id_col",  "--id-col",  str, default="Name", render=render_value,
                             description="Column holding the sequence identifier."),
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

    # -----------------------------------------------------------------------
    # count — high-level @element
    # -----------------------------------------------------------------------

    @element
    def count(
        self,
        sample: NextGenSampleElement,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
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
        sample : NextGenSampleElement
            Input sample.  Both single-end and paired-end samples are handled
            automatically.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None
            Output directory; defaults to
            ``results/counts/<version>/<sample.root>``.
        filename : Path | str | None
            Output file name; defaults to the tag's ``default_output``.
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
        fastq_r1 = sample.fastq_r1
        fastq_r2 = sample.fastq_r2 if sample.fastq_r2 else None

        tag = from_prior(
            sample.tag,
            tag,
            stage=Stage.QUANT,
            method=Method.MMFQCOUNT,
            state=State.COUNT,
            ext="tsv",
        )

        output_dir = Path(outdir or self.default_output_dir(sample.root))
        out_filename = filename or tag.default_output
        output_tsv = output_dir / out_filename

        runner = self.run_count(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            output_tsv=output_tsv,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(tag, "count")
        determinants = self.signature_determinants(params, subroutine="count")
        inputs = (fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,)

        return TableElement(
            key,
            runner,
            tag=tag,
            tsv=output_tsv,
            artifacts={"tsv": output_tsv},
            determinants=determinants,
            inputs=inputs,
            pres=(sample,),
            name=name,
        )

    # -----------------------------------------------------------------------
    # count — low-level @subroutine
    # -----------------------------------------------------------------------

    @subroutine
    def run_count(
        self,
        fastq_r1: Path | str,
        fastq_r2: Path | str | None,
        output_tsv: Path | str,
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
        output_tsv : Path | str
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
        fastq_r1 = Path(fastq_r1).absolute()
        output_tsv = Path(output_tsv).absolute()

        arguments = [
            "count",
            "--r1", str(fastq_r1),
            "--output", str(output_tsv),
        ]
        if fastq_r2 is not None:
            arguments += ["--r2", str(Path(fastq_r2).absolute())]
        in_paths = [fastq_r1] + ([Path(fastq_r2).absolute()] if fastq_r2 else [])
        out_paths = [output_tsv]

        return (
            arguments,
            "count",
            in_paths,
            out_paths,
            None,   # no piped output
            None,   # no pre-hook
            None,   # no post-hook
        )

    # -----------------------------------------------------------------------
    # match — high-level @element
    # -----------------------------------------------------------------------

    @element
    def match(
        self,
        counts: Element,
        predefined: "Element | Path | str",
        *,
        seq_col: str = "Sequence",
        r2_col: str | None = None,
        id_col: str = "Name",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
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
        predefined : Element | Path | str
            Either a pipeline Element whose ``"tsv"`` (or ``"path"``) artifact
            points to the predefined-sequences TSV, or a plain path.
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
        outdir : Path | str | None
            Output directory; defaults to the same directory as the counts TSV.
        filename : Path | str | None
            Base name for the matched output file.
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

        counts_tsv = Path(counts.tsv).absolute()

        # Resolve predefined path
        if isinstance(predefined, Element):
            pred_path = Path(
                predefined.artifacts.get("tsv")
                or predefined.artifacts.get("path")
                or next(iter(predefined.artifacts.values()))
            ).absolute()
            pred_pres: tuple = (predefined,)
        else:
            pred_path = Path(predefined).absolute()
            pred_pres = ()

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
            state=State.ANNOTATE,
            ext="tsv",
        )

        output_dir = Path(outdir or counts_tsv.parent)
        out_filename = filename or tag.default_output
        matched_tsv = output_dir / out_filename
        unmatched_tsv = output_dir / (
            out_filename.replace(".tsv", "_unmatched.tsv")
            if ".tsv" in str(out_filename)
            else f"{out_filename}_unmatched.tsv"
        )

        runner = self.run_match(
            counts_tsv=counts_tsv,
            predefined_tsv=pred_path,
            matched_tsv=matched_tsv,
            unmatched_tsv=unmatched_tsv,
            params=merged_params,
            cfg=cfg,
        )

        key, name = self.build_element_name(
            tag, "match", seq_col=seq_col, id_col=id_col
        )
        determinants = self.signature_determinants(merged_params, subroutine="match")

        return TableElement(
            key,
            runner,
            tag=tag,
            artifacts={"matched": matched_tsv, "unmatched": unmatched_tsv},
            determinants=determinants,
            inputs=(counts_tsv, pred_path),
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
            "--counts",     str(counts_tsv),
            "--predefined", str(predefined_tsv),
            "--output",     str(matched_tsv),
            "--unmatched",  str(unmatched_tsv),
        ]

        # Append column flags from params
        params = params or Params()
        cli_extras = self.to_cli(params, subroutine="match")
        arguments += cli_extras

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
        score: Callable[[int, int], float] | None = None,
        seq_col: str = "R2",
        column_prefix: str = "Count",
        exclude_score: float | None | Callable = 0.0,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> TableElement:
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
        seq_col : str, optional
            Column name in the input TSVs that contains the sequence 
            identifiers. Default is "R2".
        column_prefix : str, optional
            Prefix for the count columns in the input TSV (e.g. "Count (sample)")
        exclude_score : float | None | Callable, optional
            If set, sequences with a score equal to this will be excluded from 
            the filtered score TSV. If a callable is provided, it will be called 
            with the score and should return True if the sequence should be 
            excluded.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None, optional
            Optional output directory override.
        filename : Path | str | None, optional
            Optional filename override. 
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
        def default_score(count1: int, count2: int) -> float:
            if count2 == 0:
                return float('inf') if count1 > 0 else 0.0
            return np.log2(count1 / count2)
        if score is None:
            score = default_score

        tag = from_prior(
            compare.tag,
            tag,
            stage=Stage.QUANT,
            method=Method.MMFQCOUNT,
            state=State.SCORE,
            ext="tsv",
        )
        output_dir = Path(outdir or compare.file.parent).absolute()
        score_filename = filename or tag.default_output
        score_tsv = output_dir / score_filename
        filtered_score_tsv = score_tsv.with_suffix(".filtered.tsv")
        pres = (compare, against)
        inputs = (compare.file, against.file)
        determinants = [function_hash(score), seq_col, column_prefix, exclude_score]
        # Resolve predefined path

        runner = self.compare_counts(
            compare=compare.tag.root,
            against=against.tag.root,
            compare_tsv=compare.file,
            against_tsv=against.file,
            score_tsv=score_tsv,
            filtered_score_tsv=filtered_score_tsv,
            seq_col=seq_col,
            score_func=score,
            column_prefix=column_prefix,
            exclude_score=exclude_score,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(
            tag, "compare", param_str="_against=" + against.tag.root,
        )

        return TableElement(
            key,
            runner,
            tag=tag,
            artifacts={"tsv": score_tsv, "filtered": filtered_score_tsv},
            tsv=score_tsv,
            determinants=determinants,
            inputs=inputs,
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
        compare_tsv: Path | str,
        against_tsv: Path | str,
        score_tsv: Path | str,
        filtered_score_tsv: Path | str,
        *,
        seq_col: str = "R2",
        score_func: Callable[[int, int], float],
        column_prefix: str = "Count",
        exclude_score: float | None | Callable = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:
        """Compare two count TSVs and assign a score to each sequence."""

        compare_path = Path(compare_tsv)
        against_path = Path(against_tsv)
        score_path = Path(score_tsv)
        filtered_path = Path(filtered_score_tsv)

        def _find_count_column(df):
            if column_prefix in df.columns:
                return column_prefix
            matches = [c for c in df.columns if c.lower().startswith(column_prefix.lower())]
            if len(matches) > 0:
                return matches
            else:
                raise ValueError(
                    f"Prefix not in count columns for prefix {column_prefix!r}: {matches}"
                )
            
        def rename_column(name: str) -> Callable[[str], str]:
            def _rename(column_name: str) -> str:
                return f"{column_name} ({name})"
            return _rename

        def _runner() -> None:

            parents(score_path, filtered_path)
            compare_df = pd.read_csv(compare_path, sep="\t")
            against_df = pd.read_csv(against_path, sep="\t")
            read_cols_main = [c for c in compare_df.columns if c.startswith(column_prefix)]
            read_cols_other = [c for c in against_df.columns if c.startswith(column_prefix)]
            read_cols_common = set(read_cols_main) & set(read_cols_other)
            # rename Read Count columns
            rename_compare = rename_column(compare)
            rename_against = rename_column(against)
            compare_df = compare_df.rename(columns={c: rename_compare(c) for c in read_cols_main})
            against_df = against_df.rename(columns={c: rename_against(c) for c in read_cols_other})
            additional_columns = [rename_against(c) for c in read_cols_common]
            against_df = against_df[[seq_col] + additional_columns]

            # ḿerge on sequence
            merged = compare_df.merge(
                against_df,
                on=f"{seq_col}",
                how="left",
            )
            # sequences missing in against_df get NaN → treat as 0
            merged[additional_columns] = merged[additional_columns].fillna(0)
            # intersection of comparable count columns
            score_columns = []
            for col in read_cols_common:
                col_compare = rename_compare(col)
                col_against = rename_against(col)

                score_col = f"Score {col.replace(column_prefix, '').strip()} ({compare} vs {against})"
                score_columns.append(score_col)
                merged[score_col] = merged.apply(
                    lambda row: (
                        score_func(row[col_compare], row[col_against])
                        if pd.notna(row[col_compare]) and pd.notna(row[col_against])
                        else pd.NA
                    ),
                    axis=1
                )
            merged.to_csv(score_path, sep="\t", index=False)

            # Keep only rows where at least one count column (compare or against) is > 0.
            # Rows where every count is 0 (sequence invisible in both conditions) are discarded.
            all_count_cols = (
                [rename_compare(c) for c in read_cols_common]
                + [rename_against(c) for c in read_cols_common]
            )
            if all_count_cols:
                filtered = merged.loc[(merged[all_count_cols] > 0).any(axis=1)]
            else:
                filtered = merged

            filtered.to_csv(filtered_path, sep="\t", index=False)

        callspec = CallSpec(
            "compare_counts",
            kwargs={
                "compare": compare_tsv, "against": against_tsv, 
                "compare_tsv": compare_tsv, "against_tsv": against_tsv, 
                "score_tsv": score_tsv, "filtered_score_tsv": filtered_score_tsv, 
                "seq_col": seq_col, 
                "column_prefix": column_prefix},
        ).render()
        return Runnable(
            _runner,
            display=callspec,
        )

    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------

    def count_and_match(
        self,
        sample: NextGenSampleElement,
        predefined: "Element | Path | str",
        *,
        seq_col: str = "Sequence",
        r2_col: str | None = None,
        id_col: str = "Name",
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        count_params: Params | None = None,
        match_params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[Element, Element]:
        """Count reads and immediately match against predefined sequences.

        Convenience wrapper that chains :meth:`count` → :meth:`match`.

        Parameters
        ----------
        sample : NextGenSampleElement
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
        outdir : Path | str | None
            Output directory.
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
        count_el = self.count(
            sample,
            tag=tag,
            outdir=outdir,
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
            outdir=outdir,
            params=match_params,
            cfg=cfg,
        )
        return match_el, count_el

    def samplecount(
            self, 
            samples: Mapping[str, NextGenSampleElement], 
            *,
            tag: PartialElementTag | ElementTag | None = None,
            outdir: Path | str | None = None,
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
        samples : Mapping[str, NextGenSampleElement]
            A mapping of sample names to NextGenSampleElements to be counted.
        tag : PartialElementTag | ElementTag | None, optional
            Optional tag override applied to all count elements, by default None.
        outdir : Path | str | None, optional
            Output directory for all count elements, by default None. If None, 
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
            output_dir = Path(outdir) / sample_name if outdir else self.default_output_dir(sample)
            counted = self.count(
                sample, 
                tag=tag,
                outdir=output_dir,
                params=params
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
