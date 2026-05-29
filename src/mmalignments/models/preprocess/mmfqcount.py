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
from typing import Mapping

from mmalignments.models.elements import (
    Element,
    NextGenSampleElement,
    element,
)
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
from ..externals import External, ExternalRunConfig, Runnable, SubroutineIn, subroutine

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameter definitions (inlined — no JSON file needed for own Rust tools)
# ---------------------------------------------------------------------------

def _build_param_registry() -> ParamRegistry:
    """Construct the ParamRegistry for mmfqcount from the known stable CLI.

    Parameters are read directly from the tool's --help output structure, which
    is stable because we maintain the Rust source.  This avoids the JSON
    indirection used for third-party tools like GATK.
    """

    def _spec(
        name: str,
        flag: str,
        dtype: type,
        *,
        required: bool = False,
        default=None,
        affects_output: bool = True,
        description: str = "",
    ) -> ParamSpec:
        return ParamSpec(
            name=name,
            flag=flag,
            dtype=dtype,
            required=required,
            default=default,
            affects_output=affects_output,
            description=description,
            render=render_value,
        )

    # ── count subcommand ──────────────────────────────────────────────────
    count_specs = {
        "trim_start": _spec(
            "trim_start",
            "--trim-start",
            str,
            description=(
                "Trim read from the first occurrence of this k-mer (inclusive). "
                "Reads without the k-mer are discarded."
            ),
        ),
        "trim_stop": _spec(
            "trim_stop",
            "--trim-stop",
            str,
            description=(
                "Trim read up to (exclusive) the last occurrence of this k-mer."
            ),
        ),
        "trim_length": _spec(
            "trim_length",
            "--trim-length",
            int,
            description="Keep at most this many bases after adapter trimming.",
        ),
    }

    # ── match subcommand ──────────────────────────────────────────────────
    match_specs = {
        "seq_col": _spec(
            "seq_col",
            "--seq-col",
            str,
            default="Sequence",
            affects_output=True,
            description=(
                "Column in the predefined TSV holding the R1 sequence to match."
            ),
        ),
        "r2_col": _spec(
            "r2_col",
            "--r2-col",
            str,
            description=(
                "Column in the predefined TSV holding the R2 sequence "
                "(paired mode). Omit for single-end matching."
            ),
        ),
        "id_col": _spec(
            "id_col",
            "--id-col",
            str,
            default="Name",
            affects_output=True,
            description=(
                "Column in the predefined TSV holding the sequence identifier."
            ),
        ),
    }

    count_set = ParamSet(count_specs, "mmfqcount", "count")
    match_set = ParamSet(match_specs, "mmfqcount", "match")
    default_set = ParamSet({}, "mmfqcount", "default")

    return ParamRegistry(
        default=default_set,
        by_subcommand={"count": count_set, "match": match_set},
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
        print(sample.name, "artifacts:", sample.artifacts)
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

        return Element(
            key,
            runner,
            tag=tag,
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

        # Append trimming flags from params
        params = params or Params()
        cli_extras = self.to_cli(params, subroutine="count")
        arguments += cli_extras

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

        return Element(
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


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _build_param_registry_as_mapping() -> dict:
    """Return a plain dict so External.__init__ can accept it.

    The actual ParamRegistry is set directly in MmFqCount.__init__ after the
    super().__init__() call.
    """
    return {}
