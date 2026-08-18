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
from typing import Literal, Mapping

from mmalignments.models.artifacts import ArtifactSet, FastqArtifact, OutputSpec
from mmalignments.models.elements import Element, NextGenSample, element
from mmalignments.models.parameters import (
    ParamRegistry,
    Params,
    ParamSet,
    ParamSpec,
    ToolThreadSpec,
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

from ..externals import (
    External,
    ExternalRunConfig,
    SubroutineIn,
    subroutine,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameter definitions (inlined — no JSON file needed for own Rust tools)
# ---------------------------------------------------------------------------


def _build_param_registry() -> ParamRegistry:
    """Parameter registry for mmfqcount.

    Usage: ./NGmerge {-1 <file> -2 <file> -o <file>}  [optional arguments]
    Required arguments:
    -1  <file>       Input FASTQ file with reads from forward direction
    -2  <file>       Input FASTQ file with reads from reverse direction
    -o  <file>       Output FASTQ file(s):
                    - in 'stitch' mode (def.), the file of merged reads
                    - in 'adapter-removal' mode (-a), the output files
                        will be <file>_1.fastq and <file>_2.fastq
    Alignment parameters:
    -m  <int>        Minimum overlap of the paired-end reads (def. 20)
    -p  <float>      Mismatches to allow in the overlapped region
                        (a fraction of the overlap length; def. 0.10)
    -a               Use 'adapter-removal' mode (also sets -d option)
    -d               Option to check for dovetailing (with 3' overhangs)
    -e  <int>        Minimum overlap of dovetailed alignments (def. 50)
    -s               Option to produce shortest stitched read
    I/O options:
    -l  <file>       Log file for stitching results of each read pair
    -f  <file>       FASTQ files for reads that failed stitching
                        (output as <file>_1.fastq and <file>_2.fastq)
    -c  <file>       Log file for dovetailed reads (adapter sequences)
    -j  <file>       Log file for formatted alignments of merged reads
    -z/-y            Option to gzip (-z) or not (-y) FASTQ output(s)
    -i               Option to produce interleaved FASTQ output(s)
    -w  <file>       Use given error profile for merged qual scores
    -g               Use 'fastq-join' method for merged qual scores
    -q  <int>        FASTQ quality offset (def. 33)
    -u  <int>        Maximum input quality score (0-based; def. 40)
    -t  <char>       Delimiter for headers of paired reads (def. ' ')
    -n  <int>        Number of threads to use (def. 1)
    -v               Option to print status updates/counts to stderr
    """
    _specs: dict[str, ParamSpec] = {
        "_thread_spec": ToolThreadSpec(
            flag="-n", max_threads=None, fraction=1.0, multi=True
        ),
        "r1": ParamSpec(
            "r1",
            "-1",
            str,
            render=render_value,
            description="Input FASTQ file with reads from forward direction.",
        ),
        "r2": ParamSpec(
            "r2",
            "-2",
            str,
            render=render_value,
            description="Input FASTQ file with reads from reverse direction.",
        ),
        "o": ParamSpec(
            "o", "-o", str, render=render_value, description="Output FASTQ file(s)."
        ),
        "m": ParamSpec(
            "m",
            "-m",
            int,
            render=render_value,
            description="Minimum overlap of the paired-end reads.",
        ),
        "p": ParamSpec(
            "p",
            "-p",
            float,
            render=render_value,
            description="Mismatches to allow in the overlapped region.",
        ),
        "a": ParamSpec(
            "a",
            "-a",
            bool,
            render=render_value,
            description="Use 'adapter-removal' mode.",
        ),
        "d": ParamSpec(
            "d",
            "-d",
            bool,
            render=render_value,
            description="Option to check for dovetailing.",
        ),
        "e": ParamSpec(
            "e",
            "-e",
            int,
            render=render_value,
            description="Minimum overlap of dovetailed alignments.",
        ),
        "s": ParamSpec(
            "s",
            "-s",
            bool,
            render=render_value,
            description="Option to produce shortest stitched read.",
        ),
        "l": ParamSpec(
            "l",
            "-l",
            str,
            render=render_value,
            description="Log file for stitching results of each read pair.",
        ),
        "f": ParamSpec(
            "f",
            "-f",
            str,
            render=render_value,
            description="FASTQ files for reads that failed stitching.",
        ),
        "c": ParamSpec(
            "c",
            "-c",
            str,
            render=render_value,
            description="Log file for dovetailed reads (adapter sequences).",
        ),
        "j": ParamSpec(
            "j",
            "-j",
            str,
            render=render_value,
            description="Log file for formatted alignments of merged reads.",
        ),
        "z": ParamSpec(
            "z",
            "-z",
            bool,
            render=render_value,
            description="Option to gzip FASTQ output(s).",
        ),
        "y": ParamSpec(
            "y",
            "-y",
            bool,
            render=render_value,
            description="Option to not gzip FASTQ output(s).",
        ),
        "i": ParamSpec(
            "i",
            "-i",
            bool,
            render=render_value,
            description="Option to produce interleaved FASTQ output(s).",
        ),
        "w": ParamSpec(
            "w",
            "-w",
            str,
            render=render_value,
            description="Use given error profile for merged qual scores.",
        ),
        "g": ParamSpec(
            "g",
            "-g",
            bool,
            render=render_value,
            description="Use 'fastq-join' method for merged qual scores.",
        ),
        "q": ParamSpec(
            "q", "-q", int, render=render_value, description="FASTQ quality offset."
        ),
        "u": ParamSpec(
            "u",
            "-u",
            int,
            render=render_value,
            description="Maximum input quality score.",
        ),
        "t": ParamSpec(
            "t",
            "-t",
            str,
            render=render_value,
            description="Delimiter for headers of paired reads.",
        ),
        "n": ParamSpec(
            "n", "-n", int, render=render_value, description="Number of threads to use."
        ),
    }
    return ParamRegistry(
        default=ParamSet(_specs, "mmfqcount", "default"),
    )


# ---------------------------------------------------------------------------
# MmFqCount class
# ---------------------------------------------------------------------------


class NGmerge(External):
    """Wrapper for the ``ngmerge`` Rust CLI tool.

    ``ngmerge`` merges two FASTQ files of paired-end reads, optionally with adapter
    trimming and quality filtering.
    """

    def __init__(
        self,
        name: str = "ngmerge",
        primary_binary: str = "NGmerge",
        version: str | None = None,
        source: str = "https://github.com/jsh58/NGmerge.git",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialise the NGmerge wrapper.

        Parameters
        ----------
        name : str
            Logical tool name, default ``"ngmerge"``.
        primary_binary : str
            Executable name on ``$PATH``, default ``"ngmerge"``.
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
        """Return the ngmerge version string.

        Runs ``ngmerge --version`` and extracts the semver token.  The
        output format is stable because we own the Rust source
        (clap auto-generates ``ngmerge X.Y.Z``).

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
        return Path("results") / self.version_name / sample_name

    def default_output_spec(
        self, sample_name: str, compression: Literal["Raw", "Gzip"] = "Gzip"
    ) -> OutputSpec:
        """Return the default output directory for a given sample."""
        ext = "fq" if compression == "Raw" else "fq.gz"
        return OutputSpec(
            stem=sample_name, outdir=self.default_output_dir(sample_name), ext=ext
        )

    # -----------------------------------------------------------------------
    # count — high-level @element
    # -----------------------------------------------------------------------

    @element
    def merge(
        self,
        sample: NextGenSample,
        *,
        mode: str = "stitch",
        compression: Literal["Raw", "Gzip"] = "Gzip",
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Merge paired-end reads in a sample.

        Runs ``ngmerge`` on the FASTQ file(s) of *sample* and writes
        a merged FASTQ file.

        From the ngmerge documentation:
        In the default stitch mode, NGmerge combines paired-end reads that overlap into
        a single read that spans the full length of the original DNA fragment. The ends
        of the merged read are defined by the 5' ends of the original reads. Reads that
        fail the stitching process (due to a lack of sufficient overlap, or excessive
        sequencing errors) are placed into secondary output files, if the user requires
        them.

        The alternative adapter-removal mode returns the original reads as pairs,
        removing the 3' overhangs of those reads whose valid stitched alignment has this
        characteristic (Fig. 1B). Reads whose alignments do not have such overhangs
        (or do not align at all) will also be printed to the output files, unmodified.

        Parameters
        ----------
        sample : NextGenSample
            Input sample.  Both single-end and paired-end samples are handled
            automatically.
        mode : str
            Mode for merging reads. Options are "stitch" (default) or "adapter-removal".
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec: OutputSpec | None
            Optional output specification.  If *None*, a default is derived from the
            sample name.
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
        params = Params().update(params)
        if not hasattr(sample, "r1"):
            raise ValueError("Sample element must have 'r1' attribute.")
        if not hasattr(sample, "r2"):
            raise ValueError(
                "Sample element must have 'r2' attribute (can be None)."
            )  # noqa: E501

        fastq_r1 = sample.r1
        fastq_r2 = sample.r2
        tag = from_prior(
            sample.tag,
            tag,
            stage=Stage.PREP,
            method=Method.NGMERGE,
            state=State.MERGED,
        )
        spec = self.default_output_spec(sample.root).merge(outspec)
        output_fastq = spec.path()

        runner = self.run_ngmerge(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            output_fastq=output_fastq,
            mode=mode,
            compression=compression,
            params=params,
            cfg=cfg,
        )

        key, name = self.build_element_name(tag)
        determinants = (mode,) + self.signature_determinants(params)
        inputs = (fastq_r1, fastq_r2) if fastq_r2 is not None else (fastq_r1,)
        artifacts = ArtifactSet(FastqArtifact(output_fastq), primary_name="fastq")
        if "l" in params:
            artifacts = artifacts.with_extra("log", params["l"])
        source = Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=inputs,
            artifacts=artifacts,
            pres=sample.pres,
            name=name,
        )
        return source

    # -----------------------------------------------------------------------
    # count — low-level @subroutine
    # -----------------------------------------------------------------------

    @subroutine
    def run_ngmerge(
        self,
        fastq_r1: Path,
        fastq_r2: Path,
        output_fastq: Path | str,
        mode: str = "stitch",
        compression: Literal["Raw", "Gzip"] = "Gzip",
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``ngmerge``.

        Parameters
        ----------
        fastq_r1 : Path | str
            R1 FASTQ file (plain or gzip).
        fastq_r2 : Path | str
            R2 FASTQ file (plain or gzip), or *None* for single-end.
        output_fastq : Path | str
            Destination path for the merged FASTQ.
        mode : str
            Mode for merging reads. Options are "stitch" (default) or "adapter-removal".
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
        output_fastq = Path(output_fastq).resolve()

        arguments = [
            "-1",
            str(fastq_r1),
            "-2",
            str(fastq_r2),
            "-o",
            str(output_fastq),
        ]
        if compression == "Gzip":
            arguments.append("-z")
        if mode == "adapter-removal":
            arguments += ["-a", "-d"]

        return (
            arguments,
            "",
            (fastq_r1, fastq_r2),
            (output_fastq,),
            None,  # no piped output
            None,  # no pre-hook
            None,  # no post-hook
        )


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _build_param_registry_as_mapping() -> dict:
    """Return a plain dict so External.__init__ can accept it.

    The actual ParamRegistry is set directly in MmFqCount.__init__ after the
    super().__init__() call.
    """
    return {}
