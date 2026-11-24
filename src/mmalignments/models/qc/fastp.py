"""Module contains a fastp interface for FASTQ quality control and filtering."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Mapping, cast

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import (
    Element,
    FastqArtifact,
    NextGenSample,
    element,
)
from mmalignments.models.overlay import (
    CfgType,
    OutSpec,
    OutType,
    ParType,
    TagType,
)
from mmalignments.models.tags import Method, Stage, State

from ..externals import External, ExternalRunConfig, subroutine, SubroutineIn
from ..parameters import ParamSet

logger = logging.getLogger(__name__)


class FastP(External):
    """fastp interface for FASTQ quality control and filtering.

    fastp performs quality control, adapter trimming, quality filtering,
    and per-read quality pruning for FASTQ files. It generates HTML/JSON
    reports and filtered FASTQ files.

    Examples
    --------
    Run fastp on paired-end reads::

        fastp = FastP()
        from mmalignments.models.externals import ExternalRunConfig
        runner = fastp.run_fastp(
                fastq_r1="sample_R1.fastq.gz",
                fastq_r2="sample_R2.fastq.gz",
                out_r1="cleaned_R1.fastq.gz",
                out_r2="cleaned_R2.fastq.gz",
                json_out="fastp.json",
                html_out="fastp.html",
                cfg=ExternalRunConfig(threads=8)
        )
        runner()  # Execute

    Run fastp on single-end reads::

        fastp = FastP()
        from mmalignments.models.externals import ExternalRunConfig
        runner = fastp.run_fastp(
                fastq_r1="sample.fastq.gz",
                fastq_r2=None,
                out_r1="cleaned.fastq.gz",
                out_r2=None,
                json_out="fastp.json",
                html_out="fastp.html",
                cfg=ExternalRunConfig(threads=4)
        )
        runner()

    Using high-level Element interface::

        fastp = FastP()
        from mmalignments.models.data import Sample
        from mmalignments.models.externals import ExternalRunConfig
        sample = Sample(
            name="my_sample",
            pairing="paired",
            fastq_r1_path="raw_R1.fastq.gz",
            fastq_r2_path="raw_R2.fastq.gz"
        )
        fastp_elem = fastp.qc(sample, cfg=ExternalRunConfig(threads=8))
        fastp_elem.run()
    """

    def __init__(
        self,
        name: str = "fastp",
        primary_binary: str = "fastp",
        version: str | None = None,
        source: str = "https://github.com/OpenGene/fastp",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """Initialize fastp wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "fastp").
        primary_binary : str
            Path to fastp executable (default: "fastp").
        version : str | None
            Version string override.
        source : str
            URL/source for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | str | Path | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        parameters_file = Path(__file__).parent / "fastp.json"
        parameters = parameters or parameters_file
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters,
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get fastp version string.

        Parameters
        ----------
        fallback : str | None
                Value to return if version cannot be determined.

        Returns
        -------
        str | None
                Version string (e.g., "0.23.2") or fallback if not found.
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
                check=False,
                timeout=5,
            )
            output = (cp.stdout or cp.stderr or "").strip()
            if output:
                # fastp prints: "fastp 0.23.2"
                for token in output.split():
                    if token[0].isdigit() and "." in token:
                        return token
        except Exception:
            pass

        return fallback

    @subroutine
    def run_fastp(
        self,
        fastq_r1: Path | str,
        fastq_r2: Path | str | None,
        out_r1: Path | str,
        out_r2: Path | str | None,
        json_out: Path | str,
        html_out: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run fastp quality control and filtering (low-level).

        Creates a zero-argument callable that runs fastp on input FASTQ files,
        performs quality filtering, adapter trimming, and generates QC reports.

        Parameters
        ----------
        fastq_r1 : Path | str
            Input R1 FASTQ file (forward reads or single-end).
        fastq_r2 : Path | str | None
            Input R2 FASTQ file (reverse reads) for paired-end data.
        out_r1 : Path | str
            Output R1 FASTQ file (filtered/trimmed).
        out_r2 : Path | str | None
            Output R2 FASTQ file (filtered/trimmed) for paired-end data.
        json_out : Path | str
            Output JSON report file.
        html_out : Path | str
            Output HTML report file.
        par : ParType | None
            Additional fastp parameters (e.g., {"--qualified_quality_phred": 20}).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes fastp.

        Examples
        --------
        >>> fastp = FastP()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> runner = fastp.run_fastp(
        ...     fastq_r1="sample_R1.fastq.gz",
        ...     fastq_r2="sample_R2.fastq.gz",
        ...     out_r1="cleaned_R1.fastq.gz",
        ...     out_r2="cleaned_R2.fastq.gz",
        ...     json_out="fastp.json",
        ...     html_out="fastp.html",
        ...     cfg=ExternalRunConfig(threads=8)
        ... )
        >>> runner()
        """
        fastq_r1 = Path(fastq_r1).absolute()
        out_r1 = Path(out_r1).absolute()
        json_out = Path(json_out).absolute()
        html_out = Path(html_out).absolute()

        # Ensure output directories exist
        out_r1.parent.mkdir(parents=True, exist_ok=True)
        json_out.parent.mkdir(parents=True, exist_ok=True)
        html_out.parent.mkdir(parents=True, exist_ok=True)

        # Build command
        arguments = [
            "-i",
            str(fastq_r1),
            "-o",
            str(out_r1),
            "-j",
            str(json_out),
            "-h",
            str(html_out),
        ]

        # Add paired-end parameters if provided
        in_paths = [fastq_r1]
        out_paths = [out_r1, json_out, html_out]
        if fastq_r2 is not None:
            fastq_r2 = Path(fastq_r2).absolute()
            if out_r2 is None:
                raise ValueError("out_r2 must be provided when fastq_r2 is given")
            out_r2 = Path(out_r2).absolute()
            out_r2.parent.mkdir(parents=True, exist_ok=True)
            arguments.extend(
                [
                    "-I",
                    str(fastq_r2),
                    "-O",
                    str(out_r2),
                ]
            )
            in_paths.append(fastq_r2)
            out_paths.append(out_r2)
        subcommand = None
        return arguments, subcommand, in_paths, out_paths, None, None, None

    @element
    def qc(
        self,
        sample: NextGenSample,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Run fastp quality control and filtering (high-level).

        Creates an Element that runs fastp on input FASTQ files from a Sample
        or Element, performs quality filtering and adapter trimming, and
        generates filtered FASTQ files with QC reports.

        Parameters
        ----------
        sample : Sample | Element
            Sample or Element containing FASTQ file paths and metadata.
        tag : TagType | None
            Tag override for the output Element.
        out : OutType | None
            Output specification override. The primary file is the cleaned R1
            FASTQ (default ext ``R1.fastq.gz``); R2, JSON and HTML names are
            derived from it. Defaults to ``<sample.result_dir>/qc/fastp``.
        par : ParType | None
            Additional fastp parameters.
        cfg : CfgType | None
            Configuration for running the fastp command (e.g., working
            directory, threads).

        Returns
        -------
        Element
            An Element that executes fastp when run. Artifacts include:
            - fastq_r1: Cleaned R1 FASTQ file
            - fastq_r2: Cleaned R2 FASTQ file (if paired-end)
            - json: JSON report
            - html: HTML report

        Examples
        --------
        >>> fastp = FastP()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> fastp_elem = fastp.qc(sample, cfg=ExternalRunConfig(threads=8))
        >>> fastp_elem.run()
        >>> fastp_elem, sample
        """
        resolved_tag = sample.tag.bump(
            stage=Stage.PREP,
            method=Method.FASTP,
            state=State.TRIM,
        ).resolve(tag)

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=sample.result_dir / "qc" / "fastp", ext="R1.fastq.gz"),
            out,
        )
        resolved_par = par or ParType(thread=10)
        resolved_cfg = cfg or ExternalRunConfig(
            threads=cast(int, resolved_par.get("thread", 10))
        )

        # Define output paths
        out_r1 = resolved_out.file
        json_out = out_r1.with_name(out_r1.name.replace(".fastq.gz", ".json"))
        html_out = out_r1.with_name(out_r1.name.replace(".fastq.gz", ".html"))
        # Build artifacts
        out_r2 = (
            out_r1.with_name(out_r1.name.replace(".R1.", ".R2.")) if sample.r2 else None
        )
        input_files = [sample.r1]
        extras = {"json": json_out, "html": html_out, "fastq_r1": out_r1}
        if sample.r2 and out_r2:
            extras["fastq_r2"] = out_r2
            input_files.append(sample.r2)
        trimmed_fastq_source = FastqArtifact(out_r1, out_r2)
        artifacts = ArtifactSet(trimmed_fastq_source, primary_name="source", **extras)

        # Build runner
        runner = self.run_fastp(
            fastq_r1=sample.r1,
            fastq_r2=sample.r2 if sample.r2 else None,
            out_r1=out_r1,
            out_r2=out_r2,
            json_out=json_out,
            html_out=html_out,
            par=par,
            cfg=resolved_cfg,
        )

        # Build prerequisites list
        pres = [sample] if isinstance(sample, Element) else []
        # determinants
        determinants = self.signature_determinants(par)
        key = Element.generate_key(resolved_tag, self.version_name, "qc")
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(input_files),
            artifacts=artifacts,
            pres=tuple(pres),
        )
