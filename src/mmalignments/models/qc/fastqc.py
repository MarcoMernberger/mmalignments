"""Module contains a FastQC interface for quality control of sequencing data."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.elements import (
    Element,
    NextGenSampleElement,
    element,
    sample_fastqs,
)
from mmalignments.models.tags import ElementTag, Method, Stage, State, merge_tag
from mmalignments.models.tags import PartialElementTag as Tag

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet

logger = logging.getLogger(__name__)


class FastQC(External):
    """FastQC interface for quality control of sequencing reads.

    Provides methods for running FastQC on FASTQ or BAM files to generate
    quality control reports. Returns zero-argument callables that can be
    executed independently or chained.

    FastQC generates HTML reports and zip archives containing quality metrics
    for sequencing data.

    Examples
    --------
    Run FastQC on a single BAM file::

        fastqc = FastQC()
        qc_runner = fastqc.run_qc(
            input_file="aligned.bam",
            output_dir="qc_results"
        )
        qc_runner()  # Execute

    Run FastQC on multiple files::

        fastqc = FastQC()
        qc_runner = fastqc.run_qc(
            input_file=["sample1.bam", "sample2.bam"],
            output_dir="qc_results"
        )
        qc_runner()

    With custom parameters::

        fastqc = FastQC(parameters={"--threads": 4})
        qc_runner = fastqc.run_qc(
            input_file="data.fastq.gz",
            output_dir="qc",
            extract=True  # Extract zip files
        )
        qc_runner()
    """

    def __init__(
        self,
        name: str = "fastqc",
        primary_binary: str = "fastqc",
        version: str | None = None,
        source: str = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """Initialize FastQC wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "fastqc").
        primary_binary : str
            Path to FastQC executable (default: "fastqc").
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
        parameters_file = Path(__file__).parent / "fastqc.json"
        parameters = parameters or parameters_file
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters,
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get FastQC version string.

        Parameters
        ----------
        fallback : str | None
            Value to return if version cannot be determined.

        Returns
        -------
        str | None
            Version string (e.g., "0.12.1") or fallback if not found.
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
                check=True,
                timeout=10,
            )
            output = cp.stdout.strip()
            if output:
                # FastQC output: "FastQC v0.12.1"
                parts = output.split()
                for part in parts:
                    if part.startswith("v"):
                        return part[1:]  # Remove 'v' prefix
                    elif part[0].isdigit():
                        return part
        except Exception:
            pass

        return fallback

    ####################################################################################
    # Helper
    ####################################################################################

    def get_base(self, path: Path) -> str:
        base = None
        if path.suffix == ".gz":
            base = path.with_suffix("").stem
        elif path.suffix == ".fastq" or path.suffix == ".fq":
            base = path.stem
        else:
            raise NotImplementedError(f"Unexpected FASTQ filename: {path}")
        return base

    ####################################################################################
    # Elements
    ####################################################################################

    @subroutine
    def run_fastqc(
        self,
        input_files: list[Path | str],
        output_dir: Path | str,
        params: Params | None,
        cfg: ExternalRunConfig | None = ExternalRunConfig(threads=10),
    ) -> tuple[
        list[str],
        list[str],
        str | None,
        Callable[[], CompletedProcess] | None,
        Callable[[], None] | None,
    ]:
        """Run FastQC on FASTQ files (low-level).

        Creates a zero-argument callable that runs FastQC quality control
        on FASTQ files and generates HTML/ZIP reports.

        Parameters
        ----------
        input_files : list[Path | str]
            List of input FASTQ files to analyze.
        output_dir : Path | str
            Output directory for FastQC reports.
        params : Params | None
            Additional FastQC parameters (e.g., Params(threads=4, extract=True)).
        cfg : ExternalRunConfig | None
            External run configuration (e.g., threads, working directory).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes FastQC.

        Examples
        --------
        >>> fastqc = FastQC()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> runner = fastqc.run_fastqc(
        ...     input_files=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        ...     output_dir="qc_results/fastqc",
        ...     params=Params(threads=4, extract=True),
        ...     cfg=ExternalRunConfig(threads=4),
        ... )
        >>> runner()
        """
        input_files = [Path(f).absolute() for f in input_files]
        # Build command: fastqc [options] -o output_dir input_file(s)
        arguments = ["-o", str(output_dir)]
        # Add input files
        for f in input_files:
            arguments.append(str(f))
        return arguments, [output_dir], None, None, None

    @element
    def qc(
        self,
        sample: NextGenSampleElement,
        label: str,
        *,
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Run FastQC quality control (high-level).

        Creates an Element that runs FastQC on input FASTQ files from a Sample
        or Element and generates HTML/ZIP QC reports.

        Parameters
        ----------
        sample : Sample | Element
            Sample or Element containing FASTQ file paths and metadata.
        label : str
            Label for this QC run (e.g., "raw" or "cleaned").
            Used in Element naming and key generation.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        outdir : Path | str | None
            Directory for the sorted BAM file. If not provided, defaults to
            the same directory as the input BAM file.
        filename : Path | str | None
            Filename override. If not provided, defaults to ``tag.default_output``.
        params : Params | None
            Additional FastQC parameters.
        cfg : ExternalRunConfig | None
            Configuration for running the FastQC command (e.g., working
            directory, threads).

        Returns
        -------
        Element
            An Element that executes FastQC when run. Artifacts include:
            - fastqc_dir: Output directory
            - fastqc_html_r1: HTML report for R1
            - fastqc_zip_r1: ZIP archive for R1
            - fastqc_html_r2: HTML report for R2 (if paired-end)
            - fastqc_zip_r2: ZIP archive for R2 (if paired-end)
            - fastq_r1: Input R1 FASTQ (for traceability)
            - fastq_r2: Input R2 FASTQ (for traceability, if paired-end)

        Examples
        --------
        >>> fastqc = FastQC()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> sample = Sample(name="my_sample", pairing="paired",
        ...                 fastq_r1_path="sample_R1.fastq.gz",
        ...                 fastq_r2_path="sample_R2.fastq.gz")
        >>> fastqc_elem = fastqc.qc(
        ...     sample=sample,
        ...     outdir="results/qc/my_sample/fastqc_raw",
        ...     label="raw",
        ...     cfg=ExternalRunConfig(threads=4),
        ... )
        >>> fastqc_elem.run()
        """
        params = params or Params(extract=False)
        cfg = cfg or ExternalRunConfig()
        default_tag = ElementTag(
            root=sample.root,
            level=sample.tag.level + 1,
            stage=Stage.QC,
            method=Method.FASTQC,
            state=State.REPORT,
            omics=sample.tag.omics,
            param=label,
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag

        # Prepare output folder
        outdir = outdir or sample.fastq_r1.parent / "qc" / "fastqc"
        # Extract FASTQ paths from Sample or Element
        fastq_r1, fastq_r2, _, _ = sample_fastqs(sample)
        # Collect input files
        input_files = [fastq_r1, fastq_r2] if fastq_r2 else [fastq_r1]

        # Build runner
        runner = self.run_fastqc(
            input_files=input_files,
            output_dir=outdir,
            params=params,
            cfg=cfg,
        )
        # Build artifacts
        # FastQC names output files based on input filename
        r1_base = self.get_base(fastq_r1)
        artifacts = {
            "html_r1": outdir / f"{r1_base}_fastqc.html",
            "zip_r1": outdir / f"{r1_base}_fastqc.zip",
        }
        if sample.fastq_r2:
            r2_base = self.get_base(fastq_r2)
            artifacts["html_r2"] = outdir / f"{r2_base}_fastqc.html"
            artifacts["zip_r2"] = outdir / f"{r2_base}_fastqc.zip"

        # Build prerequisites list
        pres = [sample] if isinstance(sample, Element) else []
        determinants = self.signature_determinants(params)
        key = f"{tag.default_name}_{label}_{self.version_name}"
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=input_files,
            artifacts=artifacts,
            pres=pres,
        )

    @subroutine
    def run_qc(
        self,
        input_file: Path | str | list[Path | str],
        output_dir: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], None]:
        """Run FastQC on one or more sequencing files.

        Creates a zero-argument callable that runs FastQC quality control
        on FASTQ or BAM files and generates HTML reports.

        Parameters
        ----------
        input_file : Path | str | list[Path | str]
            Input file(s) to analyze (FASTQ, BAM, or SAM). Can be a single
            file or list of files.
        output_dir : Path | str
            Output directory for FastQC reports.
        params : Params | None
            Additional FastQC parameters (e.g., Params(extract=True)).
        cfg : ExternalRunConfig | None
            External run configuration (e.g., threads, working directory).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes FastQC.

        Examples
        --------
        >>> fastqc = FastQC()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> qc = fastqc.run_qc(
        ...     "aligned.bam", "qc_results", cfg=ExternalRunConfig(threads=4)
        ... )
        >>> qc()
        """
        params = params or Params(threads=10, extract=False)
        cfg = cfg or ExternalRunConfig(threads=Params.get("--threads", 1))
        # Handle single file or list of files
        if isinstance(input_file, (list, tuple)):
            input_files = [Path(f).absolute() for f in input_file]
        else:
            input_files = [Path(input_file).absolute()]

        # Build command: fastqc [options] -o output_dir input_file(s)
        arguments = []
        # Add output directory
        arguments.extend(["-o", str(output_dir)])
        # Add input files
        for f in input_files:
            arguments.append(str(f))
        return arguments, [output_dir], None, None, None
