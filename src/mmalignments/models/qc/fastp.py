"""Module contains a fastp interface for FASTQ quality control and filtering."""

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
        primary_binary: str = "./code/fastp",  # replace with "fastp" if fastp in PATH
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
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
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
        params : Params | None
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
        paths = [out_r1, json_out, html_out]
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
            paths.append(out_r2)
        return arguments, paths, None, None, None

    @element
    def qc(
        self,
        sample: NextGenSampleElement,
        *,
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Run fastp quality control and filtering (high-level).

        Creates an Element that runs fastp on input FASTQ files from a Sample
        or Element, performs quality filtering and adapter trimming, and
        generates filtered FASTQ files with QC reports.

        Parameters
        ----------
        sample : Sample | Element
            Sample or Element containing FASTQ file paths and metadata.
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
            Additional fastp parameters.
        cfg : ExternalRunConfig | None
            Configuration for running the fastp command (e.g., working
            directory, threads).

        Returns
        -------
        Element
            An Element that executes fastp when run. Artifacts include:
            - fastq_r1: Cleaned R1 FASTQ file
            - fastq_r2: Cleaned R2 FASTQ file (if paired-end)
            - fastp_json: JSON report
            - fastp_html: HTML report
            - sample_name: Sample name string

        Examples
        --------
        >>> fastp = FastP()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> sample = Sample(name="my_sample", pairing="paired",
        ...                 fastq_r1_path="raw_R1.fastq.gz",
        ...                 fastq_r2_path="raw_R2.fastq.gz")
        >>> fastp_elem = fastp.qc(sample, cfg=ExternalRunConfig(threads=8))
        >>> fastp_elem.run()
        """
        default_tag = ElementTag(
            root=sample.root,
            level=sample.tag.level + 1,
            stage=Stage.PREP,
            method=Method.FASTP,
            state=State.TRIM,
            omics=sample.tag.omics,
            ext=".R1.fastq.gz",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag
        fastq_r1, fastq_r2, _, _ = sample_fastqs(sample)

        outdir = outdir or sample.result_dir / "qc" / "fastp"
        filename = filename or tag.default_output
        params = params or Params(thread=10)
        cfg = cfg or ExternalRunConfig(threads=params.get("thread", 10))
        # Extract FASTQ paths from Sample or Element

        # Define output paths
        out_r1 = outdir / filename
        json_out = out_r1.with_name(out_r1.name.replace(".fastq.gz", ".json"))
        html_out = out_r1.with_name(out_r1.name.replace(".fastq.gz", ".html"))
        # Build artifacts
        out_r2 = (
            out_r1.with_name(out_r1.name.replace(".R1.", ".R2.")) if fastq_r2 else None
        )
        input_files = [fastq_r1]
        artifacts = {
            "fastq_r1": out_r1,
            "json": json_out,
            "html": html_out,
        }
        if fastq_r2:
            artifacts["fastq_r2"] = out_r2
            input_files.append(fastq_r2)

        # Build runner
        runner = self.run_fastp(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            out_r1=out_r1,
            out_r2=out_r2,
            json_out=json_out,
            html_out=html_out,
            params=params,
            cfg=cfg,
        )

        # Build prerequisites list
        pres = [sample] if isinstance(sample, Element) else []
        # determinants
        determinants = self.signature_determinants(params)
        key = f"{tag.default_name}_fastp_qc_{self.version_name}"
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=input_files,
            artifacts=artifacts,
            pres=pres,
        )
