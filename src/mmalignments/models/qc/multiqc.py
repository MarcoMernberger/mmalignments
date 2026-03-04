"""Module contains a MultiQC interface for aggregating QC reports."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Callable, Mapping

from mmalignments.models.externals import External, ExternalRunConfig, subroutine
from mmalignments.models.parameters import Params, ParamSet
from mmalignments.models.tasks import Element, element

logger = logging.getLogger(__name__)


class MultiQC(External):
    """MultiQC interface for aggregating quality control reports.

    Provides methods for running MultiQC to aggregate and summarize
    quality control reports from various tools (FastQC, GATK, samtools, etc.)
    into a single interactive HTML report.

    MultiQC automatically searches for log files and QC reports in the
    specified directory and generates a comprehensive report.

    Examples
    --------
    Run MultiQC on a directory containing QC reports::

    report = multiqc.run(
        analysis_dir="qc_results",
        output_dir="multiqc_report"
    )
    report()


    Run MultiQC with custom report name::

    multiqc = MultiQC()
    report = multiqc.run(
        analysis_dir=["qc_results", "other_qc"],
        output_dir="reports",
        report_name="project_qc",
        force=True
    )
    report()  # execute

    Run MultiQC with specific modules::

    multiqc = MultiQC(parameters={"--module": "fastqc"})
    report_runner = multiqc.run(
            analysis_dir="qc_results",
            output_dir="reports"
    )
    report_runner()


    """

    def __init__(
        self,
        name: str = "multiqc",
        primary_binary: str = "multiqc",
        version: str | None = None,
        source: str = "https://multiqc.info/",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialize MultiQC wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "multiqc").
        primary_binary : str
            Path to MultiQC executable (default: "multiqc").
        version : Optional[str]
            Version string override.
        source : str
            URL/source for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.

        Returns
        -------
        None
        """
        parameters_file = Path(__file__).parent / "multiqc.json"
        parameters = parameters or parameters_file
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get MultiQC version string.

        Parameters
        ----------
        fallback : str | None
                Value to return if version cannot be determined.

        Returns
        -------
        str | None
                Version string (e.g., "1.21") or fallback if not found.
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
                # MultiQC output: "multiqc, version 1.21"
                parts = output.split()
                for i, part in enumerate(parts):
                    if part.lower() == "version" and i + 1 < len(parts):
                        return parts[i + 1]
                    elif part[0].isdigit() and "." in part:
                        return part
        except Exception:
            pass

        return fallback

    @subroutine
    def run_multiqc(
        self,
        analysis_dir: Path | str,
        output_dir: Path | str,
        *,
        report_name: str = "multiqc_report.html",
        force: bool = True,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], None]:
        """Run MultiQC to aggregate QC reports (low-level).

        Creates a zero-argument callable that runs MultiQC on a directory
        containing QC reports and generates a summary HTML report.

        Parameters
        ----------
        analysis_dir : Path | str
                Directory to search for QC reports.
                MultiQC will recursively search this directory.
        output_dir : Path | str
                Output directory for MultiQC report.
        report_name : str
                Name for the report HTML file (default: "multiqc_report.html").
        force : bool
                If True, overwrite existing reports (default: True).
        params : Params | None
                Additional MultiQC parameters (e.g., {"--module": "fastqc"}).
        cfg : ExternalRunConfig | None
                Optional configuration for running the subprocess (e.g. working
                directory, environment variables).

        Returns
        -------
        Callable[[], None]
                Zero-argument callable that executes MultiQC.

        Examples
        --------
        >>> multiqc = MultiQC()
        >>> runner = multiqc.run_multiqc(
        ...     analysis_dir="results/qc",
        ...     output_dir="results/qc/multiqc",
        ...     report_name="project_qc.html"
        ... )
        >>> runner()
        """
        analysis_dir = Path(analysis_dir).absolute()
        # Build command: multiqc [options] -o output_dir -n report_name analysis_dir
        arguments = [
            "-o",
            str(output_dir),
            "-n",
            report_name,
            "-i",
            report_name,  # set report title to match report name
        ]

        # Add force flag if requested
        # if force:
        #     arguments.append("--force")

        # Add analysis directory
        arguments.append(str(analysis_dir))
        return arguments, [output_dir], None, None, None

    @element
    def aggregate(
        self,
        mapped_element_or_analysis_dir: Element | Path | str,
        qc_elements: list[Element],
        *,
        output_dir: Path | str | None = None,
        report_name: str = "pre_alignment_multiqc.html",
        force: bool = True,
        pres: list[Element] | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Aggregate QC reports from multiple Elements using MultiQC (high-level).

        Creates an Element that runs MultiQC to aggregate quality control
        reports from various tools into a single interactive HTML report.

        Parameters
        ----------
        mapped_element_or_analysis_dir : Element | Path | str
            An Element that produces the mapped BAM file. This is used to infer the
            default QC root directory if not provided in cfg or as an argument.
            Alternatively, a Path or string can be provided directly as analysis_dir,
            directory to search for QC reports. MultiQC will recursively search this
            directory.
        qc_elements : list[Element]
            List of QC Elements whose outputs should be aggregated.
            Their input files are extracted for signature tracking.
        output_dir : Path | str | None
            Output directory for MultiQC report. If not provided, defaults
            to "{qc_root}/multiqc/".
        report_name : str
            Name for the report HTML file (default: "pre_alignment_multiqc.html").
        force : bool
            If True, overwrite existing reports (default: True).
        pres : list[Element] | None
            Additional prerequisite Elements (default: None, uses qc_elements).
        params : Params | None
            Additional MultiQC parameters.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
                An Element that executes MultiQC when run. Artifacts include:
                - multiqc_html: Path to HTML report
                - multiqc_data_dir: Path to multiqc_data directory

        Examples
        --------
        >>> multiqc = MultiQC()
        >>> qc_elems = [fastp_elem, fastqc_raw_elem, fastqc_cleaned_elem]
        >>> report_elem = multiqc.aggregate(
        ...     qc_root="results/qc",
        ...     qc_elements=qc_elems,
        ...     output_dir="results/qc/multiqc",
        ...     report_name="sample_qc.html"
        ... )
        >>> report_elem.run()
        """
        params = params or Params(force=force)
        qc_root = cfg.qc_root if cfg and hasattr(cfg, "qc_root") else None
        if isinstance(mapped_element_or_analysis_dir, str | Path):
            analysis_dir = Path(mapped_element_or_analysis_dir).absolute()
        else:
            analysis_dir = analysis_dir or qc_root
        if analysis_dir is None:
            raise ValueError(
                "Analysis directory must be provided either through "
                "mapped_element_or_analysis_dir or cfg.qc_root"
            )
        # Default output directory
        if output_dir is None:
            output_dir = cfg.qc_root / "multiqc"
        output_dir = Path(output_dir).absolute()
        # Collect input files from QC elements
        input_files = []
        for elem in qc_elements:
            input_files.extend(elem.output_files)

        # Build runner
        runner = self.run_multiqc(
            analysis_dir=analysis_dir,
            output_dir=output_dir,
            report_name=report_name,
            force=force,
            params=params,
            cfg=cfg,
        )

        # MultiQC output files
        artifacts = {
            "multiqc_html": output_dir / report_name,
            "multiqc_data_dir": output_dir / "multiqc_data",
        }

        # Use provided pres or default to qc_elements
        if pres is None:
            pres = qc_elements

        determinants = self.signature_determinants(params)
        return Element(
            name=f"{self.name}_aggregate",
            key=f"{self.name}_aggregate_{len(qc_elements)}_elements",
            run=runner,
            determinants=determinants,
            inputs=input_files,
            artifacts=artifacts,
            pres=pres,
        )
