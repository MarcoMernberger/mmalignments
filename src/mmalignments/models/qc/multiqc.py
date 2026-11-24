"""Module contains a MultiQC interface for aggregating QC reports."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import Element, element
from mmalignments.models.externals import (
    External,
    SubroutineIn,
    subroutine,
)
from mmalignments.models.overlay import (
    CfgType,
    OutSpec,
    OutType,
    ParType,
    TagType,
)
from mmalignments.models.parameters import ParamSet
from mmalignments.models.tags import (
    Method,
    Stage,
    State,
)

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
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """Initialize MultiQC wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "multiqc").
        primary_binary : str
            Path to MultiQC executable (default: "multiqc").
        version : str | None
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

    ####################################################################################
    # Elements and Subroutines
    ####################################################################################

    @subroutine
    def run_multiqc(
        self,
        analysis_dir: Path | str,
        output_dir: Path | str,
        *,
        report_name: str = "multiqc_report.html",
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
                Additional MultiQC parameters (e.g., {"--module": "fastqc"}).
        cfg : CfgType | None
                Optional configuration for running the subprocess (e.g. working
                directory, environment variables).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes MultiQC.

        Examples
        --------
        >>> multiqc = MultiQC()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> runner = multiqc.run_multiqc(
        ...     analysis_dir="results/qc",
        ...     output_dir="results/qc/multiqc",
        ...     report_name="project_qc.html",
        ...     cfg=ExternalRunConfig(threads=2),
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

        # Add analysis directory
        arguments.append(str(analysis_dir))
        subcommand = None
        return arguments, subcommand, [analysis_dir], [output_dir], None, None, None

    @element
    def aggregate(
        self,
        analysis_dir: Path | str,
        qc_elements: Mapping[str, Element],
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Aggregate QC reports from multiple Elements using MultiQC (high-level).

        Creates an Element that runs MultiQC to aggregate quality control
        reports from various tools into a single interactive HTML report.

        Parameters
        ----------
        analysis_dir : Path | str
            Directory to search for QC reports. MultiQC will recursively search this
            directory.
        qc_elements : list[Element]
            List of QC Elements whose outputs should be aggregated.
            Their input files are extracted for signature tracking.
        tag : TagType | None
            Tag override for the output Element.
        out : OutType | None
            Output specification override (folder, file name of the HTML
            report). Defaults to ``<analysis_dir>/multiqc_report`` and the
            tag's default name.
        par : ParType | None
            Additional MultiQC parameters.
        cfg : CfgType | None
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
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> qc_elems = [fastp_elem, fastqc_raw_elem, fastqc_cleaned_elem]
        >>> report_elem = multiqc.aggregate(
        ...     analysis_dir="results/qc",
        ...     qc_elements=qc_elems,
        ...     out=OutSpec(folder="results/qc/multiqc", stem="sample_qc", ext="html"),
        ...     cfg=ExternalRunConfig(threads=2),
        ... )
        >>> report_elem.run()
        """
        par = par or ParType(force=True)  # If True, overwrite existing reports
        prior = next(iter(qc_elements.values())).tag
        resolved_tag = prior.bump(
            level=(
                max([el.tag.level for el in qc_elements.values()]) + 1
                if qc_elements
                else 1
            ),
            stage=Stage.QC,
            method=Method.MULTIQC,
            state=State.REPORT,
        ).resolve(tag)

        if not analysis_dir:
            raise ValueError(
                "Analysis directory with qc reports must be provided to search."
            )
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=Path(analysis_dir) / "multiqc_report", ext="html"), out
        )
        outdir = resolved_out.folder
        report_file = resolved_out.file
        # Collect input files from QC elements
        input_files = []
        for elem in qc_elements.values():
            if elem.files:
                input_files.extend(elem.files)
        input_files = tuple(input_files)
        # Build runner
        runner = self.run_multiqc(
            analysis_dir=analysis_dir,
            output_dir=outdir,
            report_name=report_file.name,
            par=par,
            cfg=cfg,
        )

        # MultiQC output files
        artifacts = ArtifactSet(
            report_file,
            primary_name="multiqc_html",
            multiqc_data_dir=report_file.parent / (report_file.stem + "_data"),
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "aggregate",
            n_elements=len(qc_elements),
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=input_files,
            artifacts=artifacts,
            pres=tuple(qc_elements.values()),
        )
