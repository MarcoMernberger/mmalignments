"""Module contains a GATK interface for variant calling and preprocessing."""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.data import Genome
from mmalignments.models.parameters import Params, ParamSet
from mmalignments.models.tasks import Element, MappedElement, element

from ..externals import External, ExternalRunConfig, subroutine

logger = logging.getLogger(__name__)


class GATK(External):
    """GATK (Genome Analysis Toolkit) interface.

    Provides wrapper methods for GATK tools used in variant calling pipelines.
    Each method returns a zero-argument callable that can be chained as
    post-processing steps or executed independently.

    GATK is invoked as: gatk <ToolName> <toolArgs>

    Examples
    --------
    Mark duplicates in a BAM file::

        gatk = GATK(primary_binary="code/gatk-4.6.2.0/gatk")
        mark_dups_runner = gatk.mark_duplicates(
                input_bam="sorted.bam",
                output_bam="dedup.bam",
                metrics_file="metrics.txt"
        )
        mark_dups_runner()  # Execute

    Run Mutect2 for somatic variant calling::

        mutect_runner = gatk.mutect2(
                reference="genome.fa",
                input_bams=["tumor.bam", "normal.bam"],
                tumor_sample="TUMOR",
                normal_sample="NORMAL",
                output_vcf="tumor.vcf.gz",
                intervals="targets.bed"
        )
        mutect_runner()

    Chain multiple GATK steps::

        mark_dups = gatk.mark_duplicates("in.bam", "dedup.bam", "metrics.txt")
        learn_orient = gatk.learn_read_orientation("f1r2.tar.gz", "model.tar.gz")
        # Execute in sequence
        mark_dups()
        learn_orient()
    """

    def __init__(
        self,
        name: str = "gatk",
        primary_binary: str = "/project/code/gatk-4.6.2.0/gatk",
        version: str | None = None,
        source: str = "https://github.com/broadinstitute/gatk",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """
        Initialize GATK wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "gatk").
        primary_binary : str
            Path to GATK executable (default: "gatk").
            Can be absolute path like "code/gatk-4.6.2.0/gatk".
        version : str | None
            Version string override (e.g., "4.6.2.0").
        source : str
            URL/source for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | str | Path | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
            If a file path or string is provided, it will be loaded from JSON
            and converted to ParamSet.
        """
        parameters_file = Path(__file__).parent / "gatk.json"
        parameters = parameters or parameters_file

        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str = None) -> str | None:
        """Get GATK version string.

        Parameters
        ----------
        fallback : Optional[str]
                Value to return if version cannot be determined.

        Returns
        -------
        Optional[str]
                Version string (e.g., "4.6.2.0") or fallback if not found.
        """
        if self._version:
            return self._version

        if not self.primary_binary or not self.ensure_binary():
            return fallback

        try:
            # GATK version is typically extracted from the jar path or --version
            # Many GATK installations fail on --version without java, so we
            # try to extract from the binary path or jar name
            binary_path = Path(self.primary_binary)
            if "gatk-4." in str(binary_path.parent):
                # Extract from path like "gatk-4.6.2.0/"
                version_str = str(binary_path.parent.name).replace("gatk-", "")
                return version_str

            # Try running --version (may fail without proper java setup)
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            if cp.returncode == 0 and cp.stdout:
                # Parse version from output
                for line in cp.stdout.splitlines():
                    if "GATK" in line or "gatk" in line:
                        parts = line.split()
                        for part in parts:
                            if part[0].isdigit() and "." in part:
                                return part
        except Exception:
            pass

        return fallback

    ###########################################################################
    # Mark duplicates
    ###########################################################################

    @subroutine
    def mark_duplicates(
        self,
        input_bam: Path | str,
        output_bam: Path | str,
        metrics_file: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Mark duplicate reads in a BAM file using GATK MarkDuplicates.

        Creates a zero-argument callable that marks PCR/optical duplicates
        and writes the output BAM with duplicates flagged.

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_bam : Path | str
            Output BAM file with duplicates marked.
        metrics_file : Path | str
            Output metrics file containing duplication statistics.
        params : Params | None
            Additional GATK parameters (e.g., Params(remove_duplicates=True)).
        cfg : ExternalRunConfig | None
            Configuration for running the command.


        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes MarkDuplicates.

        Examples
        --------
        >>> gatk = GATK()
        >>> mark_dups = gatk.mark_duplicates("in.bam", "dedup.bam", "metrics.txt")
        >>> mark_dups()
        """
        # Build command: gatk MarkDuplicates -I input.bam -O output.bam -M metrics.txt
        arguments = [
            "MarkDuplicates",
            "-I",
            self.abs(input_bam),
            "-O",
            self.abs(output_bam),
            "-M",
            self.abs(metrics_file),
        ]
        return arguments, [input_bam, output_bam, metrics_file], None, None, None

    @element
    def mark(
        self,
        mapped: MappedElement,
        *,
        output_bam: Path | str = None,
        metrics_file: Path | str = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """
        High-level function for mark duplicate. Marks duplicate reads in a BAM
        file using GATK MarkDuplicates.

        Creates an Element object that can be run to invoke the low-level
        mark_duplicates function.

        Parameters
        ----------
        mapped : Element
            Element that contains the input BAM file.
        output_bam : Path | str
            Output BAM file with duplicates marked.
        metrics_file : Path | str
            Output metrics file containing duplication statistics.
        params : Params | None
            Additional GATK parameters (e.g., Params(remove_duplicates=True)).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that can be run and will produce the marked duplicates
            BAM and metrics file.

        Examples
        --------
        >>> gatk = GATK()
        >>> marked = gatk.mark(mapped, "dedup.bam", "metrics.txt")
        """

        output_bam = output_bam or mapped.bam.parent / f"{mapped.bam.stem}_marked.bam"
        metrics_file = (
            metrics_file or mapped.bam.parent / f"{mapped.bam.stem}_markdup_metrics.txt"
        )

        runner = self.mark_duplicates(
            input_bam=mapped.bam,
            output_bam=output_bam,
            metrics_file=metrics_file,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        return Element(
            name=f"{mapped.bam.stem}_marked",
            key=f"{mapped.bam.stem}_mark_duplicates_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[mapped.bam],
            artifacts={"bam": output_bam, "metrics": metrics_file},
            pres=[mapped],
        )

    ###########################################################################
    # Learn Read Orientation Model
    ###########################################################################

    @subroutine
    def learn_read_orientation(
        self,
        f1r2_tar_gz: Path | str | list[Path | str],
        output_model: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Learn read orientation model for filtering artifacts.

        Creates a zero-argument callable that learns orientation bias from
        F1R2 metrics collected during Mutect2 runs (useful for FFPE samples).

        Parameters
        ----------
        f1r2_tar_gz : Path | str | list[Path | str]
            Input F1R2 metrics file(s) in .tar.gz format.
            Can be a single file or list of files.
        output_model : Path | str
            Output orientation model file (.tar.gz).
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes LearnReadOrientationModel.

        Examples
        --------
        >>> gatk = GATK()
        >>> learn = gatk.learn_read_orientation("f1r2.tar.gz", "model.tar.gz")
        >>> learn()
        """
        if isinstance(f1r2_tar_gz, (list, tuple)):
            f1r2_files = [Path(f) for f in f1r2_tar_gz]
        else:
            f1r2_files = [Path(f1r2_tar_gz)]

        # Build command: gatk LearnReadOrientationModel -I f1r2.tar.gz -O model.tar.gz
        arguments = ["LearnReadOrientationModel"]
        for f in f1r2_files:
            arguments.extend(["-I", self.strabs(f)])
        arguments.extend(["-O", self.strabs(output_model)])
        return arguments, [*f1r2_files, output_model], None, None, None

    @element
    def readorientation(
        self,
        mutect_element: Element,
        *,
        output_model_file: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for learning read orientation model.

        Creates an Element that learns orientation bias from F1R2 metrics
        collected during Mutect2 runs (useful for FFPE samples).

        Parameters
        ----------
        mutect_element : Element
            Element containing the F1R2 metrics file from a Mutect2 run.
        output_model_file : Path | str | None
            Output orientation model file (.tar.gz). If not provided, defaults
            to "{bam_stem}_orientation.tar.gz" next to the input BAM.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that learns a read orientation model when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> orientation = gatk.readorientation(mutect_element)
        """
        output_model_file = (
            output_model_file
            or mutect_element.bam.parent
            / f"{mutect_element.bam.stem}_orientation.tar.gz"
        )
        runner = self.learn_read_orientation(
            f1r2_tar_gz=mutect_element.f1r2,
            output_model=output_model_file,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        return Element(
            name=f"{mutect_element.bam.stem}.lrom.{self.name}",
            key=f"{mutect_element.bam.stem}_learn_read_orientation_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[mutect_element.f1r2],
            artifacts={"orientation": output_model_file},
            pres=[mutect_element],
        )

    ###########################################################################
    # Mutect2
    ###########################################################################

    @subroutine
    def mutect2_bam(
        self,
        reference: Path | str,
        input_bams_tumor: list[Path | str],
        output_vcf: Path | str,
        input_bams_normal: list[Path | str] | None = None,
        tumor_sample: list[str] | None = None,
        normal_sample: list[str] | None = None,
        targets_padded_bed: Path | str | None = None,
        germline_resource: Path | str | None = None,
        panel_of_normals: Path | str | None = None,
        f1r2_tar_gz: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Run Mutect2 for somatic variant calling.

        Creates a zero-argument callable that runs GATK Mutect2 to call
        somatic mutations from tumor/normal BAM pairs.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bams_tumor : list[Path | str]
            List of tumor input BAM files.
        output_vcf : Path | str
            Output VCF file (will be gzipped).
        input_bams_normal : list[Path | str] | None
            List of normal input BAM files.
        tumor_sample : list[str] | None
            Tumor sample name(s) (required for paired calling).
        normal_sample : list[str] | None
            Normal sample name(s) (required for paired calling).
        targets_padded_bed : Path | str | None
            BED file or interval list restricting calling regions.
        germline_resource : Path | str | None
            Germline resource VCF for filtering.
        panel_of_normals : Path | str | None
            Panel of normals VCF.
        f1r2_tar_gz : Path | str | None
            Output file for F1R2 metrics (for orientation model).
        params: Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes Mutect2.

        Examples
        --------
        >>> gatk = GATK()
        >>> mutect = gatk.mutect2_bam(
        ...     reference="genome.fa",
        ...     input_bams_tumor=["tumor.bam"],
        ...     output_vcf="variants.vcf.gz",
        ...     tumor_sample=["TUMOR"],
        ...     normal_sample=["NORMAL"],
        ...     input_bams_normal=["normal.bam"],
        ... )
        >>> mutect()
        """
        # Build command: gatk Mutect2 -R ref.fa -I tumor.bam -I normal.bam ...
        arguments = ["Mutect2", "-R", self.strabs(reference)]

        for bam in input_bams_tumor:
            arguments.extend(["-I", self.strabs(bam)])
        if input_bams_normal:
            for bam in input_bams_normal:
                arguments.extend(["-I", self.strabs(bam)])

        if tumor_sample:
            for s in tumor_sample:
                arguments.extend(["-tumor", s])
        if normal_sample:
            for s in normal_sample:
                arguments.extend(["-normal", s])

        if targets_padded_bed:
            arguments.extend(["-L", self.strabs(targets_padded_bed)])
        if germline_resource:
            arguments.extend(["--germline-resource", self.strabs(germline_resource)])
        if panel_of_normals:
            arguments.extend(["--panel-of-normals", self.strabs(panel_of_normals)])
        if f1r2_tar_gz:
            arguments.extend(["--f1r2-tar-gz", self.strabs(f1r2_tar_gz)])

        arguments.extend(["-O", self.strabs(output_vcf)])

        all_paths = [reference, *input_bams_tumor, output_vcf]
        if input_bams_normal:
            all_paths.extend(input_bams_normal)
        if targets_padded_bed:
            all_paths.append(targets_padded_bed)
        if germline_resource:
            all_paths.append(germline_resource)
        if panel_of_normals:
            all_paths.append(panel_of_normals)
        if f1r2_tar_gz:
            all_paths.append(f1r2_tar_gz)

        return arguments, all_paths, None, None, None

    @element
    def mutect2(
        self,
        marked_tumor: MappedElement | list[MappedElement],
        reference: Genome,
        *,
        marked_normal: MappedElement | list[MappedElement] | None = None,
        targets_padded: Element | None = None,
        germline_resource: Element | None = None,
        panel_of_normals: Element | None = None,
        f1r2_tar_gz: Path | str | None = None,
        output_vcf: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for somatic variant calling with Mutect2.

        Creates an Element that runs GATK Mutect2 on the BAM file(s) from the
        given MappedElement(s).

        Parameters
        ----------
        marked_tumor : MappedElement | list[MappedElement]
            MappedElement(s) representing the tumor sample(s).
        reference : Genome
            Reference genome object containing the FASTA file.
        marked_normal : MappedElement | list[MappedElement] | None
            MappedElement(s) representing the normal sample(s).
        targets_padded : Element | None
            Element with a BED file or interval list restricting calling regions.
        germline_resource : Element | None
            Element with a germline resource VCF for filtering.
        panel_of_normals : Element | None
            Element with a panel of normals VCF.
        f1r2_tar_gz : Path | str | None
            Output file for F1R2 metrics (for orientation model).
        output_vcf : Path | str | None
            Output VCF file. If not provided, defaults to
            "{tumor_bam_stem}_mutect2.vcf.gz" next to the tumor BAM.
        params: Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes Mutect2 when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> mutect_element = gatk.mutect2(
        ...     mapped_tumor, genome, marked_normal=mapped_normal
        ... )
        """
        if isinstance(marked_tumor, MappedElement):
            marked_tumor = [marked_tumor]
        if isinstance(marked_normal, MappedElement):
            marked_normal = [marked_normal]

        input_bams_tumor = [m.bam for m in marked_tumor]
        input_bams_normal = [m.bam for m in marked_normal] if marked_normal else []
        tumor_sample_names = [m.name for m in marked_tumor]
        normal_sample_names = [m.name for m in marked_normal] if marked_normal else []

        output = (
            output_vcf
            or marked_tumor[0].bam.parent / f"{marked_tumor[0].bam.stem}_mutect2.vcf.gz"
        )
        f1r2 = Path(
            f1r2_tar_gz or output.parent / f"{Path(output).stem}_f1r2.tar.gz"
        ).absolute()

        runner = self.mutect2_bam(
            reference=reference.fasta,
            input_bams_tumor=input_bams_tumor,
            output_vcf=output,
            input_bams_normal=input_bams_normal or None,
            tumor_sample=tumor_sample_names,
            normal_sample=normal_sample_names or None,
            targets_padded_bed=targets_padded.bed if targets_padded else None,
            germline_resource=(
                germline_resource.snp
                if isinstance(germline_resource, Element)
                else germline_resource
            ),
            panel_of_normals=(
                panel_of_normals.pon
                if isinstance(panel_of_normals, Element)
                else panel_of_normals
            ),
            f1r2_tar_gz=f1r2,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        name = (
            f"{self.name}_mutect2_{tumor_sample_names[0]}"
            if tumor_sample_names
            else f"{self.name}_mutect2"
        )
        pres = [*marked_tumor, *(marked_normal or [])]
        if targets_padded:
            pres.append(targets_padded)
        if isinstance(germline_resource, Element):
            pres.append(germline_resource)
        if isinstance(panel_of_normals, Element):
            pres.append(panel_of_normals)
        return Element(
            name=name,
            key=(
                f"{self.name}_mutect2_{','.join(tumor_sample_names)}"
                if tumor_sample_names
                else f"{self.name}_mutect2"
            ),
            run=runner,
            determinants=determinants,
            inputs=input_bams_tumor + input_bams_normal,
            artifacts={"vcf": Path(output).absolute(), "f1r2": f1r2},
            pres=pres,
        )

    ###########################################################################
    # Pileup Summaries
    ###########################################################################

    @subroutine
    def get_pileup_summaries(
        self,
        input_bam: Path | str,
        variant_sites: Path | str,
        output_table: Path | str,
        *,
        intervals: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Get pileup summaries for contamination estimation.

        Creates a zero-argument callable that calculates allele counts at
        common variant sites for contamination estimation.

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        variant_sites : Path | str
            VCF of common variant sites (e.g., gnomAD).
        output_table : Path | str
            Output pileup table.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes GetPileupSummaries.

        Examples
        --------
        >>> gatk = GATK()
        >>> pileup = gatk.get_pileup_summaries(
        ...     input_bam="tumor.bam",
        ...     variant_sites="common_variants.vcf.gz",
        ...     output_table="tumor_pileup.table",
        ... )
        >>> pileup()
        """
        # gatk GetPileupSummaries -I input.bam -V sites.vcf -O output.table
        arguments = [
            "GetPileupSummaries",
            "-I",
            self.strabs(input_bam),
            "-V",
            self.strabs(variant_sites),
            "-O",
            self.strabs(output_table),
        ]
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])

        all_paths = [input_bam, variant_sites, output_table]
        if intervals:
            all_paths.append(intervals)
        return arguments, all_paths, None, None, None

    @element
    def pilesum(
        self,
        mapped: MappedElement,
        known_variants: Element | Path | str,
        targets_padded: Element | Path | str,
        *,
        output_table: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for pileup summaries for contamination estimation.

        Creates an Element that calculates allele counts at common variant
        sites for contamination estimation.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file and name.
        known_variants : Element | Path | str
            Element or file with VCF of common variant sites (e.g., gnomAD).
        targets_padded : Element | Path | str
            Element or file with a BED file or interval list restricting regions.
        output_table : Path | str | None
            Output pileup table. If not provided, defaults to
            "{bam_stem}_pileup.table" next to the input BAM.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes GetPileupSummaries when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> pileup_element = gatk.pilesum(mapped, known_variants, targets_padded)
        """
        input_bam = mapped.bam
        output_table = (
            output_table or input_bam.parent / f"{input_bam.stem}_pileup.table"
        )
        variant_sites = (
            known_variants.snp
            if isinstance(known_variants, Element)
            else known_variants
        )
        intervals = (
            targets_padded.bed
            if isinstance(targets_padded, Element)
            else targets_padded
        )

        runner = self.get_pileup_summaries(
            input_bam=input_bam,
            variant_sites=variant_sites,
            output_table=output_table,
            intervals=intervals,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        pres = [mapped]
        if isinstance(known_variants, Element):
            pres.append(known_variants)
        if isinstance(targets_padded, Element):
            pres.append(targets_padded)
        return Element(
            name=f"{self.name}_pileup_summaries_{mapped.bam.stem}",
            key=f"{self.name}_pileup_summaries_{mapped.bam.stem}",
            run=runner,
            determinants=determinants,
            inputs=[input_bam, variant_sites, intervals],
            artifacts={"pileup": Path(output_table).absolute()},
            pres=pres,
        )

    ###########################################################################
    # Contamination
    ###########################################################################

    @subroutine
    def calculate_contamination(
        self,
        tumor_pileup: Path | str,
        output_table: Path | str,
        *,
        normal_pileup: Path | str | None = None,
        output_segments: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Calculate contamination from pileup summaries.

        Creates a zero-argument callable that estimates cross-sample
        contamination from pileup tables.

        Parameters
        ----------
        tumor_pileup : Path | str
            Tumor pileup table from GetPileupSummaries.
        output_table : Path | str
            Output contamination table.
        normal_pileup : Path | str | None
            Normal pileup table (optional).
        output_segments : Path | str | None
            Output tumor segmentation table (optional).
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CalculateContamination.

        Examples
        --------
        >>> gatk = GATK()
        >>> contamination = gatk.calculate_contamination(
        ...     tumor_pileup="tumor_pileup.table",
        ...     output_table="contamination.table",
        ... )
        >>> contamination()
        """
        # gatk CalculateContamination -I tumor_pileup -O contamination.table
        arguments = [
            "CalculateContamination",
            "-I",
            self.strabs(tumor_pileup),
            "-O",
            self.strabs(output_table),
        ]
        if normal_pileup:
            arguments.extend(["-matched", self.strabs(normal_pileup)])
        if output_segments:
            arguments.extend(["--tumor-segmentation", self.strabs(output_segments)])

        all_paths = [tumor_pileup, output_table]
        if normal_pileup:
            all_paths.append(normal_pileup)
        if output_segments:
            all_paths.append(output_segments)
        return arguments, all_paths, None, None, None

    @element
    def contamination(
        self,
        tumor_pileup: Element,
        *,
        normal_pileup: Element | None = None,
        contamination_table: Path | str | None = None,
        output_segments: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for contamination estimation.

        Creates an Element that estimates cross-sample contamination from
        pileup tables.

        Parameters
        ----------
        tumor_pileup : Element
            Element containing the tumor pileup table from GetPileupSummaries.
        normal_pileup : Element | None
            Element containing the normal pileup table (optional).
        contamination_table : Path | str | None
            Output contamination table. If not provided, defaults to
            "{pileup_stem}.contamination.table" next to the pileup table.
        output_segments : Path | str | None
            Output tumor segmentation table (optional).
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes CalculateContamination when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> contamination_element = gatk.contamination(tumor_pileup)
        """
        tumor_table = tumor_pileup.pileup
        normal_table = normal_pileup.pileup if normal_pileup else None
        output_table = (
            contamination_table
            or (
                tumor_table.parent / f"{tumor_table.stem}.contamination.table"
            ).absolute()
        )

        runner = self.calculate_contamination(
            tumor_pileup=tumor_table,
            output_table=output_table,
            normal_pileup=normal_table,
            output_segments=output_segments,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        key = f"{self.name}_calculate_contamination_{tumor_pileup.key}"
        if normal_pileup:
            key += f"_vs_{normal_pileup.key}"

        pres = [tumor_pileup]
        inputs = [tumor_table]
        if normal_pileup:
            pres.append(normal_pileup)
            inputs.append(normal_table)
        return Element(
            name=f"{tumor_pileup.name}_contamination_{self.name}",
            key=key,
            run=runner,
            determinants=determinants,
            inputs=inputs,
            artifacts={"contamination": Path(output_table).absolute()},
            pres=pres,
        )

    ###########################################################################
    # Filter Mutect Calls
    ###########################################################################

    @subroutine
    def filter_mutect_calls(
        self,
        reference: Path | str,
        input_vcf: Path | str,
        output_vcf: Path | str,
        *,
        contamination_table: Path | str | None = None,
        orientation: Path | str | None = None,
        tumor_segmentation: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Filter Mutect2 calls using contamination and orientation models.

        Creates a zero-argument callable that applies GATK filtering to
        Mutect2 VCF output.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_vcf : Path | str
            Unfiltered Mutect2 VCF.
        output_vcf : Path | str
            Filtered output VCF.
        contamination_table : Path | str | None
            Contamination table from CalculateContamination.
        orientation : Path | str | None
            Read orientation model from LearnReadOrientationModel.
        tumor_segmentation : Path | str | None
            Tumor segmentation table from CalculateContamination.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes FilterMutectCalls.

        Examples
        --------
        >>> gatk = GATK()
        >>> filter_calls = gatk.filter_mutect_calls(
        ...     reference="genome.fa",
        ...     input_vcf="unfiltered.vcf.gz",
        ...     output_vcf="filtered.vcf.gz",
        ... )
        >>> filter_calls()
        """
        # Build command: gatk FilterMutectCalls -R ref.fa -V input.vcf -O output.vcf
        arguments = [
            "FilterMutectCalls",
            "-R",
            self.strabs(reference),
            "-V",
            self.strabs(input_vcf),
            "-O",
            self.strabs(output_vcf),
        ]
        if contamination_table:
            arguments.extend(
                ["--contamination-table", self.strabs(contamination_table)]
            )
        if orientation:
            arguments.extend(["--ob-priors", self.strabs(orientation)])
        if tumor_segmentation:
            arguments.extend(["--tumor-segmentation", self.strabs(tumor_segmentation)])

        all_paths = [reference, input_vcf, output_vcf]
        if contamination_table:
            all_paths.append(contamination_table)
        if orientation:
            all_paths.append(orientation)
        if tumor_segmentation:
            all_paths.append(tumor_segmentation)
        return arguments, all_paths, None, None, None

    @element
    def filter(
        self,
        mutect_element: Element,
        *,
        contamination_table: Element | None = None,
        orientation: Element | None = None,
        tumor_segmentation: Element | None = None,
        output_vcf: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for filtering Mutect2 calls.

        Creates an Element that applies GATK filtering to Mutect2 VCF output.

        Parameters
        ----------
        mutect_element : Element
            Element containing the unfiltered Mutect2 VCF.
        contamination_table : Element | None
            Element containing the contamination table from CalculateContamination.
        orientation : Element | None
            Element containing the read orientation model from LearnReadOrientationModel
        tumor_segmentation : Element | None
            Element containing the tumor segmentation table from CalculateContamination.
        output_vcf : Path | str | None
            Output filtered VCF file. If not provided, defaults to
            "{input_vcf_stem}_filtered.vcf.gz" next to the input VCF.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes FilterMutectCalls when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> filter_element = gatk.filter(mutect_element)
        """
        input_vcf = mutect_element.vcf
        output_vcf = Path(
            output_vcf or input_vcf.parent / f"{input_vcf.stem}_filtered.vcf.gz"
        ).absolute()

        runner = self.filter_mutect_calls(
            reference=mutect_element.reference.fasta,
            input_vcf=input_vcf,
            output_vcf=output_vcf,
            contamination_table=(
                contamination_table.contamination if contamination_table else None
            ),
            orientation=orientation.orientation if orientation else None,
            tumor_segmentation=(
                tumor_segmentation.segments if tumor_segmentation else None
            ),
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        key = f"{mutect_element.key}_filtered_{self.name}"
        pres = [mutect_element]
        inputs = [input_vcf]
        if contamination_table:
            key += f"_{contamination_table.key}"
            pres.append(contamination_table)
            inputs.append(contamination_table.contamination)
        if orientation:
            key += f"_{orientation.key}"
            pres.append(orientation)
            inputs.append(orientation.orientation)
        if tumor_segmentation:
            key += f"_{tumor_segmentation.key}"
            pres.append(tumor_segmentation)
            inputs.append(tumor_segmentation.segments)
        return Element(
            name=f"{mutect_element.name}_filtered",
            key=key,
            run=runner,
            determinants=determinants,
            inputs=inputs,
            artifacts={"vcf": output_vcf},
            pres=pres,
        )

    ###########################################################################
    # BQSR
    ###########################################################################

    @subroutine
    def basecalibrate_bam(
        self,
        reference: Path | str,
        input_bams: list[Path | str],
        output_table: Path | str,
        known_sites: Path | str | list[Path | str] | None = None,
        intervals: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Run GATK BaseRecalibrator on BAM file(s).

        Creates a zero-argument callable that computes a BQSR recalibration
        table from one or more input BAMs.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bams : list[Path | str]
            List of input BAM files.
        output_table : Path | str
            Output recalibration table file.
        known_sites : Path | str | list[Path | str] | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes BaseRecalibrator.

        Examples
        --------
        >>> gatk = GATK()
        >>> recal = gatk.basecalibrate_bam(
        ...     reference="genome.fa",
        ...     input_bams=["sample.bam"],
        ...     output_table="recal.table",
        ... )
        >>> recal()
        """
        # gatk BaseRecalibrator -R ref.fa -I input.bam --known-sites sites.vcf -O recal.table  # noqa: E501
        arguments = ["BaseRecalibrator", "-R", self.strabs(reference)]
        for bam in input_bams:
            arguments.extend(["-I", self.strabs(bam)])

        all_paths = [reference, *input_bams, output_table]

        if known_sites:
            if isinstance(known_sites, (list, tuple)):
                for ks in known_sites:
                    arguments.extend(["--known-sites", self.strabs(ks)])
                    all_paths.append(ks)
            else:
                arguments.extend(["--known-sites", self.strabs(known_sites)])
                all_paths.append(known_sites)
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])
            all_paths.append(intervals)

        arguments.extend(["-O", self.strabs(output_table)])
        return arguments, all_paths, None, None, None

    @element
    def baserecalibrate(
        self,
        mapped: MappedElement | list[MappedElement],
        reference: Genome,
        *,
        output_table: Path | str | None = None,
        known_sites: Path | str | list[Path | str] | None = None,
        intervals: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for GATK BaseRecalibrator.

        Creates an Element that computes a BQSR recalibration table from one
        or more input BAMs.

        Parameters
        ----------
        mapped : MappedElement | list[MappedElement]
            MappedElement(s) containing the BAM file(s) to recalibrate.
        reference : Genome
            Reference genome object containing the FASTA file.
        output_table : Path | str | None
            Output recalibration table file. If not provided, defaults to
            "{bam_stem}.bsqr.table" next to the first input BAM.
        known_sites : Path | str | list[Path | str] | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes BaseRecalibrator when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> recal_element = gatk.baserecalibrate(mapped, genome)
        """
        mapped_list = list(mapped) if isinstance(mapped, (list, tuple)) else [mapped]
        input_bams = [m.bam for m in mapped_list]
        first_bam = input_bams[0]
        output_table = (
            Path(output_table).absolute()
            if output_table
            else first_bam.parent / f"{first_bam.stem}.bsqr.table"
        )
        runner = self.basecalibrate_bam(
            reference=reference.fasta if hasattr(reference, "fasta") else reference,
            input_bams=input_bams,
            output_table=output_table,
            known_sites=known_sites,
            intervals=intervals,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        name = f"{self.name}_basecalibrate_{first_bam.stem}"
        return Element(
            name=name,
            key=name,
            run=runner,
            determinants=determinants,
            inputs=input_bams,
            artifacts={"bsqr_table": Path(output_table).absolute()},
            pres=mapped_list,
        )

    @subroutine
    def apply_bqsr_on_bam(
        self,
        reference: Path | str,
        input_bam: Path | str,
        recal_table: Path | str,
        output_bam: Path | str,
        *,
        intervals: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Run GATK ApplyBQSR on a BAM file.

        Creates a zero-argument callable that applies a pre-computed BQSR
        recalibration table to a BAM file.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bam : Path | str
            Input BAM file.
        recal_table : Path | str
            Recalibration table from BaseRecalibrator.
        output_bam : Path | str
            Output recalibrated BAM file.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes ApplyBQSR.

        Examples
        --------
        >>> gatk = GATK()
        >>> apply = gatk.apply_bqsr_on_bam(
        ...     reference="genome.fa",
        ...     input_bam="sample.bam",
        ...     recal_table="recal.table",
        ...     output_bam="recalibrated.bam",
        ... )
        >>> apply()
        """
        # gatk ApplyBQSR -R ref.fa -I input.bam --bqsr-recal-file recal.table -O output.bam  # noqa: E501
        arguments = [
            "ApplyBQSR",
            "-R",
            self.strabs(reference),
            "-I",
            self.strabs(input_bam),
            "--bqsr-recal-file",
            self.strabs(recal_table),
            "-O",
            self.strabs(output_bam),
        ]
        print(intervals)
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])

        all_paths = [reference, input_bam, recal_table, output_bam]
        if intervals:
            all_paths.append(intervals)
        return arguments, all_paths, None, None, None

    @element
    def bsqr(
        self,
        element: MappedElement,
        reference: Genome,
        *,
        output_bam: Path | str | None = None,
        known_sites: Path | str | list[Path | str] | None = None,
        intervals: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for Base Quality Score Recalibration (BQSR).

        Creates an Element that first computes the BQSR recalibration table
        and then applies it to the input BAM.

        Parameters
        ----------
        element : MappedElement
            Input MappedElement containing the BAM file to recalibrate.
        reference : Genome
            Reference genome object containing the FASTA file.
        output_bam : Path | str | None
            Output recalibrated BAM file. If not provided, defaults to
            "{bam_stem}.bsqr.bam" next to the input BAM.
        known_sites : Path | str | list[Path | str] | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes BQSR when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> bsqr_element = gatk.bsqr(mapped, genome)
        """
        output = Path(output_bam or element.bam.parent / f"{element.bam.stem}.bsqr.bam")
        output_table = output.with_suffix(".table")
        if isinstance(intervals, Element):
            intervals = intervals.bed
        if isinstance(known_sites, Element):
            known_sites = known_sites.snp
        print("intervlas", intervals)

        recalibration = self.baserecalibrate(
            mapped=element,
            reference=reference,
            output_table=output_table,
            known_sites=known_sites,
            intervals=intervals,
            params=params,
            cfg=cfg,
        )
        bsqr_runner = self.apply_bqsr_on_bam(
            reference=reference.fasta,
            input_bam=element.bam,
            recal_table=recalibration.bsqr_table,
            output_bam=output,
            intervals=intervals,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        known_sites_names = (
            "_".join(
                [
                    Path(ks).stem
                    for ks in (
                        known_sites
                        if isinstance(known_sites, (list, tuple))
                        else [known_sites]
                    )
                ]
            )
            if known_sites
            else "no_known_sites"
        )
        interval_name = Path(intervals).stem if intervals else "no_intervals"
        return Element(
            name=f"{self.name}_apply_bqsr_{element.bam.stem}",
            key=f"{self.name}_apply_bqsr_{element.bam.stem}_{known_sites_names}_{interval_name}",  # noqa: E501
            run=bsqr_runner,
            determinants=determinants,
            inputs=[element.bam, recalibration.bsqr_table],
            artifacts={"bam": output.absolute()},
            pres=[element, recalibration],
        )

    ###########################################################################
    # Panel of Normals
    ###########################################################################

    @subroutine
    def create_panel_of_normals(
        self,
        input_vcfs: list[Path | str],
        output_pon: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Create a Panel of Normals using GATK CreateSomaticPanelOfNormals.

        Creates a zero-argument callable that builds a panel of normals VCF
        from multiple normal sample VCFs.

        Parameters
        ----------
        input_vcfs : list[Path | str]
            List of normal sample VCF files from Mutect2 tumor-only calls.
        output_pon : Path | str
            Output VCF file for the panel of normals.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CreateSomaticPanelOfNormals.

        Examples
        --------
        >>> gatk = GATK()
        >>> pon = gatk.create_panel_of_normals(
        ...     input_vcfs=["normal1.vcf.gz", "normal2.vcf.gz"],
        ...     output_pon="pon.vcf.gz",
        ... )
        >>> pon()
        """
        # gatk CreateSomaticPanelOfNormals -vcfs n1.vcf -vcfs n2.vcf -O pon.vcf.gz
        arguments = ["CreateSomaticPanelOfNormals"]
        for vcf in input_vcfs:
            arguments.extend(["-vcfs", self.strabs(vcf)])
        arguments.extend(["-O", self.strabs(output_pon)])
        return arguments, [*input_vcfs, output_pon], None, None, None

    @element
    def pons(
        self,
        all_normals: list[Element],
        reference: Genome,
        *,
        output_pon: Path | str | None = None,
        params_mutect2: Params | None = None,
        params_pon: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for building a Panel of Normals.

        Creates an Element that first calls each normal through Mutect2 in
        tumor-only mode, then combines the results into a panel of normals.

        Parameters
        ----------
        all_normals : list[Element]
            List of normal Elements. Each must have a 'bam' artifact and a
            'name' attribute.
        reference : Genome
            Reference genome to use for the analysis.
        output_pon : Path | str | None
            Output file for the panel of normals. If not provided, defaults to
            "{version}_{reference}.pon.vcf.gz" in the cache directory.
        params_mutect2 : Params | None
            Additional GATK parameters for the Mutect2 step.
        params_pon : Params | None
            Additional GATK parameters for CreateSomaticPanelOfNormals.
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            Element representing the panel of normals.

        Raises
        ------
        ValueError
            If any normal Element is missing a 'bam' artifact.

        Examples
        --------
        >>> gatk = GATK()
        >>> pon_element = gatk.pons(all_normals, genome)
        """
        output = Path(
            output_pon or self.cache_dir / f"{self.name}_{reference.name}.pon.vcf.gz"
        )
        pres_pon = []
        for normal in all_normals:
            if "bam" not in normal.artifacts:
                raise ValueError(
                    f"Each normal Element must have a 'bam' artifact. Missing in {normal.name}"  # noqa: E501
                )
            output_vcf_normal = output.parent / f"{normal.name}.mutect2.vcf.gz"
            called = self.mutect2(
                marked_tumor=normal,
                reference=reference,
                output_vcf=output_vcf_normal,
                params=params_mutect2,
                cfg=cfg,
            )
            pres_pon.append(called)

        input_vcfs = [e.vcf for e in pres_pon]
        runner = self.create_panel_of_normals(
            input_vcfs=input_vcfs,
            output_pon=output,
            params=params_pon,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params_pon)
        name = f"{self.name}_{reference.name}_build_pon"
        return Element(
            name=name,
            key=f"{name}_from_{','.join([n.name for n in all_normals])}",
            run=runner,
            determinants=determinants,
            inputs=input_vcfs,
            artifacts={"pon": output.absolute()},
            pres=pres_pon,
        )

    ###########################################################################
    # Sequence Dictionary
    ###########################################################################

    @subroutine
    def create_sequence_dictionary(
        self,
        reference: Path | str,
        output_dict: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Create a reference sequence dictionary using GATK CreateSequenceDictionary.

        Creates a zero-argument callable that generates a .dict file from a
        reference FASTA file.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        output_dict : Path | str
            Output .dict file path.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CreateSequenceDictionary.

        Examples
        --------
        >>> gatk = GATK()
        >>> create_dict = gatk.create_sequence_dictionary("genome.fa", "genome.dict")
        >>> create_dict()
        """
        # Build command: gatk CreateSequenceDictionary -R ref.fa -O ref.dict
        arguments = [
            "CreateSequenceDictionary",
            "-R",
            self.strabs(reference),
            "-O",
            self.strabs(output_dict),
        ]
        return arguments, [reference, output_dict], None, None, None

    @element
    def sequence_dict(
        self,
        reference: Genome,
        *,
        output_dict: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function to create a reference sequence dictionary.

        Creates an Element that generates a .dict file from a reference FASTA
        using GATK CreateSequenceDictionary.

        Parameters
        ----------
        reference : Genome
            Reference genome object containing the FASTA file.
        output_dict : Path | str | None
            Output .dict file path. If not provided, defaults to placing the
            .dict file next to the FASTA file.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that creates the sequence dictionary when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> dict_element = gatk.sequence_dict(genome)
        """
        output_dict = (
            Path(output_dict).absolute()
            if output_dict
            else reference.fasta.with_suffix(".dict")
        )
        print("output_dict", output_dict)
        runner = self.create_sequence_dictionary(
            reference=reference.fasta,
            output_dict=output_dict,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        name = f"{reference.name}_sequence_dict"
        return Element(
            name=name,
            key=f"{name}_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[reference.fasta],
            artifacts={"dict": output_dict.absolute()},
            pres=[],
        )

    ###########################################################################
    # BED to Interval List
    ###########################################################################

    @subroutine
    def bed_to_interval_list(
        self,
        input_bed: Path | str,
        output_interval_list: Path | str,
        sequence_dict: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Convert a BED file to Picard interval_list format.

        Creates a zero-argument callable that converts a BED file to the
        interval_list format required by Picard tools like CollectHsMetrics.

        Parameters
        ----------
        input_bed : Path | str
            Input BED file to convert.
        output_interval_list : Path | str
            Output interval_list file.
        sequence_dict : Path | str
            Reference sequence dictionary (.dict file).
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes BedToIntervalList.

        Examples
        --------
        >>> gatk = GATK()
        >>> converter = gatk.bed_to_interval_list(
        ...     input_bed="targets.bed",
        ...     output_interval_list="targets.interval_list",
        ...     sequence_dict="genome.dict",
        ... )
        >>> converter()
        """
        # Build command: gatk BedToIntervalList -I bed -O interval_list -SD dict
        arguments = [
            "BedToIntervalList",
            "-I",
            self.strabs(input_bed),
            "-O",
            self.strabs(output_interval_list),
            "-SD",
            self.strabs(sequence_dict),
        ]
        return (
            arguments,
            [input_bed, output_interval_list, sequence_dict],
            None,
            None,
            None,
        )

    @element
    def bed2interval(
        self,
        element: Element | str | Path,
        sequence_dict: Element,
        *,
        output_list: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function to convert a BED file to Picard interval_list format.

        Creates an Element that converts a BED file to the interval_list
        format required by Picard tools like CollectHsMetrics.

        Parameters
        ----------
        element : Element | str | Path
            BED file or Element containing the BED file to convert.
        sequence_dict : Element
            Element containing the reference sequence dictionary (.dict file).
        output_list : Path | str | None
            Output interval_list file. If not provided, defaults to
            "{bed_stem}.interval_list" next to the BED file.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that converts the BED to interval_list when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> interval_element = gatk.bed2interval(bed_element, dict_element)
        """
        input_bed = element.bed if isinstance(element, Element) else Path(element)
        output_interval_list = Path(
            output_list
            or Path(input_bed).parent / f"{Path(input_bed).stem}.interval_list"
        ).absolute()

        runner = self.bed_to_interval_list(
            input_bed=input_bed,
            output_interval_list=output_interval_list,
            sequence_dict=sequence_dict.dict,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        name = f"{Path(input_bed).stem}_to_interval_list"
        pres = [element] if isinstance(element, Element) else []
        pres.append(sequence_dict)
        return Element(
            name=name,
            key=f"{name}_{sequence_dict.key}_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[input_bed],
            artifacts={"intervals": output_interval_list},
            pres=pres,
        )

    ###########################################################################
    # Alignment Summary Metrics
    ###########################################################################

    @subroutine
    def alignment_summary_metrics(
        self,
        reference: Path | str,
        input_bam: Path | str,
        output_metrics: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Collect alignment summary metrics using Picard CollectAlignmentSummaryMetrics

        Creates a zero-argument callable that collects comprehensive alignment
        metrics including mapping rates, quality scores, and mismatch rates.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bam : Path | str
            Input BAM file.
        output_metrics : Path | str
            Output metrics file.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CollectAlignmentSummaryMetrics.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics = gatk.alignment_summary_metrics(
        ...     reference="genome.fa",
        ...     input_bam="sample.bam",
        ...     output_metrics="qc/alignment_summary.txt",
        ... )
        >>> metrics()
        """
        # gatk CollectAlignmentSummaryMetrics -R ref.fa -I input.bam -O metrics.txt
        arguments = [
            "CollectAlignmentSummaryMetrics",
            "-R",
            self.strabs(reference),
            "-I",
            self.strabs(input_bam),
            "-O",
            self.strabs(output_metrics),
        ]
        return arguments, [reference, input_bam, output_metrics], None, None, None

    @element
    def alignment_summary(
        self,
        mapped: MappedElement,
        reference: Genome,
        *,
        output_metrics: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for alignment summary metrics.

        Creates an Element that collects alignment summary metrics including
        mapping rates, quality scores, and mismatch rates.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        reference : Genome
            Reference genome object.
        output_metrics : Path | str | None
            Output metrics file. If not provided, defaults to
            "{bam_stem}_alignment_summary.txt" in a qc/picard subdirectory.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes CollectAlignmentSummaryMetrics when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics_elem = gatk.alignment_summary(mapped, genome)
        """
        input_bam = mapped.bam
        qc_folder = (
            Path(output_metrics).parent
            if output_metrics
            else input_bam.parent / "qc" / "picard"
        )
        output_metrics = Path(
            output_metrics or qc_folder / f"{input_bam.stem}_alignment_summary.txt"
        ).absolute()

        runner = self.alignment_summary_metrics(
            reference=reference.fasta,
            input_bam=input_bam,
            output_metrics=output_metrics,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        return Element(
            name=f"{self.name}_alignment_summary_{input_bam.stem}",
            key=f"{input_bam.stem}_alignment_summary_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[input_bam, reference.fasta],
            artifacts={"metrics": output_metrics},
            pres=[mapped],
        )

    ###########################################################################
    # Insert Size Metrics
    ###########################################################################

    @subroutine
    def insert_size_metrics(
        self,
        input_bam: Path | str,
        output_metrics: Path | str,
        output_histogram: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Collect insert size metrics using Picard CollectInsertSizeMetrics.

        Creates a zero-argument callable that collects insert size distribution
        metrics for paired-end sequencing data.

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_metrics : Path | str
            Output metrics file.
        output_histogram : Path | str
            Output histogram PDF file.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CollectInsertSizeMetrics.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics = gatk.insert_size_metrics(
        ...     input_bam="sample.bam",
        ...     output_metrics="qc/insert_size.txt",
        ...     output_histogram="qc/insert_size.pdf",
        ... )
        >>> metrics()
        """
        # gatk CollectInsertSizeMetrics -I input.bam -O metrics.txt -H histogram.pdf
        arguments = [
            "CollectInsertSizeMetrics",
            "-I",
            self.strabs(input_bam),
            "-O",
            self.strabs(output_metrics),
            "-H",
            self.strabs(output_histogram),
        ]
        return (
            arguments,
            [input_bam, output_metrics, output_histogram],
            None,
            None,
            None,
        )

    @element
    def insertsize(
        self,
        mapped: MappedElement,
        *,
        output_metrics: Path | str | None = None,
        output_histogram: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for insert size metrics.

        Creates an Element that collects insert size distribution metrics
        for paired-end sequencing data.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        output_metrics : Path | str | None
            Output metrics file. If not provided, defaults to
            "{bam_stem}_insert_size.txt" in a qc/picard subdirectory.
        output_histogram : Path | str | None
            Output histogram PDF. If not provided, defaults to
            "{bam_stem}_insert_size.pdf" in a qc/picard subdirectory.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes CollectInsertSizeMetrics when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics_elem = gatk.insertsize(mapped)
        """
        input_bam = mapped.bam
        qc_folder = (
            Path(output_metrics).parent
            if output_metrics
            else input_bam.parent / "qc" / "picard"
        )
        output_metrics = Path(
            output_metrics or qc_folder / f"{input_bam.stem}_insert_size.txt"
        ).absolute()
        output_histogram = Path(
            output_histogram or qc_folder / f"{input_bam.stem}_insert_size.pdf"
        ).absolute()

        runner = self.insert_size_metrics(
            input_bam=input_bam,
            output_metrics=output_metrics,
            output_histogram=output_histogram,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        return Element(
            name=f"{input_bam.stem}_insert_size_{self.name}",
            key=f"{input_bam.stem}_insert_size_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[input_bam],
            artifacts={"metrics": output_metrics, "histogram": output_histogram},
            pres=[mapped],
        )

    ###########################################################################
    # HS Metrics
    ###########################################################################

    @subroutine
    def hs_metrics(
        self,
        reference: Path | str,
        input_bam: Path | str,
        output_metrics: Path | str,
        bait_intervals: Path | str,
        target_intervals: Path | str,
        per_target_coverage: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Collect hybrid selection (HS) metrics using Picard CollectHsMetrics.

        Creates a zero-argument callable that collects metrics for targeted
        sequencing including on-target rate, coverage uniformity, and fold
        enrichment.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bam : Path | str
            Input BAM file.
        output_metrics : Path | str
            Output metrics file.
        bait_intervals : Path | str
            Interval list of bait regions.
        target_intervals : Path | str
            Interval list of target regions.
        per_target_coverage : Path | str | None
            Optional output file for per-target coverage.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CollectHsMetrics.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics = gatk.hs_metrics(
        ...     reference="genome.fa",
        ...     input_bam="sample.bam",
        ...     output_metrics="qc/hs_metrics.txt",
        ...     bait_intervals="targets.interval_list",
        ...     target_intervals="targets.interval_list",
        ... )
        >>> metrics()
        """
        # gatk CollectHsMetrics -R ref.fa -I bam -O metrics --BAIT_INTERVALS baits --TARGET_INTERVALS targets  # noqa: E501
        arguments = [
            "CollectHsMetrics",
            "-R",
            self.strabs(reference),
            "-I",
            self.strabs(input_bam),
            "-O",
            self.strabs(output_metrics),
            "--BAIT_INTERVALS",
            self.strabs(bait_intervals),
            "--TARGET_INTERVALS",
            self.strabs(target_intervals),
        ]
        all_paths = [
            reference,
            input_bam,
            output_metrics,
            bait_intervals,
            target_intervals,
        ]

        if per_target_coverage:
            arguments.extend(
                ["--PER_TARGET_COVERAGE", self.strabs(per_target_coverage)]
            )
            all_paths.append(per_target_coverage)

        return arguments, all_paths, None, None, None

    @element
    def hs(
        self,
        mapped: MappedElement,
        reference: Genome,
        baits: Element | Path | str,
        targets: Element | Path | str,
        *,
        output_metrics: Path | str | None = None,
        per_target_coverage: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function for hybrid selection metrics.

        Creates an Element that collects metrics for targeted sequencing
        including on-target rate, coverage uniformity, and fold enrichment.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        reference : Genome
            Reference genome object.
        baits : Element | Path | str
            Element or file with interval list of bait regions.
        targets : Element | Path | str
            Element or file with interval list of target regions.
        output_metrics : Path | str | None
            Output metrics file. If not provided, defaults to
            "{bam_stem}_hs_metrics.txt" in a qc/picard subdirectory.
        per_target_coverage : Path | str | None
            Optional output file for per-target coverage.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes CollectHsMetrics when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> metrics_elem = gatk.hs(mapped, genome, bait_element, target_element)
        """
        input_bam = mapped.bam
        qc_dir = (
            Path(output_metrics).parent
            if output_metrics
            else input_bam.parent / "qc" / "picard"
        )
        output_metrics = Path(
            output_metrics or qc_dir / f"{input_bam.stem}_hs_metrics.txt"
        ).absolute()

        bait_path = baits.intervals if isinstance(baits, Element) else baits
        target_path = targets.intervals if isinstance(targets, Element) else targets

        runner = self.hs_metrics(
            reference=reference.fasta,
            input_bam=input_bam,
            output_metrics=output_metrics,
            bait_intervals=bait_path,
            target_intervals=target_path,
            per_target_coverage=per_target_coverage,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        pres = [mapped]
        if isinstance(baits, Element):
            pres.append(baits)
        if isinstance(targets, Element) and targets != baits:
            pres.append(targets)

        artifacts = {"metrics": output_metrics}
        if per_target_coverage:
            artifacts["per_target_coverage"] = Path(per_target_coverage).absolute()

        return Element(
            name=f"{input_bam.stem}_hs_metrics_{self.name}",
            key=f"{input_bam.stem}_hs_metrics_{self.name}",
            run=runner,
            determinants=determinants,
            inputs=[input_bam, reference.fasta, bait_path, target_path],
            artifacts=artifacts,
            pres=pres,
        )

    ###########################################################################
    # Callable Loci
    ###########################################################################

    @subroutine
    def callable_loci(
        self,
        reference_fasta: Path | str,
        input_bam: Path | str,
        output_bed: Path | str,
        output_summary: Path | str,
        target_bed: Path | str,
        min_depth: int = 10,
        min_mapping_quality: int = 20,
        min_base_quality: int = 20,
        output_json: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Collect callable loci regions using GATK CallableLoci.

        Creates a zero-argument callable that identifies regions in the
        categories CALLABLE/NO_COVERAGE/LOW_COVERAGE etc. Used to determine
        whether a SNP is in a callable region to filter artifacts.

        Parameters
        ----------
        reference_fasta : Path | str
            Reference genome FASTA file.
        input_bam : Path | str
            Input BAM file.
        output_bed : Path | str
            Output BED file with callable loci categories.
        output_summary : Path | str
            Output summary file with callable loci statistics.
        target_bed : Path | str
            BED file with target regions to evaluate.
        min_depth : int
            Minimum depth to consider a locus callable (default: 10).
        min_mapping_quality : int
            Minimum mapping quality to consider a read (default: 20).
        min_base_quality : int
            Minimum base quality to consider a read (default: 20).
        output_json : Path | str | None
            Optional output JSON file with callable bases summary.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes CallableLoci.

        Examples
        --------
        >>> gatk = GATK()
        >>> loci = gatk.callable_loci(
        ...     reference_fasta="genome.fa",
        ...     input_bam="sample.bam",
        ...     output_bed="callable_loci.bed",
        ...     output_summary="callable_loci_summary.txt",
        ...     target_bed="targets.bed",
        ... )
        >>> loci()
        """
        # gatk CallableLoci -R ref -I bam -L target -summary summary -o bed
        arguments = [
            "CallableLoci",
            "-R",
            self.strabs(reference_fasta),
            "-I",
            self.strabs(input_bam),
            "-L",
            self.strabs(target_bed),
            "--minDepth",
            str(min_depth),
            "--minMappingQuality",
            str(min_mapping_quality),
            "--minBaseQuality",
            str(min_base_quality),
            "-summary",
            self.strabs(output_summary),
            "-o",
            self.strabs(output_bed),
        ]
        all_paths = [reference_fasta, input_bam, target_bed, output_bed, output_summary]
        if output_json:
            all_paths.append(output_json)

        def _post():
            callable_bases, callable_mb = self.extract_callable_mb(output_summary)
            ret = {
                "method": "CallableLoci",
                "min_dp": str(min_depth),
                "min_mq": str(min_mapping_quality),
                "min_bq": str(min_base_quality),
                "callable_bases": str(callable_bases),
                "callable_mb": str(callable_mb),
            }
            if output_json:
                with open(output_json, "w") as f:
                    json.dump(ret, f, indent=2)
            return ret

        return arguments, all_paths, None, None, _post

    @element
    def loci(
        self,
        mapped: MappedElement,
        reference: Genome,
        *,
        output_bed: Path | str | None = None,
        output_summary: Path | str | None = None,
        targets: Element | Path | str | None = None,
        min_depth: int = 10,
        min_mapping_quality: int = 20,
        min_base_quality: int = 20,
        output_json: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """High-level function to collect callable loci regions.

        Creates an Element that identifies CALLABLE/NO_COVERAGE/LOW_COVERAGE
        regions to determine whether SNPs are in callable regions.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        reference : Genome
            Reference genome object.
        output_bed : Path | str | None
            Output BED file with callable loci. If not provided, defaults to
            "{bam_stem}.callable_loci.bed" next to the input BAM.
        output_summary : Path | str | None
            Output summary file. If not provided, defaults to
            "{bam_stem}.callable_loci_summary.txt" next to the input BAM.
        targets : Element | Path | str | None
            Optional Element or file with target regions to evaluate.
        min_depth : int
            Minimum depth to consider a locus callable (default: 10).
        min_mapping_quality : int
            Minimum mapping quality to consider a read (default: 20).
        min_base_quality : int
            Minimum base quality to consider a read (default: 20).
        output_json : Path | str | None
            Optional output JSON file with callable bases summary.
        params : Params | None
            Additional GATK parameters (e.g., Params()).
        cfg : ExternalRunConfig | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes CallableLoci when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> loci_element = gatk.loci(mapped, genome)
        """
        input_bam = mapped.bam
        output_summary = Path(
            output_summary
            or input_bam.parent / f"{input_bam.stem}.callable_loci_summary.txt"
        ).absolute()
        output_bed = Path(
            output_bed or input_bam.parent / f"{input_bam.stem}.callable_loci.bed"
        ).absolute()
        output_json = output_json or output_summary.with_suffix(".json")
        target_bed = targets.bed if isinstance(targets, Element) else targets

        inputs = [input_bam, reference.fasta]
        pres = [mapped]
        if isinstance(targets, Element):
            pres.append(targets)
            inputs.append(targets.bed)

        runner = self.callable_loci(
            reference_fasta=reference.fasta,
            input_bam=input_bam,
            output_bed=output_bed,
            output_summary=output_summary,
            target_bed=target_bed,
            min_depth=min_depth,
            min_mapping_quality=min_mapping_quality,
            min_base_quality=min_base_quality,
            output_json=output_json,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)
        name = f"{input_bam.stem}_callable_loci_{self.name}"
        key = f"{input_bam.stem}_callable_loci_{self.name}_{min_base_quality}_{min_mapping_quality}_{min_depth}"  # noqa: E501
        if target_bed:
            key += f"_{Path(target_bed).stem}"

        return Element(
            name=name,
            key=key,
            run=runner,
            determinants=determinants,
            inputs=inputs,
            artifacts={
                "summary": output_summary,
                "bed": output_bed,
                "callable": output_json,
            },
            pres=pres,
        )

    @staticmethod
    def extract_callable_mb(summary_file: Path) -> tuple[int, float]:
        callable_bases = 0

        with open(summary_file) as f:
            for line in f:
                if line.startswith("CALLABLE"):
                    callable_bases = int(line.split()[1])
                    break

        callable_mb = callable_bases / 1_000_000
        return callable_bases, callable_mb
