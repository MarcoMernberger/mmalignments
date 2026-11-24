"""Module contains a GATK interface for variant calling and preprocessing."""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping, Sequence

import pandas as pd

from mmalignments.models.aligners.samtools import (
    Samtools,  # type: ignore[import]
)
from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.genome import (
    EnsemblGenome as Genome,  # type: ignore[import]
)

from mmalignments.models.elements import (
    Element,
    FileSource,
    MappedElement,
    element,
)
from mmalignments.models.parameters import ParamSet
from mmalignments.models.overlay import (
    CfgType,
    ElementTag,
    OutSpec,
    FileSpec,
    OutType,
    ParType,
    SpecValues,
    TagType,
    pick_spec,
)
from mmalignments.models.tags import (
    Method,
    Omics,
    Stage,
    State,
)
from mmalignments.models.specs import Runnable, CallSpec
from mmalignments.services.io import read_frame, read_vcf, write_frame
from ..externals import External, SubroutineIn, subroutine

logger = logging.getLogger(__name__)


def _candidate_allele_metrics(
    vcf: pd.DataFrame,
    sample: str,
    wanted: set[tuple[str, int, str, str]],
    vcf_path: Path,
) -> dict[tuple[str, int, str, str], tuple[bool, int, int, int, float, float]]:
    """Read allele depths and the reported AF and FRAG_VAF for each match."""
    observations = {}
    columns = ["CHROM", "POS", "REF", "ALT", "FILTER", "FORMAT", sample]
    for chrom, pos, ref, alt, filters, fmt, values in vcf[columns].itertuples(
        index=False, name=None
    ):
        alt_alleles = alt.split(",")
        keys = [(chrom, pos, ref, allele) for allele in alt_alleles]
        if wanted.isdisjoint(keys):
            continue
        fields = dict(zip(fmt.split(":"), values.split(":")))
        ad = fields.get("AD", ".")
        if ad in {"", "."}:
            depths = [0] * (len(alt_alleles) + 1)
        else:
            depths = [0 if value == "." else int(value) for value in ad.split(",")]
            if len(depths) != len(alt_alleles) + 1 or min(depths) < 0:
                raise ValueError(f"Invalid AD {ad!r} at {chrom}:{pos} in {vcf_path}")
        ref_count = depths[0]
        total_alt = sum(depths[1:])
        fractions = []
        for field in ("AF", "FRAG_VAF"):
            raw = fields.get(field, ".")
            values = raw.split(",") if raw not in {"", "."} else ["."] * len(keys)
            if len(values) != len(keys):
                raise ValueError(
                    f"Invalid {field} {raw!r} at {chrom}:{pos} in {vcf_path}"
                )
            fractions.append(
                [
                    float(value) if value not in {"", "."} else float("nan")
                    for value in values
                ]
            )
        for key, alt_count, af, frag_vaf in zip(keys, depths[1:], *fractions):
            if key not in wanted:
                continue
            if key in observations:
                raise ValueError(f"Duplicate candidate allele {key} in {vcf_path}")
            observations[key] = (
                filters == "PASS",
                ref_count,
                alt_count,
                total_alt - alt_count,
                af,
                frag_vaf,
            )
    return observations


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

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get GATK version string.

        Parameters
        ----------
        fallback : str | None
            Value to return if version cannot be determined.

        Returns
        -------
        str | None
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
    # Helper
    ###########################################################################
    @staticmethod
    def check_mapped(element: Element):
        if not hasattr(element, "bam"):
            raise ValueError(
                f"Each normal Element must have a 'bam' artifact. Missing in {element.name}"  # noqa: E501
            )

    @staticmethod
    def index_by_sam(bam_file: str | Path) -> Callable[[], CompletedProcess]:
        """
        Convenience method for calling samtools to index a bam.

        Use as post-callable for GATK steps that produce BAM outputs, since
        GATK does not automatically index BAM files. This will create a
        zero-argument callable that can be chained as a post-processing step to
        index the BAM file using samtools.

        Parameters
        ----------
        bam_file : str | Path
            Path to the BAM file to be indexed.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that indexes the BAM file using samtools.
        """
        return Samtools().index_bam(bam_file)

    ###########################################################################
    # Mark duplicates
    ###########################################################################

    @subroutine
    def mark_duplicates(
        self,
        input_bam: Path | str,
        output_bam: Path | str,
        metrics_file: Path | str,
        *,
        index_off: bool = False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        index_off : bool
            If True, do not create BAM index after marking duplicates. Default
            is False.
        par : ParType | None
            Additional GATK parameters (e.g., Par(remove_duplicates=True)).
        cfg : CfgType | None
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
        subcommand = "MarkDuplicates"
        arguments = [
            subcommand,
            "-I",
            self.abs(input_bam),
            "-O",
            self.abs(output_bam),
            "-M",
            self.abs(metrics_file),
        ]
        post = self.index_by_sam(output_bam) if not index_off else None
        return (
            arguments,
            subcommand,
            [input_bam],
            [output_bam, metrics_file],
            None,
            None,
            post,
        )

    @element
    def mark(
        self,
        mapped: MappedElement,
        *,
        index_off: bool = False,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        High-level function for mark duplicate. Marks duplicate reads in a BAM
        file using GATK MarkDuplicates.

        Creates an Element object that can be run to invoke the low-level
        mark_duplicates function.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement that contains the input BAM file.
        index_off : bool
            If True, do not create BAM index after marking duplicates.
            Default is False.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BAM and extension ``bam``. The metrics file is written next
            to the BAM as ``<bam>.metrics.txt``.
        par : ParType | None
            Additional GATK parameters (e.g., Par(remove_duplicates=True)).
        cfg : CfgType | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that can be run and will produce the marked duplicates
            BAM and metrics file.

        Examples
        --------
        >>> gatk = GATK()
        >>> marked = gatk.mark(mapped)
        """

        resolved_tag = mapped.tag.bump(method=Method.GATK, state=State.DEDUP).resolve(
            tag
        )
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=mapped.bam.parent, ext="bam"), out
        )
        output_bam = resolved_out.file
        metrics_file = output_bam.with_suffix(".metrics.txt")

        runner = self.mark_duplicates(
            input_bam=mapped.bam,
            output_bam=output_bam,
            metrics_file=metrics_file,
            index_off=index_off,
            par=par,
            cfg=cfg,
        )
        artifacts = ArtifactSet(output_bam, primary_name="bam").with_extra(
            "metrics", metrics_file
        )
        if not index_off:
            artifacts = artifacts.with_extra(
                "bai", output_bam.with_suffix(output_bam.suffix + ".bai")
            )

        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag, self.version_name, "MarkDuplicates", input=mapped.bam.stem
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(mapped.bam,),
            pres=(mapped,),
        )

    ###########################################################################
    # Learn Read Orientation Model
    ###########################################################################

    @subroutine
    def learn_read_orientation(
        self,
        f1r2_tar_gz: Path | str | list[Path | str],
        output_model: Path | str,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "LearnReadOrientationModel"
        arguments = [subcommand]
        for f in f1r2_files:
            arguments.extend(["-I", self.strabs(f)])
        arguments.extend(["-O", self.strabs(output_model)])
        return arguments, subcommand, [*f1r2_files], [output_model], None, None, None

    @element
    def readorientation(
        self,
        mutect_element: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function for learning read orientation model.

        Creates an Element that learns orientation bias from F1R2 metrics
        collected during Mutect2 runs (useful for FFPE samples).

        Parameters
        ----------
        mutect_element : Element
            Element containing the F1R2 metrics file from a Mutect2 run.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            F1R2 file and extension ``tar.gz``.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        resolved_tag = mutect_element.tag.bump(
            method=Method.GATK, state=State.MODEL
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=mutect_element.f1r2.parent, ext="tar.gz"), out
        )
        output_model_file = resolved_out.file
        runner = self.learn_read_orientation(
            f1r2_tar_gz=mutect_element.f1r2,
            output_model=output_model_file,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "LearnReadOrientationModel",
            f1r2=f"{mutect_element.f1r2.stem}",
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(mutect_element.f1r2,),
            artifacts=ArtifactSet(output_model_file, primary_name="orientation"),
            pres=(mutect_element,),
        )

    ###########################################################################
    # Mutect2
    ###########################################################################

    @subroutine
    def mutect2_bam(
        self,
        reference: Path | str,
        input_bams_tumor: Mapping[str, Path | str],
        output_vcf: Path | str,
        *,
        input_bams_normal: Mapping[str, Path | str] | None = None,
        targets_padded_bed: Path | str | None = None,
        germline_resource: Path | str | None = None,
        panel_of_normals: Path | str | None = None,
        f1r2_tar_gz: Path | str | None = None,
        index_off: bool = False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run Mutect2 for somatic variant calling.

        Creates a zero-argument callable that runs GATK Mutect2 to call
        somatic mutations from tumor/normal BAM pairs.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        input_bams_tumor : Mapping[str, Path | str]
            Mapping of tumor sample name -> tumor BAM path.
        output_vcf : Path | str
            Output VCF file (will be gzipped).
        input_bams_normal : Mapping[str, Path | str] | None
            Mapping of normal sample name -> normal BAM path (optional).
        targets_padded_bed : Path | str | None
            BED file or interval list restricting calling regions.
        germline_resource : Path | str | None
            Germline resource VCF for filtering.
        panel_of_normals : Path | str | None
            Panel of normals VCF.
        f1r2_tar_gz : Path | str | None
            Output file for F1R2 metrics (for orientation model).
        index_off : bool
            If True, do not create an index for the output VCF file. By default,
            an index will be created.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "Mutect2"
        arguments = [subcommand, "-R", self.strabs(reference)]
        in_paths = [reference]
        names = []

        for name, bam in input_bams_tumor.items():
            arguments.extend(["-I", self.strabs(bam)])
            names.extend(["-tumor", name])
            in_paths.append(bam)
        if input_bams_normal:
            for name, bam in input_bams_normal.items():
                arguments.extend(["-I", self.strabs(bam)])
                names.extend(["-normal", name])
                in_paths.append(bam)

        if targets_padded_bed:
            arguments.extend(["-L", self.strabs(targets_padded_bed)])
            in_paths.append(targets_padded_bed)
        if germline_resource:
            arguments.extend(["--germline-resource", self.strabs(germline_resource)])
            in_paths.append(germline_resource)
        if panel_of_normals:
            arguments.extend(["--panel-of-normals", self.strabs(panel_of_normals)])
            in_paths.append(panel_of_normals)
        if f1r2_tar_gz:
            arguments.extend(["--f1r2-tar-gz", self.strabs(f1r2_tar_gz)])
            in_paths.append(f1r2_tar_gz)

        output = self.strabs(output_vcf)
        arguments.extend(["-O", output])
        # by default also calculate an index
        post = self.index_feature_file(output) if not index_off else None
        return arguments, subcommand, in_paths, [output_vcf], None, None, post

    @element
    def mutect2(
        self,
        marked_tumor: MappedElement | Mapping[str, MappedElement],
        reference: Genome,
        *,
        marked_normal: MappedElement | Mapping[str, MappedElement] | None = None,
        bam_indices: Element | Mapping[str, Element] | None = None,
        targets_padded: Element | None = None,
        germline_resource: Element | None = None,
        panel_of_normals: Element | None = None,
        f1r2_tar_gz: bool = True,
        index_off: bool = False,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function for somatic variant calling with Mutect2.

        Creates an Element that runs GATK Mutect2 on the BAM file(s) from the
        given MappedElement(s).

        Parameters
        ----------
        marked_tumor : MappedElement | Mapping[str, MappedElement]
            MappedElement(s) representing the tumor sample(s).
        reference : Genome
            Reference genome object containing the FASTA file.
        marked_normal : MappedElement | Mapping[str, MappedElement]
            MappedElement(s) representing the normal sample(s).
        bam_indices : Element | Mapping[str, Element] | None
            Optional BAM-index Element(s), keyed by sample name when multiple
            tumor/normal BAMs are supplied. These prerequisites ensure indexes
            are generated before Mutect2 runs.
        targets_padded : Element | None
            Element with a BED file or interval list restricting calling regions.
        germline_resource : Element | None
            Element with a germline resource VCF for filtering.
        panel_of_normals : Element | None
            Element with a panel of normals VCF.
        f1r2_tar_gz : bool
            If to output a file for F1R2 metrics (for orientation model).
        index_off : bool
            If True, do not create an index for the output VCF file. By default,
            an index will be created.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            tumor BAM file and extension ``vcf.gz``.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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

        def validate_input(marked_tumor):
            if marked_tumor is None:
                return marked_tumor
            if isinstance(marked_tumor, Element):
                marked_tumor = {marked_tumor.tag.root: marked_tumor}
            elif isinstance(marked_tumor, Mapping):
                for element in marked_tumor.values():
                    GATK.check_mapped(element)
            else:
                raise ValueError(
                    f"marked_tumor argument must be an Element with a bam file or a mapping thereof, was {type(marked_tumor)}"  # noqa: E501
                )
            return marked_tumor

        marked_tumor = validate_input(marked_tumor)
        marked_normal = validate_input(marked_normal)
        root = next(iter(marked_tumor))
        input_bams_tumor = {name: item.bam for name, item in marked_tumor.items()}
        pres = tuple(marked_tumor.values())
        inputs = tuple(input_bams_tumor.values())
        input_bams_normal = None
        suffix = ""
        if marked_normal:
            input_bams_normal = {name: item.bam for name, item in marked_normal.items()}
            pres += tuple(marked_normal.values())
            inputs += tuple(input_bams_normal.values())
            suffix = f"_against_{','.join([m.tag.default_name for m in marked_normal.values()])}"  # noqa: E501
        if bam_indices is None:
            indexes: dict[str, Element] = {}
        elif isinstance(bam_indices, Element):
            indexes = {bam_indices.tag.root: bam_indices}
        elif isinstance(bam_indices, Mapping):
            if not all(isinstance(index, Element) for index in bam_indices.values()):
                raise TypeError("Every BAM index must be an Element.")
            indexes = dict(bam_indices)
        else:
            raise TypeError(
                "bam_indices must be an Element, a mapping of Elements, or None."
            )
        required_index_names = set(input_bams_tumor)
        if input_bams_normal:
            required_index_names.update(input_bams_normal)
        if bam_indices is not None:
            missing_indexes = sorted(required_index_names - set(indexes))
            if missing_indexes:
                raise ValueError(
                    "BAM index Elements are missing for samples: " f"{missing_indexes}."
                )
            selected_indexes = tuple(
                indexes[name] for name in sorted(required_index_names)
            )
            mismatched_indexes = [
                (name, indexes[name].tag.root)
                for name in sorted(required_index_names)
                if indexes[name].tag.root != name
            ]
            if mismatched_indexes:
                raise ValueError(
                    "BAM index tags must match their sample names: "
                    f"{mismatched_indexes}."
                )
            pres += selected_indexes
            inputs += tuple(index.file for index in selected_indexes)
        resolved_tag = ElementTag(
            root=root,
            level=marked_tumor[root].tag.level + 1 if marked_tumor else 1,
            stage=Stage.CALL,
            method=Method.GATK,
            state=State.RAW,
            omics=marked_tumor[root].tag.omics if marked_tumor else None,
        ).resolve(tag)

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=marked_tumor[root].bam.parent, ext="vcf.gz"), out
        )
        output = resolved_out.file
        output_index = output.with_suffix(output.suffix + ".tbi")
        f1r2 = output.with_suffix(".f1r2.tar.gz") if f1r2_tar_gz else None
        runner = self.mutect2_bam(
            reference=reference.fasta,
            input_bams_tumor=input_bams_tumor,
            output_vcf=output,
            input_bams_normal=input_bams_normal,
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
            index_off=index_off,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par, "Mutect2")
        tumor_str = ",".join([m.tag.default_name for m in marked_tumor.values()])
        key = Element.generate_key(
            resolved_tag, self.version_name, "Mutect2", on=f"{tumor_str}{suffix}"
        )
        pres += tuple(
            x
            for x in [targets_padded, germline_resource, panel_of_normals]
            if isinstance(x, Element)
        )
        artifacts = ArtifactSet(output, primary_name="vcf")
        artifacts = artifacts.with_extra(
            "stats", output.with_suffix(output.suffix + ".stats")
        )
        if f1r2:
            artifacts = artifacts.with_extra("f1r2", f1r2)
        if not index_off:
            artifacts = artifacts.with_extra("tbi", output_index)
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            artifacts=artifacts,
            inputs=inputs,
            pres=pres,
        )

    ###########################################################################
    # Mutation Count
    ###########################################################################

    @element
    def mutationcount(
        self,
        mutation_table: Element,
        mutect_elements: Sequence[Element],
        *,
        block: str | int,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Build mutation counts and AF and FRAG_VAF matrices from VCFs."""

        from mmalignments.services.io import safe_name

        resolved_tag = mutation_table.tag.bump(
            stage=Stage.QUANT, method=Method.GATK, state=State.COUNT, omics=Omics.DNA
        ).resolve(tag)
        stem = f"{safe_name(block)}"
        matrix_stem = f"{stem}_vaf_matrix"
        counts_stem = FileSpec(f"{stem}_mutation_counts", "tsv")

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(
                folder=mutation_table.file.parent,
                stem=matrix_stem,
                ext="tsv",
                additional_output={
                    "counts": counts_stem,
                    "frag_vaf_matrix": FileSpec(f"{stem}_frag_vaf_matrix", "tsv"),
                },
            ),
            out,
        )
        determinants = ()
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            f"mutationcount_{safe_name(block)}",
        )
        pres = tuple([mutation_table] + list(mutect_elements))
        artifacts = ArtifactSet.from_outspec(resolved_out)

        runner = self.count_vaf_from_vcf(
            [el.file for el in mutect_elements],
            mutation_table.file,
            matrix_path=artifacts.primary.resolve(),
            counts_path=artifacts["counts"].resolve(),
            frag_matrix_path=artifacts["frag_vaf_matrix"].resolve(),
        )

        return Element(
            key,
            runner,
            tag=resolved_tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(mutation_table.file,),
            pres=pres,
        )

    def count_vaf_from_vcf(
        self,
        mutect_paths: Sequence[Path],
        mutation_table: Path,
        matrix_path: Path,
        counts_path: Path,
        *,
        frag_matrix_path: Path | None = None,
    ) -> Runnable:
        """Write candidate-by-sample VAF and long-form counts from Mutect2 VCFs.

        Supply one VCF per sample. Candidates require block, mutation_id,
        chrom, pos, ref, and alt columns. Matching uses chrom, pos, ref, and
        each individual ALT allele. Optional called_<sample> columns provide
        the call flags. Otherwise, matching PASS records are marked as called.

        Counts and depth come from FORMAT/AD. The two matrices contain
        FORMAT/AF and FORMAT/FRAG_VAF, respectively. The counts table includes
        both fields as af and frag_vaf, with vaf and vaf_percent derived from AF.
        Absent alleles receive zeros. Missing fields on matching records are
        written as empty values.

        frag_matrix_path defaults to a sibling of matrix_path with a
        frag_vaf_matrix suffix in the filename stem.
        """
        if frag_matrix_path is None:
            prefix = matrix_path.stem.removesuffix("vaf_matrix")
            if prefix == matrix_path.stem:
                prefix += "_"
            frag_matrix_path = matrix_path.with_name(
                f"{prefix}frag_vaf_matrix{matrix_path.suffix}"
            )
        output_paths = (matrix_path, counts_path, frag_matrix_path)
        if len({path.resolve() for path in output_paths}) != len(output_paths):
            raise ValueError(
                "AF matrix, fragment VAF matrix, and counts paths must differ."
            )

        def call():
            allele_columns = ["chrom", "pos", "ref", "alt"]
            variant_columns = ["block", "mutation_id", *allele_columns]
            depth_columns = [
                "called_by_mutect",
                "ref_count",
                "alt_count",
                "other_count",
            ]
            fraction_columns = ["af", "frag_vaf"]
            count_columns = [
                "sample",
                *depth_columns,
                "depth",
                "vaf",
                "vaf_percent",
                *fraction_columns,
            ]
            read_options = (
                {} if mutation_table.suffix == ".parquet" else {"dtype": {"chrom": str}}
            )
            candidates = read_frame(mutation_table, **read_options).reset_index(
                drop=True
            )
            missing = set(variant_columns) - set(candidates.columns)
            if missing:
                raise ValueError(
                    f"Candidate table {mutation_table} is missing columns: "
                    f"{sorted(missing)}"
                )
            candidate_keys = list(
                candidates[allele_columns]
                .astype({"chrom": str, "pos": int})
                .itertuples(index=False, name=None)
            )
            wanted = set(candidate_keys)
            matrix = candidates[variant_columns].copy()
            frag_matrix = matrix.copy()
            counts = []

            for vcf_path in mutect_paths:
                vcf = read_vcf(vcf_path)
                if len(vcf.columns) != 10:
                    raise ValueError(f"Expected one Mutect2 sample in {vcf_path}")
                sample = vcf.columns[9]
                if sample in matrix.columns:
                    raise ValueError(
                        f"Duplicate or reserved sample name {sample!r} in {vcf_path}"
                    )
                observations = _candidate_allele_metrics(vcf, sample, wanted, vcf_path)
                sample_counts = candidates[variant_columns].copy()
                sample_counts["sample"] = sample
                called_column = f"called_{sample}"
                rows = [
                    observations.get(key, (False, 0, 0, 0, 0.0, 0.0))
                    for key in candidate_keys
                ]
                metrics = pd.DataFrame(rows, columns=depth_columns + fraction_columns)
                if called_column in candidates:
                    metrics["called_by_mutect"] = candidates[called_column].fillna(
                        False
                    )
                metrics["depth"] = (
                    metrics["ref_count"] + metrics["alt_count"] + metrics["other_count"]
                )
                metrics["vaf"] = metrics["af"]
                metrics["vaf_percent"] = metrics["af"] * 100
                matrix[sample] = metrics["af"]
                frag_matrix[sample] = metrics["frag_vaf"]
                frag_matrix[sample] = metrics["frag_vaf percent"] = (
                    metrics["frag_vaf"] * 100
                )
                counts.append(
                    sample_counts.join(metrics)[variant_columns + count_columns]
                )
            count_table = (
                pd.concat(counts, ignore_index=True)
                if counts
                else pd.DataFrame(columns=variant_columns + count_columns)
            )
            write_frame(matrix, matrix_path, index=False)
            write_frame(frag_matrix, frag_matrix_path, index=False)
            write_frame(count_table, counts_path, index=False)

        return Runnable(
            call,
            display=CallSpec(
                path=(
                    "gatk",
                    "count_vaf_from_vcf",
                ),
            ).render(),
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
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "GetPileupSummaries"
        arguments = [
            subcommand,
            "-I",
            self.strabs(input_bam),
            "-V",
            self.strabs(variant_sites),
            "-O",
            self.strabs(output_table),
        ]
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])

        in_paths = [input_bam, variant_sites]
        if intervals:
            in_paths.append(intervals)
        return arguments, subcommand, in_paths, [output_table], None, None, None

    @element
    def pilesum(
        self,
        mapped: MappedElement,
        known_variants: Element | Path | str,
        targets_padded: Element | Path | str,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BAM and extension ``table``.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        input_bam = Path(mapped.bam)
        resolved_tag = mapped.tag.bump(
            stage=Stage.CALL, method=Method.GATK, state=State.PILE
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bam.parent, ext="table"), out
        )
        output_table = resolved_out.file
        variant_sites = (
            known_variants.vcf
            if isinstance(known_variants, Element)
            else Path(known_variants)
        )
        intervals = (
            targets_padded.bed
            if isinstance(targets_padded, Element)
            else Path(targets_padded)
        )

        runner = self.get_pileup_summaries(
            input_bam=input_bam,
            variant_sites=variant_sites,
            output_table=output_table,
            intervals=intervals,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        pres = (mapped,)
        if isinstance(known_variants, Element):
            pres += (known_variants,)
        if isinstance(targets_padded, Element):
            pres += (targets_padded,)
        key = Element.generate_key(
            resolved_tag, self.version_name, "GetPileupSummaries", input=mapped.bam.stem
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_bam, variant_sites, intervals),
            artifacts=ArtifactSet(output_table.absolute(), primary_name="pileup"),
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "CalculateContamination"
        arguments = [
            subcommand,
            "-I",
            self.strabs(tumor_pileup),
            "-O",
            self.strabs(output_table),
        ]
        if normal_pileup:
            arguments.extend(["-matched", self.strabs(normal_pileup)])
        if output_segments:
            arguments.extend(["--tumor-segmentation", self.strabs(output_segments)])

        in_paths = [tumor_pileup]
        if normal_pileup:
            in_paths.append(normal_pileup)
        out_paths = [output_table]
        if output_segments:
            out_paths.append(output_segments)
        return arguments, subcommand, in_paths, out_paths, None, None, None

    @element
    def contamination(
        self,
        tumor_pileup: Element,
        *,
        normal_pileup: Element | None = None,
        output_segments: Path | str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            tumor pileup table and extension ``table``.
        output_segments : Path | str | None
            Output tumor segmentation table (optional).
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        resolved_tag = tumor_pileup.tag.bump(
            root=f"{tumor_pileup.name}.contamination",
            stage=Stage.CALL,
            method=Method.GATK,
            state=State.MODEL,
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=tumor_table.parent, ext="table"), out
        )
        output_table = resolved_out.file.absolute()

        runner = self.calculate_contamination(
            tumor_pileup=tumor_table,
            output_table=output_table,
            normal_pileup=normal_table,
            output_segments=output_segments,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        pres = [tumor_pileup]
        inputs = [tumor_table]
        if normal_pileup:
            pres.append(normal_pileup)
            inputs.append(normal_table)
        params_str = tumor_pileup.key
        if normal_pileup:
            params_str += f"_vs_{normal_pileup.key}"
        key = Element.generate_key(
            resolved_tag, self.version_name, "CalculateContamination", params=params_str
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(inputs),
            artifacts=ArtifactSet(output_table, primary_name="contamination"),
            pres=tuple(pres),
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
        contamination: Path | str | None = None,
        orientation: Path | str | None = None,
        tumor_segmentation: Path | str | None = None,
        stats: Path | str | None = None,
        index_off: bool = False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        contamination : Path | str | None
            Contamination table from CalculateContamination.
        orientation : Path | str | None
            Read orientation model from LearnReadOrientationModel.
        tumor_segmentation : Path | str | None
            Tumor segmentation table from CalculateContamination.
        stats : Path | str | None
            Stats file from Mutect2.
        index_off : bool
            If True, do not create an index for the output VCF file. By default,
            an index will be created.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "FilterMutectCalls"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference),
            "-V",
            self.strabs(input_vcf),
            "-O",
            self.strabs(output_vcf),
        ]
        if contamination:
            arguments.extend(["--contamination-table", self.strabs(contamination)])
        if orientation:
            arguments.extend(["--ob-priors", self.strabs(orientation)])
        if tumor_segmentation:
            arguments.extend(["--tumor-segmentation", self.strabs(tumor_segmentation)])
        if stats:
            arguments.extend(["--stats", self.strabs(stats)])
        in_paths = [reference, input_vcf]
        if contamination:
            in_paths.append(contamination)
        if orientation:
            in_paths.append(orientation)
        if tumor_segmentation:
            in_paths.append(tumor_segmentation)
        post = self.index_feature_file(output_vcf, par, cfg) if not index_off else None
        return arguments, subcommand, in_paths, [output_vcf], None, None, post

    @element
    def filter(
        self,
        mutect_element: Element,
        reference: Genome,
        *,
        contamination: Element | None = None,
        orientation: Element | None = None,
        tumor_segmentation: Element | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        index_off: bool = False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function for filtering Mutect2 calls.

        Creates an Element that applies GATK filtering to Mutect2 VCF output.

        Parameters
        ----------
        mutect_element : Element
            Element containing the unfiltered Mutect2 VCF.
        contamination : Element | None
            Element containing the contamination table from CalculateContamination.
        orientation : Element | None
            Element containing the read orientation model from LearnReadOrientationModel
        tumor_segmentation : Element | None
            Element containing the tumor segmentation table from CalculateContamination.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input VCF and extension ``vcf.gz``.
        index_off : bool
            If True, do not create an index for the output VCF file. By default,
            an index will be created.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        stats = mutect_element.artifacts["stats"].resolve()
        resolved_tag = mutect_element.tag.bump(state=State.FILTER).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_vcf.parent, ext="vcf.gz"), out
        )
        output_vcf = resolved_out.file.absolute()

        runner = self.filter_mutect_calls(
            reference=reference.fasta,
            input_vcf=input_vcf,
            output_vcf=output_vcf,
            contamination=(contamination.contamination if contamination else None),
            orientation=orientation.orientation if orientation else None,
            tumor_segmentation=(
                tumor_segmentation.segments if tumor_segmentation else None
            ),
            stats=stats,
            index_off=index_off,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        pres = [mutect_element]
        inputs = [input_vcf]
        if contamination:
            pres.append(contamination)
            inputs.append(contamination.contamination)
        if orientation:
            pres.append(orientation)
            inputs.append(orientation.orientation)
        if tumor_segmentation:
            pres.append(tumor_segmentation)
            inputs.append(tumor_segmentation.segments)

        con = contamination.key if contamination else None
        ori = orientation.key if orientation else None
        seg = tumor_segmentation.key if tumor_segmentation else None
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "FilterMutectCalls",
            on=mutect_element.key,
            contamination=con,
            orientation=ori,
            segmentation=seg,
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(inputs),
            artifacts=ArtifactSet(output_vcf, primary_name="vcf"),
            pres=tuple(pres),
        )

    ###########################################################################
    # BQSR
    ###########################################################################

    @subroutine
    def basecalibrate_bam(
        self,
        reference: Path | str,
        sequence_dictionary: Path | str,
        input_bams: list[Path | str],
        output_table: Path | str,
        known_sites: Path | str | list[Path | str] | None = None,
        intervals: Path | str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run GATK BaseRecalibrator on BAM file(s).

        Creates a zero-argument callable that computes a BQSR recalibration
        table from one or more input BAMs.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        sequence_dictionary: Path | str
            Reference sequence dictionary file (e.g., .dict from Picard).
        input_bams : list[Path | str]
            List of input BAM files.
        output_table : Path | str
            Output recalibration table file.
        known_sites : Path | str | list[Path | str] | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "BaseRecalibrator"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference),
            "--sequence-dictionary",
            self.strabs(sequence_dictionary),
        ]
        for bam in input_bams:
            arguments.extend(["-I", self.strabs(bam)])

        in_paths = [reference, *input_bams]
        if sequence_dictionary:
            in_paths.append(sequence_dictionary)

        if known_sites:
            if isinstance(known_sites, (list, tuple)):
                for ks in known_sites:
                    arguments.extend(["--known-sites", self.strabs(ks)])
                    in_paths.append(ks)
            else:
                arguments.extend(["--known-sites", self.strabs(known_sites)])
                in_paths.append(known_sites)
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])
            in_paths.append(intervals)

        arguments.extend(["-O", self.strabs(output_table)])
        return arguments, subcommand, in_paths, [output_table], None, None, None

    @element
    def modelbsqr(
        self,
        mapped: MappedElement | list[MappedElement],
        reference: Genome,
        refdict_element: Element | None = None,
        *,
        known_sites: Path | str | list[Path | str] | Element | None = None,
        intervals: Path | str | Element | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
        pre: list[Element] | None = None,
    ) -> Element:
        """High-level function for GATK BaseRecalibrator.

        Creates an Element that computes a BQSR recalibration model and outputs
        a corresponding table from one or more input BAMs.

        Parameters
        ----------
        mapped : MappedElement | list[MappedElement]
            MappedElement(s) containing the BAM file(s) to recalibrate.
        reference : Genome
            Reference genome object containing the FASTA file.
        refdict_element : Element | None
            Optional Element containing the reference sequence dictionary. If not
            provided, the sequence dictionary will be inferred from the reference
            FASTA file.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            first input BAM and extension ``table``.
        known_sites : Path | str | list[Path | str] | Element | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | Element | None
            BED file or interval list restricting regions.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that executes BaseRecalibrator when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> recal_element = gatk.modelbsqr(mapped, genome)
        """
        mapped_list = list(mapped) if isinstance(mapped, (list, tuple)) else [mapped]
        input_bams = [m.bam for m in mapped_list]
        first_bam = input_bams[0]
        resolved_tag = (
            mapped_list[0]
            .tag.bump(
                stage=Stage.CALL,
                method=Method.GATK,
                state=State.MODEL,
                flag="bsqr",
            )
            .resolve(tag)
        )
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=first_bam.parent, ext="table"), out
        )
        output_table = resolved_out.file

        if isinstance(intervals, Element) and intervals.files:
            intervals = next(iter(intervals.files))

        known_sites = (
            list(known_sites.artifacts.values())
            if isinstance(known_sites, Element)
            else known_sites
        )
        runner = self.basecalibrate_bam(
            reference=reference.fasta if hasattr(reference, "fasta") else reference,
            sequence_dictionary=(
                refdict_element.dict
                if isinstance(refdict_element, Element)
                else refdict_element
            ),
            input_bams=tuple(input_bams),
            output_table=output_table,
            known_sites=known_sites,
            intervals=intervals,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag, self.version_name, "BaseRecalibrator", first=first_bam.stem
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(input_bams),
            artifacts=ArtifactSet(output_table.absolute(), primary_name="bsqr_table"),
            pres=tuple(mapped_list),
        )

    @subroutine
    def apply_bqsr_on_bam(
        self,
        reference: Path | str,
        input_bam: Path | str,
        recal_table: Path | str,
        output_bam: Path | str,
        *,
        index_off: bool = False,
        intervals: Path | str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        index_off : bool
            If True, do not create BAM index after recalibration. Default is
            False.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "ApplyBQSR"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference),
            "-I",
            self.strabs(input_bam),
            "--bqsr-recal-file",
            self.strabs(recal_table),
            "-O",
            self.strabs(output_bam),
        ]
        if intervals:
            arguments.extend(["-L", self.strabs(intervals)])

        in_paths = [reference, input_bam, recal_table]
        if intervals:
            in_paths.append(intervals)

        post = None if index_off else GATK.index_by_sam(output_bam)
        return arguments, subcommand, in_paths, [output_bam], None, None, post

    @element
    def applybsqr(
        self,
        mapped: MappedElement,
        bsqrmodel: Element,
        reference: Genome,
        *,
        intervals: Path | str | Element | None = None,
        index_off: bool = False,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
        pre: list[Element] | None = None,
    ) -> MappedElement:
        """High-level function for Base Quality Score Recalibration (BQSR).

        Creates an Element that applies a pre-computed BQSR recalibration
        table to the input BAM.

        Parameters
        ----------
        mapped : MappedElement
            Input MappedElement containing the BAM file to recalibrate.
        bsqrmodel : Element
            Element with the BSQR model table.
        reference : Genome
            Reference genome object containing the FASTA file.
        intervals : Path | str | None
            BED file or interval list restricting regions.
        index_off : bool
            If True, do not create BAM index after recalibration. Default is
            False.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BAM and extension ``bam``.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        pres = [mapped, bsqrmodel] + (pre or [])
        resolved_tag = mapped.tag.bump(
            method=Method.GATK, state=State.RECAL, flag="bsqr"
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=mapped.bam.parent, ext="bam"), out
        )
        output = resolved_out.file
        if isinstance(intervals, Element):
            intervals = intervals.bed

        bsqr_runner = self.apply_bqsr_on_bam(
            reference=reference.fasta,
            input_bam=mapped.bam,
            recal_table=bsqrmodel.bsqr_table,
            output_bam=output,
            intervals=intervals,
            index_off=index_off,
            par=par,
            cfg=cfg,
        )
        artifacts = ArtifactSet(output.absolute(), primary_name="bam")
        if not index_off:
            artifacts = artifacts.with_extra(
                "bai", output.with_suffix(output.suffix + ".bai")
            )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "ApplyBQSR",
            input=f"{mapped.bam.stem}",
            bsqr=f"{bsqrmodel.key}",
        )
        return MappedElement(
            key,
            bsqr_runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(mapped.bam, bsqrmodel.bsqr_table),
            artifacts=artifacts,
            pres=tuple(pres),
        )

    def bsqr(
        self,
        mapped: MappedElement,
        reference: Genome,
        *,
        refdict_element: Element | None = None,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        known_sites: Path | str | list[Path | str] | Element | None = None,
        intervals: Path | str | Element | None = None,
        index_off: bool = False,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
        pre: list[Element] | None = None,
    ) -> MappedElement:
        """High-level function for Base Quality Score Recalibration (BQSR).

        Creates an Element that first computes the BQSR recalibration table
        and then applies it to the input BAM.

        Parameters
        ----------
        mapped : MappedElement
            Input MappedElement containing the BAM file to recalibrate.
        reference : Genome
            Reference genome object containing the FASTA file.
        refdict_element : Element | None
            Optional Element containing the reference sequence dictionary. If
            not provided, the sequence dictionary will be inferred from the
            reference FASTA file.
        tag : TagType | SpecValues[TagType] | None
            Tag override, either for both steps or per step using
            ``SpecValues``. For per-step overrides, declare
            ``valid_steps=("modelbsqr", "applybsqr")`` and put values in the
            ``values`` mapping.
        out : OutType | SpecValues[OutType] | None
            Output specification override, either for both steps or per step
            using ``SpecValues``. Use per-step values if a non-default ``out``
            is needed, since the model table and the BAM have different
            extensions.
        known_sites : Path | str | list[Path | str] | Element | None
            Known variant sites for BQSR. Can be a single VCF or a list.
        intervals : Path | str | Element | None
            BED file or interval list restricting regions.
        index_off : bool
            If True, do not create BAM index after recalibration. Default is
            False.
        par : ParType | SpecValues[ParType] | None
            Additional GATK parameters (e.g., Par()), either for both steps or
            per step using ``SpecValues``.
        cfg : CfgType | SpecValues[CfgType] | None
            Configuration for running the command, either for both steps or per
            step using ``SpecValues``.

        Returns
        -------
        Element
            An Element that executes BQSR when run.

        Examples
        --------
        >>> gatk = GATK()
        >>> bsqr_element = gatk.bsqr(mapped, genome)
        """
        model = self.modelbsqr(  # type: ignore[reportCallIssue]
            mapped=mapped,
            reference=reference,
            refdict_element=refdict_element,
            tag=pick_spec(tag, "modelbsqr"),
            out=pick_spec(out, "modelbsqr"),
            known_sites=known_sites,
            intervals=intervals,
            par=pick_spec(par, "modelbsqr"),
            cfg=pick_spec(cfg, "modelbsqr"),
            pre=pre,
        )
        recalibrated = self.applybsqr(  # type: ignore[reportCallIssue]
            mapped=mapped,
            bsqrmodel=model,
            reference=reference,
            intervals=intervals,
            index_off=index_off,
            tag=pick_spec(tag, "applybsqr"),
            out=pick_spec(out, "applybsqr"),
            par=pick_spec(par, "applybsqr"),
            cfg=pick_spec(cfg, "applybsqr"),
            pre=pre,
        )
        return recalibrated

    ###########################################################################
    # Panel of Normals
    ###########################################################################

    @subroutine
    def create_panel_of_normals(
        self,
        input_vcfs: list[Path | str],
        output_pon: Path | str,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Create a Panel of Normals using GATK CreateSomaticPanelOfNormals.

        Creates a zero-argument callable that builds a panel of normals VCF
        from multiple normal sample VCFs.

        Parameters
        ----------
        input_vcfs : list[Path | str]
            List of normal sample VCF files from Mutect2 tumor-only calls.
        output_pon : Path | str
            Output VCF file for the panel of normals.
        par : ParType | None
            Additional GATK parameters (e.g., Par()).
        cfg : CfgType | None
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
        subcommand = "CreateSomaticPanelOfNormals"
        arguments = [subcommand]
        input_vcfs = sorted(input_vcfs, key=str)
        for vcf in input_vcfs:
            arguments.extend(["-V", self.strabs(vcf)])
        arguments.extend(["-O", self.strabs(output_pon)])
        return arguments, subcommand, [*input_vcfs], [output_pon], None, None, None

    @element
    def pons(
        self,
        all_normals: list[MappedElement],
        reference: Genome,
        *,
        index_off: bool = False,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
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
        index_off : bool
            If True, do not create BAM index after recalibration. Default is
            False.
        tag : TagType | None
            Partial or full Element tag override for the panel of normals.
        out : OutType | None
            Output specification override for the panel of normals. Defaults to
            ``cache/PoNs`` and extension ``vcf.gz``. The per-normal Mutect2
            calls are written to the same folder.
        par : ParType | SpecValues[ParType] | None
            Additional GATK parameters. For per-step overrides, declare
            ``valid_steps=("mutect2", "pon")`` and put values in the
            ``values`` mapping.
        cfg : CfgType | SpecValues[CfgType] | None
            Configuration for running the commands, optionally overridden
            per step using the same ``valid_steps`` and ``values`` pattern.

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
        all_normals = sorted(
            all_normals, key=lambda e: e.name
        )  # we need sorting for the name
        resolved_tag = ElementTag(
            root=f"{reference.name}",
            level=all_normals[0].tag.level + 1,
            stage=Stage.PREP,
            method=Method.GATK,
            state=State.MODEL,
            omics=all_normals[0].tag.omics,
            flag="PoN",
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=Path("cache/PoNs"), ext="vcf.gz"), out
        )
        output = resolved_out.file
        par_mutect2 = pick_spec(par, "mutect2")
        par_pon = pick_spec(par, "pon")
        pres_pon = []
        for normal in all_normals:
            GATK.check_mapped(normal)
            called = self.mutect2(
                marked_tumor=normal,
                reference=reference,
                index_off=index_off,
                out=OutSpec(
                    folder=output.parent,
                    stem=f"{normal.name}.mutect2",
                    ext="vcf.gz",
                ),
                par=par_mutect2,
                cfg=pick_spec(cfg, "mutect2"),
            )
            pres_pon.append(called)
        pres_pon = sorted(pres_pon, key=lambda e: e.name)
        input_vcfs = [e.vcf for e in pres_pon]
        runner = self.create_panel_of_normals(
            input_vcfs=input_vcfs,
            output_pon=output,
            par=par_pon,
            cfg=pick_spec(cfg, "pon"),
        )
        determinants = self.signature_determinants(par_pon)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "CreateSomaticPanelOfNormals",
            ref=f"{reference.name}",
            fromf=f"{','.join(sorted([n.name for n in all_normals]))}",
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(input_vcfs),
            artifacts=ArtifactSet(output.absolute(), primary_name="pon"),
            pres=tuple(pres_pon),
        )

    ###########################################################################
    # Sequence Dictionary
    ###########################################################################

    @subroutine
    def create_sequence_dictionary(
        self,
        reference: Path | str,
        output_dict: Path | str,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Create a reference sequence dictionary using GATK CreateSequenceDictionary.

        Creates a zero-argument callable that generates a .dict file from a
        reference FASTA file.

        Parameters
        ----------
        reference : Path | str
            Reference genome FASTA file.
        output_dict : Path | str
            Output .dict file path.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        subcommand = "CreateSequenceDictionary"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference),
            "-O",
            self.strabs(output_dict),
        ]
        return arguments, subcommand, [reference], [output_dict], None, None, None

    @element
    def seqdict(
        self,
        reference: Genome,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function to create a reference sequence dictionary.

        Creates an Element that generates a .dict file from a reference FASTA
        using GATK CreateSequenceDictionary.

        Parameters
        ----------
        reference : Genome
            Reference genome object containing the FASTA file.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the FASTA folder, ext
            ``dict``.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        resolved_tag = ElementTag(
            root=reference.name,
            level=1,
            stage=Stage.PREP,
            method=Method.GATK,
            state=State.RAW,
            omics=Omics.DNA,
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=reference.fasta.parent, ext="dict"), out
        )
        output_dict = resolved_out.file
        runner = self.create_sequence_dictionary(
            reference=reference.fasta,
            output_dict=output_dict,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "CreateSequenceDictionary",
            ref=reference.name,
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(reference.fasta,),
            artifacts=ArtifactSet(output_dict.absolute(), primary_name="dict"),
            pres=(),
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes BedToIntervalList.

        Examples
        --------
        >>> gatk = GATK()
        >>> converter = gatk.bed_to_interval_list(
               ...     sequence_dict="genome.dict",
        ... )
        >>> converter()
        """
        # Build command: gatk BedToIntervalList -I bed -O interval_list -SD dict
        subcommand = "BedToIntervalList"
        arguments = [
            subcommand,
            "-I",
            self.strabs(input_bed),
            "-O",
            self.strabs(output_interval_list),
            "-SD",
            self.strabs(sequence_dict),
        ]
        return (
            arguments,
            subcommand,
            [input_bed, sequence_dict],
            [output_interval_list],
            None,
            None,
            None,
        )

    @element
    def bed2interval(
        self,
        bed_element: Element | FileSource,
        sequence_dict: Element | FileSource,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function to convert a BED file to Picard interval_list format.

        Creates an Element that converts a BED file to the interval_list
        format required by Picard tools like CollectHsMetrics.

        Parameters
        ----------
        bed_element : Element | FileSource
            BED file or Element containing the BED file to convert.
        sequence_dict : Element | FileSource
            Element containing the reference sequence dictionary (.dict file).
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the BED folder, ext
            ``intervals``.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        input_bed = bed_element.artifacts.primary.resolve()
        root = bed_element.tag.root
        resolved_tag = ElementTag(
            root=root,
            level=bed_element.tag.level + 1,
            stage=Stage.PREP,
            method=Method.GATK,
            state=State.RAW,
            omics=bed_element.tag.omics,
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bed.parent, stem=input_bed.stem, ext="intervals"),
            out,
        )
        output_interval_list = resolved_out.file

        runner = self.bed_to_interval_list(
            input_bed=input_bed,
            output_interval_list=output_interval_list,
            sequence_dict=sequence_dict.dict,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        pres = [bed_element]
        pres.append(sequence_dict)
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "BedToIntervalList",
            sdict=f"{sequence_dict.key}",
            input=f"{input_bed.stem}",
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_bed,),
            artifacts=ArtifactSet(output_interval_list, primary_name="intervals"),
            pres=tuple(pres),
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        subcommand = "CollectAlignmentSummaryMetrics"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference),
            "-I",
            self.strabs(input_bam),
            "-O",
            self.strabs(output_metrics),
        ]
        return (
            arguments,
            subcommand,
            [reference, input_bam],
            [output_metrics],
            None,
            None,
            None,
        )

    @element
    def alignmetrics(
        self,
        mapped: MappedElement,
        reference: Genome,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``<bam folder>/qc/picard``,
            ext ``txt``.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        resolved_tag = mapped.tag.bump(
            method=Method.PICARD, state=State.STAT, flag="alignsummary"
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bam.parent / "qc" / "picard", ext="txt"), out
        )
        output_metrics = resolved_out.file

        runner = self.alignment_summary_metrics(
            reference=reference.fasta,
            input_bam=input_bam,
            output_metrics=output_metrics,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag, self.version_name, "CollectAlignmentSummaryMetrics"
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_bam, reference.fasta),
            artifacts=ArtifactSet(output_metrics, primary_name="metrics"),
            pres=(mapped,),
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        subcommand = "CollectInsertSizeMetrics"
        arguments = [
            subcommand,
            "-I",
            self.strabs(input_bam),
            "-O",
            self.strabs(output_metrics),
            "-H",
            self.strabs(output_histogram),
        ]
        return (
            arguments,
            subcommand,
            [input_bam],
            [output_metrics, output_histogram],
            None,
            None,
            None,
        )

    @element
    def insertmetrics(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
        output_histogram_name: Path | str | None = None,
    ) -> Element:
        """High-level function for insert size metrics.

        Creates an Element that collects insert size distribution metrics
        for paired-end sequencing data.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``<bam folder>/qc/picard``,
            ext ``txt``.
        output_histogram_name : str | None
            Output histogram PDF name. If not provided, the metrics file with
            suffix ``.hist.pdf`` is used.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        resolved_tag = mapped.tag.bump(
            stage=Stage.QC, method=Method.PICARD, state=State.STAT, flag="insize"
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bam.parent / "qc" / "picard", ext="txt"), out
        )
        output_metrics = resolved_out.file
        output_histogram = (
            output_metrics.parent / output_histogram_name
            if output_histogram_name
            else output_metrics.with_suffix(".hist.pdf")
        )

        runner = self.insert_size_metrics(
            input_bam=input_bam,
            output_metrics=output_metrics,
            output_histogram=output_histogram,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag, self.version_name, "CollectInsertSizeMetrics"
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_bam,),
            artifacts=ArtifactSet(
                output_metrics, primary_name="metrics", histogram=output_histogram
            ),
            pres=(mapped,),
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        subcommand = "CollectHsMetrics"
        arguments = [
            subcommand,
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
        in_paths = [
            reference,
            input_bam,
            bait_intervals,
            target_intervals,
        ]

        if per_target_coverage:
            arguments.extend(
                ["--PER_TARGET_COVERAGE", self.strabs(per_target_coverage)]
            )
            in_paths.append(per_target_coverage)

        return arguments, subcommand, in_paths, [output_metrics], None, None, None

    @element
    def hs(
        self,
        mapped: MappedElement,
        reference: Genome,
        baits: Element | Path | str,
        targets: Element | Path | str,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        per_target_coverage: Path | str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BAM and extension ``txt``.
        per_target_coverage : Path | str | None
            Optional output file for per-target coverage.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        resolved_tag = mapped.tag.bump(
            stage=Stage.QC, method=Method.PICARD, state=State.STAT, flag="hs"
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bam.parent / "qc" / "picard", ext="txt"), out
        )
        output_metrics = resolved_out.file.absolute()

        bait_path = Path(
            baits.intervals if isinstance(baits, Element) else baits
        ).absolute()
        target_path = Path(
            targets.intervals if isinstance(targets, Element) else targets
        ).absolute()

        runner = self.hs_metrics(
            reference=reference.fasta,
            input_bam=input_bam,
            output_metrics=output_metrics,
            bait_intervals=bait_path,
            target_intervals=target_path,
            per_target_coverage=per_target_coverage,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)
        pres: list[Element] = [mapped]
        if isinstance(baits, Element):
            pres.append(baits)
        if isinstance(targets, Element) and targets != baits:
            pres.append(targets)

        artifacts = ArtifactSet(output_metrics, primary_name="metrics")
        if per_target_coverage:
            artifacts = artifacts.with_extra(
                "per_target_coverage", Path(per_target_coverage).absolute()
            )

        key = Element.generate_key(
            resolved_tag, self.version_name, "CollectHsMetrics", input=input_bam.stem
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_bam, reference.fasta, bait_path, target_path),
            artifacts=artifacts,
            pres=tuple(pres),
        )

    ###########################################################################
    # Index Feature File
    ###########################################################################

    @element
    def indexfeature(
        self,
        element: Element | str | Path,
        filetype: str = "vcf",
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """High-level function to index a feature file (e.g., VCF).

        Creates an Element that indexes a feature file using GATK IndexFeatureFile.

        Parameters
        ----------
        element : Element | str | Path
            Element or file to index.
        filetype : str
            Type of the feature file (e.g., "vcf", "bed", "intervals").
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the input folder and stem
            with ext ``idx``.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
            Configuration for running the command.

        Returns
        -------
        Element
            An Element that indexes the feature file when run.
        """
        if isinstance(element, Element):
            input_path = getattr(element, filetype, None)
            if not isinstance(input_path, Path):
                raise ValueError(
                    f"Element {element} does not have a '{filetype}' artifact of type Path: {input_path}"  # noqa: E501
                )
            resolved_tag = element.tag.bump(
                stage=Stage.PREP, method=Method.GATK, state=State.INDEX
            ).resolve(tag)
            pres = (element,)
        else:
            input_path = Path(element)
            pres = ()
            resolved_tag = ElementTag(
                root=input_path.stem,
                level=1,
                stage=Stage.PREP,
                method=Method.GATK,
                state=State.INDEX,
                omics=Omics.DNA,
            ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_path.parent, stem=input_path.stem, ext="idx"), out
        )
        output_index = resolved_out.file
        runner = self.index_feature_file(input_path, par, cfg)
        determinants = self.signature_determinants(par)
        key = Element.generate_key(
            resolved_tag, self.version_name, "IndexFeatureFile", input=input_path.stem
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_path,),
            artifacts=ArtifactSet(output_index, primary_name="idx"),
            pres=pres,
        )

    @subroutine
    def index_feature_file(
        self,
        input_path: str | Path,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """
        Create a callable that indexes a feature file using GATK IndexFeatureFile.

        Parameters
        ----------
        input_path : str | Path
            The file to be indexed.
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
            Configuration for running the command.

        Returns
        -------
        Callable[[], CompletedProcess]
            A zero-argument callable that executes the indexing command.
        """

        subcommand = "IndexFeatureFile"
        arguments = [
            subcommand,
            "-I",
            self.strabs(input_path),
        ]
        return arguments, subcommand, [input_path], [], None, None, None

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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional GATK parameters.
        cfg : CfgType | None
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
        subcommand = "CallableLoci"
        arguments = [
            subcommand,
            "-R",
            self.strabs(reference_fasta),
            "-I",
            self.strabs(input_bam),
            "-L",
            self.strabs(target_bed),
            # "--min-depth",
            # str(min_depth),
            # "--min-mapping-quality",
            # str(min_mapping_quality),
            # "--min-base-quality",
            # str(min_base_quality),
            "-summary",
            self.strabs(output_summary),
            "-O",
            self.strabs(output_bed),
        ]
        output_summary = Path(output_summary)
        in_paths = [reference_fasta, input_bam, target_bed]
        out_paths = [output_bed, output_summary]
        if output_json:
            out_paths.append(output_json)

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

        return arguments, subcommand, in_paths, out_paths, None, None, _post

    @element
    def loci(
        self,
        mapped: MappedElement,
        reference: Genome,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        output_name_bed: str | None = None,
        output_name_json: str | None = None,
        targets: Element | Path | str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override for the summary file. Defaults to the
            BAM folder, ext ``txt``.
        output_name_bed : str | None
            Output BED file name. Defaults to the summary name with ``.bed``.
        output_name_json : str | None
            Output JSON file name. Defaults to the summary name with ``.json``.
        targets : Element | Path | str | None
            Optional Element or file with target regions to evaluate.
        par : ParType | None
            Additional GATK parameters. Defaults to min_depth=5,
            min_mapping_quality=20 and min_base_quality=20.
        cfg : CfgType | None
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
        par = par or ParType(
            min_depth=5,
            min_mapping_quality=20,
            min_base_quality=20,
        )
        input_bam = mapped.bam
        resolved_tag = mapped.tag.bump(
            stage=Stage.CALL, method=Method.GATK, state=State.LOCI
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bam.parent, ext="txt"), out
        )
        output_summary = resolved_out.file
        output_bed = output_summary.parent / (
            output_name_bed or output_summary.with_suffix(".bed").name
        )
        output_json = output_summary.parent / (
            output_name_json or output_summary.with_suffix(".json").name
        )
        target_bed = targets.bed if isinstance(targets, Element) else targets

        inputs = [input_bam, reference.fasta]
        pres: list[Element] = [mapped]
        if isinstance(targets, Element):
            pres.append(targets)
            inputs.append(targets.bed)

        runner = self.callable_loci(
            reference_fasta=reference.fasta,
            input_bam=input_bam,
            output_bed=output_bed,
            output_summary=output_summary,
            target_bed=target_bed,
            # min_depth=min_depth,
            # min_mapping_quality=min_mapping_quality,
            # min_base_quality=min_base_quality,
            output_json=output_json,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par, subroutine="CallableLoci")
        targ = {Path(target_bed).stem} if target_bed else None
        key = Element.generate_key(
            resolved_tag, self.version_name, "CallableLoci", targets=targ
        )  # noqa: E501
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(inputs),
            artifacts=ArtifactSet(
                output_summary,
                primary_name="summary",
                bed=output_bed,
                callable=output_json,
            ),
            pres=tuple(pres),
        )

    @staticmethod
    def extract_callable_mb(summary_file: Path) -> tuple[int, float]:
        callable_bases = 0

        with open(summary_file) as f:
            for line in f:
                if "CALLABLE" in line:
                    callable_bases = int(line.split("CALLABLE")[1].strip())
                    break

        callable_mb = callable_bases / 1000000
        return callable_bases, callable_mb
