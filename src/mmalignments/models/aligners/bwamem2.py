"""Module contains an aligner interface for BWA-MEM2."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.data import Genome, Sample, sample_fastqs
from mmalignments.models.tags import ElementTag, Method, Omics, Stage, State
from mmalignments.models.tasks import Element, MappedElement, element
from mmalignments.services.io import parents

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet
from .samtools import Samtools

logger = logging.getLogger(__name__)


class BWAMem2(External):
    """BWA-MEM2 aligner interface.

    BWA-MEM2 is a faster reimplementation of BWA-MEM with identical output.
    This class provides methods for building genome indices and aligning
    sequencing reads (FASTQ) to a reference genome.

    Examples
    --------
    Build a genome index::

        aligner = BWAMem2()
        runner = aligner.index(
                fasta_file="genome.fa",
                output_prefix="genome_index/bwa_mem2_index"
        )

    Align paired-end reads::

        aligner = BWAMem2(parameters={"-t": 8})
        runner = aligner.align(
                index_prefix="genome_index/bwa_mem2_index",
                fastq_r1="sample_R1.fastq.gz",
                fastq_r2="sample_R2.fastq.gz",
                output_sam="aligned.sam"
        )
    """

    def __init__(
        self,
        name: str = "bwa-mem2",
        primary_binary: str = "bwa-mem2",
        version: str | None = None,
        folder: Path | None = None,
        source: str | None = None,
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
        genome_name: str | None = None,
    ) -> None:
        """Initialize BWA-MEM2 aligner.

        Parameters
        ----------
        name : str
            Tool name (default: "bwa-mem2").
        primary_binary : str
            Binary executable name (default: "bwa-mem2").
        version : Optional[str]
            Version string override.
        folder: Path | None
            Optional path to a folder where the tool should write its output.
        source : str | None
            URL/source for the tool.
        genome_name: str | None
            Name of the genome being aligned against.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        if source is None:
            source = "https://github.com/bwa-mem2/bwa-mem2"
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            folder=folder,
            source=source,
            parameters=parameters or {},
        )
        self.index_extensions = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        self.index_prefix = f"{self.name}_index"

        if folder is None:
            folder = Path("results") / "alignments" / self.version_name
            if genome_name is not None:
                folder = folder / genome_name
        self._folder = folder

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get BWA-MEM2 version string.

        Parameters
        ----------
        fallback : Optional[str]
            Value to return if version cannot be determined (default: None).

        Returns
        -------
        Optional[str]
            Version string (e.g., "2.2.1") or fallback if not found.
        """
        if not self.primary_binary or not self.ensure_binary():
            return fallback

        try:
            cp = subprocess.run(
                [self.primary_binary, "version"],
                capture_output=True,
                text=True,
                check=True,
            )
            output = cp.stdout.strip() or cp.stderr.strip()
            if output:
                return output
        except subprocess.CalledProcessError:
            return fallback

    ###########################################################################
    # Helpers
    ###########################################################################

    def _configure(self, output_bam: Path) -> ExternalRunConfig:
        cfg = ExternalRunConfig()
        cfg.cwd = output_bam.parent
        return cfg

    def default_aligned_dir(self, sample_name: str, reference_name) -> Path:
        """Return directory path for aligned BAM files for a given sample.

        Parameters
        ----------
        sample_name : str
                Name of the sample.

        Returns
        -------
        Path
                Directory path where aligned BAM files for the sample are stored.
        """
        return (
            Path("results")
            / "aligned"
            / f"{self.version_name}"
            / f"{reference_name}"
            / sample_name
        )

    def default_index_dir(self, genome: Genome) -> Path:
        """Return directory path for BWA-MEM2 index files for a given genome.

        Parameters
        ----------
        genome : Genome
            The genome object.

        Returns
        -------
        Path
            Directory path where BWA-MEM2 index files for the genome are stored.
        """
        return Path(genome.base) / self.name / self.version_name / genome.name / "index"

    def index_filenames_for_prefix(self, prefix: Path | str) -> list[Path]:
        """Return list of expected index files for a given prefix.

        Parameters
        ----------
        prefix : Path | str
                Index prefix path.

        Returns
        -------
        list[Path]
                List of expected index file paths.
        """
        prefix = Path(prefix)
        return [Path(str(prefix) + ext) for ext in self.index_extensions]

    def index_exists(self, prefix: Path | str) -> bool:
        """Check if all required index files exist.

        Parameters
        ----------
        prefix : Path | str
                Index prefix path.

        Returns
        -------
        bool
                True if all index files exist; False otherwise.
        """
        return all(f.exists() for f in self.index_filenames_for_prefix(prefix))

    def rg(self, sample_name: str) -> str:
        """
        Generate a read group (RG) string for a given sample.

        Parameters
        ----------
        sample_name : str
            Name of the sample.

        Returns
        -------
        str
            Read group string formatted for BWA-MEM2.
        """
        return f"@RG\\tID:{sample_name}_rg1\\tSM:{sample_name}\\tLB:{sample_name}_251215\\tPL:ILLUMINA\\tPU:unit1"  # noqa: E501

    ###########################################################################
    # Indexing (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def index(
        self,
        genome: Genome,
        *,
        output_dir: Path | None = None,
        prefix: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """
        Create an Element representing the BWA-MEM2 index for a given reference
        genome.

        This is a high level function that wraps the low-level
        ``index_reference`` method. It checks paths, constructs the output
        prefix, and defines the expected index files as artifacts. The returned
        Element can be used as a prerequisite for alignment steps.
        Also takes care of determinants and signatures for caching and reruns.

        Parameters
        ----------
        genome : Genome
            The reference genome for which the index is being created.
        output_prefix : Path | None
            Optional prefix for the output index files.
        params : Params | None
            Additional parameters for indexing call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            The Element representing the created index.
        """
        fasta_file = genome.fasta
        output_dir = output_dir or self.default_index_dir(genome)
        output_prefix = output_dir / self.index_prefix
        runner = self.index_reference(
            fasta_file=fasta_file, output_prefix=output_prefix, params=params, cfg=cfg
        )
        tag = ElementTag(
            root=genome.name,
            level=0,
            omics=Omics.DNA,
            stage=Stage.PREP,
            method=Method.BWAMEM2,
            state=State.INDEX,
            ext=None,
        )
        index_files = self.index_filenames_for_prefix(output_prefix)
        key = f"{tag.default_name}_index_{self.version_name}"
        artifacts = {output_file.suffix[1:]: output_file for output_file in index_files}
        artifacts.update(
            {
                "index_prefix": output_prefix,
            }
        )
        determinants = self.signature_determinants(params, subroutine="index")

        result = Element(
            key,
            runner,
            tat=tag,
            name=None,  # f"{genome.name} index ({self.version_name})",
            artifacts=artifacts,
            determinants=determinants,
            inputs=[fasta_file],
            store_attributes={"multi_core": False},
        )
        return result

    @subroutine
    def index_reference(
        self,
        fasta_file: Path | str,
        output_prefix: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """
        Build BWA-MEM2 index from a reference FASTA file.

        Creates index files required for alignment. The index files will be
        created with the specified prefix (e.g., prefix.0123, prefix.amb,
        prefix.ann, prefix.bwt.2bit.64, prefix.pac).

        Parameters
        ----------
        fasta_file : Path | str
                Path to the reference genome FASTA file.
        output_prefix : Path | str
                Prefix for index files (e.g., "genome_index/bwa_mem2_index").
        params : Params | None
            Additional parameters for indexing call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Callable[[], subprocess.CompletedProcess]
                A zero-argument callable that builds the index when invoked.
                The callable is also executed immediately once.

        Examples
        --------
        >>> aligner = BWAMem2()
        >>> runner = aligner.index("genome.fa", "index/genome")
        >>> # Index is already built; can re-run with: runner()
        """
        fasta_file = Path(fasta_file).absolute()
        output_prefix = Path(output_prefix).absolute()

        # Ensure output directory exists
        parents(output_prefix)

        # bwa-mem2 index -p <prefix> <fasta>
        return (
            ["index", "-p", str(output_prefix), str(fasta_file)],
            [output_prefix],
            None,
            None,
            None,
        )

    ###########################################################################
    # Mapping (High-level and Ĺow-level wrapper)
    ###########################################################################

    @element
    def align(
        self,
        sample: Element | Sample,
        index: Element,
        *,
        output_bam: Path | str | None = None,
        read_group: str | None = None,
        post: list[Callable[[], CompletedProcess]] | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """
        Align FASTQ reads to a reference genome using BWA-MEM2.

        Performs sequence alignment of single-end or paired-end reads to an
        indexed reference genome. Parameters can be overridden for this
        specific alignment via ``parameter``.

        Parameters
        ----------
        sample : Element | Sample
            Sample or element containing FASTQ file paths and a sample name.
        index : Element
            Index element containing the BWA-MEM2 index prefix and metadata.
        output_bam : Path | str | None
            Path to output BAM file.
        read_group : str | None
            Read group header line (e.g., "@RG\\tID:sample1\\tSM:sample1").
        post : list[Callable[[], CompletedProcess]] | None
            Optional list of parameterless callables that will be executed
            sequentially after the alignment has completed successfully.
            Each callable must be a zero-argument function (like the
            return value of this ``align`` method).
        params : Params | None
            Additional parameters for mem call (alignment).
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Examples
        --------
        Single-end alignment::

            aligner = BWAMem2(parameters={"-t": 8})
            runner = aligner.align(
                    index_prefix="genome_index/bwa",
                    fastq_r1="sample.fastq.gz",
                    output_sam="aligned.sam"
            )

        Paired-end alignment with read group::

            runner = aligner.align(
                    index_prefix="genome_index/bwa",
                    fastq_r1="sample_R1.fastq.gz",
                    fastq_r2="sample_R2.fastq.gz",
                    output_sam="aligned.sam",
                    read_group="@RG\\tID:lane1\\tSM:sample1\\tPL:ILLUMINA"
            )

        With post-processing callbacks::

            def sam_to_bam():
                import subprocess
                subprocess.run(
                    ["samtools", "view", "-b", "-o", "aligned.bam", "aligned.sam"],
                    check=True
                )

        # pass a list of zero-argument callables; they will be run in
        # the order given after the alignment finished successfully
        runner = aligner.align(
            index_prefix="genome_index/bwa",
            fastq_r1="sample_R1.fastq.gz",
            fastq_r2="sample_R2.fastq.gz",
            output_sam="aligned.sam",
            post=[sam_to_bam]
        )
        """
        index_prefix = index.index_prefix.absolute()
        fastq_r1, fastq_r2, sample_name, rg = sample_fastqs(sample)
        output_bam = (
            output_bam
            or self.default_aligned_dir(sample_name, index.name)
            / f"{sample_name}_aligned.bam"
        ).resolve()

        rg = read_group or self.rg(sample_name)

        runner = self.align_fastq(
            index_prefix=index_prefix,
            output_bam=output_bam,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            read_group=rg,
            post=post,
            params=params,
            cfg=cfg,
        )
        name = f"{sample_name}_mapped_to_{index.name}_by_{self.name}"
        key = f"{sample_name}_mapped_to_{index.index_prefix}_by_{self.version_name}_rg={rg}"  # noqa: E501
        determinants = self.signature_determinants(params, subroutine="mem")

        result = MappedElement(
            name,
            key,
            runner,
            determinants=determinants,
            inputs=[fastq_r1] + ([fastq_r2] if fastq_r2 else []),
            artifacts={"bam": output_bam},
            pres=[index],
            store_attributes={"multi_core": True},
        )
        return result

    @subroutine
    def align_fastq(
        self,
        index_prefix: Path | str,
        output_bam: Path | str,
        fastq_r1: Path | str,
        *,
        fastq_r2: Path | str | None = None,
        read_group: str | None = None,
        post: list[Callable[[], CompletedProcess]] | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Align FASTQ reads to a reference genome using BWA-MEM2.

        Performs sequence alignment of single-end or paired-end reads to an
        indexed reference genome. Parameters can be overridden for this
        specific alignment via ``parameter``.

        Parameters
        ----------
        index_prefix : Path | str
            Path to the BWA-MEM2 index prefix (without file extensions).
        fastq_r1 : Path | str
            Path to R1 (forward) FASTQ file or single-end FASTQ.
        output_bam : Path | str
            Path to output BAM file.
        fastq_r2 : Path | str | None
            Path to R2 (reverse) FASTQ file for paired-end data.
        read_group : str | None
            Read group header line (e.g., "@RG\\tID:sample1\\tSM:sample1").
        post : list[Callable[[], CompletedProcess]] | None
            Optional list of parameterless callables that will be executed
            sequentially after the alignment has completed successfully.
            Each callable must be a zero-argument function (like the
            return value of this ``align`` method).
        params : Params | None
            Additional parameters for mem call (alignment).
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], subprocess.CompletedProcess]
            A zero-argument callable that performs the alignment when invoked.
            The callable is also executed immediately once.

        Examples
        --------
        Single-end alignment::

            aligner = BWAMem2(parameters={"-t": 8})
            runner = aligner.align(
                    index_prefix="genome_index/bwa",
                    fastq_r1="sample.fastq.gz",
                    output_sam="aligned.sam"
            )

        Paired-end alignment with read group::

            runner = aligner.align(
                    index_prefix="genome_index/bwa",
                    fastq_r1="sample_R1.fastq.gz",
                    fastq_r2="sample_R2.fastq.gz",
                    output_sam="aligned.sam",
                    read_group="@RG\\tID:lane1\\tSM:sample1\\tPL:ILLUMINA"
            )

        With post-processing callbacks::

            def sam_to_bam():
                import subprocess
                subprocess.run(
                    ["samtools", "view", "-b", "-o", "aligned.bam", "aligned.sam"],
                    check=True
                )

            # pass a list of zero-argument callables; they will be run in
            # the order given after the alignment finished successfully
            runner = aligner.align(
                index_prefix="genome_index/bwa",
                fastq_r1="sample_R1.fastq.gz",
                fastq_r2="sample_R2.fastq.gz",
                output_sam="aligned.sam",
                post=[sam_to_bam]
            )
        """
        index_prefix = Path(index_prefix).absolute()
        fastq_r1 = Path(fastq_r1).absolute()
        output_bam = Path(output_bam).absolute()
        # Build positional arguments: mem <index> <fastq_r1> [fastq_r2]
        arguments = ["mem"]
        # Add read group if provided
        if read_group:
            arguments.extend(["-R", read_group])
        # Add index and FASTQ files as positional arguments
        arguments.extend(
            [
                str(index_prefix),
                str(fastq_r1),
            ]
        )
        if fastq_r2:
            fastq_r2 = Path(fastq_r2).absolute()
            arguments.append(str(fastq_r2))

        # Output to BAM file
        arguments.extend(["-o", str(output_bam)])
        return arguments, [output_bam], None, None, None

    ###########################################################################
    # Con   venience methods for common workflows
    ###########################################################################

    def alignsort(
        self,
        sample: Sample,
        genome: Genome,
        indexer: Element,
        *,
        output_bam: Path | str | None = None,
        read_group: str | None = None,
        parameters: Mapping[str, Params] | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[Element, Element]:
        """
        Align FASTQ reads to a reference genome and sort the resulting BAM file.

        This is a high-level convenience method that performs both alignment
        and sorting in one step. It invokes the ``align`` method to perform the
        alignment, and then uses Samtools to sort the resulting BAM file.
        The sorted BAM file and its index will be the final outputs, wrapped in
        output Elements.
        Parameters
        ----------
        sample : Sample
            Sample containing FASTQ file paths and a sample name.
        genome : Genome
            Genome reference for alignment.
        indexer : Element
            Indexer element for the reference genome.
        output_bam : Optional[Path | str]
            Optional path to the output BAM file. If not provided, a default
            path will be used based on the sample and genome names.
        read_group : Optional[str], optional
            Read group information for the alignment, by default None
        parameters : Mapping[str, Params] | None, optional
            Additional parameters for mem call (alignment) and subsequent
            sorting.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        tuple[Element, Element]
            A tuple containing:
            - The Element representing the BAM index file.
            - The Element representing the sorted BAM file.
            First return value is the final Element.
        """
        parameters = parameters or {}
        output_bam = (
            output_bam
            or self.default_aligned_dir(sample.name, genome.name)
            / f"{sample.name}_aligned.bam"
        )
        output_bam = Path(output_bam).resolve()
        mapped = self.align(
            sample=sample,
            index=indexer,
            output_bam=output_bam,
            read_group=read_group,
            params=parameters.get("mem", None),
            cfg=cfg,
        )

        st = Samtools()
        sorted_bam = output_bam.parent / f"{sample.name}_sorted.bam"
        mapped_sorted = st.sort(
            mapped=mapped,
            output_bam_file=sorted_bam,
            params=parameters.get("sort", Params(threads=10)),
            cfg=ExternalRunConfig(cwd=output_bam.parent, threads=10),
        )
        mapped_sorted_index = st.index(
            mapped=mapped_sorted,
            params=parameters.get("index", None),
            cfg=cfg,
        )
        return mapped_sorted, mapped_sorted_index
