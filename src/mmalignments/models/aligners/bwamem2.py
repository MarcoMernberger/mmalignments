"""Module contains an aligner interface for BWA-MEM2."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.data import Genome
from mmalignments.models.elements import (
    Element,
    MappedElement,
    NextGenSample,
    element,
    #    sample_fastqs,
)
from mmalignments.models.overlay import (
    CfgType,
    ElementTag,
    OutSpec,
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
from mmalignments.services.io import parents

from ..externals import External, ExternalRunConfig, SubroutineIn, subroutine
from ..parameters import ParamSet
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
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
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
        parameters_file = Path(__file__).parent / f"{Path(__file__).stem}.json"
        parameters = parameters or parameters_file
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
        tag: TagType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
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
        par : ParType | None
            Additional parameters for indexing call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            The Element representing the created index.
        """
        fasta_file = genome.fasta
        output_dir = output_dir or self.default_index_dir(genome)
        output_prefix = str(output_dir / self.index_prefix)
        runner = self.index_reference(
            fasta_file=fasta_file, output_prefix=output_prefix, par=par, cfg=cfg
        )
        resolved_tag = ElementTag(
            root=genome.name,
            level=1,
            omics=Omics.DNA,
            stage=Stage.PREP,
            method=Method.BWAMEM2,
            state=State.INDEX,
        ).resolve(tag)
        index_files = self.index_filenames_for_prefix(output_prefix)
        artifacts = ArtifactSet(
            output_prefix,
            primary_name="index_prefix",
            **{output_file.suffix[1:]: output_file for output_file in index_files},
        )
        determinants = self.signature_determinants(par, subroutine="index")
        key = Element.generate_key(resolved_tag, self.version_name, "index")

        result = Element(
            key,
            runner,
            tag=resolved_tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(fasta_file,),
        )
        return result

    @subroutine
    def index_reference(
        self,
        fasta_file: Path | str,
        output_prefix: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional parameters for indexing call.
        cfg : CfgType | None
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
        subcommand = "index"
        ret = (
            [subcommand, "-p", str(output_prefix), str(fasta_file)],
            subcommand,
            [fasta_file],
            [output_prefix],
            None,
            None,
            None,
        )
        print(ret)
        return ret

    ###########################################################################
    # Mapping (High-level and Ĺow-level wrapper)
    ###########################################################################

    @element
    def align(
        self,
        sample: NextGenSample,
        index: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        read_group: str | None = None,
        index_off: bool = True,
        post: list[Callable[[], CompletedProcess]] | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> MappedElement:
        """
        Align FASTQ reads to a reference genome using BWA-MEM2.

        Performs sequence alignment of single-end or paired-end reads to an
        indexed reference genome. Parameters can be overridden for this
        specific alignment via ``parameter``.

        Parameters
        ----------
        sample : NextGenSample
            Element containing FASTQ file paths and a sample name.
        index : Element
            Index element containing the BWA-MEM2 index prefix and metadata.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the aligned directory of
            the sample, ext ``bam``.
        read_group : str | None
            Read group header line (e.g., "@RG\\tID:sample1\\tSM:sample1").
        index_off : bool
            If True, skip BAM indexing after alignment. By default, the
            resulting BAM file will be indexed using Samtools, and the index
            file will be included in the post-processing steps.
        post : list[Callable[[], CompletedProcess]] | None
            Optional list of parameterless callables that will be executed
            sequentially after the alignment has completed successfully.
            Each callable must be a zero-argument function (like the
            return value of this ``align`` method).
        par : ParType | None
            Additional parameters for mem call (alignment).
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            The Element representing the aligned BAM file, with metadata and
            post-processing steps defined.

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
        index_prefix = Path(index.index_prefix).resolve()
        fastq_r1 = sample.r1
        fastq_r2 = sample.r2
        sample_name = sample.root
        rg = sample.read_group
        resolved_tag = sample.tag.bump(
            stage=Stage.ALIGN,
            method=Method.BWAMEM2,
            state=State.MAP,
            omics=Omics.DNA,
            flag=index.tag.root,
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(
                folder=self.default_aligned_dir(sample.root, index.tag.root),
                ext="bam",
            ),
            out,
        )
        output_bam = resolved_out.file
        rg = read_group or self.rg(sample_name)

        runner = self.align_fastq(
            index_prefix=index_prefix,
            output_bam=output_bam,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            read_group=rg,
            index_off=index_off,
            post=post,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par, subroutine="mem")
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "mem",
            index=f"{index_prefix.name}",
            rg=f"{rg}",
        )
        result = MappedElement(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,),
            artifacts=ArtifactSet(output_bam, primary_name="bam"),
            pres=(index,),
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
        index_off: bool = True,
        post: Callable[[], CompletedProcess] | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        index_off : bool
            If True, skip BAM indexing after alignment. By default, the
            resulting BAM file will be indexed using Samtools, and the index
            file will be included in the post-processing steps.
        post : list[Callable[[], CompletedProcess]] | None
            Optional list of parameterless callables that will be executed
            sequentially after the alignment has completed successfully.
            Each callable must be a zero-argument function (like the
            return value of this ``align`` method).
        par : ParType | None
            Additional parameters for mem call (alignment).
        cfg : CfgType | None
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
        index_prefix = Path(index_prefix)
        fastq_r1 = Path(fastq_r1)
        output_bam = Path(output_bam)

        # Build positional arguments: mem <index> <fastq_r1> [fastq_r2]
        arguments = ["mem"]

        # Add read group if provided
        if read_group:
            arguments.extend(["-R", read_group])
        # Add index and FASTQ files as positional arguments

        arguments.extend(
            [
                self.strabs(index_prefix),
                self.strabs(fastq_r1),
            ]
        )
        in_paths = [index_prefix, fastq_r1]
        if fastq_r2:
            arguments.append(self.strabs(fastq_r2))
            in_paths.append(Path(fastq_r2))

        output_bam = self.strabs(output_bam)
        # Output to BAM file
        arguments.extend(
            [
                "-o",
                output_bam,
            ]
        )

        if not index_off:
            compressed_bam_file = Path(str(output_bam) + ".bam")
            Samtools().bgfz_compress(
                bam_file=output_bam,
                compressed_bam_file=compressed_bam_file,
            )
            sam_index = Samtools().index_bam(
                bam_file=compressed_bam_file,
                bai_file=Path(str(output_bam) + ".bai"),
            )
            # calculate the index if needed. post however has priority
            post = post or sam_index
        subcommand = "mem"
        return arguments, subcommand, in_paths, [output_bam], None, None, post

    ###########################################################################
    # Con   venience methods for common workflows
    ###########################################################################

    def alignsort(
        self,
        sample: NextGenSample,
        indexer: Element,
        *,
        read_group: str | None = None,
        compress: bool = True,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> tuple[Element, ...]:
        """
        Align FASTQ reads to a reference genome, then sort and index the BAM.

        This convenience method chains alignment, sorting and indexing. If
        ``compress`` is True, the sorted BAM is first BGZF-compressed via
        samtools view before indexing.

        The sorted (or compressed-sorted) BAM and its index are returned as
        Elements; the unsorted mapped BAM is also returned.
        Parameters
        ----------
        sample : NextGenSample
            NextGenSample containing FASTQ file paths and a sample name.
        indexer : Element
            Indexer element for the reference genome.
        read_group : Optional[str], optional
            Read group information for the alignment, by default None
        compress : bool
            If True (default), compress the sorted BAM using
            ``samtools view -b`` before running ``samtools index``.
        tag : TagType | SpecValues[TagType] | None
            Tag overrides, optionally per step
            (``mem``, ``sort``, ``compress``, ``index``).
        out : OutType | SpecValues[OutType] | None
            Output specification overrides, optionally per step (``mem``,
            ``sort``, ``compress``, ``index``).
        par : ParType | SpecValues[ParType] | None
            Additional parameters, optionally per step (``mem``, ``sort``,
            ``compress``, ``index``).
        cfg : CfgType | SpecValues[CfgType] | None
            Run configuration, optionally per step (``mem``, ``sort``,
            ``compress``, ``index``).

        Notes
        -----
        A plain ``tag``, ``out``, ``par`` or ``cfg`` value is shared by all
        applicable steps. To override individual steps, pass a ``SpecValues``
        with ``valid_steps=("mem", "sort", "compress", "index")`` and put
        overrides in ``values``; for example::

            par = SpecValues(
                default=None,
                valid_steps=("mem", "sort", "compress", "index"),
                values={"mem": mem_value},
            )

        Unspecified steps use ``default``. Override names and lookups are
        validated against ``valid_steps``.

        Returns
        -------
        tuple[Element, Element, Element]
            ``(final_bam, final_bam_index, mapped_bam)`` where ``final_bam``
            is sorted BAM when ``compress=False`` and compressed BAM when
            ``compress=True``. The first Element is the terminal output of
            the workflow.
        """
        mapped = self.align(
            sample=sample,
            index=indexer,
            tag=pick_spec(tag, "mem"),
            out=pick_spec(out, "mem"),
            read_group=read_group,
            par=pick_spec(par, "mem"),
            cfg=pick_spec(cfg, "mem"),
        )

        st = Samtools()
        mapped_sorted = st.sort(
            mapped=mapped,
            tag=pick_spec(tag, "sort"),
            out=pick_spec(out, "sort"),
            par=pick_spec(par, "sort"),
            cfg=pick_spec(cfg, "sort"),
        )
        mapped_final = mapped_sorted
        if compress:
            mapped_final = st.compress(
                mapped=mapped_sorted,
                tag=pick_spec(tag, "compress"),
                out=pick_spec(out, "compress"),
                par=pick_spec(par, "compress"),
                cfg=pick_spec(cfg, "compress"),
            )
        mapped_sorted_index = st.index(
            mapped=mapped_final,
            tag=pick_spec(tag, "index"),
            out=pick_spec(out, "index"),
            par=pick_spec(par, "index"),
            cfg=pick_spec(cfg, "index"),
        )
        return mapped_final, mapped_sorted_index, mapped
