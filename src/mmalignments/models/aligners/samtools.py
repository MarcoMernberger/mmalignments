"""Module contains a samtools interface for post-alignment processing."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Callable, Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.data import Genome
from mmalignments.models.elements import (
    Element,
    MappedElement,
    element,
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
    Stage,
    State,
)
from mmalignments.services.io import parents

from ..externals import External, SubroutineIn, subroutine
from ..parameters import ParamSet

logger = logging.getLogger(__name__)


class Samtools(External):
    """Samtools interface for BAM file processing.

    Provides methods for sorting and indexing BAM files using samtools.
    Returns zero-argument callables that can be chained as post-processing
    steps after alignment.

    Examples
    --------
    Sort and index a BAM file::

        st = Samtools()
        sort_runner = st.sort(input_bam="aligned.bam", output_bam="sorted.bam")
        index_runner = st.index(bam_file="sorted.bam")

    Chain with BWA-MEM2 alignment::

        from mmalignments.models.aligners.bwamem2 import BWAMem2

        aligner = BWAMem2()
        st = Samtools()
        output_bam = Path("aligned.bam")
        sorted_bam = Path("sorted.bam")

        # Create post-processing callables
        sort_callable = st.sort(input_bam=output_bam, output_bam=sorted_bam)
        index_callable = st.index(bam_file=sorted_bam)

        # Pass as post list to align
        runner = aligner.align(
                index_prefix="genome_index/bwa",
                fastq_r1="sample_R1.fastq.gz",
                fastq_r2="sample_R2.fastq.gz",
                output_bam=output_bam,
                post=[sort_callable, index_callable]
        )
    """

    def __init__(
        self,
        name: str = "samtools",
        primary_binary: str = "samtools",
        version: str | None = None,
        source: str = "https://github.com/samtools/samtools",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """Initialize Samtools wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "samtools").
        primary_binary : str
            Binary executable name (default: "samtools").
        version : str | None
            Version string override.
        source : str
            URL/source for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        parameters_file = Path(__file__).parent / f"{primary_binary}.json"
        parameters = parameters or parameters_file
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Get samtools version string.

        Parameters
        ----------
        fallback : Optional[str]
                Value to return if version cannot be determined (default: None).

        Returns
        -------
        Optional[str]
                Version string (e.g., "1.18") or fallback if not found.
        """
        if not self.primary_binary or not self.ensure_binary():
            return fallback

        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            output = cp.stdout.strip()
            if output:
                first_line = output.splitlines()[0]
                parts = first_line.split()
                if len(parts) >= 2 and parts[0] == "samtools":
                    return parts[1]
        except subprocess.CalledProcessError:
            return fallback

        return fallback

    ###########################################################################
    # Sort (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def sort(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> MappedElement:
        """Sort a BAM file using samtools sort.

        Creates a zero-argument callable that sorts the input BAM file and
        writes the sorted output to the specified path. The callable can be
        passed to BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file to be sorted.
        tag : TagType | None
            Partial or full Element tag override. Defaults are derived from the
            input Element's tag.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BAM and extension ``bam``.
        par : ParType | None
            Parameters for sorting call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        MappedElement
            New MappedElement representing the sorted BAM file, with the
            sort command as its run callable and the original mapped element
            as a pre-requisite.

        Examples
        --------
        >>> st = Samtools()
        >>> sort_runner = st.sort(mapped_element)
        >>> sort_runner()  # Execute the sort
        """
        input_bam = mapped.bam
        tag = mapped.tag.bump(method=Method.SAMTOOLS, state=State.SORT).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(folder=input_bam.parent, ext="bam"), out
        )
        sorted_bam = out.file

        runner = self.sort_bam(
            input_bam=input_bam,
            output_bam=sorted_bam,
            cfg=cfg,
            par=par,
        )
        key = Element.generate_key(tag, self.version_name, "sort")
        return MappedElement(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="sort"),
            inputs=(input_bam,),
            artifacts=ArtifactSet(sorted_bam, primary_name="bam"),
            pres=(mapped,),
        )

    @subroutine
    def sort_bam(
        self,
        input_bam: Path | str,
        output_bam: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Sort a BAM file using samtools sort.

        Creates a zero-argument callable that sorts the input BAM file and
        writes the sorted output to the specified path. The callable can be
        passed to BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        input_bam : Path | str
            Path to the unsorted BAM file.
        output_bam : Path | str
            Path for the sorted output BAM file.
        par : ParType | None
            Parameters for sorting call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that performs the sort when invoked.

        Examples
        --------
        >>> st = Samtools()
        >>> sort_runner = st.sort_bam("aligned.bam", "sorted.bam")
        >>> sort_runner()  # Execute the sort
        """
        # Build command: samtools sort -o <output> <input>
        subcommand = "sort"
        arguments = [subcommand, "-o", self.strabs(output_bam), self.strabs(input_bam)]
        return arguments, subcommand, [input_bam], [output_bam], None, None, None

    ###########################################################################
    # Compress (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def compress(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> MappedElement:
        """Compress a BAM file using samtools view -b.

        Converts an uncompressed BAM into a BGZF-compressed BAM suitable for
        ``samtools index``.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file to compress.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``<bam>.compressed.bam``
            next to the input BAM.
        par : ParType | None
            Parameters for compression call.
        cfg : CfgType | None
            Optional configuration for running the subprocess.

        Returns
        -------
        MappedElement
            New MappedElement representing the compressed BAM.
        """
        input_bam = mapped.bam
        tag = mapped.tag.bump(
            method=Method.SAMTOOLS, state=State.TRANSFORMED, flag="compress"
        ).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(
                folder=input_bam.parent, stem=f"{input_bam.stem}.compressed", ext="bam"
            ),
            out,
        )
        compressed_bam = out.file

        runner = self.bgfz_compress(
            input_bam=input_bam,
            output_bam=compressed_bam,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(tag, self.version_name, "compress")
        return MappedElement(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="view"),
            inputs=(input_bam,),
            artifacts=ArtifactSet(compressed_bam, primary_name="bam"),
            pres=(mapped,),
        )

    @subroutine
    def bgfz_compress(
        self,
        input_bam: Path | str,
        output_bam: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Compress a BAM file to BGZF BAM using samtools view.

        Parameters
        ----------
        input_bam : Path | str
            Path to the BAM file to compress.
        output_bam : Path | str
            Path for the compressed BAM output.
        par : ParType | None
            Parameters for compression call.
        cfg : CfgType | None
            Optional configuration for running the subprocess.

        Returns
        -------
        SubroutineIn
            Subroutine tuple for ``samtools view -b <input> -o <output>``.
        """
        arguments = [
            "view",
            "-b",
            self.strabs(input_bam),
            "-o",
            self.strabs(output_bam),
        ]
        subcommand = "view"
        return arguments, subcommand, [input_bam], [output_bam], None, None, None

    ###########################################################################
    # Index (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def index(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create a BAM index (.bai) using samtools index.

        Creates a zero-argument callable that indexes the BAM file, producing
        a .bai file alongside the input. The callable can be passed to
        BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file to index.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``<bam>.bai`` next to
            the input BAM.
        par : ParType | None
            Parameters for indexing call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            Zero-argument callable that performs the indexing when invoked.

        Examples
        --------
        >>> st = Samtools()
        >>> index_runner = st.index(mapped_element)
        >>> index_runner()  # Execute the index
        """
        bam_file = mapped.bam
        tag = mapped.tag.bump(method=Method.SAMTOOLS, state=State.INDEX).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(folder=bam_file.parent, stem=bam_file.name, ext="bai"), out
        )
        bai_file = out.file

        runner = self.index_bam(
            bam_file=bam_file,
            bai_file=bai_file,
            par=par,
            cfg=cfg,
        )

        key = Element.generate_key(tag, self.version_name, "index")
        return Element(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="index"),
            inputs=(bam_file,),
            artifacts=ArtifactSet(bai_file, primary_name="bai"),
            pres=(mapped,),
        )

    @subroutine
    def index_bam(
        self,
        bam_file: Path | str,
        *,
        bai_file: Path | str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Create a BAM index (.bai) using samtools index.

        Creates a zero-argument callable that indexes the BAM file, producing
        a .bai file alongside the input. The callable can be passed to
        BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        bam_file : Path | str
            Path to the BAM file to index.
        par : ParType | None
            Parameters for indexing call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], CompletedProcess]
                Zero-argument callable that performs the indexing when invoked.

        Examples
        --------
        >>> st = Samtools()
        >>> index_runner = st.index_bam("sorted.bam")
        >>> index_runner()  # Execute the index
        """
        bam_file = Path(bam_file).resolve()
        bai_file = bai_file or bam_file.with_suffix(bam_file.suffix + ".bai")
        # Build command: samtools index <bam_file>
        arguments = ["index", self.strabs(bam_file), "-o", self.strabs(bai_file)]
        subcommand = "index"
        return arguments, subcommand, [bam_file], [bai_file], None, None, None

    ###########################################################################
    # Faidx (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def faidx(
        self,
        genome: Genome,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        cut2: bool = True,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        Run faidx on a Genome-like object and return an Element.

        Parameters
        ----------
        genome : Genome
            Genome-like object containing the reference FASTA file.
        tag : TagType | None
            Partial or full Element tag override. The default tag is derived
            from the genome name.
        out : OutType | None
            Output specification override. Defaults to the FASTA file's parent
            directory and extension ``genome``.
        cut2 : bool, optional
            Whether to cut the .fai output to create a .genome file
            (default: True).
            If True, a .genome file will be created with the same base name as
            the output_file, containing the first two columns of the .fai
            (chromosome name and length).
        par : ParType | None
            Parameters for the faidx call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            Element representing the faidx operation.
        """
        fasta = (genome.fasta).absolute()
        tag = ElementTag(
            root=genome.name,
            level=1,
            stage=Stage.PREP,
            method=Method.SAMTOOLS,
            state=State.INDEX,
            omics=None,
        ).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(folder=fasta.parent, ext="genome"), out
        )
        output = out.file

        runner = self.faidx_fasta(
            fasta=fasta, output_file=output, cut2=cut2, par=par, cfg=cfg
        )
        artifacts = ArtifactSet(output, primary_name="fai")
        if cut2:
            artifacts = artifacts.with_extra("faicut", output.with_suffix(".faicut"))
        key = Element.generate_key(tag, self.version_name, "faidx", cut2=cut2)
        return Element(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="faidx"),
            inputs=(fasta,),
            artifacts=artifacts,
            pres=(),
        )

    @subroutine
    def faidx_fasta(
        self,
        fasta: Path | str,
        output_file: Path | str,
        cut2: bool = True,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Create a fasta index (.fai) with samtools faidx and write a two-column genome
        file.

        Parameters
        ----------
        fasta : Path | str
            Path to the reference FASTA file.
        output_file : Path | str
            Output two-column genome file (chr\tlength).
        cut2 : bool
            Whether to cut the .fai output to create a .genome file (default: True).
            If True, a .genome file will be created with the same base name as
            the output_file, containing the first two columns of the .fai
            (chromosome name and length).
        cwd : Optional[Path]
            Working directory for the samtools command (default: parent of fasta).
        parameters : Any
            Additional parameter overrides.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that performs the actions when invoked.
        """

        def _faidx_and_cut_callable(genome_out: Path) -> Callable:
            fai_path = Path(str(fasta) + ".fai")
            """Create genome file from the .fai after indexing."""

            def __post_call():
                # Read .fai and write first two columns to genome_out
                with (
                    fai_path.open("r", encoding="utf-8") as fin,
                    genome_out.open("w", encoding="utf-8") as fout,
                ):
                    for line in fin:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 2:
                            fout.write(f"{parts[0]}\t{parts[1]}\n")

            return __post_call

        fasta = Path(fasta).absolute()
        output_file = Path(output_file).absolute()

        # Ensure output directory exists
        parents(output_file.parent)

        # Build samtools faidx command
        arguments = ["faidx", str(fasta)]

        call_afterwards = None
        if cut2:
            call_afterwards = _faidx_and_cut_callable(
                output_file.with_suffix(".faicut")
            )
        subcommand = "faidx"
        return (
            arguments,
            subcommand,
            [fasta],
            [output_file],
            None,
            None,
            call_afterwards,
        )

    ###########################################################################
    # Stats (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def stats(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> "Element":
        """Run samtools stats from a MappedElement.

        Creates an Element that runs samtools stats to collect comprehensive
        alignment statistics.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``qc/samtools/``
            relative to the BAM file's directory and extension ``txt``.
        par : ParType | None
            Parameters for the stats call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            An Element that executes samtools stats when run.

        Examples
        --------
        >>> st = Samtools()
        >>> stats_elem = st.stats(mapped)
        >>> stats_elem.run()
        """
        input_bam = mapped.bam
        tag = mapped.tag.bump(
            stage=Stage.QC, method=Method.SAMTOOLS, state=State.STAT
        ).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(folder=input_bam.parent / "qc" / "samtools", ext="txt"), out
        )
        output_file = out.file

        runner = self.stats_bam(
            input_bam=input_bam,
            output_file=output_file,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(tag, self.version_name, "stats")
        return Element(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="stats"),
            inputs=(input_bam,),
            artifacts=ArtifactSet(output_file, primary_name="stats"),
            pres=(mapped,),
        )

    @subroutine
    def stats_bam(
        self,
        input_bam: Path | str,
        output_file: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run samtools stats on a BAM file.

        Creates a zero-argument callable that runs samtools stats to
        collect comprehensive alignment statistics.

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_file : Path | str
            Output file for stats results.
        par : ParType | None
            Parameters for the stats call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).
        Returns
        -------
        Callable[[], CompletedProcess]
                Zero-argument callable that executes samtools stats.

        Examples
        --------
        >>> st = Samtools()
        >>> stats = st.stats_bam("sample.bam", "qc/stats.txt")
        >>> stats()
        """
        # samtools stats input.bam > output.txt
        arguments = ["stats", self.strabs(input_bam)]
        subcommand = "stats"
        return (
            arguments,
            subcommand,
            [input_bam],
            [output_file],
            output_file,
            None,
            None,
        )

    ###########################################################################
    # Flagstat (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def flagstat(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> "Element":
        """Run samtools flagstat from a MappedElement.

        Creates an Element that runs samtools flagstat to collect basic
        alignment statistics.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to ``qc/samtools/``
            relative to the BAM file's directory and extension ``json``
            (``json`` or ``tsv`` are supported).
        par : ParType | None
            Parameters for the flagstat call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            An Element that executes samtools flagstat when run.

        Examples
        --------
        >>> st = Samtools()
        >>> flagstat_elem = st.flagstat(mapped)
        >>> flagstat_elem.run()
        """
        input_bam = mapped.bam
        tag = mapped.tag.bump(
            stage=Stage.QC, method=Method.SAMTOOLS, state=State.STAT, flag="flagstat"
        ).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(
            OutSpec(folder=input_bam.parent / "qc" / "samtools", ext="json"), out
        )
        output_file = out.file
        runner = self.flagstat_bam(
            input_bam=input_bam,
            output_file=output_file,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(tag, self.version_name, "flagstat")

        return Element(
            key,
            runner,
            tag=tag,
            determinants=self.signature_determinants(par, subroutine="flagstat"),
            inputs=(input_bam,),
            artifacts=ArtifactSet(output_file, primary_name="flagstat"),
            pres=(mapped,),
        )

    @subroutine
    def flagstat_bam(
        self,
        input_bam: Path | str,
        output_file: Path | str,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run samtools flagstat on a BAM file.

        Creates an Element that runs samtools flagstat to
        collect basic alignment statistics (mapped/unmapped, properly paired, etc.).

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_file : Path | str
            Output file for flagstat results.
        par : ParType | None
            Parameters for the flagstat call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes samtools flagstat.

        Examples
        --------
        >>> st = Samtools()
        >>> flagstat = st.flagstat_bam("sample.bam", "qc/flagstat.txt")
        >>> flagstat()
        """
        suffix = Path(output_file).suffix.lower()
        if suffix not in [".json", ".tsv"]:
            raise ValueError(
                f"Output file extension {suffix} not supported for samtools stats. "
                "Use .json or .tsv."
            )

        # samtools flagstat input.bam > output.txt
        arguments = ["flagstat", str(input_bam), "-O", suffix[1:]]
        subcommand = "flagstat"
        return (
            arguments,
            subcommand,
            [input_bam],
            [output_file],
            output_file,
            None,
            None,
        )

    ###########################################################################
    # Convenience methods for common workflows
    ###########################################################################

    def compress_index(
        self,
        mapped: MappedElement,
        *,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> tuple[MappedElement, Element]:
        """Compress a BAM and then create its index.

        Parameters
        ----------
        mapped : MappedElement
            Input mapped BAM element.
        tag : TagType | SpecValues[TagType] | None
            Tag overrides, optionally per step (``compress``, ``index``).
        out : OutType | SpecValues[OutType] | None
            Output overrides, optionally per step (``compress``, ``index``).
        par : ParType | SpecValues[ParType] | None
            Parameter overrides, optionally per step (``compress``, ``index``).
        cfg : CfgType | SpecValues[CfgType] | None
            Run configuration overrides, optionally per step
            (``compress``, ``index``).

        Returns
        -------
        tuple[MappedElement, Element]
            ``(compressed_bam, compressed_bam_index)``.
        """
        compressed = self.compress(
            mapped=mapped,
            tag=pick_spec(tag, "compress"),
            out=pick_spec(out, "compress"),
            par=pick_spec(par, "compress"),
            cfg=pick_spec(cfg, "compress"),
        )
        compressed_index = self.index(
            mapped=compressed,
            tag=pick_spec(tag, "index"),
            out=pick_spec(out, "index"),
            par=pick_spec(par, "index"),
            cfg=pick_spec(cfg, "index"),
        )
        return compressed, compressed_index
