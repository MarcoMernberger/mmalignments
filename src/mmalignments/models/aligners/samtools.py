"""Module contains a samtools interface for post-alignment processing."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.data import Genome
from mmalignments.models.elements import Element, MappedElement, element
from mmalignments.models.tags import ElementTag, Method, Stage, State, merge_tag
from mmalignments.models.tags import PartialElementTag, from_prior
from mmalignments.services.io import parents

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet

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
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> MappedElement:
        """Sort a BAM file using samtools sort.

        Creates a zero-argument callable that sorts the input BAM file and
        writes the sorted output to the specified path. The callable can be
        passed to BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file to be sorted.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        outdir : Path | str | None
            Directory for the sorted BAM file. If not provided, defaults to
            the same directory as the input BAM file.
        filename : Path | str | None
            Filename override. If not provided, defaults to
            ``tag.default_output``.
        params : Params | None
            Parameters for sorting call.
        cfg : ExternalRunConfig | None
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
        tag = from_prior(mapped.tag, tag, method=Method.SAMTOOLS, state=State.SORT)
        sorted_bam = Path(outdir or input_bam.parent) / (filename or tag.default_output)

        # Create the runner by calling External.run
        runner = self.sort_bam(
            input_bam=input_bam,
            output_bam=sorted_bam,
            cfg=cfg,
            params=params,
        )
        element = MappedElement(
            name=f"{mapped.name}.sorted",
            key=f"{mapped.name}.({self.name}.sort).{self.version_name}",
            run=runner,
            tag=tag,
            determinants=self.signature_determinants(params, subroutine="sort"),
            inputs=[input_bam],
            artifacts={"bam": sorted_bam},
            pres=[mapped],
        )
        return element

    @subroutine
    def sort_bam(
        self,
        input_bam: Path | str,
        output_bam: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
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
        params : Params | None
            Parameters for sorting call.
        cfg : ExternalRunConfig | None
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
        arguments = ["sort", "-o", self.strabs(output_bam), self.strabs(input_bam)]
        return arguments, [output_bam], None, None, None

    ###########################################################################
    # Index (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def index(
        self,
        mapped: MappedElement,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Create a BAM index (.bai) using samtools index.

        Creates a zero-argument callable that indexes the BAM file, producing
        a .bai file alongside the input. The callable can be passed to
        BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file to index.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        params : Params | None
            Parameters for indexing call.
        cfg : ExternalRunConfig | None
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
        tag = from_prior(
            mapped.tag,
            tag=tag,
            level=mapped.tag.level,
            method=Method.SAMTOOLS,
            state=State.INDEX,
            ext="bai",
        )
        # must always be the same as the sorted file
        bai_file = bam_file.with_suffix(bam_file.suffix + ".bai")

        runner = self.index_bam(
            bam_file=bam_file,
            bai_file=bai_file,
            params=params,
            cfg=cfg,
        )

        key = f"{tag.default_name}_index_{self.version_name}"
        element = Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=self.signature_determinants(params, subroutine="index"),
            inputs=[bam_file],
            artifacts={tag.ext: bai_file},
            pres=[mapped],
        )
        return element

    @subroutine
    def index_bam(
        self,
        bam_file: Path | str,
        *,
        bai_file: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Create a BAM index (.bai) using samtools index.

        Creates a zero-argument callable that indexes the BAM file, producing
        a .bai file alongside the input. The callable can be passed to
        BWAMem2.align() as part of the `post` list.

        Parameters
        ----------
        bam_file : Path | str
            Path to the BAM file to index.
        params : Params | None
            Parameters for indexing call.
        cfg : ExternalRunConfig | None
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
        bam_file = Path(bam_file).absolute()
        bai_file = bai_file or bam_file.with_suffix(bam_file.suffix + ".bai")
        # Build command: samtools index <bam_file>
        arguments = ["index", self.strabs(bam_file), "-o", self.strabs(bai_file)]
        return arguments, [bai_file], None, None, None

    ###########################################################################
    # Faidx (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def faidx(
        self,
        genome: Genome,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        cut2: bool = True,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """
        Run faidx on a Genome-like object and return an Element.

        Parameters
        ----------
        genome : Genome
            Genome-like object containing the reference FASTA file.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the genome name.
        outdir : Path | str | None
            Directory for the output genome file. If not provided, defaults to
            the FASTA file's parent directory.
        filename : Path | str | None
            Filename override. If not provided, defaults to
            ``tag.default_output``.
        cut2 : bool, optional
            Whether to cut the .fai output to create a .genome file
            (default: True).
            If True, a .genome file will be created with the same base name as
            the output_file, containing the first two columns of the .fai
            (chromosome name and length).
        params : Params | None
            Parameters for the faidx call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            Element representing the faidx operation.
        """
        fasta = (genome.fasta).absolute()
        default_tag = ElementTag(
            root=genome.name,
            level=1,
            stage=Stage.PREP,
            method=Method.SAMTOOLS,
            state=State.INDEX,
            omics=None,
            ext="genome",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag
        output = Path(outdir or fasta.parent) / (filename or tag.default_output)

        runner = self.faidx_fasta(
            fasta=fasta, output_file=output, cut2=cut2, params=params, cfg=cfg
        )
        artifacts = {"fai": output}
        if cut2:
            artifacts["faicut"] = output.with_suffix(".faicut")
        element = Element(
            name=f"{genome.name}.faidx",
            key=f"{self.name}_faidx_{self.version_name}_{cut2}",
            run=runner,
            tag=tag,
            determinants=self.signature_determinants(params, subroutine="faidx"),
            inputs=[fasta],
            artifacts=artifacts,
            pres=[],
        )
        return element

    @subroutine
    def faidx_fasta(
        self,
        fasta: Path | str,
        output_file: Path | str,
        cut2: bool = True,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
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
            fai_path = fasta.with_suffix(fasta.suffix + ".fai")
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
        return arguments, [output_file], None, None, {"post:": call_afterwards}

    ###########################################################################
    # Stats (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def stats(
        self,
        mapped: MappedElement,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> "Element":
        """Run samtools stats from a MappedElement.

        Creates an Element that runs samtools stats to collect comprehensive
        alignment statistics.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        outdir : Path | str | None
            Directory for the output stats file. If not provided, defaults to
            "qc/samtools/" relative to the BAM file's directory.
        filename : Path | str | None
            Filename override. If not provided, defaults to
            ``tag.default_output``.
        params : Params | None
            Parameters for the stats call.
        cfg : ExternalRunConfig | None
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
        tag = from_prior(
            mapped.tag,
            tag,
            level=1,
            stage=Stage.QC,
            method=Method.SAMTOOLS,
            state=State.RAW,
            ext="txt",  # may be tsv or json
        )
        default_outdir = input_bam.parent / "qc" / "samtools"
        output_file = Path(outdir or default_outdir) / (filename or tag.default_output)

        runner = self.stats_bam(
            input_bam=input_bam,
            output_file=output_file,
            params=params,
            cfg=cfg,
        )

        return Element(
            name=f"{input_bam.stem}_stats_{self.name}",
            key=f"{input_bam.stem}_stats_{self.version_name}",
            run=runner,
            tag=tag,
            determinants=self.signature_determinants(params, subroutine="stats"),
            inputs=[input_bam],
            artifacts={"stats": output_file},
            pres=[mapped],
        )

    @subroutine
    def stats_bam(
        self,
        input_bam: Path | str,
        output_file: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Run samtools stats on a BAM file.

        Creates a zero-argument callable that runs samtools stats to
        collect comprehensive alignment statistics.

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_file : Path | str
            Output file for stats results.
        params : Params | None
            Parameters for the stats call.
        cfg : ExternalRunConfig | None
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
        return arguments, [output_file], output_file, None, None

    ###########################################################################
    # Flagstat (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def flagstat(
        self,
        mapped: MappedElement,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> "Element":
        """Run samtools flagstat from a MappedElement.

        Creates an Element that runs samtools flagstat to collect basic
        alignment statistics.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        outdir : Path | str | None
            Directory for the output flagstat file. If not provided, defaults
            to "qc/samtools/" relative to the BAM file's directory.
        filename : Path | str | None
            Filename override. If not provided, defaults to
            ``tag.default_output``.
        params : Params | None
            Parameters for the flagstat call.
        cfg : ExternalRunConfig | None
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
        tag = from_prior(
            mapped.tag,
            tag,
            stage=Stage.QC,
            method=Method.SAMTOOLS,
            state=State.STAT,
            ext="json",
        )
        outdir = Path(outdir or input_bam.parent / "qc" / "samtools")
        filename = filename or tag.default_output
        output_file = outdir / filename
        # generate runner
        runner = self.flagstat_bam(
            input_bam=input_bam,
            output_file=output_file,
            params=params,
            cfg=cfg,
        )

        return Element(
            key=f"{tag.default_name}_{mapped.name}_flagstat_{self.version_name}",
            run=runner,
            tag=tag,
            determinants=self.signature_determinants(params, subroutine="flagstat"),
            inputs=[input_bam],
            artifacts={"flagstat": output_file},
            pres=[mapped],
        )

    @subroutine
    def flagstat_bam(
        self,
        input_bam: Path | str,
        output_file: Path | str,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Run samtools flagstat on a BAM file.

        Creates an Element that runs samtools flagstat to
        collect basic alignment statistics (mapped/unmapped, properly paired, etc.).

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_file : Path | str
            Output file for flagstat results.
        params : Params | None
            Parameters for the flagstat call.
        cfg : ExternalRunConfig | None
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
        return arguments, [output_file], output_file, None, None
