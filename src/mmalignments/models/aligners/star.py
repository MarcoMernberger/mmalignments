"""
Module contains an aligner interface for STAR  (Spliced Transcripts Alignment to a
Reference).
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.data import Genome
from mmalignments.models.elements import (
    Element,
    MappedElement,
    NextGenSample,
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
    Omics,
    Stage,
    State,
)
from mmalignments.services.io import parents

from ..externals import External, SubroutineIn, subroutine
from ..parameters import ParamSet
from .samtools import Samtools

logger = logging.getLogger(__name__)

# Files produced inside a STAR genome directory
_GENOME_FILES = [
    "Genome",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrNameLength.txt",
    "chrName.txt",
    "chrStart.txt",
    "genomeParameters.txt",
]

# Default output files produced per alignment run (relative to outFileNamePrefix)
_ALIGN_OUTPUTS = {
    "bam": "Aligned.sortedByCoord.out.bam",
    "log_final": "Log.final.out",
    "log": "Log.out",
    "sj": "SJ.out.tab",
}

_ALIGN_OUTPUTS_UNSORTED = {
    "bam": "Aligned.out.bam",
    "log_final": "Log.final.out",
    "log": "Log.out",
    "sj": "SJ.out.tab",
}


class STAR(External):
    """STAR RNA-seq aligner interface.

    STAR (Spliced Transcripts Alignment to a Reference) is an ultrafast
    universal RNA-seq aligner. It supports two main run modes:

    * ``genomeGenerate`` — build a genome index from a FASTA + optional GTF
    * ``alignReads``     — align FASTQ reads against the index

    Additional run modes (``inputAlignmentsFromBAM``, ``liftOver``,
    ``soloCellFiltering``) are exposed via lower-level subroutines.

    Examples
    --------
    Build a genome index::

        star = STAR()
        idx = star.index(genome, gtf=Path("annotation.gtf"))

    Align paired-end reads and sort by coordinate::

        mapped = star.align(sample, idx)

    With 2-pass mapping and gene counts::

        mapped = star.align(
            sample, idx,
            par=ParType(twopassMode="Basic", quantMode="GeneCounts"),
        )

    STARsolo (single-cell)::

        mapped = star.align(
            sample, idx,
            par=ParType(
                soloType="CB_UMI_Simple",
                soloCBwhitelist="barcodes.txt",
                soloCBlen=16,
                soloUMIlen=12,
            ),
        )
    """

    def __init__(
        self,
        name: str = "STAR",
        primary_binary: str = "STAR",
        version: str | None = None,
        folder: Path | None = None,
        source: str | None = None,
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
        genome_name: str | None = None,
    ) -> None:
        """Initialise the STAR wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: ``"STAR"``).
        primary_binary : str
            Executable name/path (default: ``"STAR"``).
        version : str | None
            Version string override.  When ``None`` the version is queried
            from the binary.
        folder : Path | None
            Default output root.  Derived automatically if ``None``.
        source : str | None
            URL / citation for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | str | Path | None
            Parameter overrides.  Defaults to the bundled ``star.json``.
        genome_name : str | None
            Name of the reference genome (used to build default paths).
        """
        if source is None:
            source = "https://github.com/alexdobin/STAR"
        parameters = parameters or Path(__file__).parent / "star.json"
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            folder=folder,
            source=source,
            parameters=parameters,
        )
        if folder is None:
            folder = Path("results") / "alignments" / self.version_name
            if genome_name is not None:
                folder = folder / genome_name
        self._folder = folder

    # ------------------------------------------------------------------
    # Version detection
    # ------------------------------------------------------------------

    def get_version(self, fallback: str | None = None) -> str | None:
        """Return the STAR version string (e.g. ``"2.7.11b"``).

        Parameters
        ----------
        fallback : str | None
            Returned when the version cannot be determined.
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
            )
            output = (cp.stdout or cp.stderr or "").strip()
            m = re.search(r"(\d+\.\d+[\.\w]*)", output)
            if m:
                return m.group(1)
        except Exception:
            pass
        return fallback

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    def default_index_dir(self, genome: Genome) -> Path:
        """Return the default genome index directory for *genome*."""
        return Path(genome.base) / self.name / self.version_name / genome.name / "index"

    def default_aligned_dir(self, sample_name: str, genome_name: str) -> Path:
        """Return the default output directory for an alignment run."""
        return (
            Path("results") / "aligned" / self.version_name / genome_name / sample_name
        )

    def genome_files(self, genome_dir: Path) -> list[Path]:
        """Return the list of expected genome index files in *genome_dir*."""
        return [genome_dir / f for f in _GENOME_FILES]

    def genome_exists(self, genome_dir: Path) -> bool:
        """Return ``True`` when all expected genome index files are present."""
        return all(f.exists() for f in self.genome_files(genome_dir))

    def rg(self, sample_name: str) -> str:
        """Build a minimal read-group string for *sample_name*."""
        return f"ID:{sample_name}_rg1 SM:{sample_name} LB:{sample_name} PL:ILLUMINA"

    # ------------------------------------------------------------------
    # Genome index — high-level Element
    # ------------------------------------------------------------------

    @element
    def index(
        self,
        genome: Genome,
        *,
        gtf: Path | str | None = None,
        output_dir: Path | str | None = None,
        sjdb_overhang: int = 100,
        tag: TagType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Build a STAR genome index.

        Wraps ``--runMode genomeGenerate``.  The resulting :class:`Element`
        exposes ``element.genome_dir`` (a ``Path``) as its primary artifact
        and can be passed directly to :meth:`align`.

        Parameters
        ----------
        genome : Genome
            Reference genome (provides ``.fasta`` and ``.name``).
        gtf : Path | str | None
            Optional GTF annotation file.  When provided, splice junctions
            are included in the index.
        output_dir : Path | str | None
            Directory for genome index files.  Auto-derived when ``None``.
        sjdb_overhang : int
            ``--sjdbOverhang`` value (ideally read_length - 1).
        tag : TagType | None
            Tag override for naming / provenance.
        par : ParType | None
            Additional ``genomeGenerate`` parameters.
        cfg : CfgType | None
            Optional subprocess configuration (threads, env, …).

        Returns
        -------
        Element
            Artifacts: ``{"genome_dir": Path, "gtf": Path | None}``.
        """
        fasta_file = genome.fasta
        output_dir = Path(output_dir or self.default_index_dir(genome))
        gtf_path = Path(gtf) if gtf else None

        runner = self.generate_genome(
            fasta_file=fasta_file,
            genome_dir=output_dir,
            gtf=gtf_path,
            sjdb_overhang=sjdb_overhang,
            par=par,
            cfg=cfg,
        )

        resolved_tag = ElementTag(
            root=genome.name,
            level=1,
            omics=Omics.RNA,
            stage=Stage.PREP,
            method=Method.STAR,
            state=State.INDEX,
        ).resolve(tag)

        artifacts = ArtifactSet(output_dir, primary_name="genome_dir")
        if gtf_path:
            artifacts = artifacts.with_extra("gtf", gtf_path)

        in_paths = [fasta_file] + ([gtf_path] if gtf_path else [])
        determinants = self.signature_determinants(par, subroutine="genomeGenerate")
        key = Element.generate_key(resolved_tag, self.version_name, "genomeGenerate")

        return Element(
            key,
            runner,
            tag=resolved_tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=tuple(in_paths),
        )

    # ------------------------------------------------------------------
    # Genome index — low-level subroutine
    # ------------------------------------------------------------------

    @subroutine
    def generate_genome(
        self,
        fasta_file: Path | str,
        genome_dir: Path | str,
        *,
        gtf: Path | str | None = None,
        sjdb_overhang: int = 100,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run ``STAR --runMode genomeGenerate``.

        Parameters
        ----------
        fasta_file : Path | str
            Reference FASTA (uncompressed).
        genome_dir : Path | str
            Directory that will contain the generated index files.
        gtf : Path | str | None
            Optional GTF annotation file.
        sjdb_overhang : int
            ``--sjdbOverhang`` value.
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        SubroutineIn
        """
        fasta_file = Path(fasta_file).absolute()
        genome_dir = Path(genome_dir).absolute()
        parents(genome_dir / "_placeholder")  # ensure dir exists

        arguments = [
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(genome_dir),
            "--genomeFastaFiles",
            str(fasta_file),
            "--sjdbOverhang",
            str(sjdb_overhang),
        ]

        in_paths: list[Path] = [fasta_file]
        if gtf:
            gtf = Path(gtf).absolute()
            arguments += ["--sjdbGTFfile", str(gtf)]
            in_paths.append(gtf)

        subcommand = "genomeGenerate"
        out_paths: list[Path] = self.genome_files(genome_dir)
        return arguments, subcommand, in_paths, out_paths, None, None, None

    # ------------------------------------------------------------------
    # Align reads — high-level Element
    # ------------------------------------------------------------------

    @element
    def align(
        self,
        sample: NextGenSample,
        index: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        read_group: str | None = None,
        sorted_bam: bool = True,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> MappedElement:
        """Align FASTQ reads with STAR (``--runMode alignReads``).

        Parameters
        ----------
        sample : NextGenSample
            Element carrying R1 (and optional R2) FASTQ paths.
        index : Element
            Genome index element produced by :meth:`index`.  Must expose
            ``element.genome_dir``.
        tag : TagType | None
            Tag override.
        out : OutType | None
            Output specification override.  STAR writes all files with fixed
            names, so only ``out.folder`` is used; it is auto-derived when not
            given.
        read_group : str | None
            ``--outSAMattrRGline`` value.  Auto-generated when ``None``.
        sorted_bam : bool
            When ``True`` (default), ``--outSAMtype BAM SortedByCoordinate``
            is used.  Set to ``False`` for unsorted BAM output.
        par : ParType | None
            Additional ``alignReads`` parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        MappedElement
            Artifacts: ``{"bam": Path, "log_final": Path, "log": Path,
            "sj": Path}``.  When ``sorted_bam=True`` the BAM is coordinate-
            sorted.
        """
        genome_dir = Path(index.genome_dir).absolute()

        resolved_tag = sample.tag.bump(
            stage=Stage.ALIGN,
            method=Method.STAR,
            state=State.MAP,
            omics=Omics.RNA,
            flag=index.tag.root,
        ).resolve(tag)

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(
                folder=self.default_aligned_dir(sample.root, index.tag.root),
                ext="bam",
            ),
            out,
        )
        outdir = Path(resolved_out.folder)

        runner = self.align_fastq(
            genome_dir=genome_dir,
            fastq_r1=sample.r1,
            fastq_r2=sample.r2,
            outdir=outdir,
            read_group=sample.read_group,
            sorted_bam=sorted_bam,
            par=par,
            cfg=cfg,
        )

        output_map = _ALIGN_OUTPUTS if sorted_bam else _ALIGN_OUTPUTS_UNSORTED
        extras = {k: outdir / v for k, v in output_map.items() if k != "bam"}
        artifacts = ArtifactSet(
            outdir / output_map["bam"], primary_name="bam", **extras
        )

        determinants = self.signature_determinants(par, subroutine="alignReads")
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "alignReads",
            genome=genome_dir.name,
            rg=sample.read_group,
        )

        return MappedElement(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(sample.r1, sample.r2) if sample.r2 else (sample.r1,),
            artifacts=artifacts,
            pres=(index,),
        )

    # ------------------------------------------------------------------
    # Align reads — low-level subroutine
    # ------------------------------------------------------------------

    @subroutine
    def align_fastq(
        self,
        genome_dir: Path | str,
        fastq_r1: Path | str,
        outdir: Path | str,
        *,
        fastq_r2: Path | str | None = None,
        read_group: str | None = None,
        sorted_bam: bool = True,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Low-level STAR alignment subroutine.

        Runs ``STAR --runMode alignReads`` with the supplied inputs and
        writes all output files into *outdir*.

        Parameters
        ----------
        genome_dir : Path | str
            Path to the STAR genome index directory.
        fastq_r1 : Path | str
            Path to R1 FASTQ (or single-end FASTQ).  Can be gzip-compressed.
        outdir : Path | str
            Directory for all output files (STAR's ``--outFileNamePrefix``
            is set to ``outdir/``).
        fastq_r2 : Path | str | None
            Path to R2 FASTQ for paired-end data.
        read_group : str | None
            ``--outSAMattrRGline`` string.
        sorted_bam : bool
            ``True`` → ``BAM SortedByCoordinate``; ``False`` → ``BAM Unsorted``.
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        SubroutineIn
        """
        genome_dir = Path(genome_dir).absolute()
        fastq_r1 = Path(fastq_r1).absolute()
        outdir = Path(outdir).absolute()
        parents(outdir / "_placeholder")

        out_prefix = str(outdir) + "/"

        # Detect gzip and add readFilesCommand if not already in params
        read_files_cmd: list[str] = []
        if str(fastq_r1).endswith(".gz"):
            read_files_cmd = ["--readFilesCommand", "zcat"]

        bam_sort_type = "BAM SortedByCoordinate" if sorted_bam else "BAM Unsorted"

        arguments = [
            "--runMode",
            "alignReads",
            "--genomeDir",
            str(genome_dir),
            "--outFileNamePrefix",
            out_prefix,
            "--outSAMtype",
            bam_sort_type,
            "--readFilesIn",
            str(fastq_r1),
        ]
        if fastq_r2:
            arguments.append(str(Path(fastq_r2).absolute()))

        if read_files_cmd:
            arguments += read_files_cmd

        if read_group:
            arguments += ["--outSAMattrRGline", read_group]

        in_paths: list[Path] = [fastq_r1]
        if fastq_r2:
            in_paths.append(Path(fastq_r2).absolute())

        output_map = _ALIGN_OUTPUTS if sorted_bam else _ALIGN_OUTPUTS_UNSORTED
        out_paths = [outdir / v for v in output_map.values()]

        subcommand = "alignReads"
        return arguments, subcommand, in_paths, out_paths, None, None, None

    # ------------------------------------------------------------------
    # inputAlignmentsFromBAM — high-level Element
    # ------------------------------------------------------------------

    @element
    def wig_from_bam(
        self,
        bam: Element,
        *,
        wig_type: str = "bedGraph",
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Generate a wiggle/bedGraph from an existing sorted BAM.

        Uses ``--runMode inputAlignmentsFromBAM`` with ``--outWigType``.

        Parameters
        ----------
        bam : Element
            BAM element (must expose ``element.bam``).
        wig_type : str
            ``--outWigType`` value, e.g. ``"bedGraph"`` or
            ``"bedGraph read1_5p"``.
        tag : TagType | None
            Tag override.
        out : OutType | None
            Output specification override.  Only ``out.folder`` is used
            (STAR fixes the file names).
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        Element
            Artifacts: ``{"wig_plus": Path, "wig_minus": Path}``.
        """
        bam_file = Path(bam.bam).absolute()

        resolved_tag = bam.tag.bump(
            stage=Stage.ALIGN,
            method=Method.STAR,
            state=State.STAT,
        ).resolve(tag)

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=bam_file.parent / "wig", ext="bedGraph"), out
        )
        outdir = Path(resolved_out.folder).absolute()

        runner = self.bam_to_wig(
            bam_file=bam_file,
            outdir=outdir,
            wig_type=wig_type,
            par=par,
            cfg=cfg,
        )

        # STAR produces Stranded files by default
        artifacts = ArtifactSet(
            outdir / "Signal.Unique.str1.out.bg",
            primary_name="wig_plus",
            wig_minus=outdir / "Signal.Unique.str2.out.bg",
        )

        determinants = self.signature_determinants(par, subroutine="alignReads")
        key = Element.generate_key(
            resolved_tag, self.version_name, "inputAlignmentsFromBAM"
        )

        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(bam_file,),
            artifacts=artifacts,
            pres=(bam,),
        )

    @subroutine
    def bam_to_wig(
        self,
        bam_file: Path | str,
        outdir: Path | str,
        *,
        wig_type: str = "bedGraph",
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run ``STAR --runMode inputAlignmentsFromBAM``.

        Parameters
        ----------
        bam_file : Path | str
            Sorted, indexed BAM file.
        outdir : Path | str
            Output directory.
        wig_type : str
            ``--outWigType`` value.
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        SubroutineIn
        """
        bam_file = Path(bam_file).absolute()
        outdir = Path(outdir).absolute()
        parents(outdir / "_placeholder")

        arguments = [
            "--runMode",
            "inputAlignmentsFromBAM",
            "--inputBAMfile",
            str(bam_file),
            "--outFileNamePrefix",
            str(outdir) + "/",
            "--outWigType",
            wig_type,
        ]

        out_paths = [
            outdir / "Signal.Unique.str1.out.bg",
            outdir / "Signal.Unique.str2.out.bg",
        ]

        subcommand = "inputAlignmentsFromBAM"
        return arguments, subcommand, [bam_file], out_paths, None, None, None

    # ------------------------------------------------------------------
    # liftOver — high-level Element
    # ------------------------------------------------------------------

    @element
    def liftover(
        self,
        gtf: Path | str,
        chain_files: list[Path | str],
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Lift-over a GTF file between genome assemblies (``--runMode liftOver``).

        Parameters
        ----------
        gtf : Path | str
            GTF file to lift over.
        chain_files : list[Path | str]
            One or more chain files.
        tag : TagType | None
            Tag override.
        out : OutType | None
            Output specification override.  Only ``out.folder`` is used
            (STAR fixes the file name).
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        Element
            Artifacts: ``{"gtf": Path}`` — lifted-over GTF.
        """
        gtf = Path(gtf).absolute()
        chain_paths = [Path(c).absolute() for c in chain_files]

        resolved_tag = ElementTag(
            root=gtf.stem,
            level=1,
            omics=Omics.RNA,
            stage=Stage.PREP,
            method=Method.STAR,
            state=State.ANNOTATE,
        ).resolve(tag)

        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=gtf.parent / "liftover", ext="gtf"), out
        )
        outdir = Path(resolved_out.folder).absolute()

        runner = self.run_liftover(
            gtf=gtf,
            chain_files=chain_paths,
            outdir=outdir,
            par=par,
            cfg=cfg,
        )

        out_gtf = outdir / (gtf.stem + ".liftOver.gtf")

        in_paths = [gtf] + chain_paths
        key = Element.generate_key(resolved_tag, self.version_name, "liftOver")

        return Element(
            key,
            runner,
            tag=resolved_tag,
            inputs=tuple(in_paths),
            artifacts=ArtifactSet(out_gtf, primary_name="gtf"),
        )

    @subroutine
    def run_liftover(
        self,
        gtf: Path | str,
        chain_files: list[Path],
        outdir: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run ``STAR --runMode liftOver``.

        Parameters
        ----------
        gtf : Path | str
            Source GTF file.
        chain_files : list[Path]
            Chain files for liftover.
        outdir : Path | str
            Output directory.
        par : ParType | None
            Additional parameters.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        SubroutineIn
        """
        gtf = Path(gtf).absolute()
        outdir = Path(outdir).absolute()
        parents(outdir / "_placeholder")

        arguments = [
            "--runMode",
            "liftOver",
            "--sjdbGTFfile",
            str(gtf),
            "--genomeChainFiles",
            *[str(c) for c in chain_files],
            "--outFileNamePrefix",
            str(outdir) + "/",
        ]

        out_gtf = outdir / (gtf.stem + ".liftOver.gtf")
        in_paths = [gtf] + list(chain_files)

        subcommand = "liftOver"
        return arguments, subcommand, in_paths, [out_gtf], None, None, None

    # ------------------------------------------------------------------
    # Convenience: index + align
    # ------------------------------------------------------------------

    def alignsort(
        self,
        sample: NextGenSample,
        index: Element,
        *,
        read_group: str | None = None,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> tuple[Element, Element, Element]:
        """Align, then sort + index the BAM using Samtools.

        This is a convenience wrapper that chains :meth:`align` →
        ``Samtools.sort`` → ``Samtools.index``.

        Parameters
        ----------
        sample : NextGenSample
            Input sample.
        index : Element
            STAR genome index element.
        read_group : str | None
            Read group string.
        tag, out, par, cfg
            Overrides, optionally per step (``align``, ``sort``, ``index``)
            via ``SpecValues``.

        Returns
        -------
        tuple[Element, Element, Element]
            ``(sorted_bam, sorted_bam_index, raw_mapped)`` — first entry is
            the terminal Element of the chain.
        """
        mapped = self.align(
            sample=sample,
            index=index,
            tag=pick_spec(tag, "align"),
            out=pick_spec(out, "align"),
            read_group=read_group,
            par=pick_spec(par, "align"),
            cfg=pick_spec(cfg, "align"),
        )

        st = Samtools()
        mapped_sorted = st.sort(
            mapped=mapped,
            tag=pick_spec(tag, "sort"),
            out=pick_spec(out, "sort"),
            par=pick_spec(par, "sort"),
            cfg=pick_spec(cfg, "sort"),
        )
        mapped_sorted_index = st.index(
            mapped=mapped_sorted,
            tag=pick_spec(tag, "index"),
            out=pick_spec(out, "index"),
            par=pick_spec(par, "index"),
            cfg=pick_spec(cfg, "index"),
        )

        return mapped_sorted, mapped_sorted_index, mapped
