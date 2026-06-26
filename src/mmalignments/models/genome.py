"""Genome reference models: file-based and Ensembl-backed.

Design overview
---------------
``GenomeBase``
    Abstract base class providing the standard genome interface:
    ``.name``, ``.fasta``, ``.gtf``, ``.gff3``, ``.chromosome_lengths``,
    ``.get_sequence()``.  Elements for common preparation steps (faidx, …)
    are returned as factory methods so they plug directly into the pipeline.

``LocalGenome``
    Wraps already-present FASTA / GTF files.  No network access needed.
    Use this when you have a custom reference or a pre-downloaded genome.

``EnsemblGenome``
    Lazy Ensembl genome: files live in a shared prebuild directory
    (``/machine/…``) that is shared across containers so that each genome
    only needs to be downloaded once.  It offers Element-returning methods
    for every download / preparation step so that everything integrates with
    the pipeline's skip/cache logic.

Compatibility
-------------
Both ``LocalGenome`` and ``EnsemblGenome`` expose the same attributes as the
legacy ``Genome`` dataclass (``name``, ``fasta``, ``base``) so existing
GATK / BWAMem2 / STAR tool wrappers continue to work without modification.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
import urllib.request
from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from mmalignments.models.elements import (
    Element,
    ElementTag,
    Runnable,
    element
)
from mmalignments.models.tags import Method, Omics, Stage, State
from mmalignments.services.environment import get_variable, hostname, prebuilt_path
from mmalignments.services.io import parents
from mmalignments.services.logging import current_call_to_string

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Shared prebuild root – identical to the path that mbf uses.
# Override via ``EnsemblGenome(base_dir=…)``.
# ---------------------------------------------------------------------------
_DEFAULT_PREBUILD_ROOT = Path(prebuilt_path()) / "ppg2" / hostname() / "ensembl"

# ---------------------------------------------------------------------------
# Known Ensembl assembly names (avoids a REST call in the common case).
# Extend as needed.
# ---------------------------------------------------------------------------
_KNOWN_ASSEMBLIES: dict[str, str] = {
    "Homo_sapiens": "GRCh38",
    "Mus_musculus": "GRCm39",
    "Rattus_norvegicus": "mRatBN7.2",
    "Danio_rerio": "GRCz11",
    "Drosophila_melanogaster": "BDGP6.46",
    "Caenorhabditis_elegans": "WBcel235",
    "Saccharomyces_cerevisiae": "R64-1-1",
    "Ustilago_maydis": "KM212",
}


###############################################################################
# GenomeBase
###############################################################################

class GenomeBase(ABC):
    """Abstract base for all genome references.

    Subclasses must implement the private path resolvers; everything else
    (chromosome lengths, sequence queries, Elements) is implemented here.
    """

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------


    @property
    @abstractmethod
    def name(self) -> str:
        """Unique genome identifier, e.g. ``Homo_sapiens_113``."""

    @property
    @abstractmethod
    def base(self) -> Path:
        """Root directory where genome files are stored."""

    @property
    @abstractmethod
    def fasta(self) -> Path:
        """Path to the (uncompressed) genome FASTA."""

    @property
    def gtf(self) -> Path | None:
        """Path to the gene annotation GTF file, or *None* if not available."""
        return None

    @property
    def gff3(self) -> Path | None:
        """Path to the gene annotation GFF3 file, or *None* if not available."""
        return None

    @property
    def default_prebuild_path(self) -> Path:
        """Default path for prebuilt genome files."""
        return prebuilt_path / "ppg2" / hostname() / "ensembl"

    # ------------------------------------------------------------------
    # Derived properties (lazy, no network / subprocess required)
    # ------------------------------------------------------------------

    @property
    def fai(self) -> Path:
        """Expected path of the samtools FASTA index (.fai)."""
        return Path(str(self.fasta) + ".fai")

    @property
    def genome_file(self) -> Path:
        """Two-column BEDtools genome file (chr\\tlength)."""
        return self.fasta.with_suffix(".genome")

    @cached_property
    def chromosome_lengths(self) -> dict[str, int]:
        """Chromosome → length mapping, read lazily from the ``.fai`` index.

        Raises ``FileNotFoundError`` when the index does not exist yet.
        Call :meth:`faidx` first or run the element returned by that method.
        """
        if not self.fai.exists():
            raise FileNotFoundError(
                f".fai index not found: {self.fai}\n"
                "Run the element returned by genome.faidx() first."
            )
        lengths: dict[str, int] = {}
        with self.fai.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
        return lengths

    def get_sequence(self, chrom: str, start: int, stop: int) -> str:
        """Fetch a genomic sub-sequence (0-based, half-open).

        Requires ``pysam`` and a valid FASTA index.

        Parameters
        ----------
        chrom : str
            Chromosome / contig name.
        start : int
            0-based start position.
        stop : int
            0-based exclusive end position.

        Returns
        -------
        str
            DNA sequence (uppercase).
        """
        try:
            import pysam  # type: ignore[import]
        except ImportError as exc:
            raise ImportError("pysam is required for get_sequence()") from exc
        with pysam.FastaFile(str(self.fasta)) as fa:
            return fa.fetch(chrom, start, stop).upper()

    # ------------------------------------------------------------------
    # Compatibility with the legacy Genome dataclass
    # ------------------------------------------------------------------

    def as_data(self):
        """Return a lightweight ``data.Genome`` dataclass for tool wrappers.

        This is a bridge for code that was written against the old
        ``data.Genome`` dataclass, e.g.::

            gatk.seqdict(genome.as_data())
        """
        from mmalignments.models.data import Genome as _LegacyGenome

        # Build a minimal Genome pointing at our fasta / base
        obj = object.__new__(_LegacyGenome)
        object.__setattr__(obj, "species", getattr(self, "species", self.name))
        object.__setattr__(obj, "revision", getattr(self, "revision", 0))
        object.__setattr__(obj, "prebuild_prefix", "")
        object.__setattr__(obj, "genetic_code", None)
        # Monkey-patch the properties to point at our paths
        cls_override = type(
            "_BridgeGenome",
            (_LegacyGenome,),
            {
                "fasta": property(lambda _self, _fasta=self.fasta: _fasta),
                "base": property(lambda _self, _base=self.base: _base),
                "name": property(lambda _self, _name=self.name: _name),
            },
        )
        bridge = cls_override.__new__(cls_override)
        object.__setattr__(bridge, "species", getattr(self, "species", self.name))
        object.__setattr__(bridge, "revision", getattr(self, "revision", 0))
        object.__setattr__(bridge, "prebuild_prefix", "")
        object.__setattr__(bridge, "genetic_code", None)
        return bridge

    # ------------------------------------------------------------------
    # Pipeline Elements
    # ------------------------------------------------------------------

    def faidx(
        self,
        *,
        outdir: Path | None = None,
        tag: ElementTag | None = None,
    ) -> Element:
        """Return an Element that creates a samtools FASTA index (``.fai``).

        The Element writes ``{fasta}.fai`` and also produces a two-column
        BEDtools ``.genome`` file next to it.

        Parameters
        ----------
        outdir : Path | None
            Ignored (output is always co-located with the FASTA). Present for
            API symmetry.
        tag : ElementTag | None
            Override the auto-generated tag.
        """
        from mmalignments.models.aligners.samtools import Samtools

        st = Samtools()
        return st.faidx(self, tag=tag, outdir=outdir)

    def download_fasta(self) -> Element:  # noqa: D102
        raise NotImplementedError(
            f"{type(self).__name__} does not support on-demand FASTA download. "
            "Use EnsemblGenome instead."
        )

    def download_gtf(self) -> Element:  # noqa: D102
        raise NotImplementedError(
            f"{type(self).__name__} does not support on-demand GTF download. "
            "Use EnsemblGenome instead."
        )

    def download_gff3(self) -> Element:  # noqa: D102
        raise NotImplementedError(
            f"{type(self).__name__} does not support on-demand GFF3 download. "
            "Use EnsemblGenome instead."
        )

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(name={self.name!r}, fasta={self.fasta})"
        )


###############################################################################
# LocalGenome
###############################################################################

class LocalGenome(GenomeBase):
    """Genome reference backed by local, already-present files.

    Examples
    --------
    ::

        genome = LocalGenome(
            name="hg38_custom",
            fasta=Path("/data/references/hg38/genome.fa"),
            gtf=Path("/data/references/hg38/genes.gtf"),
        )
        faidx_el = genome.faidx()
    """

    def __init__(
        self,
        name: str,
        fasta: Path | str,
        *,
        gtf: Path | str | None = None,
        gff3: Path | str | None = None,
        base: Path | str | None = None,
    ) -> None:
        self._name = name
        self._fasta = Path(fasta)
        self._gtf = Path(gtf) if gtf else None
        self._gff3 = Path(gff3) if gff3 else None
        self._base = Path(base) if base else self._fasta.parent

    @property
    def name(self) -> str:
        return self._name

    @property
    def base(self) -> Path:
        return self._base

    @property
    def fasta(self) -> Path:
        return self._fasta

    @property
    def gtf(self) -> Path | None:
        return self._gtf

    @property
    def gff3(self) -> Path | None:
        return self._gff3


###############################################################################
# EnsemblGenome
###############################################################################

_ENSEMBL_FTP = "https://ftp.ensembl.org/pub"


class EnsemblGenome(GenomeBase):
    """Ensembl genome with lazy, on-demand downloading.

    Files are stored in a shared prebuild directory so that the genome
    is downloaded only once, even across multiple analysis containers.

    The class mirrors the interface of the legacy ``data.Genome`` dataclass
    so that existing tool wrappers (GATK, BWAMem2, STAR, …) work unchanged.

    Parameters
    ----------
    species : str
        Ensembl species string, capitalised like ``"Homo_sapiens"``.
    revision : int
        Ensembl release number, e.g. ``113``.
    base_dir : Path | None
        Root of the shared prebuild store.  Defaults to
        ``/machine/ffs/prebuild/externals/ppg2/clara/ensembl``.
    assembly : str | None
        Assembly name (e.g. ``"GRCh38"``).  Looked up automatically for
        well-known species; must be supplied for exotic genomes.
    toplevel : bool
        When *True*, download the ``toplevel`` FASTA instead of
        ``primary_assembly``.  Some species (e.g. yeast) only have toplevel.

    Examples
    --------
    ::

        genome = EnsemblGenome("Homo_sapiens", 113)

        # Elements for downloads / preparation
        fasta_el   = genome.download_fasta()
        gtf_el     = genome.download_gtf()
        faidx_el   = genome.faidx()

        # Sequence queries (after faidx is done)
        seq = genome.get_sequence("chr1", 0, 100)

        # Pass to tools that accept the legacy Genome dataclass
        bwa.align(sample, bwa.index(genome))
    """

    def __init__(
        self,
        species: str,
        revision: int,
        *,
        base_dir: Path | str | None = None,
        assembly: str | None = None,
        toplevel: bool = False,
    ) -> None:
        if not re.match(r"^[A-Z][a-z]+_[a-z]+$", species):
            raise ValueError(
                f"Species must be formatted like 'Homo_sapiens', got {species!r}"
            )
        self.species = species
        self.revision = int(revision)
        self.toplevel = toplevel

        self._base_dir = Path(base_dir) if base_dir else _DEFAULT_PREBUILD_ROOT
        self._assembly = assembly
        self._fai_cache: dict[str, int] | None = None

    # ------------------------------------------------------------------
    # Assembly name resolution
    # ------------------------------------------------------------------

    @cached_property
    def assembly(self) -> str:
        """Ensembl assembly name, resolved from known list or REST API."""
        if self._assembly:
            return self._assembly
        if self.species in _KNOWN_ASSEMBLIES:
            return _KNOWN_ASSEMBLIES[self.species]
        return self._fetch_assembly_from_rest()

    def _fetch_assembly_from_rest(self) -> str:
        """Ask the Ensembl REST API for the assembly name (requires network)."""
        url = (
            f"https://rest.ensembl.org/info/assembly/{self.species}"
            f"?content-type=application/json"
        )
        logger.info("Fetching assembly name for %s from Ensembl REST …", self.species)
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310
                import json
                data = json.loads(resp.read())
                return data["assembly_name"]
        except Exception as exc:
            raise RuntimeError(
                f"Could not determine assembly name for {self.species} from "
                f"Ensembl REST API.\n"
                f"  URL: {url}\n"
                f"  Error: {exc}\n"
                "Pass assembly= explicitly to EnsemblGenome(…)."
            ) from exc

    # ------------------------------------------------------------------
    # GenomeBase interface
    # ------------------------------------------------------------------

    @property
    def name(self) -> str:
        return f"{self.species}_{self.revision}"

    @property
    def base(self) -> Path:
        return self._base_dir / self.species / self.name

    @property
    def fasta(self) -> Path:
        return self.base / "genome.fasta"

    @property
    def gtf(self) -> Path | None:
        p = self.base / "genes.gtf"
        return p  # return path even if not yet downloaded

    @property
    def gff3(self) -> Path | None:
        p = self.base / "genes.gff3"
        return p

    # ------------------------------------------------------------------
    # FTP URL builders
    # ------------------------------------------------------------------

    def _species_lower(self) -> str:
        return self.species.lower()

    def _fasta_url(self) -> str:
        suffix = "toplevel" if self.toplevel else "primary_assembly"
        filename = f"{self.species}.{self.assembly}.dna.{suffix}.fa.gz"
        return (
            f"{_ENSEMBL_FTP}/release-{self.revision}/fasta"
            f"/{self._species_lower()}/dna/{filename}"
        )

    def _gtf_url(self) -> str:
        filename = f"{self.species}.{self.assembly}.{self.revision}.gtf.gz"
        return (
            f"{_ENSEMBL_FTP}/release-{self.revision}/gtf"
            f"/{self._species_lower()}/{filename}"
        )

    def _gff3_url(self) -> str:
        filename = f"{self.species}.{self.assembly}.{self.revision}.gff3.gz"
        return (
            f"{_ENSEMBL_FTP}/release-{self.revision}/gff3"
            f"/{self._species_lower()}/{filename}"
        )

    # ------------------------------------------------------------------
    # Download helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _download_and_gunzip(url: str, dest: Path) -> None:
        """Download *url*, gunzip, and write to *dest*."""
        parents(dest)
        logger.info("Downloading %s → %s", url, dest)
        tmp = dest.with_suffix(".downloading")
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:  # noqa: S310
                with gzip.GzipFile(fileobj=resp) as gz:  # type: ignore[arg-type]
                    with open(tmp, "wb") as out:
                        shutil.copyfileobj(gz, out)
            tmp.rename(dest)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise

    # ------------------------------------------------------------------
    # Download Elements
    # ------------------------------------------------------------------

    @element
    def download_fasta(self) -> Element:
        """Return an Element that downloads and gunzips the genome FASTA.

        The FASTA is written to ``{base}/genome.fasta``.
        The Element is skipped if the file already exists.

        Returns
        -------
        Element
        """
        dest = self.fasta
        print(self.fasta)
        raise ValueError(f"Debug: {self.fasta}")
        url = self._fasta_url()
        tag = ElementTag(
            root=self.name,
            level=0,
            stage=Stage.PREP,
            method=Method.ENSEMBL,
            state=State.DOWNLOAD,
            omics=Omics.DNA,
            ext="fasta",
        )
        key = f"{self.name}_download_fasta"

        def _run() -> None:
            self._download_and_gunzip(url, dest)

        runner = Runnable(
            _run,
            cmd=[f"download {url} → {dest}"],
            display="download_fasta",
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"fasta": dest},
            determinants=(url,),
            inputs=(),
            name=f"{self.name} download genome FASTA",
        )

    @element
    def download_gtf(self) -> Element:
        """Return an Element that downloads and gunzips the gene annotation GTF.

        The GTF is written to ``{base}/genes.gtf``.

        Returns
        -------
        Element
        """
        dest = self.gtf
        print(self.fasta)
        raise ValueError(f"Debug: {self.fasta}")
        url = self._gtf_url()
        tag = ElementTag(
            root=self.name,
            level=0,
            stage=Stage.PREP,
            method=Method.ENSEMBL,
            state=State.DOWNLOAD,
            omics=Omics.RNA,
            ext="gtf",
        )
        key = f"{self.name}_download_gtf"

        def _run() -> None:
            self._download_and_gunzip(url, dest)

        runner = Runnable(
            _run,
            cmd=[f"download {url} → {dest}"],
            display="download_gtf",
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"gtf": dest},
            determinants=(url,),
            inputs=(),
            name=f"{self.name} download GTF",
        )

    @element
    def download_gff3(self) -> Element:
        """Return an Element that downloads and gunzips the GFF3 annotation.

        The GFF3 is written to ``{base}/genes.gff3``.

        Returns
        -------
        Element
        """
        dest = self.gff3
        url = self._gff3_url()
        tag = ElementTag(
            root=self.name,
            level=0,
            stage=Stage.PREP,
            method=Method.ENSEMBL,
            state=State.DOWNLOAD,
            omics=Omics.RNA,
            ext="gff3",
        )
        key = f"{self.name}_download_gff3"

        def _run() -> None:
            self._download_and_gunzip(url, dest)

        runner = Runnable(
            _run,
            cmd=[f"download {url} → {dest}"],
            display="download_gff3",
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"gff3": dest},
            determinants=(url,),
            inputs=(),
            name=f"{self.name} download GFF3",
        )

    
    def prepare(
        self,
        *,
        gtf: bool = True,
        gff3: bool = False,
    ) -> list[Element]:
        """Convenience: return all download + faidx Elements in order.

        Parameters
        ----------
        gtf : bool
            Include GTF download (default: True).
        gff3 : bool
            Include GFF3 download (default: False).

        Returns
        -------
        list[Element]
            ``[download_fasta, download_gtf?, download_gff3?, faidx]``
        """
        fasta_el = self.download_fasta()
        elements: list[Element] = [fasta_el]
        if gtf:
            elements.append(self.download_gtf())
        if gff3:
            elements.append(self.download_gff3())
        elements.append(self.faidx())
        return elements

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        asm = self._assembly or _KNOWN_ASSEMBLIES.get(self.species, "?")
        return (
            f"EnsemblGenome(species={self.species!r}, revision={self.revision}, "
            f"assembly={asm!r}, base={self.base})"
        )
