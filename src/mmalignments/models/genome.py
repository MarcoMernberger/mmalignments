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
import json
import logging
import re
import shutil
import sys
import urllib.request
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Callable, Sequence

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import Element, element
from mmalignments.models.overlay import ElementTag
from mmalignments.models.specs import Runnable
from mmalignments.models.tags import Method, Omics, Stage, State
from mmalignments.services.environment import hostname, prebuilt_path

logger = logging.getLogger(__name__)


def _default_prebuild_root() -> Path:
    """Shared prebuild root, identical to the one the legacy mbf system used."""
    root = prebuilt_path()
    if not root:
        raise RuntimeError(
            "MBF_EXTERNAL_PREBUILD_PATH is not set. Set it or pass base_dir= "
            "to EnsemblGenome."
        )
    return Path(root) / "ppg2" / hostname() / "ensembl"

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
        return _default_prebuild_root()

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

# Probed in this order when looking for the server that hosts a species.
_ENSEMBL_SERVERS = (
    "https://ftp.ensembl.org/pub/release-{rev}/",
    "https://ftp.ensemblgenomes.org/pub/release-{rev}/fungi/",
    "https://ftp.ensemblgenomes.org/pub/release-{rev}/metazoa/",
    "https://ftp.ensemblgenomes.org/pub/release-{rev}/plants/",
    "https://ftp.ensemblgenomes.org/pub/release-{rev}/protists/",
)

# Tables the legacy mbf EnsemblGenome downloaded from the core database dump.
_SQL_TABLES = (
    "gene",
    "transcript",
    "translation",
    "stable_id_event",
    "external_db",
    "object_xref",
    "xref",
    "alt_allele",
    "seq_region",
    "seq_region_attrib",
    "attrib_type",
    "external_synonym",
)


class _ElementPath(type(Path())):  # type: ignore[misc]
    """A normal :class:`~pathlib.Path` that can also be *called*.

    ``genome.fasta`` behaves like the plain path to ``genome.fasta`` (so every
    tool wrapper that does ``genome.fasta.parent`` keeps working), while
    ``genome.fasta()`` returns the Element that produces that file.
    """

    _element_factory: Callable[[], Element]

    def with_segments(self, *segments):  # derived paths (.parent, ...) are plain
        return Path(*segments)

    def __call__(self) -> Element:
        return self._element_factory()


@dataclass(frozen=True)
class _FileSpec:
    """Where a file lives in the prebuild store and how to fetch it.

    Attributes
    ----------
    group
        Sub-folder below ``<root>/<Species>/``. These are exactly the job names
        of the legacy system, e.g. ``dna`` or ``sql/core/gene``.
    filename
        Name of the file inside the group folder.
    remote_dir
        Directory on the server (relative to the release URL) to list.
    patterns
        Regular expressions (tried in order) selecting the entry in the listing.
        The first pattern with exactly one hit wins.
    remote_name
        Maps the selected listing entry to the path of the file to download.
        Defaults to the entry itself.
    gunzip
        Decompress while downloading.
    verifiable
        True if a ``.url`` file next to the file records the release it came
        from. Files computed locally can only be trusted when they were
        produced by this genome (history hash or ``by_element``).
    """

    group: str
    filename: str
    omics: Omics
    remote_dir: str = ""
    patterns: tuple[str, ...] = ()
    remote_name: Callable[[str], str] | None = None
    gunzip: bool = False
    verifiable: bool = True


def _file_specs(species: str) -> dict[str, _FileSpec]:
    low = species.lower()
    specs = {
        "server": _FileSpec("server", "url.txt", Omics.DNA),
        "fasta": _FileSpec(
            "dna",
            "genome.fasta",
            Omics.DNA,
            f"fasta/{low}/dna/",
            (
                rf"{species}\..+\.dna\.primary_assembly\.fa\.gz",
                rf"{species}\..+\.dna\.toplevel\.fa\.gz",
            ),
            gunzip=True,
        ),
        "gtf": _FileSpec(
            "gtf",
            "genes.gtf",
            Omics.RNA,
            f"gtf/{low}/",
            (rf"{species}\..+\.\d+\.gtf\.gz",),
            gunzip=True,
        ),
        "gff3": _FileSpec(
            "gff3",
            "genes.gff3",
            Omics.RNA,
            f"gff3/{low}/",
            (rf"{species}\..+\.\d+\.gff3\.gz",),
            gunzip=True,
        ),
        "cdna": _FileSpec(
            "cdna",
            "cdna.fasta",
            Omics.RNA,
            f"fasta/{low}/cdna/",
            (rf"{species}\..+\.cdna\.all\.fa\.gz",),
            gunzip=True,
        ),
        "pep": _FileSpec(
            "pep",
            "pep.fasta",
            Omics.PROTEIN,
            f"fasta/{low}/pep/",
            (rf"{species}\..+\.pep\.all\.fa\.gz",),
            gunzip=True,
        ),
        "references": _FileSpec(
            "chromosomes_and_contigs",
            "references.txt",
            Omics.DNA,
            verifiable=False,
        ),
        "sql_definitions": _FileSpec(
            "sql/core/sql_def",
            "core.sql.gz",
            Omics.DNA,
            "mysql/",
            (rf"{low}_core_.+",),
            remote_name=lambda d: f"{d}/{d}.sql.gz",
        ),
    }
    for table in _SQL_TABLES:
        specs[f"sql_{table}"] = _FileSpec(
            f"sql/core/{table}",
            f"{table}.txt.gz",
            Omics.DNA,
            "mysql/",
            (rf"{low}_core_.+",),
            remote_name=lambda d, table=table: f"{d}/{table}.txt.gz",
        )
    return specs


# --------------------------------------------------------------------------
# Locating existing files
# --------------------------------------------------------------------------


def _usable(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _release_of(path: Path) -> int | None:
    """Ensembl release a downloaded file came from, if it was recorded."""
    marker = path if path.name == "url.txt" else Path(str(path) + ".url")
    try:
        match = re.search(r"release-(\d+)", marker.read_text())
    except OSError:
        return None
    return int(match.group(1)) if match else None


def _ppg2_histories() -> list[dict[str, str]]:
    """Read shared-job histories for the running script and the legacy layout."""
    filename = "SharedMultiFileGeneratingJobs.json"
    script = Path(sys.argv[0]).name
    paths = (
        Path(".ppg/per_script") / script / "history" / filename,
        Path(".ppg/history") / filename,
    )
    histories: list[dict[str, str]] = []
    for path in paths:
        try:
            data = json.loads(path.read_text())
        except FileNotFoundError:
            continue
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            logger.warning("Cannot read ppg2 shared-job history %s: %s", path, exc)
            continue
        if not isinstance(data, dict) or not all(
            isinstance(key, str) and isinstance(value, str)
            for key, value in data.items()
        ):
            logger.warning("Invalid ppg2 shared-job history format: %s", path)
            continue
        histories.append(data)
    return histories


def _history_file(history: dict[str, str], group_dir: Path, filename: str) -> Path | None:
    digest = history.get(str(group_dir))
    if not digest or not re.fullmatch(r"[0-9a-fA-F]{128}", digest):
        return None
    path = group_dir / "by_input" / digest / filename
    return path if _usable(path) else None


def _find_existing(group_dir: Path, spec: _FileSpec, revision: int) -> Path | None:
    """Find a file for ``revision`` created earlier by the legacy system.

    A shared build's URL may name an older release: ppg2 deduplicates identical
    outputs across releases. A history mapping is therefore trusted only when
    its corresponding server job records the requested release. Otherwise
    only a candidate with a matching URL marker can be selected.
    """
    by_input = group_dir / "by_input"
    species_dir = group_dir
    for _ in Path(spec.group).parts:
        species_dir = species_dir.parent
    for history in _ppg2_histories():
        recorded = _history_file(history, group_dir, spec.filename)
        server = _history_file(history, species_dir / "server", "url.txt")
        if recorded and server and _release_of(server) == revision:
            return recorded
        if recorded and _release_of(recorded) == revision:
            return recorded
    if not spec.verifiable:
        return None
    hits = [p for p in by_input.glob(f"*/{spec.filename}") if _usable(p)]
    hits.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    for hit in hits:
        if _release_of(hit) == revision:
            return hit
    return None


# --------------------------------------------------------------------------
# Network helpers (only used when a file is really missing)
# --------------------------------------------------------------------------


def _http_get(url: str, timeout: int = 120) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as resp:  # noqa: S310
        return resp.read()


def _list_directory(url: str) -> list[str]:
    """Entry names of an Ensembl FTP/HTTP directory listing."""
    html = _http_get(url).decode("utf-8", "replace")
    names = [
        href.rstrip("/")
        for href in re.findall(r'href="([^"?#]+)"', html)
        if "://" not in href and not href.startswith(("/", ".."))
    ]
    return list(dict.fromkeys(names))


def _pick_entry(names: Sequence[str], patterns: Sequence[str], url: str) -> str:
    for pattern in patterns:
        hits = [n for n in names if re.fullmatch(pattern, n)]
        if len(hits) == 1:
            return hits[0]
    raise ValueError(
        f"Could not select a unique entry in {url} with patterns {list(patterns)}. "
        f"Entries: {list(names)[:30]}"
    )


def _download(
    url: str, dest: Path, gunzip: bool, append: Sequence[Path] = ()
) -> None:
    """Download ``url`` to ``dest`` atomically and record the source in ``.url``."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".downloading")
    logger.info("Downloading %s -> %s", url, dest)
    try:
        with urllib.request.urlopen(url, timeout=300) as resp, open(tmp, "wb") as out:  # noqa: S310
            if gunzip:
                with gzip.GzipFile(fileobj=resp) as gz:  # type: ignore[arg-type]
                    shutil.copyfileobj(gz, out)
            else:
                shutil.copyfileobj(resp, out)
            for extra in append:
                out.write(Path(extra).read_bytes())
        tmp.replace(dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    Path(str(dest) + ".url").write_text(url)


def _ensure_fai(fasta: Path) -> None:
    if not Path(str(fasta) + ".fai").exists():
        import pysam

        pysam.faidx(str(fasta))


def _read_base_url(server_file: Path) -> str:
    """Release URL stored in ``url.txt`` (old files may say ``ftp://``)."""
    url = server_file.read_text().strip().replace("ftp://", "https://", 1)
    return url if url.endswith("/") else url + "/"


class EnsemblGenome(GenomeBase):
    """Ensembl genome stored in the shared prebuild store, built from Elements.

    Every file the legacy ``mbf`` ``EnsemblGenome`` produced has an Element
    here. Each file is looked up in this order:

    1. ``<root>/<Species>/<group>/by_element/<revision>/<file>`` (written by
       this class),
    2. ``<root>/<Species>/<group>/by_input/<hash>/<file>`` (written by the
       legacy system; the running script's pypipegraph2 history and its server
       job identify the release, or the file's ``.url`` provides a match),
    3. otherwise the path under ``by_element`` is used and the Element
       downloads the file when it runs.

    An existing file without reliable release metadata cannot safely be reused.
    Pypipegraph2 may deduplicate files across releases, so an older release in
    a file's URL does not necessarily mean its history mapping is wrong.

    Parameters
    ----------
    species
        Capitalised species, e.g. ``"Homo_sapiens"``.
    revision
        Ensembl release number.
    assembly
        Assembly name (informational; defaults to a known value if any).
    prebuild_prefix
        Accepted for compatibility with the legacy ``Genome`` dataclass; unused.
    base_dir
        Store root (``.../ppg2/<host>/ensembl``). Defaults to the shared
        ``MBF_EXTERNAL_PREBUILD_PATH`` location.
    toplevel
        Prefer the ``toplevel`` FASTA over ``primary_assembly``.
    extra_fastas
        Additional FASTA files appended to a freshly downloaded genome.

    Examples
    --------
    ::

        genome = EnsemblGenome("Homo_sapiens", 116, "GRCh38", "")
        fasta = genome.fasta()      # Element producing genome.fasta
        fasta.file                  # exact Path of genome.fasta
        genome.fasta                # same path, usable by every tool wrapper
        gtf = genome.gtf()          # likewise cdna, pep, gff3, references, ...
        indexer = aligner.index(genome)
    """

    # Elements are shared between instances describing the same genome so that
    # the registry never sees two different objects with one key.
    _element_memo: dict[tuple, Element] = {}

    def __init__(
        self,
        species: str,
        revision: int,
        assembly: str | None = None,
        prebuild_prefix: str = "",
        *,
        base_dir: Path | str | None = None,
        toplevel: bool = False,
        extra_fastas: Sequence[Path | str] = (),
    ) -> None:
        if not re.match(r"^[A-Z][a-z]+_[a-z]+$", species):
            raise ValueError(
                f"Species must be formatted like 'Homo_sapiens', got {species!r}"
            )
        self.species = species
        self.revision = int(revision)
        self.toplevel = toplevel
        self.prebuild_prefix = prebuild_prefix
        self.extra_fastas = tuple(Path(p) for p in extra_fastas)
        self._base_dir = Path(base_dir) if base_dir else None
        self._assembly = assembly
        self._specs = _file_specs(species)
        self._paths: dict[str, Path] = {}

    # ------------------------------------------------------------------
    # Identity and locations
    # ------------------------------------------------------------------

    @cached_property
    def assembly(self) -> str:
        """Assembly name given at construction or a known default."""
        return self._assembly or _KNOWN_ASSEMBLIES.get(self.species, "unknown")

    @property
    def name(self) -> str:
        return f"{self.species}_{self.revision}"

    @property
    def root(self) -> Path:
        """Root of the prebuild store."""
        return self._base_dir or _default_prebuild_root()

    @property
    def base(self) -> Path:
        """Species folder, e.g. ``.../ensembl/Homo_sapiens``."""
        return self.root / self.species

    @property
    def _identity(self) -> tuple:
        return (
            self.species,
            self.revision,
            self.toplevel,
            str(self.root),
            tuple(str(p) for p in self.extra_fastas),
        )

    def _location(self, item: str) -> Path:
        """Path where ``item`` exists or will be created (cached)."""
        if item not in self._paths:
            spec = self._specs[item]
            group_dir = self.base / spec.group
            fresh = group_dir / "by_element" / str(self.revision) / spec.filename
            found = None if _usable(fresh) else _find_existing(
                group_dir, spec, self.revision
            )
            self._paths[item] = found or fresh
        return self._paths[item]

    def _with_element(self, item: str) -> Path:
        path = _ElementPath(self._location(item))
        path._element_factory = lambda: self._element_for(item)
        return path

    def _element_for(self, item: str) -> Element:
        if item == "server":
            return self.server()
        if item == "references":
            return self.references_element()
        return self._download_element(item)

    # ------------------------------------------------------------------
    # File paths; calling one returns the Element that produces it
    # ------------------------------------------------------------------

    @property
    def fasta(self) -> Path:
        """Genome FASTA. ``genome.fasta()`` is the Element creating it."""
        return self._with_element("fasta")

    @property
    def gtf(self) -> Path:
        """Gene annotation GTF. ``genome.gtf()`` is the Element creating it."""
        return self._with_element("gtf")

    @property
    def gff3(self) -> Path:
        """Gene annotation GFF3. ``genome.gff3()`` is the Element creating it."""
        return self._with_element("gff3")

    @property
    def cdna(self) -> Path:
        """Transcript sequences. ``genome.cdna()`` is the Element creating it."""
        return self._with_element("cdna")

    @property
    def pep(self) -> Path:
        """Protein sequences. ``genome.pep()`` is the Element creating it."""
        return self._with_element("pep")

    @property
    def references(self) -> Path:
        """Contig headers of the FASTA (``genome.references()`` is its Element)."""
        return self._with_element("references")

    @property
    def sql_definitions(self) -> Path:
        """Core database schema (``genome.sql_definitions()`` is its Element)."""
        return self._with_element("sql_definitions")

    @property
    def server_file(self) -> Path:
        """``url.txt`` holding the release URL (``genome.server_file()``)."""
        return self._with_element("server")

    def sql_table_file(self, table: str) -> Path:
        """Path of one SQL table dump; call it to get the Element."""
        self._check_table(table)
        path = _ElementPath(self._location(f"sql_{table}"))
        path._element_factory = lambda: self.sql_table(table)
        return path

    # ------------------------------------------------------------------
    # Elements
    # ------------------------------------------------------------------

    def _tag(self, item: str, omics: Omics, state: State = State.DOWNLOAD) -> ElementTag:
        return ElementTag(
            root=self.name,
            level=0,
            stage=Stage.PREP,
            method=Method.ENSEMBL,
            state=state,
            omics=omics,
            flag=item,
        )

    def _memoized(self, item: str, build: Callable[[], Element]) -> Element:
        memo_key = (self._identity, item)
        if memo_key not in EnsemblGenome._element_memo:
            EnsemblGenome._element_memo[memo_key] = build()
        return EnsemblGenome._element_memo[memo_key]

    @element
    def server(self) -> Element:
        """Element finding the Ensembl server hosting this species (``url.txt``)."""

        def build() -> Element:
            dest = self._location("server")
            species, revision = self.species, self.revision

            def run() -> None:
                if _usable(dest):
                    return
                for template in _ENSEMBL_SERVERS:
                    base_url = template.format(rev=revision)
                    try:
                        listing = _http_get(base_url + "fasta/").decode(
                            "utf-8", "replace"
                        )
                    except OSError:
                        continue
                    if species.lower() in listing:
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        dest.write_text(base_url)
                        return
                raise ValueError(
                    f"Could not find {species} (release {revision}) on any Ensembl server"
                )

            tag = self._tag("server", Omics.DNA)
            return Element(
                Element.generate_key(tag, "ensembl", "server"),
                Runnable(run, display=f"ensembl server {self.name}"),
                tag=tag,
                artifacts=ArtifactSet(dest, primary_name="server"),
                determinants=(f"release={revision}",),
            )

        return self._memoized("server", build)

    @element
    def _download_element(self, item: str) -> Element:
        def build() -> Element:
            spec = self._specs[item]
            dest = self._location(item)
            server = self.server()
            server_file = server.file
            patterns = spec.patterns
            if item == "fasta" and self.toplevel:
                patterns = tuple(reversed(patterns))
            extras = self.extra_fastas if item == "fasta" else ()
            pres = [server]
            if item.startswith("sql_") and item != "sql_definitions":
                pres.append(self._download_element("sql_definitions"))

            def run() -> None:
                if not _usable(dest):
                    listing_url = _read_base_url(server_file) + spec.remote_dir
                    entry = _pick_entry(_list_directory(listing_url), patterns, listing_url)
                    remote = spec.remote_name(entry) if spec.remote_name else entry
                    _download(listing_url + remote, dest, spec.gunzip, extras)
                if item == "fasta":
                    _ensure_fai(dest)

            artifacts = ArtifactSet(dest, primary_name=item)
            if item == "fasta":
                artifacts = ArtifactSet(dest, primary_name="fasta").with_extra(
                    "fai", Path(str(dest) + ".fai")
                )
            tag = self._tag(item, spec.omics)
            return Element(
                Element.generate_key(
                    tag,
                    "ensembl",
                    "download",
                    suffix=item,
                    toplevel=self.toplevel or None,
                ),
                Runnable(run, display=f"ensembl download {item} {self.name}"),
                tag=tag,
                artifacts=artifacts,
                determinants=(f"release={self.revision}", *(str(p) for p in extras)),
                inputs=extras,
                pres=tuple(pres),
            )

        return self._memoized(item, build)

    @element
    def references_element(self) -> Element:
        """Element writing the FASTA headers (``references.txt``) of the genome."""

        def build() -> Element:
            dest = self._location("references")
            fasta_element = self._element_for("fasta")
            fasta_path = fasta_element.file

            def run() -> None:
                if _usable(dest):
                    return
                with open(fasta_path) as handle:
                    keys = [ln[1:].rstrip("\n") for ln in handle if ln.startswith(">")]
                dest.parent.mkdir(parents=True, exist_ok=True)
                tmp = dest.with_name(dest.name + ".tmp")
                tmp.write_text("\n".join(keys))
                tmp.replace(dest)

            tag = self._tag("references", Omics.DNA, State.GENERATED)
            return Element(
                Element.generate_key(tag, "ensembl", "references"),
                Runnable(run, display=f"ensembl references {self.name}"),
                tag=tag,
                artifacts=ArtifactSet(dest, primary_name="references"),
                inputs=(fasta_path,),
                pres=(fasta_element,),
            )

        return self._memoized("references", build)

    @staticmethod
    def _check_table(table: str) -> None:
        if table not in _SQL_TABLES:
            raise ValueError(f"Unknown SQL table {table!r}; choose from {_SQL_TABLES}")

    def sql_table(self, table: str) -> Element:
        """Element downloading one core-database table dump (``<table>.txt.gz``).

        The core schema (:attr:`sql_definitions`) is a prerequisite because the
        dumps are headerless and the column names come from it.
        """
        self._check_table(table)
        return self._download_element(f"sql_{table}")

    def sql_tables(self) -> list[Element]:
        """Elements for all SQL tables the legacy class downloaded."""
        return [self.sql_table(t) for t in _SQL_TABLES]

    def prepare(
        self,
        *,
        gtf: bool = True,
        gff3: bool = False,
        cdna: bool = False,
        pep: bool = False,
        references: bool = False,
        sql: bool = False,
    ) -> list[Element]:
        """Return the requested Elements in dependency order, ending with faidx.

        Parameters
        ----------
        gtf, gff3, cdna, pep, references, sql
            Select which optional files to include. The genome FASTA is always
            included; ``sql`` adds the schema and all SQL tables.
        """
        elements: list[Element] = [self.fasta()]
        if gtf:
            elements.append(self.gtf())
        if gff3:
            elements.append(self.gff3())
        if cdna:
            elements.append(self.cdna())
        if pep:
            elements.append(self.pep())
        if references:
            elements.append(self.references())
        if sql:
            elements.append(self.sql_definitions())
            elements.extend(self.sql_tables())
        elements.append(self.faidx())
        return elements

    # ------------------------------------------------------------------
    # Convenience readers (valid after the matching Element ran)
    # ------------------------------------------------------------------

    def true_chromosomes(self) -> list[str]:
        """Names of real chromosomes (all contigs if none are labelled)."""
        keys = Path(self._location("references")).read_text().split("\n")
        chroms = [k for k in keys if "chromosome:" in k] or keys
        return [k.split(" ", 1)[0] for k in chroms]

    def __repr__(self) -> str:
        return (
            f"EnsemblGenome(species={self.species!r}, revision={self.revision}, "
            f"assembly={self.assembly!r}, base={self.base})"
        )
