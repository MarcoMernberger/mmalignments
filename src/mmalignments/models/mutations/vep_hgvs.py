"""Offline VEP annotation with reference-checked genomic HGVS."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import tarfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, BinaryIO
from urllib.request import urlopen

from mmalignments.models.artifacts import ArtifactSet, TableArtifact
from mmalignments.models.elements import Element, FileSource, TableElement, element
from mmalignments.models.externals import External, Runnable, SubroutineIn, subroutine
from mmalignments.models.overlay import (
    CfgType,
    Ecfg,
    ElementTag,
    OutSpec,
    OutType,
    ParType,
    SpecValues,
    TagType,
    pick_spec,
)
from mmalignments.models.parameters import ParamSet, ParamSpec, render_flag
from mmalignments.models.tags import Method, Omics, Stage, State
from mmalignments.services.io import ensure, parents

if TYPE_CHECKING:
    from pandas import DataFrame

    from mmalignments.models.genome import EnsemblGenome

_ENST = re.compile(r"ENST[0-9]+(?:\.[0-9]+)?")
_BASES = re.compile(r"[ACGT]+")
_CONTIG = re.compile(r"[A-Za-z0-9_.][A-Za-z0-9_.-]*")
_REFSEQ = re.compile(r"N[CGTWZ]_[A-Z0-9]+\.[0-9]+")
_COLUMNS = (
    "vep_row",
    "vep_variant_id",
    "vep_alt",
    "vep_transcript_id",
    "vep_consequence",
    "hgvs_g",
    "hgvs_c",
    "hgvs_p",
    "spdi",
)


class VEPAnnotationError(ValueError):
    """Invalid input, reference, or VEP configuration or output."""


class _Reference:
    """Read an uncompressed, faidx-indexed FASTA."""

    def __init__(self, stream: BinaryIO, index: Path) -> None:
        self.stream = stream
        self.contigs: dict[str, tuple[int, int, int, int]] = {}
        self.block: tuple[str, int, str] = ("", 0, "")
        for line in index.read_text().splitlines():
            fields = line.split("\t")
            try:
                name = fields[0]
                length, offset, bases, width = map(int, fields[1:5])
                if (
                    not _CONTIG.fullmatch(name)
                    or name in self.contigs
                    or length < 1
                    or offset < 1
                    or not 0 < bases <= width
                ):
                    raise ValueError
            except ValueError as exc:
                raise VEPAnnotationError(
                    f"Invalid FASTA index entry in {index}: {line!r}"
                ) from exc
            self.contigs[name] = length, offset, bases, width
        if not self.contigs:
            raise VEPAnnotationError(f"Empty FASTA index: {index}")

    def length(self, chrom: str) -> int:
        """Return the contig length."""
        return self.contigs[chrom][0]

    def resolve(self, chrom: str) -> str:
        """Resolve exact names and conventional chr/MT aliases."""
        if chrom in self.contigs:
            return chrom
        bare = chrom.removeprefix("chr")
        aliases = (
            {"M", "MT", "chrM", "chrMT"}
            if bare in {"M", "MT"}
            else {bare, f"chr{bare}"}
        )
        matches = aliases & self.contigs.keys()
        if len(matches) != 1:
            raise VEPAnnotationError(
                f"Chromosome {chrom!r} is missing or ambiguous in the FASTA"
            )
        return matches.pop()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch zero-based half-open reference sequence."""
        length, offset, bases, width = self.contigs[chrom]
        if not 0 <= start <= end <= length:
            raise VEPAnnotationError(
                f"Interval {chrom}:{start + 1}-{end} is outside the FASTA"
            )
        if start == end:
            return ""
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases
        self.stream.seek(first)
        seq = self.stream.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"")
        if len(seq) != end - start or not re.fullmatch(rb"[A-Za-z]+", seq):
            raise VEPAnnotationError(
                f"FASTA and index disagree at {chrom}:{start + 1}-{end}"
            )
        return seq.decode("ascii").upper()

    def base(self, chrom: str, pos: int) -> str:
        """Read one base, reusing a block for repeat traversal."""
        start = pos // 4096 * 4096
        if self.block[:2] != (chrom, start):
            self.block = (
                chrom,
                start,
                self.fetch(chrom, start, min(start + 4096, self.length(chrom))),
            )
        return self.block[2][pos - start]


@dataclass(frozen=True)
class _Edit:
    start: int
    ref: str
    alt: str

    def shifted(self, fasta: _Reference, chrom: str, direction: int) -> _Edit:
        """Shift a pure indel and rotate its allele to preserve the haplotype."""
        if self.ref and self.alt:
            return self
        motif = self.ref or self.alt
        start = self.start
        length = fasta.length(chrom)
        while True:
            rotation = (start - self.start) % len(motif)
            if direction == 1:
                adjacent = start + len(self.ref)
                if adjacent == length or fasta.base(chrom, adjacent) != motif[rotation]:
                    break
            elif start == 0 or fasta.base(chrom, start - 1) != motif[rotation - 1]:
                break
            start += direction
        rotation = (start - self.start) % len(motif)
        motif = motif[rotation:] + motif[:rotation]
        return _Edit(start, motif if self.ref else "", motif if self.alt else "")

    def vcf(self, fasta: _Reference, chrom: str) -> tuple[int, str, str]:
        """Pad an edit as a minimal VCF record."""
        if self.ref and self.alt:
            return self.start + 1, self.ref, self.alt
        if self.start:
            anchor = fasta.base(chrom, self.start - 1)
            return self.start, anchor + self.ref, anchor + self.alt
        if len(self.ref) == fasta.length(chrom):
            raise VEPAnnotationError(
                f"Cannot represent a whole-contig deletion in VCF: {chrom}"
            )
        anchor = fasta.base(chrom, len(self.ref))
        return 1, self.ref + anchor, self.alt + anchor

    def hgvs(self, fasta: _Reference, chrom: str) -> str | None:
        """Format a 3′-shifted edit on a linear genomic reference."""
        start = self.start
        end = start + len(self.ref)
        interval = str(end) if start + 1 == end else f"{start + 1}_{end}"
        mitochondrial = chrom in {"M", "MT", "chrM", "chrMT", "NC_012920.1"}
        if len(self.ref) == len(self.alt) == 1:
            prefix = "m" if mitochondrial else "g"
            return f"{prefix}.{interval}{self.ref}>{self.alt}"
        if mitochondrial:
            return None
        if not self.alt:
            return f"g.{interval}del"
        if not self.ref:
            dup_start = start - len(self.alt)
            if dup_start >= 0 and fasta.fetch(chrom, dup_start, start) == self.alt:
                interval = (
                    str(start) if len(self.alt) == 1 else f"{dup_start + 1}_{start}"
                )
                return f"g.{interval}dup"
            if 0 < start < fasta.length(chrom):
                return f"g.{start}_{start + 1}ins{self.alt}"
            return None
        if self.ref.translate(str.maketrans("ACGT", "TGCA"))[::-1] == self.alt:
            return f"g.{interval}inv"
        return f"g.{interval}delins{self.alt}"


@dataclass(frozen=True)
class _Variant:
    row: int
    id: str
    chrom: str
    vcf: tuple[int, str, str]
    alt: str
    transcript: str | None
    genomic: str | None


def _trim(start: int, ref: str, alt: str) -> _Edit:
    """Remove common context while keeping the variant in one edit."""
    end = 0
    limit = min(len(ref), len(alt))
    while end < limit and ref[-end - 1] == alt[-end - 1]:
        end += 1
    if end:
        ref, alt = ref[:-end], alt[:-end]
    begin = 0
    limit = min(len(ref), len(alt))
    while begin < limit and ref[begin] == alt[begin]:
        begin += 1
    return _Edit(start + begin, ref[begin:], alt[begin:])


def _variants(frame: DataFrame, fasta: _Reference) -> list[_Variant]:
    """Validate input and split alternatives into stable biallelic records."""
    from pandas import isna

    if not frame.columns.is_unique or set(_COLUMNS) & set(frame.columns):
        raise VEPAnnotationError(
            "Input columns must be unique and cannot use annotation column names"
        )
    chrom_col = "chromosome" if "chromosome" in frame else "chrom"
    pos_col = "position" if "position" in frame else "pos"
    required = [chrom_col, pos_col, "ref", "alt"]
    if not set(required) <= set(frame.columns):
        raise VEPAnnotationError(
            "Input requires chromosome/chrom, position/pos, ref and alt"
        )
    has_transcript = "transcript_id" in frame
    columns = required + (["transcript_id"] if has_transcript else [])
    result = []
    for row, values in enumerate(frame[columns].itertuples(index=False, name=None)):
        chrom, pos, ref, alt = map(str, values[:4])
        ref = ref.upper()
        alts = alt.upper().split(",")
        if (
            not _CONTIG.fullmatch(chrom)
            or not re.fullmatch(r"[0-9]+", pos)
            or int(pos) < 1
            or not _BASES.fullmatch(ref)
            or any(not _BASES.fullmatch(a) or a == ref for a in alts)
        ):
            raise VEPAnnotationError(
                f"Invalid chromosome, 1-based position, REF or ALT in input row {row}"
            )
        transcript = None
        if has_transcript and not isna(values[4]) and values[4] != "":
            transcript = str(values[4])
            if not _ENST.fullmatch(transcript):
                raise VEPAnnotationError(
                    f"Invalid ENST transcript in input row {row}: {transcript}"
                )
        chrom = fasta.resolve(chrom)
        start = int(pos) - 1
        actual = fasta.fetch(chrom, start, start + len(ref))
        if actual != ref:
            raise VEPAnnotationError(
                f"Reference mismatch in input row {row} at {chrom}:{pos}: REF {ref}, FASTA {actual}"
            )
        for allele, alt in enumerate(alts):
            edit = _trim(start, ref, alt)
            left = edit.shifted(fasta, chrom, -1)
            right = left.shifted(fasta, chrom, 1)
            result.append(
                _Variant(
                    row,
                    f"row{row}_alt{allele}",
                    chrom,
                    left.vcf(fasta, chrom),
                    alt,
                    transcript,
                    right.hgvs(fasta, chrom),
                )
            )
    return result


def _write_vcf(path: Path, variants: list[_Variant], fasta: _Reference) -> None:
    """Write sorted, biallelic VCF records with the original row/allele IDs."""
    parents(path)
    with path.open("w") as stream:
        stream.write("##fileformat=VCFv4.2\n")
        for chrom in sorted({v.chrom for v in variants}):
            stream.write(f"##contig=<ID={chrom},length={fasta.length(chrom)}>\n")
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
        for variant in sorted(variants, key=lambda v: (v.chrom, v.vcf[0], v.id)):
            pos, ref, alt = variant.vcf
            writer.writerow((variant.chrom, pos, variant.id, ref, alt, ".", ".", "."))


def _read_vep(path: Path, ids: set[str]) -> dict[str, list[dict[str, str]]]:
    """Read VEP tab output, including releases/configurations using Extra."""
    with path.open() as stream:
        for line in stream:
            if line.startswith("#Uploaded_variation\t"):
                header = line[1:].rstrip("\r\n").split("\t")
                break
        else:
            raise VEPAnnotationError("Missing VEP tabular header")
        if len(set(header)) != len(header) or not {
            "Feature",
            "Feature_type",
            "Consequence",
        } <= set(header):
            raise VEPAnnotationError("Invalid VEP tabular header")
        result: dict[str, list[dict[str, str]]] = defaultdict(list)
        for fields in csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE):
            if not fields:
                continue
            if len(fields) != len(header) or fields[0] not in ids:
                raise VEPAnnotationError(
                    f"Malformed VEP row or unknown variant ID: {fields[0]!r}"
                )
            record = dict(zip(header, fields, strict=True))
            for field in record.get("Extra", "").split(";"):
                key, sep, value = field.partition("=")
                if sep:
                    record.setdefault(key, value)
            result[fields[0]].append(record)
        return result


def _value(record: dict[str, str], key: str) -> str | None:
    """Convert VEP placeholders to null."""
    value = record.get(key)
    return value if value not in {None, "", "-", "."} else None


def _select(records: list[dict[str, str]], policy: str) -> list[dict[str, str]]:
    """Apply transcript policy, preferring MANE Select when available."""
    transcripts = [r for r in records if r["Feature_type"] == "Transcript"]
    if policy == "all":
        return transcripts
    if policy == "mane":
        return [r for r in transcripts if _value(r, "MANE_SELECT")] or transcripts
    if policy == "canonical":
        return [r for r in transcripts if r.get("CANONICAL") == "YES"]
    if "." in policy:
        return [r for r in transcripts if r["Feature"] == policy]
    return [r for r in transcripts if r["Feature"].partition(".")[0] == policy]


def _annotate(
    frame: DataFrame,
    variants: list[_Variant],
    annotations: dict[str, list[dict[str, str]]],
    policy: str,
) -> DataFrame:
    """Expand rows by ALT and selected transcript while retaining the input index."""
    accessions: dict[str, set[str]] = defaultdict(set)
    for variant in variants:
        for record in annotations.get(variant.id, []):
            value = _value(record, "HGVSg")
            accession = value.partition(":")[0] if value else ""
            if _REFSEQ.fullmatch(accession):
                accessions[variant.chrom].add(accession)
    names = {
        chrom: next(iter(ids)) for chrom, ids in accessions.items() if len(ids) == 1
    }
    rows = []
    fields = []
    for variant in variants:
        records = annotations.get(variant.id, [])
        accession = names.get(variant.chrom, variant.chrom)
        genomic = f"{accession}:{variant.genomic}" if variant.genomic else None
        spdi = next(
            (value for record in records if (value := _value(record, "SPDI"))), None
        )
        for record in _select(records, variant.transcript or policy) or [{}]:
            rows.append(variant.row)
            fields.append(
                (
                    variant.row,
                    variant.id,
                    variant.alt,
                    _value(record, "Feature"),
                    _value(record, "Consequence"),
                    genomic,
                    _value(record, "HGVSc"),
                    _value(record, "HGVSp"),
                    _value(record, "SPDI") or spdi,
                )
            )
    result = frame.iloc[rows].copy()
    for index, column in enumerate(_COLUMNS):
        result[column] = [row[index] for row in fields]
    return result


def _probe(tool: External, args: list[str]) -> subprocess.CompletedProcess[str]:
    """Probe VEP without relying on help text to enumerate supported flags."""
    try:
        return subprocess.run(
            tool.build_cmd(args), capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise VEPAnnotationError(
            f"Cannot start VEP executable {tool.primary_binary}: {exc}"
        ) from exc


def _table(source: Element | FileSource) -> TableArtifact:
    """Resolve a source's primary table through the artifact interface."""
    table = source.primary
    if isinstance(table, TableArtifact):
        return table
    if isinstance(table, Path) and table.suffix in {".parquet", ".tsv", ".csv"}:
        return TableArtifact(table)
    raise VEPAnnotationError("VEPAnnotator requires a primary table artifact")


def _cache_metadata(info: Path, assembly: str) -> str:
    """Read and validate VEP cache metadata."""
    if not info.is_file():
        raise VEPAnnotationError(f"Missing VEP cache metadata: {info}")
    text = info.read_text()
    match = re.search(r"^assembly\t([^\r\n]+)", text, re.M)
    if not text.strip() or (match and match[1].partition(".")[0] != assembly):
        raise VEPAnnotationError(f"Invalid or incompatible VEP cache metadata: {info}")
    return text


class VEPAnnotator(External):
    """Build an offline VEP DAG for short human variants.

    The workflow prepares reference-checked input, runs VEP, and reports its
    transcript annotations alongside independently normalized genomic HGVS.
    """

    def __init__(
        self,
        *,
        genome: EnsemblGenome,
        cache_dir: Path | str = Path("cache/vep"),
        executable: Path | str = "vep",
        predict_p: bool = False,
        folder: Path | None = None,
    ) -> None:
        """Configure the external VEP wrapper.

        Parameters
        ----------
        genome : EnsemblGenome
            Human Ensembl reference with an uncompressed FASTA and adjacent .fai.
        cache_dir : Path or str
            Root of the local human Ensembl cache for the genome revision.
        executable : Path or str
            Local VEP executable.
        predict_p : bool
            Request VEP's predicted protein HGVS notation.
        folder : Path or None
            Default output directory, falling back to cache/vep.
        """
        if genome.species.lower() != "homo_sapiens":
            raise VEPAnnotationError("VEPAnnotator requires a human Ensembl genome")
        assembly = genome.assembly
        if not re.fullmatch(r"GRCh[0-9]+", assembly):
            raise VEPAnnotationError(
                f"Expected a human GRCh assembly, got {assembly!r}"
            )
        self.reference = genome.fasta.resolve()
        self.index = genome.fai.resolve()
        self.cache = Path(cache_dir).resolve()
        self.release = str(genome.revision)
        if not self.reference.is_file() or not self.index.is_file():
            raise VEPAnnotationError(
                f"Provide a local FASTA and its .fai index: {self.reference}, {self.index}"
            )
        binary = shutil.which(str(executable))
        if binary is None:
            raise VEPAnnotationError(f"VEP executable not found: {executable}")
        self.binary = Path(binary).absolute()
        super().__init__(
            name="vep",
            primary_binary=str(self.binary),
            version="probe",
            folder=folder or Path("cache/vep"),
            parameters=ParamSet(
                {"refseq": ParamSpec("refseq", "--refseq", bool, render=render_flag)}
            ),
        )
        self.assembly = assembly
        help_result = _probe(self, ["--help"])
        version_data = help_result.stdout + help_result.stderr
        match = re.search(r"ensembl-vep\s*:\s*(\d+)([^\s]*)", version_data)
        if help_result.returncode or not match:
            raise VEPAnnotationError(
                f"Cannot determine Ensembl VEP version: {version_data.strip()}"
            )
        self._version = match[1] + match[2]
        self.component_versions = (
            version_data.partition("Versions:")[2].partition("\n\n")[0].strip()
            or self.version
        )
        flags = [
            "--hgvs",
            "--hgvsg",
            "--shift_hgvs",
            "1",
            "--transcript_version",
            "--canonical",
            "--check_ref",
            "--no_escape",
        ]
        if predict_p:
            flags.append("--hgvsp_use_prediction")
        if assembly == "GRCh38":
            flags.append("--mane")
        for optional in ("--hgvsg_use_accession", "--spdi"):
            probe = _probe(self, ["--help", optional])
            if probe.returncode == 0:
                flags.append(optional)
            elif not re.search(
                rf"unknown option:\s*{optional[2:]}\b", probe.stderr, re.I
            ):
                raise VEPAnnotationError(
                    f"VEP option probe failed: {probe.stderr.strip()}"
                )
        probe = _probe(self, ["--help", *flags])
        if probe.returncode:
            raise VEPAnnotationError(
                f"VEP {self.version} lacks required HGVS options: {probe.stderr.strip()}"
            )
        self.flags = tuple(flags)

    @element(auto_outputs=False)
    def download_cache(
        self,
        *,
        refseq: bool = False,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create the Ensembl VEP cache installation element.

        Parameters
        ----------
        refseq : bool
            Install the RefSeq transcript cache when enabled.
        tag : TagType or None
            Optional tag override shared by callers of this cache.
        out : OutType or None
            Output spec for the metadata file at VEP's required cache path.
        par : ParType or None
            Reserved for API symmetry.
        cfg : CfgType or None
            Reserved for API symmetry.

        Returns
        -------
        Element
            Deferred installation with cache metadata as its primary artifact.
        """
        species = "homo_sapiens_refseq" if refseq else "homo_sapiens"
        cache_info = (
            self.cache / species / f"{self.release}_{self.assembly}" / "info.txt"
        )
        tag = ElementTag(
            root=f"{species}_{self.release}_{self.assembly}",
            level=0,
            stage=Stage.PREP,
            method=Method.ENSEMBL,
            state=State.DOWNLOAD,
            omics=Omics.DNA,
            flag="vep-cache",
        ).resolve(tag)
        spec = OutSpec("info", cache_info.parent, ext="txt").patch(out)
        info = spec.file.resolve()
        if info != cache_info:
            raise VEPAnnotationError(f"VEP requires cache metadata at {cache_info}")
        branch = "grch37/" if self.assembly == "GRCh37" else ""
        url = (
            f"https://ftp.ensembl.org/pub/{branch}release-{self.release}"
            f"/variation/indexed_vep_cache/{species}_vep_{self.release}_{self.assembly}.tar.gz"
        )
        key, _ = self.build_element_name(tag, "download_cache")
        return Element(
            key,
            self._download_cache(url, self.cache, info, self.assembly),
            tag,
            ArtifactSet(info, primary_name="cache_info"),
            determinants=(url, str(info)),
        )

    @staticmethod
    def _download_cache(url: str, cache: Path, info: Path, assembly: str) -> Runnable:
        """Build a runner that streams and installs an indexed cache archive."""

        def run() -> None:
            if info.is_file():
                _cache_metadata(info, assembly)
                return
            ensure(cache)
            with TemporaryDirectory(prefix=".vep-cache-", dir=cache) as folder:
                work = Path(folder)
                with urlopen(url, timeout=300) as response:  # noqa: S310
                    with tarfile.open(fileobj=response, mode="r|gz") as archive:
                        archive.extractall(work, filter="data")
                extracted_info = work / info.relative_to(cache)
                _cache_metadata(extracted_info, assembly)
                ensure(info.parent.parent)
                try:
                    extracted_info.parent.rename(info.parent)
                except OSError:
                    if not info.is_file():
                        raise
                    _cache_metadata(info, assembly)

        return Runnable(
            run, cmd=["download", url, str(info)], display="Download VEP cache"
        )

    @element
    def prepare(
        self,
        source: Element | FileSource,
        vcf: Path | None = None,
        manifest: Path | None = None,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create the biallelic VCF and row mapping element.

        Parameters
        ----------
        source : Element or FileSource
            Variant table producer or file source.
        vcf : Path or None
            Destination for the reference-checked VCF.
        manifest : Path or None
            Destination for the variant mapping and reference checksum.

        tag : TagType or None
            Optional override for the derived tag.
        out : OutType or None
            Optional primary output specification.
        par : ParType or None
            Reserved for API symmetry.
        cfg : CfgType or None
            Reserved for API symmetry.

        Returns
        -------
        Element
            Deferred input preparation.
        """
        table = _table(source)
        tag = source.tag.bump(flag="vep-input").resolve(tag)
        if vcf is None or out is not None:
            vcf = OutSpec(tag.default_name, self.folder, ext="vcf").patch(out).file
        vcf = vcf.resolve()
        manifest = (manifest or vcf.with_suffix(".json")).resolve()
        key, _ = self.build_element_name(tag, "prepare")
        return Element(
            key,
            self._prepare_input(
                table.resolve(), self.reference, self.index, vcf, manifest
            ),
            tag,
            ArtifactSet(vcf, primary_name="vcf", manifest=manifest),
            pres=(source,),
            inputs=(self.reference, self.index),
            determinants=(str(vcf), str(manifest)),
        )

    @staticmethod
    def _prepare_input(
        table: Path, reference: Path, index: Path, vcf: Path, manifest: Path
    ) -> Runnable:
        """Build a runner that writes the VCF and row mapping."""

        def run() -> None:
            frame = TableArtifact(table).read(drop_unnamed_columns=False)
            parents(vcf, manifest)
            with TemporaryDirectory(prefix=".vep-input-", dir=vcf.parent) as folder:
                work = Path(folder)
                with reference.open("rb") as stream:
                    if stream.peek(2)[:2] == b"\x1f\x8b":
                        raise VEPAnnotationError("Reference FASTA must be uncompressed")
                    ref = _Reference(stream, index)
                    variants = _variants(frame, ref)
                    _write_vcf(work / "input.vcf", variants, ref)
                    stream.seek(0)
                    checksum = hashlib.file_digest(stream, "sha256").hexdigest()
                mapping = {
                    "reference_sha256": checksum,
                    "variants": [asdict(v) for v in variants],
                }
                (work / "manifest.json").write_text(json.dumps(mapping) + "\n")
                (work / "input.vcf").replace(vcf)
                (work / "manifest.json").replace(manifest)

        return Runnable(run, display=f"Prepare VEP input {vcf}")

    @element
    def process(
        self,
        prepared: Element,
        output: Path | None = None,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create the offline VEP execution element.

        Parameters
        ----------
        prepared : Element
            The VCF and mapping produced by prepare.
        output : Path or None
            Destination for the VEP tabular output.

        tag : TagType or None
            Optional override for the derived tag.
        out : OutType or None
            Optional primary output specification.
        par : ParType or None
            VEP options. Set ``refseq=True`` to use the RefSeq transcript cache.
        cfg : CfgType or None
            VEP subprocess configuration.

        Returns
        -------
        Element
            Deferred VEP invocation.
        """
        vcf = prepared.artifacts["vcf"]
        manifest = prepared.artifacts["manifest"]
        tag = prepared.tag.bump(flag="vep-annotation").resolve(tag)
        if output is None or out is not None:
            output = OutSpec(tag.default_name, self.folder, ext="tsv").patch(out).file
        output = output.resolve()
        key, _ = self.build_element_name(tag, "process")
        cache = self.download_cache(refseq=bool(par and par.get("refseq")))
        cache_info = cache.file
        return Element(
            key,
            self._process_vcf(
                vcf, manifest, output, cache_info, self.assembly, par, cfg
            ),
            tag,
            ArtifactSet(output, primary_name="vep_tsv"),
            pres=(prepared, cache),
            inputs=(cache_info, self.binary),
            determinants=(
                self.assembly,
                *self.flags,
                str(output),
                *self.signature_determinants(par),
                repr(cfg),
            ),
        )

    def _process_vcf(
        self,
        vcf: Path,
        manifest: Path,
        output: Path,
        cache_info: Path,
        assembly: str,
        par: ParType | None,
        cfg: CfgType | None,
    ) -> Runnable:
        """Build a runner for offline VEP and validate its reference checks."""
        run_cfg = Ecfg().resolve(cfg) if cfg is not None else None

        def run() -> None:
            _cache_metadata(cache_info, assembly)
            parents(output)
            with TemporaryDirectory(prefix=".vep-run-", dir=output.parent) as folder:
                work = Path(folder)
                tab = work / "output.tsv"
                skipped = work / "skipped.txt"
                if json.loads(manifest.read_text())["variants"]:
                    self._run_vep(vcf, tab, skipped, par=par, cfg=run_cfg)()
                    if skipped.exists():
                        for line in skipped.read_text().splitlines():
                            if "does not match Ensembl reference allele" in line:
                                raise VEPAnnotationError(
                                    f"VEP reference mismatch: {line}"
                                )
                else:
                    tab.write_text(
                        "#Uploaded_variation\tFeature\tFeature_type\tConsequence\n"
                    )
                tab.replace(output)

        return Runnable(run, display=f"Run offline VEP on {vcf}")

    @subroutine
    def _run_vep(
        self,
        vcf: Path,
        tab: Path,
        skipped: Path,
        *,
        par: ParType | None = None,
        cfg: Ecfg | None = None,
    ) -> SubroutineIn:
        """Build the offline VEP command and its file dependencies."""
        arguments = [
            "--offline",
            "--cache",
            "--species",
            "homo_sapiens",
            "--dir",
            str(tab.parent),
            "--dir_cache",
            str(self.cache),
            "--cache_version",
            self.release,
            "--assembly",
            self.assembly,
            "--fasta",
            str(self.reference),
            "--format",
            "vcf",
            "--input_file",
            str(vcf),
            "--output_file",
            str(tab),
            "--skipped_variants_file",
            str(skipped),
            "--tab",
            "--no_stats",
            *self.flags,
        ]
        return arguments, None, (vcf, self.reference), (tab, skipped), None, None, None

    @element
    def report(
        self,
        source: Element | FileSource,
        prepared: Element,
        processed: Element,
        output: Path | None = None,
        metadata: Path | None = None,
        *,
        transcript_policy: str = "all",
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> TableElement:
        """Create the annotated table and provenance element.

        Parameters
        ----------
        source : Element or FileSource
            Original table containing the columns and index to preserve.
        prepared : Element
            VCF preparation element that records the row and ALT mapping.
        processed : Element
            Element that produces raw VEP tabular output.
        output : Path or None
            Destination for the annotated Parquet table.
        metadata : Path or None
            Destination for the JSON provenance record.

        transcript_policy : str
            Transcript selection policy for this report.
        tag : TagType or None
            Optional override for the derived tag.
        out : OutType or None
            Optional primary output specification.
        par : ParType or None
            Reserved for API symmetry.
        cfg : CfgType or None
            Reserved for API symmetry.

        Returns
        -------
        TableElement
            Deferred report with the table as its primary artifact.
        """
        if transcript_policy not in {
            "all",
            "mane",
            "canonical",
        } and not _ENST.fullmatch(transcript_policy):
            raise VEPAnnotationError(
                f"Invalid transcript policy: {transcript_policy!r}"
            )
        table = _table(source)
        manifest = prepared.artifacts["manifest"]
        tab = processed.artifacts["vep_tsv"]
        tag = processed.tag.bump(flag="vep-report").resolve(tag)
        if output is None or out is not None:
            spec = OutSpec(tag.default_name, self.folder, ext="parquet").patch(out)
            if spec.ext != "parquet":
                raise VEPAnnotationError("VEPAnnotator output must be Parquet")
            output = spec.file
        output = output.resolve()
        metadata = (metadata or output.with_suffix(".vep.json")).resolve()
        key, _ = self.build_element_name(tag, "report")
        _, cache = processed.pres
        cache_info = cache.artifacts["cache_info"]
        vep_flags = (
            (*self.flags, "--refseq")
            if cache_info.parent.parent.name == "homo_sapiens_refseq"
            else self.flags
        )
        provenance: dict[str, object] = {
            "assembly": self.assembly,
            "reference_fasta": str(self.reference),
            "vep_executable": self.primary_binary,
            "vep_version": self.version,
            "vep_component_versions": self.component_versions,
            "ensembl_cache_release": int(self.release),
            "ensembl_cache_path": str(cache_info.parent),
            "vep_flags": vep_flags,
            "transcript_policy": transcript_policy,
        }
        return TableElement(
            key,
            self._write_report(
                table.resolve(),
                manifest,
                tab,
                output,
                metadata,
                cache_info,
                transcript_policy,
                provenance,
            ),
            tag,
            ArtifactSet(
                TableArtifact(output), primary_name="parquet", metadata=metadata
            ),
            pres=(source, prepared, processed),
            determinants=(transcript_policy, str(output), str(metadata)),
        )

    @staticmethod
    def _write_report(
        table: Path,
        manifest: Path,
        tab: Path,
        output: Path,
        metadata: Path,
        cache_info: Path,
        transcript_policy: str,
        provenance: dict[str, object],
    ) -> Runnable:
        """Build a runner that writes the annotation table and provenance."""

        def run() -> None:
            mapping = json.loads(manifest.read_text())
            variants = [
                _Variant(**{**record, "vcf": tuple(record["vcf"])})
                for record in mapping["variants"]
            ]
            annotations = _read_vep(tab, {v.id for v in variants})
            frame = TableArtifact(table).read(drop_unnamed_columns=False)
            result = _annotate(frame, variants, annotations, transcript_policy)
            report_metadata = {
                **provenance,
                "reference_sha256": mapping["reference_sha256"],
                "ensembl_cache_info": cache_info.read_text(),
                "normalization": {
                    "vcf": "reference-checked, minimal and left-aligned biallelic records",
                    "hgvs_g": "FASTA-verified on linear references; 3-prime shifted pure indels",
                    "hgvs_c": "VEP shift_hgvs=1; independently unverified",
                    "hgvs_p": "VEP protein HGVS output when available",
                    "spdi": "VEP output when available",
                },
            }
            parents(output, metadata)
            with TemporaryDirectory(prefix=".vep-report-", dir=output.parent) as folder:
                work = Path(folder)
                result.to_parquet(work / "result.parquet", index=True)
                (work / "metadata.json").write_text(
                    json.dumps(report_metadata, indent=2) + "\n"
                )
                (work / "metadata.json").replace(metadata)
                (work / "result.parquet").replace(output)

        return Runnable(run, display=f"Write VEP annotation table {output}")

    def annotate(
        self,
        source: Element | FileSource | Path | str,
        *,
        transcript_policy: str = "all",
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> tuple[TableElement, Element, Element]:
        """Composite function that prepares, processes, and reports VEP annotations.

        Parameters
        ----------
        source : Element or FileSource or Path or str
            Upstream table or path to a Parquet, TSV, or CSV variant table.
        transcript_policy : str
            all, mane, canonical, or an ENST identifier. An input transcript_id
            overrides this policy. A versioned ENST requires an exact match.
        tag : TagType or SpecValues or None
            Tag override shared by all steps or keyed by step name.
        out : OutType or SpecValues or None
            Shared output naming based on the report path. Use SpecValues for
            overrides keyed by prepare, process, or report.
        par : ParType or SpecValues or None
            Options shared by all steps or keyed by step name, such as
            ``Par(refseq=True)`` for process.
        cfg : CfgType or SpecValues or None
            Runtime configuration shared by all steps or keyed by step name.

        Returns
        -------
        tuple[TableElement, Element, Element]
            Report, VEP execution, and VCF preparation, in last-to-first order.
        """
        prior = (
            FileSource(TableArtifact(Path(source).resolve()))
            if isinstance(source, (str, Path))
            else source
        )
        table = _table(prior)
        default_tag = ElementTag.from_prior(prior.tag, flag="vep")
        spec = OutSpec(
            stem=default_tag.default_name, folder=self.folder, ext="parquet"
        ).resolve(pick_spec(out, "report"))
        if spec.ext != "parquet":
            raise VEPAnnotationError(
                "VEPAnnotator output must be Parquet to preserve the input index"
            )
        output = spec.file.resolve()
        metadata = output.with_suffix(".vep.json")
        vcf = output.with_suffix(".input.vcf")
        tab = output.with_suffix(".vep.tsv")
        if isinstance(out, SpecValues):
            vcf = (
                OutSpec(vcf.stem, vcf.parent, ext="vcf")
                .resolve(pick_spec(out, "prepare"))
                .file.resolve()
            )
            tab = (
                OutSpec(tab.stem, tab.parent, ext="tsv")
                .resolve(pick_spec(out, "process"))
                .file.resolve()
            )
        manifest = vcf.with_suffix(".json")
        outputs = {output, metadata, tab, vcf, manifest}
        if len(outputs) != 5 or outputs & {table.resolve(), self.reference, self.index}:
            raise VEPAnnotationError(
                "VEP outputs must be distinct from each other and the inputs"
            )
        prepared = self.prepare(
            prior,
            vcf,
            manifest,
            tag=pick_spec(tag, "prepare"),
            par=pick_spec(par, "prepare"),
            cfg=pick_spec(cfg, "prepare"),
        )
        processed = self.process(
            prepared,
            tab,
            tag=pick_spec(tag, "process"),
            par=pick_spec(par, "process"),
            cfg=pick_spec(cfg, "process"),
        )
        reported = self.report(
            prior,
            prepared,
            processed,
            output,
            metadata,
            transcript_policy=transcript_policy,
            tag=pick_spec(tag, "report"),
            par=pick_spec(par, "report"),
            cfg=pick_spec(cfg, "report"),
        )
        return reported, processed, prepared
