from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Literal

from Bio import SeqIO  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.models.artifacts import ArtifactSet, OutputSpec, TableArtifact
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FastqSource,
    FileSource,
    element,
    generate_element_key_name,
)
from mmalignments.models.externals import Runnable
from mmalignments.models.tags import (
    ElementTag,
    Method,
    Omics,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.genomic import (
    hamming_early_break,
    reverse_complement,
)
from mmalignments.services.io import read_frame, write_frame


@element
def uniqueness(
    source: Element | FileSource,
    genomic: Element | FileSource,
    *,
    max_hamming: int = 5,
    barcode_columns: list[str] = ["barcode"],
    direction: Literal["forward", "reverse", "both"] = "both",
    tag: PartialElementTag | ElementTag | None = None,
    outspec: OutputSpec | None = None,
) -> Element:
    """
    Checks for a dataframe of barcodes and their uniqueness in a genomic
    sequence, allowing for a specified number of mismatches. This we can use
    to estimate how error-tolerant our barcodes are in the context of a
    given short genomic sequence. This may be useful for amplicon
    demultiplexing.


    Parameters
    ----------
    source : Element | FileSource
        Element or FileSource containing a dataframe of barcodes to check.
    genomic : Element | FileSource
        Element or FileSource containing a genomic sequence to check against.
    max_hamming : int
        Maximum number of mismatches allowed for a barcode to be considered
        a match in the genomic sequence.
    barcode_columns : list[str]
        List of column names in the dataframe that contain the barcodes.
    direction : Literal["forward", "reverse", "both"]
        Direction to check for uniqueness. "forward" checks the barcodes as-is,
        "reverse" checks the reverse complement of the barcodes, and "both"
        checks both directions.
    tag : PartialElementTag | ElementTag | None
        Optional tag override.
    outspec : OutputSpec | None
        Optional output specification for the configuration file.
    params : Params | None
        Unused; kept for API symmetry.
    cfg : ExternalRunConfig | None
        Unused; kept for API symmetry.

    Returns
    -------
    Element
        The dataframe of barcodes and their uniqueness radius will be saved to
        the specified output file.
    """
    fasta_path = genomic.artifacts.primary.resolve()
    frame_path = source.artifacts.primary.resolve()
    tag = from_prior(
        source.tag,
        tag,
        root=source.root,
        stage=Stage.PREP,
        method=Method.CUSTOM,
        state=State.PREPROCESS,
        omics=Omics.DNA,
        param="uniqueness",
    )
    spec = OutputSpec(tag.default_name, Path("cache/barcodes"), ext="tsv").merge(
        outspec
    )
    outfile = spec.path()
    pres = source.pres if isinstance(source, FileSource) else (source,)
    pres += genomic.pres if isinstance(genomic, FileSource) else (genomic,)
    determinants = tuple(barcode_columns) + (str(max_hamming),)
    runner = uniqueness_radius(
        frame_path,
        fasta_path,
        outfile,
        max_hamming,
        barcode_columns,
        direction=direction,
    )
    key, name = generate_element_key_name(tag, "barcodes", subcommand="uniqueness")
    artifacts = ArtifactSet(TableArtifact(outfile), primary_name=spec.ext)
    return Element(
        key,
        runner,
        tag=tag,
        artifacts=artifacts,
        determinants=determinants,
        pres=pres,
        name=name,
    )


@element
def sampleadapters(
    source: Element | FastqSource,
    *,
    start_length: int = 10,
    sample_size: int = 10000,
    tag: PartialElementTag | ElementTag | None = None,
    outspec: OutputSpec | None = None,
) -> Element:
    """
    Constructs a dataframe with the most frequent starting sequences in a sample
    of reads. Input is an Element or a FileSource, that encapsulates a
    FastqSource. This may contain paired-end sequencing files or single-end.
    The first ``sample_size`` reads are used to determine the most frequent
    starting sequences of the reads up to a length of ``start_length``. The
    output is a dataframe with the most frequent starting sequences in the
    sample of reads, along with their counts and frequencies.


    Parameters
    ----------
    source : Element | FastqSource
        Element or FastqSource containing the sequencing reads.
    start_length : int
        Length of the starting sequences to consider.
    sample_size : int
        Number of reads to sample.
    tag : PartialElementTag | ElementTag | None
        Optional tag override.
    outspec : OutputSpec | None
        Optional output specification for the configuration file.
    params : Params | None
        Unused; kept for API symmetry.
    cfg : ExternalRunConfig | None
        Unused; kept for API symmetry.

    Returns
    -------
    Element
        The dataframe of barcodes and their uniqueness radius will be saved to
        the specified output file.
    """
    fastq = source.artifacts.primary
    tag = from_prior(
        source.tag,
        tag,
        root=source.root,
        stage=Stage.PREP,
        method=Method.CUSTOM,
        state=State.PREPROCESS,
        omics=Omics.DNA,
        param="sampleadapters",
    )
    r1 = fastq.r1
    r2 = fastq.r2  # may be None
    spec = OutputSpec(tag.default_name, Path("cache/barcodes"), ext="tsv").merge(
        outspec
    )
    outfile = spec.path()
    pres = source.pres if isinstance(source, FastqSource) else (source,)
    determinants = (str(start_length), str(sample_size))
    runner = sample_start_reads(
        r1,
        r2,
        outfile,
        start_length=start_length,
        sample_size=sample_size,
    )
    key, name = generate_element_key_name(tag, "barcodes", subcommand="sampleadapters")
    artifacts = ArtifactSet(TableArtifact(outfile), primary_name=spec.ext)
    return Element(
        key,
        runner,
        tag=tag,
        artifacts=artifacts,
        determinants=determinants,
        pres=pres,
        name=name,
    )


################################################################################
# low-level function
################################################################################


def _input_fastq_name(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _count_fastq_prefixes(
    path: Path,
    *,
    start_length: int,
    sample_size: int,
) -> list[tuple[str, int]]:
    counts: Counter[str] = Counter()
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            for idx, record in enumerate(SeqIO.parse(handle, "fastq")):
                if idx >= sample_size:
                    break
                prefix = str(record.seq)[:start_length].upper()
                counts[prefix] += 1
    else:
        with path.open("rt") as handle:
            for idx, record in enumerate(SeqIO.parse(handle, "fastq")):
                if idx >= sample_size:
                    break
                prefix = str(record.seq)[:start_length].upper()
                counts[prefix] += 1

    return sorted(counts.items(), key=lambda x: (-x[1], x[0]))


def _single_end_frame(
    input_name: str,
    counts: list[tuple[str, int]],
) -> DataFrame:
    rows = [
        {
            "input": input_name,
            "Sequence (R1)": seq,
            "Count (R1)": count,
        }
        for seq, count in counts
    ]
    return DataFrame(
        rows,
        columns=["input", "Sequence (R1)", "Count (R1)"],
    )


def _paired_end_frame(
    input_name: str,
    r1_counts: list[tuple[str, int]],
    r2_counts: list[tuple[str, int]],
) -> DataFrame:
    n_rows = max(len(r1_counts), len(r2_counts))
    rows = []

    for i in range(n_rows):
        seq_r1, count_r1 = r1_counts[i] if i < len(r1_counts) else (None, None)
        seq_r2, count_r2 = r2_counts[i] if i < len(r2_counts) else (None, None)
        rows.append(
            {
                "input": input_name,
                "Sequence (R1)": seq_r1,
                "Count (R1)": count_r1,
                "Sequence (R2)": seq_r2,
                "Count (R2)": count_r2,
            }
        )

    return DataFrame(
        rows,
        columns=[
            "input",
            "Sequence (R1)",
            "Count (R1)",
            "Sequence (R2)",
            "Count (R2)",
        ],
    )


def _format_hit_example(
    label: str,
    position: int,
    distance: int,
    window: str,
) -> str:
    return f"{label}@{position}:d{distance}:{window}"


def _flank_queries(
    flank: str,
    direction: Literal["forward", "reverse", "both"],
) -> dict[str, str]:
    reverse_flank = reverse_complement(flank).upper()

    if direction == "forward":
        return {"fwd": flank}

    if direction == "reverse":
        return {"rev": reverse_flank}

    if direction == "both":
        return {"fwd": flank, "rev": reverse_flank}

    raise ValueError("direction must be 'forward', 'reverse', or 'both'")


def _flank_uniqueness_radius_details(
    flank: str,
    genomic: str,
    max_distance: int | None = None,
    direction: Literal["forward", "reverse", "both"] = "both",
    *,
    example_limit: int = 3,
) -> tuple[int, list[str]]:
    flank = flank.upper()
    genomic = genomic.upper()
    queries = _flank_queries(flank, direction)

    L = len(flank)

    distance_hits = defaultdict(int)
    distance_examples: dict[int, list[str]] = defaultdict(list)

    for i in range(len(genomic) - L + 1):
        window = genomic[i : i + L]

        for label, query in queries.items():
            d = hamming_early_break(query, window, max_distance)

            # If max_distance is exceeded, ignore this window
            if max_distance is not None and d > max_distance:
                continue

            distance_hits[d] += 1

            if d > 0:
                examples = distance_examples[d]
                if len(examples) < example_limit:
                    examples.append(_format_hit_example(label, i, d, window))

    if not distance_hits:
        return -1, []

    cumulative_hits = 0

    for d in sorted(distance_hits):
        cumulative_hits += distance_hits[d]

        if cumulative_hits > 1:
            return d - 1, distance_examples.get(d, [])[:example_limit]

    radius = max(distance_hits)
    return radius, distance_examples.get(radius, [])[:example_limit]


def sample_start_reads(
    r1: Path,
    r2: Path | None,
    outfile: Path,
    *,
    start_length: int = 10,
    sample_size: int = 10000,
) -> Runnable:
    def __call() -> DataFrame:
        input_name = _input_fastq_name(r1)
        r1_counts = _count_fastq_prefixes(
            r1,
            start_length=start_length,
            sample_size=sample_size,
        )

        if r2 is None:
            frame = _single_end_frame(input_name, r1_counts)
        else:
            r2_counts = _count_fastq_prefixes(
                r2,
                start_length=start_length,
                sample_size=sample_size,
            )
            frame = _paired_end_frame(input_name, r1_counts, r2_counts)

        write_frame(frame, outfile, index=False)
        return frame

    display = CallSpec(
        path=("sample_start_reads",),
        kwargs={
            "r1": r1,
            "r2": r2,
            "outfile": outfile,
            "start_length": start_length,
            "sample_size": sample_size,
        },
    ).render()
    return Runnable(__call, display=display)


def flank_uniqueness_radius(
    flank: str,
    genomic: str,
    max_distance: int | None = None,
    direction: Literal["forward", "reverse", "both"] = "both",
) -> int:
    """
    Determine the maximum Hamming distance for which a flank remains unique.

    Parameters
    ----------
    flank : str
        The flank sequence to check for uniqueness.
    genomic : str
        The genomic sequence to check against.
    max_distance : int | None
        The maximum Hamming distance to consider. If None, all distances are
        considered.
    direction : Literal["forward", "reverse", "both"]
        The direction to check for uniqueness.

    Returns
    -------
    int
        The maximum Hamming distance for which the flank remains unique.
        Returns -1 if the flank is not unique even at exact matching.
    """

    radius, _ = _flank_uniqueness_radius_details(
        flank,
        genomic,
        max_distance=max_distance,
        direction=direction,
    )
    return radius


def barcode_check_uniqueness_radius(
    barcode_path: Path,
    genomic_path: Path,
    outfile: Path,
    max_hamming: int,
    barcode_columns: list[str],
    direction: Literal["forward", "reverse", "both"] = "both",
) -> DataFrame:
    """
    Check if the barcodes are unique in the genomic sequence, allowing for a
    specified number of mismatches. For each barcode, determine the maximum
    Hamming distance for which it remains unique (Uniqueness radius).

    Parameters
    ----------
    barcode_path : Path
        The path to the barcode file.
    genomic_path : Path
        The path to the genomic FASTA file.
    outfile : Path
        The path to the output file where results will be saved.
    max_hamming : int
        Maximum number of mismatches allowed for a barcode to be considered
        a match in the genomic sequence.
    barcode_columns : list[str]
        List of column names in the barcode file that contain the barcodes.
    direction : Literal["forward", "reverse", "both"]
        Direction to check for uniqueness. "forward" checks the barcodes as-is,
        "reverse" checks the reverse complement of the barcodes, and "both"
        checks both directions.

    Returns
    -------
    DataFrame
        A DataFrame containing the barcodes and their uniqueness radius.
    """
    barcodes = read_frame(barcode_path)
    record = SeqIO.read(genomic_path, "fasta")
    genomic_sequence = str(record.seq)

    for column in barcode_columns:
        barcodes[f"Non-exact hit examples ({column})"] = barcodes[column].apply(
            lambda x: "; ".join(
                _flank_uniqueness_radius_details(
                    x,
                    genomic_sequence,
                    max_hamming,
                    direction,
                )[1]
            )
        )
        barcodes[f"Unique radius ({column})"] = barcodes[column].apply(
            lambda x: _flank_uniqueness_radius_details(
                x,
                genomic_sequence,
                max_hamming,
                direction,
            )[0]
        )
    write_frame(barcodes, outfile)
    return barcodes


def uniqueness_radius(
    barcode_path: Path,
    genomic_path: Path,
    outfile: Path,
    max_hamming: int,
    barcode_columns: list[str],
    direction: Literal["forward", "reverse", "both"] = "both",
) -> Runnable:
    def __call():
        return barcode_check_uniqueness_radius(
            barcode_path,
            genomic_path,
            outfile,
            max_hamming,
            barcode_columns,
            direction,  # noqa: E501
        )

    display = CallSpec(
        path=("uniqueness_radius",),
        kwargs={
            "barcode_path": barcode_path,
            "genomic_path": genomic_path,
            "outfile": outfile,
            "max_hamming": max_hamming,
            "barcode_columns": barcode_columns,
        },
    ).render()
    return Runnable(__call, display=display)


def check_flank_code(
    flank: str, genomic: str, max_hamming: int
) -> tuple[dict[int, tuple[str, int]], bool]:
    """
    Find all occurrences of a flank sequence in genomic sequence allowing
    up to max_hamming mismatches.

    Returns:
        hits:
            dict mapping genomic position -> (matched sequence, hamming distance)

        unique:
            True if exactly one genomic position matches within the allowed
            mismatch threshold.
    """

    flank = flank.upper()
    genomic = genomic.upper()

    L = len(flank)

    hits = {}

    for i in range(len(genomic) - L + 1):
        window = genomic[i : i + L]

        distance = hamming_early_break(flank, window, max_distance=max_hamming)

        if distance <= max_hamming:
            hits[i] = (window, distance)

    return hits, len(hits) == 1
