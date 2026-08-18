from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Literal

import pandas as pd  # type: ignore[import]
from Bio import SeqIO  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.models.artifacts import (
    ArtifactSet,
    FileArtifact,
    FileSpec,
    OutputSpec,
    TableArtifact,
)
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FastqSource,
    FileSource,
    element,
    generate_element_key_name,
)
from mmalignments.models.externals import Runnable
from mmalignments.models.reports import write_sampleadapters_report
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
    longest_common_prefix,
    match_pattern,
    reverse_complement,
)
from mmalignments.services.io import read_frame, write_frame


@element
def uniqueness(
    source: Element | FileSource,
    genomic: Element | FileSource,
    *,
    max_hamming: int = 5,
    barcode_columns: list[str] = ["used_start_barcode", "used_end_barcode"],
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
    spec = OutputSpec(
        tag.default_name,
        Path("cache/barcodes"),
        ext="tsv",
    ).merge(outspec)
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
    key, name = generate_element_key_name(
        tag,
        "barcodes",
        subcommand="uniqueness",
    )
    artifacts = ArtifactSet(
        TableArtifact(outfile),
        primary_name=spec.ext,
        # html=FileArtifact(outfile.with_suffix(".html")),
    )
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
    barcodes: Element | FileSource | None = None,
    *,
    start_length: int = 18,
    sample_size: int = 50000,
    forward_col: list[str] | None = ["used_start_barcode"],
    reverse_col: list[str] | None = ["used_end_barcode"],
    sample_col: str = "sample",
    max_edit_distance: int = 8,
    min_count: int = 1,
    leven: bool = False,
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
    barcodes : Element | FileSource | None
        Optional barcode table used to classify top start sequences.
    start_length : int
        Length of the starting sequences to consider.
    sample_size : int
        Number of reads to sample.
    forward_col : list[str] | None
        Barcode columns containing forward/start barcodes.
    reverse_col : list[str] | None
        Barcode columns containing reverse/end barcodes.
    sample_col : str
        Column in barcode table containing sample names.
    max_edit_distance : int
        Maximum distance for assigning a barcode match.
    min_count : int
        Minimum count required before a sequence is classified.
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
    spec = OutputSpec(
        tag.default_name,
        Path("cache/barcodes"),
        ext="tsv",
    ).merge(outspec)
    outfile = spec.path()
    pres = source.pres if isinstance(source, FastqSource) else (source,)
    barcode_path = None
    if barcodes is not None:
        barcode_path = barcodes.artifacts.primary.resolve()
        pres += barcodes.pres if isinstance(barcodes, FileSource) else (barcodes,)

    determinants = (
        str(start_length),
        str(sample_size),
        str(forward_col),
        str(reverse_col),
        sample_col,
        str(max_edit_distance),
        str(min_count),
        str(leven),
    )
    runner = sample_start_reads(
        r1,
        r2,
        outfile,
        start_length=start_length,
        sample_size=sample_size,
        barcode_path=barcode_path,
        barcode_columns_r1=forward_col,
        barcode_columns_r2=reverse_col,
        barcode_sample_column=sample_col,
        max_edit_distance=max_edit_distance,
        min_count=min_count,
        leven=leven,
    )
    key, name = generate_element_key_name(
        tag,
        "barcodes",
        subcommand="sampleadapters",
    )
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
def sampleadapters_report(
    source: Element | FileSource,
    *,
    tag: PartialElementTag | ElementTag | None = None,
    outspec: OutputSpec | None = None,
) -> Element:
    """
    Builds a QC report from the TSV output of ``sampleadapters``.

    The report stays downstream of the table-producing element, so the DAG
    keeps the detailed assignments and the report as separate nodes.
    """
    table_path = source.artifacts.primary.resolve()
    tag = from_prior(
        source.tag,
        tag,
        root=source.root,
        stage=Stage.SUMMARY,
        method=Method.CUSTOM,
        state=State.REPORT,
        omics=Omics.DNA,
        param="sampleadapters_report",
    )
    spec = OutputSpec(
        tag.default_name,
        Path("cache/barcodes"),
        ext="html",
    ).merge(outspec)
    spec = spec.add_output(
        "tsv",
        FileSpec(spec.stem or tag.default_name, "tsv"),
    )
    report_path = spec.path()
    summary_path = spec.files["tsv"]

    pres = source.pres if isinstance(source, FileSource) else (source,)

    def __call() -> Path:
        return write_sampleadapters_report(
            table_path,
            report_path,
            summary_path,
        )

    display = CallSpec(
        path=("sampleadapters_report",),
        kwargs={
            "table_path": table_path,
            "report_path": report_path,
            "summary_path": summary_path,
        },
    ).render()
    key, name = generate_element_key_name(
        tag,
        "barcodes",
        subcommand="sampleadapters_report",
    )
    artifacts = ArtifactSet(
        FileArtifact(report_path),
        primary_name=spec.ext,
        tsv=TableArtifact(summary_path),
    )
    return Element(
        key,
        Runnable(__call, display=display),
        tag=tag,
        artifacts=artifacts,
        determinants=(),
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


def _count_fastq_prefix_pairs(
    r1_path: Path,
    r2_path: Path,
    *,
    start_length: int,
    sample_size: int,
) -> list[tuple[tuple[str, str], int]]:
    counts: Counter[tuple[str, str]] = Counter()

    r1_opener = gzip.open if r1_path.suffix == ".gz" else Path.open
    r2_opener = gzip.open if r2_path.suffix == ".gz" else Path.open

    with (
        r1_opener(r1_path, "rt") as r1_handle,
        r2_opener(r2_path, "rt") as r2_handle,
    ):
        r1_records = SeqIO.parse(r1_handle, "fastq")
        r2_records = SeqIO.parse(r2_handle, "fastq")

        for idx, (r1, r2) in enumerate(zip(r1_records, r2_records)):
            if idx >= sample_size:
                break

            r1_prefix = str(r1.seq[:start_length]).upper()
            r2_prefix = str(r2.seq[:start_length]).upper()

            counts[(r1_prefix, r2_prefix)] += 1

    return sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0]),
    )


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


# def _paired_end_frame(
#     input_name: str,
#     r1_counts: list[tuple[str, int]],
#     r2_counts: list[tuple[str, int]],
# ) -> DataFrame:
#     n_rows = max(len(r1_counts), len(r2_counts))
#     rows = []

#     for i in range(n_rows):
#         seq_r1, count_r1 = r1_counts[i] if i < len(r1_counts) else (None, None)
#         seq_r2, count_r2 = r2_counts[i] if i < len(r2_counts) else (None, None)
#         rows.append(
#             {
#                 "input": input_name,
#                 "Sequence (R1)": seq_r1,
#                 "Count (R1)": count_r1,
#                 "Sequence (R2)": seq_r2,
#                 "Count (R2)": count_r2,
#             }
#         )

#     return DataFrame(
#         rows,
#         columns=[
#             "input",
#             "Sequence (R1)",
#             "Count (R1)",
#             "Sequence (R2)",
#             "Count (R2)",
#         ],
#     )


def _paired_end_frame(
    input_name: str,
    r1_r2_counts: list[tuple[tuple[str, str], int]],
) -> DataFrame:
    n_rows = len(r1_r2_counts)
    rows = []

    for i in range(n_rows):
        (seq_r1, seq_r2), count = r1_r2_counts[i]
        rows.append(
            {
                "input": input_name,
                "Sequence (R1)": seq_r1,
                "Count (R1)": count,
                "Sequence (R2)": seq_r2,
                "Count (R2)": count,
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


def _levenshtein_distance(
    a: str,
    b: str,
    max_distance: int | None = None,
) -> int:
    if a == b:
        return 0

    if len(a) < len(b):
        a, b = b, a

    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        row_min = i

        for j, cb in enumerate(b, start=1):
            cost = 0 if ca == cb else 1
            ins = current[j - 1] + 1
            delete = previous[j] + 1
            sub = previous[j - 1] + cost
            value = min(ins, delete, sub)
            current.append(value)
            if value < row_min:
                row_min = value

        if max_distance is not None and row_min > max_distance:
            return max_distance + 1

        previous = current

    return previous[-1]


def _barcode_candidates(
    barcodes: DataFrame,
    sample_column: str,
    barcode_columns: list[str],
) -> list[tuple[str, str, str]]:
    candidates: list[tuple[str, str, str]] = []
    for _, row in barcodes.iterrows():
        sample = row.get(sample_column)
        if sample is None:
            continue

        sample_name = str(sample)
        for col in barcode_columns:
            value = row.get(col)
            if value is None:
                continue
            barcode = str(value).strip().upper()
            if barcode and barcode.lower() != "nan":
                candidates.append((sample_name, col, barcode))
    return candidates


def _classify_sequence(
    sequence: str,
    candidates: list[tuple[str, str, str]],
    max_edit_distance: int,
    distance_func: Callable[[str, str, int | None], int],
) -> tuple[str, int, str, str, int, str]:

    if not candidates:
        return (
            "unknown",
            -1,
            "unknown",
            "",
            0,
            "",
        )

    best_sample = "unknown"
    best_col = "unknown"
    best_barcode = ""
    best_distance = max_edit_distance + 1

    for sample_name, barcode_col, barcode in candidates:
        distance = distance_func(
            sequence,
            barcode,
            max_edit_distance,
        )

        if distance < best_distance:
            best_distance = distance
            best_sample = sample_name
            best_col = barcode_col
            best_barcode = barcode

    if best_distance > max_edit_distance:
        return (
            "unknown",
            -1,
            "unknown",
            "",
            0,
            "",
        )

    prefix_length = longest_common_prefix(
        sequence,
        best_barcode,
    )

    pattern = match_pattern(
        sequence,
        best_barcode,
    )

    return (
        best_sample,
        best_distance,
        f"{best_col}:{best_sample}:{best_barcode}",
        best_barcode,
        prefix_length,
        pattern,
    )


def _classify_row_with_threshold(
    seq: str,
    count: int,
    candidates: list[tuple[str, str, str]],
    max_edit_distance: int,
    min_count: int,
    distance_func: Callable[[str, str, int | None], int],
) -> tuple[str, int, str, str, int, str]:
    if count < min_count:
        return "unknown", -1, "unknown", "", 0, ""
    return _classify_sequence(
        seq, candidates, max_edit_distance, distance_func=distance_func
    )


def _classify_column(
    frame: DataFrame,
    sequence_column: str,
    count_column: str,
    candidates: list[tuple[str, str, str]],
    max_edit_distance: int,
    min_count: int,
    distance_func: Callable[[str, str, int | None], int],
) -> tuple:
    matches = frame.apply(
        lambda row: _classify_row_with_threshold(
            str(row[sequence_column]).upper(),
            int(row[count_column]) if not pd.isna(row[count_column]) else 0,
            candidates,
            max_edit_distance,
            min_count,
            distance_func=distance_func,
        ),
        axis=1,
    )
    sample_series = matches.apply(lambda x: x[0])
    distance_series = matches.apply(lambda x: x[1])
    best_series = matches.apply(lambda x: x[2])
    barcode_series = matches.apply(lambda x: x[3])
    prefix_series = matches.apply(lambda x: x[4])
    pattern_series = matches.apply(lambda x: x[5])

    return (
        sample_series,
        distance_series,
        best_series,
        barcode_series,
        prefix_series,
        pattern_series,
    )


def _classify_single_end_frame(
    frame: DataFrame,
    barcodes: DataFrame,
    barcode_sample_column: str,
    merged_columns: list[str],
    max_edit_distance: int,
    min_count: int,
    distance_func: Callable[[str, str, int | None], int],
) -> DataFrame:
    candidates = _barcode_candidates(
        barcodes,
        barcode_sample_column,
        merged_columns,
    )
    (
        r1_sample,
        r1_dist,
        r1_best,
        r1_barcode,
        r1_prefix,
        r1_pattern,
    ) = _classify_column(
        frame,
        "Sequence (R1)",
        "Count (R1)",
        candidates=candidates,
        max_edit_distance=max_edit_distance,
        min_count=min_count,
        distance_func=distance_func,
    )
    classified = frame.copy()
    classified["Sample (R1)"] = r1_sample
    classified["Edit distance (R1)"] = r1_dist
    classified["Best (R1)"] = r1_best
    classified["Matched barcode (R1)"] = r1_barcode
    classified["Prefix match length (R1)"] = r1_prefix
    classified["Match pattern (R1)"] = r1_pattern
    return classified


def _classify_paired_frame(
    frame: DataFrame,
    barcodes: DataFrame,
    barcode_sample_column: str,
    start_columns: list[str],
    end_columns: list[str],
    max_edit_distance: int,
    min_count: int,
    distance_func: Callable[[str, str, int | None], int],
) -> DataFrame:
    all_candidates = _barcode_candidates(
        barcodes,
        barcode_sample_column,
        list(dict.fromkeys(start_columns + end_columns)),
    )
    (
        r1_sample,
        r1_dist,
        r1_best,
        r1_barcode,
        r1_prefix,
        r1_pattern,
    ) = _classify_column(
        frame,
        "Sequence (R1)",
        "Count (R1)",
        all_candidates,
        max_edit_distance,
        min_count,
        distance_func,
    )

    (
        r2_sample,
        r2_dist,
        r2_best,
        r2_barcode,
        r2_prefix,
        r2_pattern,
    ) = _classify_column(
        frame,
        "Sequence (R2)",
        "Count (R2)",
        all_candidates,
        max_edit_distance,
        min_count,
        distance_func,
    )
    classified = frame.copy()

    classified["Sample (R1)"] = r1_sample
    classified["Edit distance (R1)"] = r1_dist
    classified["Best (R1)"] = r1_best
    classified["Matched barcode (R1)"] = r1_barcode
    classified["Prefix match length (R1)"] = r1_prefix
    classified["Match pattern (R1)"] = r1_pattern

    classified["Sample (R2)"] = r2_sample
    classified["Edit distance (R2)"] = r2_dist
    classified["Best (R2)"] = r2_best
    classified["Matched barcode (R2)"] = r2_barcode
    classified["Prefix match length (R2)"] = r2_prefix
    classified["Match pattern (R2)"] = r2_pattern

    return classified


def _format_hit_example(
    label: str,
    position: int,
    distance: int,
    window: str,
) -> str:
    return f"{label}@{position}:d{distance}:{window}"


# def _flank_queries(
#     flank: str,
#     direction: Literal["forward", "reverse", "both"],
# ) -> dict[str, str]:
#     reverse_flank = reverse_complement(flank).upper()

#     if direction == "forward":
#         return {"fwd": flank}

#     if direction == "reverse":
#         return {"rev": reverse_flank}

#     if direction == "both":
#         return {"fwd": flank, "rev": reverse_flank}

#     raise ValueError("direction must be 'forward', 'reverse', or 'both'")


# def _flank_uniqueness_radius_details(
#     flank: str,
#     genomic: str,
#     max_distance: int | None = None,
#     direction: Literal["forward", "reverse", "both"] = "both",
#     *,
#     example_limit: int = 3,
# ) -> tuple[int, list[str]]:
#     flank = flank.upper()
#     genomic = genomic.upper()
#     queries = _flank_queries(flank, direction)

#     L = len(flank)

#     distance_hits = defaultdict(int)
#     distance_examples: dict[int, list[str]] = defaultdict(list)

#     for i in range(len(genomic) - L + 1):
#         window = genomic[i : i + L]

#         for label, query in queries.items():
#             d = hamming_early_break(query, window, max_distance)

#             # If max_distance is exceeded, ignore this window
#             if max_distance is not None and d > max_distance:
#                 continue

#             distance_hits[d] += 1

#             if d > 0:
#                 examples = distance_examples[d]
#                 if len(examples) < example_limit:
#                     examples.append(_format_hit_example(label, i, d, window))

#     if not distance_hits:
#         return -1, []

#     cumulative_hits = 0

#     for d in sorted(distance_hits):
#         cumulative_hits += distance_hits[d]

#         if cumulative_hits > 1:
#             return d - 1, distance_examples.get(d, [])[:example_limit]

#     radius = max(distance_hits)
#     return radius, distance_examples.get(radius, [])[:example_limit]


def _flank_queries(
    flank: str,
    direction: Literal["forward", "reverse", "both"],
) -> tuple[tuple[str, str], ...]:
    """
    Return the query sequences for the requested direction(s).

    For ``both``, forward and reverse-complement queries are returned.
    """
    flank = flank.upper()

    if direction == "forward":
        return (("fwd", flank),)

    reverse_flank = reverse_complement(flank).upper()

    if direction == "reverse":
        return (("rev", reverse_flank),)

    if direction == "both":
        return (
            ("fwd", flank),
            ("rev", reverse_flank),
        )

    raise ValueError("direction must be 'forward', 'reverse', or 'both'")


def _flank_uniqueness_radius_details(
    flank: str,
    genomic: str,
    max_distance: int | None = None,
    direction: Literal["forward", "reverse", "both"] = "both",
    *,
    example_limit: int = 3,
) -> tuple[int, list[str]]:
    """
    Determine the uniqueness radius of a flank in a genomic sequence.

    A genomic position counts as one hit, regardless of whether the
    forward query, reverse-complement query, or both queries match there.

    For ``direction="both"``, the best (smallest) Hamming distance of
    the forward and reverse-complement query at each genomic position
    is used.

    The uniqueness radius is the largest Hamming distance for which
    at most one distinct genomic position is matched.

    If zero or one genomic positions match within ``max_distance``,
    ``max_distance`` is returned.

    Parameters
    ----------
    flank : str
        Barcode/flank sequence to search for.
    genomic : str
        Genomic sequence to search.
    max_distance : int | None
        Maximum Hamming distance to consider. This is primarily a
        performance limit, not a biological limit.
    direction : Literal["forward", "reverse", "both"]
        Which orientation(s) to search.
    example_limit : int
        Maximum number of non-exact collision examples to return.

    Returns
    -------
    tuple[int, list[str]]
        The uniqueness radius and examples of hits at the distance
        where uniqueness is lost.
    """
    flank = flank.upper()
    genomic = genomic.upper()

    length = len(flank)

    if length == 0:
        raise ValueError("flank must not be empty")

    if length > len(genomic):
        return max_distance or 0, []

    queries = _flank_queries(flank, direction)

    # Number of DISTINCT genomic positions hit at each distance.
    distance_hits: defaultdict[int, int] = defaultdict(int)

    # Examples are only needed for distances at which multiple
    # genomic positions exist.
    distance_examples: defaultdict[int, list[str]] = defaultdict(list)

    for position in range(len(genomic) - length + 1):
        window = genomic[position : position + length]

        # Find the best orientation for this genomic position.
        #
        # For "forward" / "reverse" this loop contains only one query.
        # For "both", we use the smaller distance of the two orientations.
        best_distance = None
        best_label = None

        for label, query in queries:
            distance = hamming_early_break(
                query,
                window,
                max_distance,
            )

            # This orientation cannot produce a relevant hit.
            if max_distance is not None and distance > max_distance:
                continue

            # Keep the best orientation for this genomic position.
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_label = label

                # Exact match is the theoretical minimum.
                # There is no reason to test the other orientation.
                if distance == 0:
                    break

        # No orientation matches within max_distance.
        if best_distance is None:
            continue

        distance_hits[best_distance] += 1

        if best_distance > 0:
            examples = distance_examples[best_distance]

            if len(examples) < example_limit:
                examples.append(
                    _format_hit_example(
                        best_label,
                        position,
                        best_distance,
                        window,
                    )
                )

        # Two exact hits mean radius 0. Nothing can make the radius
        # smaller, so we can stop immediately.
        if best_distance == 0 and distance_hits[0] > 1:
            return (
                0,
                distance_examples.get(0, [])[:example_limit],
            )

    # No genomic position matched within max_distance.
    #
    # This means the barcode is unique throughout the entire
    # search radius.
    if not distance_hits:
        return max_distance or 0, []

    cumulative_hits = 0

    for distance in sorted(distance_hits):
        cumulative_hits += distance_hits[distance]

        # The second DISTINCT genomic position appears at this
        # Hamming distance. Therefore the previous distance was
        # the last distance at which the barcode was unique.
        if cumulative_hits > 1:
            return (
                distance - 1,
                distance_examples[distance][:example_limit],
            )

    # There was at most one distinct genomic position within the
    # entire search radius.
    return (
        max_distance if max_distance is not None else max(distance_hits),
        distance_examples.get(
            max(distance_hits),
            [],
        )[:example_limit],
    )


def barcode_check_uniqueness_radius(
    barcode_path: Path,
    genomic_path: Path,
    outfile: Path,
    max_hamming: int,
    barcode_columns: list[str],
    direction: Literal["forward", "reverse", "both"] = "both",
) -> DataFrame:
    """
    Check barcode uniqueness in a genomic sequence.

    For each barcode, determine the largest Hamming distance for
    which at most one distinct genomic position matches.

    ``max_hamming`` acts as a computational search limit. If zero
    or one genomic positions match within this limit, the returned
    uniqueness radius is ``max_hamming``.

    Results are written to ``outfile`` and returned as a DataFrame.
    """
    if max_hamming < 0:
        raise ValueError("max_hamming must be >= 0")

    barcodes = read_frame(barcode_path)

    record = SeqIO.read(genomic_path, "fasta")
    genomic_sequence = str(record.seq).upper()

    for column in barcode_columns:
        results = barcodes[column].map(
            lambda barcode: _flank_uniqueness_radius_details(
                barcode,
                genomic_sequence,
                max_hamming,
                direction,
            )
        )

        barcodes[f"Unique radius ({column})"] = results.str[0]

        barcodes[f"Non-exact hit examples ({column})"] = results.str[1].str.join("; ")

    write_frame(barcodes, outfile)

    return barcodes


def sample_start_reads(
    r1: Path,
    r2: Path | None,
    outfile: Path,
    *,
    start_length: int = 10,
    sample_size: int = 10000,
    barcode_path: Path | None = None,
    barcode_columns_r1: list[str] | None = None,
    barcode_columns_r2: list[str] | None = None,
    barcode_sample_column: str = "sample",
    max_edit_distance: int = 5,
    min_count: int = 1,
    leven: bool = False,
) -> Runnable:
    def __call() -> DataFrame:
        distance_func = _levenshtein_distance if leven else hamming_early_break
        input_name = _input_fastq_name(r1)

        if r2 is None:
            counts = _count_fastq_prefixes(
                r1,
                start_length=start_length,
                sample_size=sample_size,
            )
            frame = _single_end_frame(input_name, counts)
        else:
            counts = _count_fastq_prefix_pairs(
                r1,
                r2,
                start_length=start_length,
                sample_size=sample_size,
            )
            frame = _paired_end_frame(input_name, counts)

        if barcode_path is not None:
            barcodes = read_frame(barcode_path)

            all_columns = []
            if barcode_columns_r1:
                all_columns.extend(barcode_columns_r1)
            if barcode_columns_r2:
                all_columns.extend(barcode_columns_r2)

            merged_columns = list(dict.fromkeys(all_columns))

            if r2 is None:
                frame = _classify_single_end_frame(
                    frame,
                    barcodes,
                    barcode_sample_column,
                    merged_columns,
                    max_edit_distance,
                    min_count,
                    distance_func=distance_func,
                )
            else:
                start_columns = list(barcode_columns_r1 or ())
                end_columns = list(barcode_columns_r2 or ())
                if not start_columns and merged_columns:
                    start_columns = merged_columns
                if not end_columns:
                    end_columns = start_columns
                frame = _classify_paired_frame(
                    frame,
                    barcodes,
                    barcode_sample_column,
                    start_columns,
                    end_columns,
                    max_edit_distance,
                    min_count,
                    distance_func=distance_func,
                )

        write_frame(frame, outfile, index=False)
        # write_sampleadapters_report(outfile, outfile.with_suffix(".html"))
        return frame

    display = CallSpec(
        path=("sample_start_reads",),
        kwargs={
            "r1": r1,
            "r2": r2,
            "outfile": outfile,
            "start_length": start_length,
            "sample_size": sample_size,
            "barcode_path": barcode_path,
            "barcode_columns_r1": barcode_columns_r1,
            "barcode_columns_r2": barcode_columns_r2,
            "barcode_sample_column": barcode_sample_column,
            "max_edit_distance": max_edit_distance,
            "min_count": min_count,
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


# def barcode_check_uniqueness_radius(
#     barcode_path: Path,
#     genomic_path: Path,
#     outfile: Path,
#     max_hamming: int,
#     barcode_columns: list[str],
#     direction: Literal["forward", "reverse", "both"] = "both",
# ) -> DataFrame:
#     """
#     Check if the barcodes are unique in the genomic sequence, allowing for a
#     specified number of mismatches. For each barcode, determine the maximum
#     Hamming distance for which it remains unique (Uniqueness radius).

#     Parameters
#     ----------
#     barcode_path : Path
#         The path to the barcode file.
#     genomic_path : Path
#         The path to the genomic FASTA file.
#     outfile : Path
#         The path to the output file where results will be saved.
#     max_hamming : int
#         Maximum number of mismatches allowed for a barcode to be considered
#         a match in the genomic sequence.
#     barcode_columns : list[str]
#         List of column names in the barcode file that contain the barcodes.
#     direction : Literal["forward", "reverse", "both"]
#         Direction to check for uniqueness. "forward" checks the barcodes as-is,
#         "reverse" checks the reverse complement of the barcodes, and "both"
#         checks both directions.

#     Returns
#     -------
#     DataFrame
#         A DataFrame containing the barcodes and their uniqueness radius.
#     """
#     barcodes = read_frame(barcode_path)
#     record = SeqIO.read(genomic_path, "fasta")
#     genomic_sequence = str(record.seq)

#     for column in barcode_columns:
#         results = barcodes[column].map(
#             lambda barcode: _flank_uniqueness_radius_details(
#                 barcode,
#                 genomic_sequence,
#                 max_hamming,
#                 direction,
#             )
#         )

#         barcodes[f"Unique radius ({column})"] = results.str[0]
#         barcodes[f"Non-exact hit examples ({column})"] = results.str[1].str.join("; ")

#     write_frame(barcodes, outfile)
#     return barcodes


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
            dict mapping genomic position -> tuple of matched sequence and
            hamming distance

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
