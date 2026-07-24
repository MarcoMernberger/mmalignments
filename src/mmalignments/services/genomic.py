import math
from pathlib import Path
from typing import Optional

import pandas as pd  # type: ignore[import]
import pyBigWig  # type: ignore[import]
from Bio.Seq import Seq  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]
from scipy.spatial.distance import hamming


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def longest_common_prefix(a: str, b: str) -> int:
    """
    Return the number of identical bases at the beginning of two sequences.
    """
    length = 0

    for x, y in zip(a, b):
        if x != y:
            break
        length += 1

    return length


def match_pattern(a: str, b: str) -> str:
    """
    Return a visual representation of matching positions.

    Example:
        ACTGAA
        ACTCAA

        |||.||
    """
    return "".join("|" if x == y else "." for x, y in zip(a, b))


def hamming_early_break(a: str, b: str, max_distance: Optional[int] = None) -> int:
    if len(b) < len(a):
        a, b = b, a

    mismatches = 0

    for x, y in zip(a, b[: len(a)]):
        if x != y:
            mismatches += 1
            if max_distance is not None and mismatches > max_distance:
                return mismatches

    return mismatches


def hamming_scipy(a: str, b: str, relative: bool = False) -> int:
    distance = hamming(a, b)
    if not relative:
        distance = distance * len(a)
    return distance


def export_sgrna_bed_and_bigwig(
    df: DataFrame,
    chrom_col: str = "Chromosome",
    start_col: str = "Start",
    stop_col: str = "Stop",
    name_col: str = "Feature",
    strand_col: str = "Direction",
    score_cols: Optional[list[str]] = None,
    folder: Path = Path("."),
    out_prefix: str = "sgrna_export",
):
    folder.mkdir(parents=True, exist_ok=True)
    bed_path = folder / f"{out_prefix}.bed"
    bw_paths = {}
    for col in [chrom_col, name_col, start_col, stop_col, strand_col] + (
        score_cols or []
    ):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame.")
        df = df.sort_values([chrom_col, start_col, stop_col])
    # -----------------------
    # BED EXPORT
    # -----------------------
    with open(bed_path, "w") as bed:
        bed.write('track name=sgRNAs description="sgRNAs" visibility=2\n')

        for _, row in df.iterrows():
            try:
                chrom = str(row[chrom_col])
                start = int(row[start_col])
                stop = int(row[stop_col])
            except Exception:
                continue

            name = row[name_col]

            strand = row[strand_col]
            if strand == 1:
                strand = "+"
            elif strand == -1:
                strand = "-"
            else:
                strand = "."

            bed.write(f"{chrom}\t{start}\t{stop}\t{name}\t0\t{strand}\n")

    # -----------------------
    # BIGWIG EXPORT (one per score column)
    # -----------------------
    if score_cols:
        # Prepare normalized, sorted coordinates once for all score columns.
        coords = df[[chrom_col, start_col, stop_col]].copy()
        coords[start_col] = pd.to_numeric(coords[start_col], errors="coerce")
        coords[stop_col] = pd.to_numeric(coords[stop_col], errors="coerce")
        coords = coords.dropna(subset=[chrom_col, start_col, stop_col])
        coords[start_col] = coords[start_col].astype(int)
        coords[stop_col] = coords[stop_col].astype(int)
        coords = coords[
            (coords[start_col] >= 0) & (coords[stop_col] > coords[start_col])
        ]

        if coords.empty:
            return str(bed_path), bw_paths

        # BigWig header expects (chromosome_name, chromosome_size) for each chromosome.
        chrom_sizes = (
            coords.groupby(chrom_col, sort=True)[stop_col].max().astype(int).add(1)
        )
        header = [(str(chrom), int(size)) for chrom, size in chrom_sizes.items()]

        for col in score_cols:
            if col not in df.columns:
                continue

            bw_path = folder / f"{out_prefix}_{col}.bw"
            bw_paths[col] = str(bw_path)

            bw = pyBigWig.open(str(bw_path), "w")

            bw.addHeader(header)

            # BigWig entries must be sorted and non-overlapping. Collapse to 1bp
            # signal per start position and aggregate duplicates.
            values = pd.to_numeric(df[col], errors="coerce")
            entries = coords.copy()
            entries["_value"] = values
            entries = entries.dropna(subset=["_value"])

            if entries.empty:
                bw.close()
                continue

            entries = (
                entries.groupby([chrom_col, start_col], as_index=False, sort=True)[
                    "_value"
                ]
                .mean()
                .sort_values([chrom_col, start_col], kind="mergesort")
            )

            chroms = []
            starts = []
            ends = []
            vals = []

            for _, row in entries.iterrows():
                chrom = str(row[chrom_col])
                start = int(row[start_col])
                val = float(row["_value"])
                if start < 0 or not math.isfinite(val):
                    continue
                chroms.append(chrom)
                starts.append(start)
                ends.append(start + 1)
                vals.append(val)

            if chroms:
                bw.addEntries(chroms, starts, ends=ends, values=vals)

            bw.close()

    return str(bed_path), bw_paths


def export_sgrna_bed(
    df,
    output_path: Path = Path("sgrnas.bed"),
    chrom_col: str = "Chromosome",
    start_col: str = "Start",
    stop_col: str = "Stop",
    name_col: str = "Feature",
    score_cols: Optional[list[str]] = None,
    strand_col: str = "Direction",
):
    out = Path(output_path)

    with open(out, "w") as f:
        f.write('track name=sgRNAs description="sgRNA tracks" visibility=2\n')

        for _, row in df.iterrows():
            try:
                start = int(row[start_col])
                stop = int(row[stop_col])
            except Exception:
                continue

            name = str(row.get(name_col, "sgRNA"))

            # score must be 0–1000 in BED convention
            if score_cols:
                if isinstance(score_cols, list) and len(score_cols) > 1:
                    score = 0.0
                    score_strings = []
                    for score_col in score_cols:
                        if score_col in df.columns:
                            try:
                                score = str(row[score_col])
                            except Exception:
                                score = "0.0"
                            score_strings.append(f"{score_col}={score}")
                        else:
                            score_strings.append(f"{score_col}=NA")
                    name = f"name={name}|{'|'.join(score_strings)}"
                elif isinstance(score_cols, list) and len(score_cols) == 1:
                    score_col = score_cols[0]
                    if score_col in df.columns:
                        try:
                            score = float(row[score_col])
                        except Exception:
                            score = 0.0
                    else:
                        score = 0.0
                    name = f"name={name}|{score_col}={score}"
            else:
                score = 0.0

            # BED score must be int 0-1000
            score = max(0, min(1000, int(score)))

            strand = row.get(strand_col, ".")
            if strand == 1:
                strand = "+"
            elif strand == -1:
                strand = "-"
            else:
                strand = "."

            f.write(
                "\t".join(
                    [
                        str(row[chrom_col]),
                        str(start),
                        str(stop),
                        name,
                        str(score),
                        strand,
                    ]
                )
                + "\n"
            )

    return str(out)
