from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.services.io import parents, read_frame, write_frame


def _format_int(value: int) -> str:
    return f"{value:,}"


def _table_block(frame: DataFrame, *, index: bool = False) -> str:
    if frame.empty:
        return "<p>No rows available.</p>"
    return frame.to_html(
        index=index,
        classes=["table", "table-striped", "table-sm"],
    )


def _summary_table(frame: DataFrame) -> DataFrame:
    rows: list[dict[str, str]] = [{"Metric": "Rows", "Value": _format_int(len(frame))}]

    for count_column in ("Count (R1)", "Count (R2)"):
        if count_column in frame.columns:
            total = int(frame[count_column].fillna(0).sum())
            rows.append(
                {
                    "Metric": f"Total {count_column}",
                    "Value": _format_int(total),
                }
            )

    return DataFrame(rows, columns=["Metric", "Value"])


def _classification_tables(frame: DataFrame) -> list[tuple[str, DataFrame]]:
    sections: list[tuple[str, DataFrame]] = []

    for read_label in ("R1", "R2"):
        sample_column = f"Sample ({read_label})"
        count_column = f"Count ({read_label})"
        distance_column = f"Edit distance ({read_label})"

        if sample_column not in frame.columns:
            continue

        classified = frame[frame[sample_column].fillna("unknown") != "unknown"].copy()
        if classified.empty:
            continue

        if count_column in classified.columns:
            counts = classified.groupby(sample_column, dropna=False)[count_column].sum()
        else:
            counts = classified.groupby(sample_column, dropna=False).size()

        summary = counts.reset_index(name="Assigned count")
        if distance_column in classified.columns:
            min_distance = classified.groupby(sample_column, dropna=False)[
                distance_column
            ].min()
            summary["Best edit distance"] = summary[sample_column].map(min_distance)

        summary = summary.sort_values(
            by=["Assigned count", sample_column],
            ascending=[False, True],
        )
        sections.append((f"Assignments ({read_label})", summary))

    return sections


def _top_rows(frame: DataFrame, *, limit: int = 25) -> DataFrame:
    sort_columns = [col for col in ("Count (R1)", "Count (R2)") if col in frame.columns]
    if not sort_columns:
        return frame.head(limit)
    return frame.sort_values(by=sort_columns, ascending=False).head(limit)


def _read_labels(frame: DataFrame) -> tuple[str, ...]:
    labels = [label for label in ("R1", "R2") if f"Count ({label})" in frame.columns]
    return tuple(labels)


def _empty_assignment_summary() -> DataFrame:
    return DataFrame(
        columns=[
            "read",
            "sample",
            "barcode_type",
            "matched_barcode",
            "assigned_count",
            "read_total",
            "share_of_read",
            "mean_edit_distance",
            "min_edit_distance",
        ]
    )


def _sampleadapters_summary_usecols(column_name: str) -> bool:
    needed_columns = {
        f"{prefix} ({read_label})"
        for read_label in ("R1", "R2")
        for prefix in (
            "Assigned",
            "Best",
            "Count",
            "Edit distance",
            "Matched barcode",
        )
    }
    return column_name in needed_columns


def _assignment_summary(frame: DataFrame) -> DataFrame:
    summaries = []

    read_labels = _read_labels(frame)

    for read_label in read_labels:
        assigned_column = f"Assigned ({read_label})"
        count_column = f"Count ({read_label})"
        best_column = f"Best ({read_label})"
        distance_column = f"Edit distance ({read_label})"
        barcode_column = f"Matched barcode ({read_label})"

        if assigned_column not in frame.columns or count_column not in frame.columns:
            continue

        assigned = frame[assigned_column]

        counts = pd.to_numeric(frame[count_column], errors="coerce")

        mask = assigned.notna() & assigned.ne("unknown") & counts.notna()

        if not mask.any():
            continue

        idx = mask[mask].index

        # Nur relevante Zeilen kopieren
        sample = assigned.loc[idx].astype(str)
        count = counts.loc[idx].astype(float)

        if best_column in frame:
            best = frame.loc[idx, best_column].fillna("").astype(str)
        else:
            best = pd.Series("", index=idx)

        if barcode_column in frame:
            barcode = frame.loc[idx, barcode_column].fillna("").astype(str)
        else:
            barcode = pd.Series("", index=idx)

        if distance_column in frame:
            distance = pd.to_numeric(
                frame.loc[idx, distance_column], errors="coerce"
            ).fillna(0)
        else:
            distance = pd.Series(0.0, index=idx)

        # Barcode-Typ schneller extrahieren
        barcode_type = best.str.partition(":")[0]
        barcode_type = barcode_type.where(
            best.str.contains(":", regex=False), "unknown"
        )

        classified = pd.DataFrame(
            {
                "read": read_label,
                "sample": sample.values,
                "barcode_type": barcode_type.values,
                "matched_barcode": barcode.values,
                "count": count.values,
                "edit_distance": distance.values,
            }
        )

        summaries.append(
            classified.groupby(
                [
                    "read",
                    "sample",
                    "barcode_type",
                    "matched_barcode",
                ],
                as_index=False,
                dropna=False,
            ).agg(
                assigned_count=("count", "sum"),
                mean_edit_distance=("edit_distance", "mean"),
                min_edit_distance=("edit_distance", "min"),
            )
        )

    if not summaries:
        return _empty_assignment_summary()

    summary = pd.concat(
        summaries,
        ignore_index=True,
    )

    summary["read_total"] = summary.groupby("read")["assigned_count"].transform("sum")

    summary["share_of_read"] = summary["assigned_count"] / summary[
        "read_total"
    ].replace(0, pd.NA)

    summary["assigned_count"] = summary["assigned_count"].round().astype(int)

    summary["read_total"] = summary["read_total"].round().astype(int)

    summary["mean_edit_distance"] = summary["mean_edit_distance"].round(3)

    summary["min_edit_distance"] = summary["min_edit_distance"].round().astype(int)

    return summary.sort_values(
        [
            "read",
            "assigned_count",
            "sample",
            "barcode_type",
        ],
        ascending=[
            True,
            False,
            True,
            True,
        ],
    )


def _read_metrics_table(frame: DataFrame) -> DataFrame:
    rows: list[dict[str, object]] = []

    for read_label in _read_labels(frame):
        count_column = f"Count ({read_label})"
        sample_column = f"Sample ({read_label})"
        distance_column = f"Edit distance ({read_label})"

        if count_column not in frame.columns:
            continue

        counts = pd.to_numeric(frame[count_column], errors="coerce").fillna(0)
        total = float(counts.sum())

        if sample_column in frame.columns:
            assigned_mask = frame[sample_column].fillna("unknown") != "unknown"
        else:
            assigned_mask = pd.Series(False, index=frame.index)

        assigned = float(counts[assigned_mask].sum())
        unknown = total - assigned

        if distance_column in frame.columns:
            distances = pd.to_numeric(frame[distance_column], errors="coerce").fillna(0)
            perfect = float(counts[(assigned_mask) & (distances == 0)].sum())
            one_mismatch = float(counts[(assigned_mask) & (distances == 1)].sum())
            multi_mismatch = float(counts[(assigned_mask) & (distances >= 2)].sum())
            weighted_distance = float(
                (counts[assigned_mask] * distances[assigned_mask]).sum()
            )
            mean_distance = weighted_distance / assigned if assigned else 0.0
        else:
            perfect = 0.0
            one_mismatch = 0.0
            multi_mismatch = 0.0
            mean_distance = 0.0

        rows.append(
            {
                "Read": read_label,
                "Total count": int(round(total)),
                "Assigned count": int(round(assigned)),
                "Unknown count": int(round(unknown)),
                "Assigned %": assigned / total if total else 0.0,
                "Perfect match %": perfect / total if total else 0.0,
                "1 mismatch %": one_mismatch / total if total else 0.0,
                ">=2 mismatches %": multi_mismatch / total if total else 0.0,
                "Mean edit distance": round(mean_distance, 3),
            }
        )

    return DataFrame(
        rows,
        columns=[
            "Read",
            "Total count",
            "Assigned count",
            "Unknown count",
            "Assigned %",
            "Perfect match %",
            "1 mismatch %",
            ">=2 mismatches %",
            "Mean edit distance",
        ],
    )


def _read_abundance_table(frame: DataFrame, read_label: str) -> DataFrame:
    summary = _assignment_summary(frame)
    summary = summary[summary["read"] == read_label].copy()
    if summary.empty:
        return DataFrame(
            columns=["sample", "barcode_type", "assigned_count", "share_of_read"]
        )

    abundance = summary.groupby(
        ["sample", "barcode_type"], dropna=False, as_index=False
    )["assigned_count"].sum()
    read_total = float(abundance["assigned_count"].sum())
    abundance["share_of_read"] = (
        abundance["assigned_count"] / read_total if read_total else 0.0
    )
    return abundance.sort_values(
        by=["assigned_count", "sample", "barcode_type"],
        ascending=[False, True, True],
    )


def _read_heatmap_table(frame: DataFrame, read_label: str) -> DataFrame:
    abundance = _read_abundance_table(frame, read_label)
    if abundance.empty:
        return DataFrame()

    heatmap = abundance.pivot(
        index="sample",
        columns="barcode_type",
        values="assigned_count",
    ).fillna(0)
    return heatmap.sort_index(axis=0).sort_index(axis=1)


def _styled_table_block(frame: DataFrame, *, bar_column: str | None = None) -> str:
    if frame.empty:
        return "<p>No rows available.</p>"

    styler = frame.style
    if bar_column and bar_column in frame.columns:
        styler = styler.bar(subset=[bar_column], color="#4c78a8")

    numeric_columns = [
        col for col in frame.columns if pd.api.types.is_numeric_dtype(frame[col])
    ]
    if numeric_columns:
        formatters: dict[str, str] = dict.fromkeys(numeric_columns, "{:,.0f}")
        if "share_of_read" in frame.columns:
            formatters["share_of_read"] = "{:.1%}"
        if "Assigned %" in frame.columns:
            formatters["Assigned %"] = "{:.1%}"
        if "Perfect match %" in frame.columns:
            formatters["Perfect match %"] = "{:.1%}"
        if "1 mismatch %" in frame.columns:
            formatters["1 mismatch %"] = "{:.1%}"
        if ">=2 mismatches %" in frame.columns:
            formatters[">=2 mismatches %"] = "{:.1%}"
        if "Mean edit distance" in frame.columns:
            formatters["Mean edit distance"] = "{:.3f}"
        styler = styler.format(formatters)  # type: ignore[arg-type]

    return styler.to_html()


def _heatmap_block(frame: DataFrame) -> str:
    if frame.empty:
        return "<p>No rows available.</p>"
    return frame.style.format("{:,.0f}").background_gradient(cmap="Blues").to_html()


def _sampleadapters_qmd(table_path: Path, frame: DataFrame) -> str:
    title = table_path.stem.replace("_", " ")
    metrics = _read_metrics_table(frame)
    assignments = _assignment_summary(frame)
    parts = [
        "---",
        f'title: "Sample Adapters QC Report: {title}"',
        "format:",
        "  html:",
        "    toc: true",
        "    code-fold: false",
        "execute:",
        "  enabled: false",
        "---",
        "",
        f"Source table: `{table_path}`",
        "",
        "## Summary",
        "",
        _table_block(metrics),
        "",
    ]

    if not assignments.empty:
        parts.extend(
            [
                "## Aggregated Assignments",
                "",
                _table_block(assignments.head(50)),
                "",
            ]
        )

    for read_label in _read_labels(frame):
        abundance = _read_abundance_table(frame, read_label)
        if abundance.empty:
            continue

        heatmap = _read_heatmap_table(frame, read_label)
        parts.extend(
            [
                f"## Barcode Abundance ({read_label})",
                "",
                _styled_table_block(abundance, bar_column="assigned_count"),
                "",
                f"### Sample by Barcode Type ({read_label})",
                "",
                _heatmap_block(heatmap),
                "",
            ]
        )

    parts.extend(
        [
            "## Top Observed Sequences",
            "",
            _table_block(_top_rows(frame)),
            "",
            "## Full Table",
            "",
            _table_block(frame),
            "",
        ]
    )

    return "\n".join(parts)


def render_quarto_html(qmd_text: str, outfile: Path) -> Path:
    outfile = outfile.resolve()
    parents(outfile)
    to_render = outfile.with_suffix(".report.qmd")
    with to_render.open("w", encoding="utf-8") as f:
        f.write(qmd_text)

        subprocess.run(
            [
                "quarto",
                "render",
                str(to_render),
                "--to",
                "html",
                "--output",
                outfile.name,
                "--output-dir",
                str(outfile.parent),
            ],
            check=True,
        )

    return outfile


def write_sampleadapters_report(
    table_path: Path,
    outfile: Path | None = None,
    summary_outfile: Path | None = None,
) -> Path:
    table_path = table_path.resolve()
    report_path = outfile if outfile is not None else table_path.with_suffix(".html")
    if summary_outfile is None:
        return report_path

    frame = read_frame(table_path, usecols=_sampleadapters_summary_usecols)
    # qmd_text = _sampleadapters_qmd(table_path, frame)
    # render_quarto_html(qmd_text, report_path)

    summary_frame = _assignment_summary(frame)
    write_frame(summary_frame, summary_outfile, index=False)

    return report_path
