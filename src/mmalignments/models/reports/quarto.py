from __future__ import annotations

import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

from pandas import DataFrame  # type: ignore[import]

from mmalignments.services.io import read_frame


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


def _sampleadapters_qmd(table_path: Path, frame: DataFrame) -> str:
    title = table_path.stem.replace("_", " ")
    parts = [
        "---",
        f'title: "Sample Adapters Report: {title}"',
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
        _table_block(_summary_table(frame)),
        "",
    ]

    for heading, summary in _classification_tables(frame):
        parts.extend(
            [
                f"## {heading}",
                "",
                _table_block(summary),
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
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(prefix="mmalignments-quarto-") as tmp_dir:
        tmp_path = Path(tmp_dir)
        qmd_path = tmp_path / "report.qmd"
        qmd_path.write_text(qmd_text, encoding="utf-8")

        subprocess.run(
            [
                "quarto",
                "render",
                str(qmd_path),
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


def write_sampleadapters_report(table_path: Path, outfile: Path | None = None) -> Path:
    table_path = table_path.resolve()
    report_path = outfile if outfile is not None else table_path.with_suffix(".html")
    frame = read_frame(table_path)
    qmd_text = _sampleadapters_qmd(table_path, frame)
    return render_quarto_html(qmd_text, report_path)
