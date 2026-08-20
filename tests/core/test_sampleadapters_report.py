from pathlib import Path

import pandas as pd

from mmalignments.models.reports import quarto


def test_assignment_summary_groups_classified_reads_only():
    frame = pd.DataFrame(
        {
            "Assigned (R1)": ["S1", "S1", "unknown", "S2"],
            "Count (R1)": [10, 5, 7, 3],
            "Best (R1)": ["fw:abc", "fw:abc", "fw:abc", "rev:def"],
            "Matched barcode (R1)": ["AAA", "AAA", "AAA", "CCC"],
            "Edit distance (R1)": [0, 1, 0, 2],
            "Assigned (R2)": ["unknown", "S3", "S3", None],
            "Count (R2)": [4, 6, 2, 1],
            "Best (R2)": ["rev:xyz", "rev:xyz", "oops", "rev:xyz"],
            "Matched barcode (R2)": ["GGG", "GGG", "TTT", "GGG"],
            "Edit distance (R2)": [0, 0, 3, 0],
        }
    )

    result = quarto._assignment_summary(frame).reset_index(drop=True)
    expected = pd.DataFrame(
        [
            {
                "read": "R1",
                "sample": "S1",
                "barcode_type": "fw",
                "matched_barcode": "AAA",
                "assigned_count": 15,
                "mean_edit_distance": 0.5,
                "min_edit_distance": 0,
                "read_total": 18,
                "share_of_read": 15 / 18,
            },
            {
                "read": "R1",
                "sample": "S2",
                "barcode_type": "rev",
                "matched_barcode": "CCC",
                "assigned_count": 3,
                "mean_edit_distance": 2.0,
                "min_edit_distance": 2,
                "read_total": 18,
                "share_of_read": 3 / 18,
            },
            {
                "read": "R2",
                "sample": "S3",
                "barcode_type": "rev",
                "matched_barcode": "GGG",
                "assigned_count": 6,
                "mean_edit_distance": 0.0,
                "min_edit_distance": 0,
                "read_total": 8,
                "share_of_read": 6 / 8,
            },
            {
                "read": "R2",
                "sample": "S3",
                "barcode_type": "unknown",
                "matched_barcode": "TTT",
                "assigned_count": 2,
                "mean_edit_distance": 3.0,
                "min_edit_distance": 3,
                "read_total": 8,
                "share_of_read": 2 / 8,
            },
        ]
    )

    pd.testing.assert_frame_equal(result, expected)


def test_write_sampleadapters_report_skips_read_without_summary(monkeypatch, tmp_path):
    called = False

    def fake_read_frame(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("read_frame should not be called")

    monkeypatch.setattr(quarto, "read_frame", fake_read_frame)

    table_path = tmp_path / "input.tsv"
    report_path = quarto.write_sampleadapters_report(table_path)

    assert report_path == table_path.with_suffix(".html")
    assert called is False


def test_write_sampleadapters_report_reads_summary_columns_only(monkeypatch, tmp_path):
    source = pd.DataFrame(
        {
            "Assigned (R1)": ["S1", "unknown"],
            "Count (R1)": [10, 3],
            "Best (R1)": ["fw:abc", "fw:abc"],
            "Matched barcode (R1)": ["AAA", "AAA"],
            "Edit distance (R1)": [0, 1],
            "Irrelevant": ["keep", "out"],
        }
    )
    captured: dict[str, object] = {}

    def fake_read_frame(path: Path, **kwargs):
        usecols = kwargs["usecols"]
        captured["selected_columns"] = [
            column for column in source.columns if usecols(column)
        ]
        return source[captured["selected_columns"]]

    def fake_write_frame(frame, path: Path, index: bool = False):
        captured["written_frame"] = frame.reset_index(drop=True)
        captured["written_path"] = path
        captured["written_index"] = index

    monkeypatch.setattr(quarto, "read_frame", fake_read_frame)
    monkeypatch.setattr(quarto, "write_frame", fake_write_frame)

    table_path = tmp_path / "input.tsv"
    summary_path = tmp_path / "summary.tsv"
    report_path = tmp_path / "report.html"

    returned_path = quarto.write_sampleadapters_report(
        table_path,
        report_path,
        summary_path,
    )

    assert returned_path == report_path
    assert captured["selected_columns"] == [
        "Assigned (R1)",
        "Count (R1)",
        "Best (R1)",
        "Matched barcode (R1)",
        "Edit distance (R1)",
    ]
    pd.testing.assert_frame_equal(
        captured["written_frame"],
        quarto._assignment_summary(source).reset_index(drop=True),
    )
    assert captured["written_path"] == summary_path
    assert captured["written_index"] is False
