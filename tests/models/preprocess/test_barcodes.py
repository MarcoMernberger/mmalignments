from __future__ import annotations

import gzip
import importlib
import io
import sys
import types
import typing
from pathlib import Path

import pandas as pd
import pytest

# Keep imports isolated from package-level issues in this repository.
PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = PROJECT_ROOT / "src"
MMALIGNMENTS_ROOT = SRC_ROOT / "mmalignments"


def _ensure_namespace_package(name: str, path: Path) -> None:
    if name in sys.modules:
        return
    module = types.ModuleType(name)
    module.__path__ = [str(path)]
    sys.modules[name] = module


_ensure_namespace_package("mmalignments", MMALIGNMENTS_ROOT)
_ensure_namespace_package("mmalignments.models", MMALIGNMENTS_ROOT / "models")
_ensure_namespace_package("mmalignments.services", MMALIGNMENTS_ROOT / "services")

# elements.py applies runtime_checkable to non-Protocol classes.
typing.runtime_checkable = lambda cls: cls

artifacts_module = importlib.import_module("mmalignments.models.artifacts")
overlay_module = importlib.import_module("mmalignments.models.overlay")
elements_module = importlib.import_module("mmalignments.models.elements")

# barcodes.py imports these from artifacts/elements, so provide fallbacks when
# they are absent in current source layout.
if not hasattr(artifacts_module, "OutputSpec"):
    artifacts_module.OutputSpec = overlay_module.OutputSpec
if not hasattr(artifacts_module, "FileSpec"):
    artifacts_module.FileSpec = overlay_module.FileSpec
if not hasattr(elements_module, "generate_element_key_name"):
    elements_module.generate_element_key_name = lambda *args, **kwargs: ("k", "n")

barcodes_module = importlib.import_module("mmalignments.models.preprocess.barcodes")


class DummyRecord:
    """Minimal record-like object exposing .seq for SeqIO parser tests."""

    def __init__(self, seq: str):
        self.seq = seq


class DummyRunnable:
    """Simple Runnable stand-in that executes the stored callback."""

    def __init__(self, fn, display=None):
        self.fn = fn
        self.display = display

    def __call__(self):
        return self.fn()


class DummyCallSpec:
    """Lightweight replacement for CallSpec.render() behavior."""

    def __init__(self, path, kwargs):
        self.path = path
        self.kwargs = kwargs

    def render(self):
        return f"{self.path}:{sorted(self.kwargs)}"


class DummyTableArtifact:
    """Path wrapper used in high-level element-constructor tests."""

    def __init__(self, path: Path):
        self.path = path

    def resolve(self) -> Path:
        return self.path


class DummyArtifactSet:
    """Simple artifact collection with primary handle."""

    def __init__(self, primary, primary_name=None, **extras):
        self.primary = primary
        self.primary_name = primary_name
        self.extras = extras


class DummyOutputSpec:
    """OutputSpec stand-in that supports merge/path/add_output used here."""

    def __init__(self, stem: str, folder: Path, ext: str = "tsv"):
        self.stem = stem
        self.folder = Path(folder)
        self.ext = ext
        self.files: dict[str, Path] = {}

    def merge(self, outspec):
        return outspec or self

    def add_output(self, name: str, filespec):
        self.files[name] = self.folder / f"{filespec.stem}.{filespec.ext}"
        return self

    def path(self) -> Path:
        return self.folder / f"{self.stem}.{self.ext}"


class DummyFileSpec:
    """FileSpec stand-in used by DummyOutputSpec.add_output()."""

    def __init__(self, stem: str, ext: str):
        self.stem = stem
        self.ext = ext


class DummyElement:
    """Capture constructor arguments for high-level element factories."""

    def __init__(self, key, runner, **kwargs):
        self.key = key
        self.runner = runner
        for k, v in kwargs.items():
            setattr(self, k, v)


def test_input_fastq_name_handles_all_suffixes():
    """Ensure FASTQ name normalization strips known FASTQ suffixes only."""
    assert barcodes_module._input_fastq_name(Path("a.fastq.gz")) == "a"
    assert barcodes_module._input_fastq_name(Path("a.fq.gz")) == "a"
    assert barcodes_module._input_fastq_name(Path("a.fastq")) == "a"
    assert barcodes_module._input_fastq_name(Path("a.fq")) == "a"
    assert barcodes_module._input_fastq_name(Path("a.custom")) == "a"


def test_count_fastq_prefixes_plain_and_gzip_branches(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Cover both file-opening branches and verify sorting/count truncation."""
    plain = tmp_path / "reads.fastq"
    gz = tmp_path / "reads.fastq.gz"
    plain.write_text("", encoding="utf-8")
    with gzip.open(gz, "wt") as handle:
        handle.write("")

    records = [DummyRecord("AACCCC"), DummyRecord("aaGGGG"), DummyRecord("AACCAA")]
    monkeypatch.setattr(barcodes_module.SeqIO, "parse", lambda _h, _fmt: iter(records))

    plain_counts = barcodes_module._count_fastq_prefixes(
        plain,
        start_length=2,
        sample_size=2,
    )
    gz_counts = barcodes_module._count_fastq_prefixes(
        gz,
        start_length=2,
        sample_size=10,
    )

    # sample_size=2 keeps only first two records for plain path.
    assert plain_counts == [("AA", 2)]
    # gzip path sees all three records, sorted by count desc then sequence asc.
    assert gz_counts == [("AA", 3)]


def test_single_end_frame_builds_expected_columns():
    """Verify single-end frame schema and row content."""
    frame = barcodes_module._single_end_frame("sample", [("AAAA", 3), ("TTTT", 1)])
    assert list(frame.columns) == ["input", "Sequence (R1)", "Count (R1)"]
    assert frame.iloc[0].to_dict()["Sequence (R1)"] == "AAAA"


def test_paired_end_frame_handles_uneven_lengths():
    """Ensure paired-end frame pads missing side with None values."""
    frame = barcodes_module._paired_end_frame("sample", [("A", 3), ("B", 2)], [("X", 1)])
    assert list(frame.columns) == [
        "input",
        "Sequence (R1)",
        "Count (R1)",
        "Sequence (R2)",
        "Count (R2)",
    ]
    assert len(frame) == 2
    assert frame.iloc[1]["Sequence (R2)"] is None


def test_levenshtein_distance_exact_swap_and_early_break():
    """Cover equality path, unequal-length swap path, and max-distance break."""
    assert barcodes_module._levenshtein_distance("AAAA", "AAAA") == 0
    assert barcodes_module._levenshtein_distance("AAA", "AA") == 1
    assert barcodes_module._levenshtein_distance("AAAA", "TTTT", max_distance=1) == 2


def test_barcode_candidates_filters_empty_none_and_nan():
    """Ensure candidate extraction keeps only usable barcode values."""
    frame = pd.DataFrame(
        [
            {"sample": "S1", "b1": "acg", "b2": None},
            {"sample": None, "b1": "TTT", "b2": "CCC"},
            {"sample": "S2", "b1": "nan", "b2": " gga "},
        ]
    )
    out = barcodes_module._barcode_candidates(frame, "sample", ["b1", "b2"])
    assert out == [("S1", "b1", "ACG"), ("S2", "b2", "GGA")]


def test_classify_sequence_handles_empty_candidates():
    """When no candidates exist, sequence classification should return unknown tuple."""
    out = barcodes_module._classify_sequence("AAAA", [], 1, lambda *_: 0)
    assert out == ("unknown", "unknown", -1, "unknown", "", 0, "")


def test_classify_sequence_picks_best_and_builds_details(
    monkeypatch: pytest.MonkeyPatch,
):
    """Verify best-match selection, threshold assignment, and detail fields."""
    monkeypatch.setattr(barcodes_module, "longest_common_prefix", lambda a, b: 2)
    monkeypatch.setattr(barcodes_module, "match_pattern", lambda a, b: "xx..")

    candidates = [("S1", "c1", "TTTT"), ("S2", "c2", "AAAT")]

    def dist(seq, barcode, maxd):
        return 1 if barcode == "AAAT" else 4

    out = barcodes_module._classify_sequence("AAAA", candidates, 2, dist)
    assert out[0] == "S2"
    assert out[1] == "S2"
    assert out[2] == 1
    assert out[3] == "c2:S2:AAAT"
    assert out[4] == "AAAT"
    assert out[5] == 2
    assert out[6] == "xx.."


def test_classify_row_with_threshold_short_circuit_and_delegate(
    monkeypatch: pytest.MonkeyPatch,
):
    """Cover min-count guard and delegated classification path."""
    low = barcodes_module._classify_row_with_threshold(
        "AAAA",
        count=0,
        candidates=[("S", "c", "AAAA")],
        max_edit_distance=1,
        min_count=1,
        distance_func=lambda *_: 0,
    )
    assert low == ("unknown", "unknown", -1, "unknown", "", 0, "")

    monkeypatch.setattr(
        barcodes_module,
        "_classify_sequence",
        lambda *args, **kwargs: ("S", "S", 0, "c:S:AAAA", "AAAA", 4, "||||"),
    )
    high = barcodes_module._classify_row_with_threshold(
        "AAAA",
        count=3,
        candidates=[("S", "c", "AAAA")],
        max_edit_distance=1,
        min_count=1,
        distance_func=lambda *_: 0,
    )
    assert high[0] == "S"


def test_classify_column_maps_tuple_fields_to_series():
    """Ensure tuple results are split into dedicated pandas Series outputs."""
    frame = pd.DataFrame(
        [
            {"Sequence (R1)": "AA", "Count (R1)": 2},
            {"Sequence (R1)": "TT", "Count (R1)": None},
        ]
    )

    def dist(*_args):
        return 0

    out = barcodes_module._classify_column(
        frame,
        "Sequence (R1)",
        "Count (R1)",
        candidates=[("S", "c", "AA")],
        max_edit_distance=1,
        min_count=1,
        distance_func=dist,
    )
    assert len(out) == 7
    assert out[0].iloc[0] in ("S", "unknown")


def test_classify_single_end_frame_adds_assignment_columns(
    monkeypatch: pytest.MonkeyPatch,
):
    """Check single-end classification enriches frame with all annotation columns."""
    frame = pd.DataFrame([{"Sequence (R1)": "AA", "Count (R1)": 2}])
    barcodes = pd.DataFrame([{"sample": "S", "b": "AA"}])
    monkeypatch.setattr(barcodes_module, "_barcode_candidates", lambda *_: [("S", "b", "AA")])
    monkeypatch.setattr(
        barcodes_module,
        "_classify_column",
        lambda *args, **kwargs: tuple(pd.Series([x]) for x in ["S", "S", 0, "b:S:AA", "AA", 2, "||"]),
    )

    out = barcodes_module._classify_single_end_frame(
        frame,
        barcodes,
        "sample",
        ["b"],
        1,
        1,
        lambda *_: 0,
    )
    assert "Assigned (R1)" in out.columns
    assert "Match pattern (R1)" in out.columns


def test_classify_paired_frame_adds_both_r1_and_r2_columns(
    monkeypatch: pytest.MonkeyPatch,
):
    """Validate paired-end classification adds R1 and R2 annotation fields."""
    frame = pd.DataFrame(
        [{"Sequence (R1)": "AA", "Count (R1)": 2, "Sequence (R2)": "TT", "Count (R2)": 2}]
    )
    barcodes = pd.DataFrame([{"sample": "S", "b1": "AA", "b2": "TT"}])
    monkeypatch.setattr(barcodes_module, "_barcode_candidates", lambda *_: [("S", "b1", "AA")])
    monkeypatch.setattr(
        barcodes_module,
        "_classify_column",
        lambda *args, **kwargs: tuple(pd.Series([x]) for x in ["S", "S", 0, "b:S:AA", "AA", 2, "||"]),
    )

    out = barcodes_module._classify_paired_frame(
        frame,
        barcodes,
        "sample",
        ["b1"],
        ["b2"],
        1,
        1,
        lambda *_: 0,
    )
    assert "Assigned (R1)" in out.columns
    assert "Assigned (R2)" in out.columns


def test_format_hit_example_and_flank_queries(monkeypatch: pytest.MonkeyPatch):
    """Cover hit formatting and query-direction fanout, including validation."""
    assert barcodes_module._format_hit_example("fwd", 4, 2, "ACGT") == "fwd@4:d2:ACGT"

    monkeypatch.setattr(barcodes_module, "reverse_complement", lambda s: "TGCA")
    assert barcodes_module._flank_queries("ACGT", "forward") == {"fwd": "ACGT"}
    assert barcodes_module._flank_queries("ACGT", "reverse") == {"rev": "TGCA"}
    assert barcodes_module._flank_queries("ACGT", "both") == {"fwd": "ACGT", "rev": "TGCA"}
    with pytest.raises(ValueError, match="direction must be"):
        barcodes_module._flank_queries("ACGT", "invalid")


def test_flank_uniqueness_radius_details_paths(monkeypatch: pytest.MonkeyPatch):
    """Cover no-hit, cumulative break, and the single-distance uniqueness branch."""
    monkeypatch.setattr(barcodes_module, "reverse_complement", lambda s: s)

    # No-hit branch via max_distance filter.
    monkeypatch.setattr(barcodes_module, "hamming_early_break", lambda q, w, m: 5)
    radius, examples = barcodes_module._flank_uniqueness_radius_details(
        "AAAA",
        "AAAATTTT",
        max_distance=1,
    )
    assert radius == -1
    assert examples == []

    # Cumulative break branch: two exact hits make uniqueness radius -1.
    monkeypatch.setattr(barcodes_module, "hamming_early_break", lambda q, w, m: 0)
    radius2, _examples2 = barcodes_module._flank_uniqueness_radius_details(
        "AA",
        "AAAA",
        max_distance=2,
    )
    assert radius2 == -1

    # Unique-radius branch: in forward-only direction, a single match bucket keeps the
    # radius at the observed distance and records the example.
    def scripted_hamming(query, window, max_distance):
        return 1 if window == "AC" else 5

    monkeypatch.setattr(barcodes_module, "hamming_early_break", scripted_hamming)
    radius3, examples3 = barcodes_module._flank_uniqueness_radius_details(
        "AA",
        "AC",
        max_distance=3,
        direction="forward",
    )
    assert radius3 == 1
    assert isinstance(examples3, list)
    assert examples3 and examples3[0].startswith("fwd@0:d1:")


def test_sample_start_reads_single_end_and_write_frame(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Validate single-end sampling pipeline and writer invocation."""
    monkeypatch.setattr(barcodes_module, "Runnable", DummyRunnable)
    monkeypatch.setattr(barcodes_module, "CallSpec", DummyCallSpec)
    monkeypatch.setattr(barcodes_module, "_input_fastq_name", lambda _p: "sample")
    monkeypatch.setattr(barcodes_module, "_count_fastq_prefixes", lambda *args, **kwargs: [("AA", 2)])
    monkeypatch.setattr(
        barcodes_module,
        "_single_end_frame",
        lambda name, counts: pd.DataFrame([{"input": name, "Sequence (R1)": "AA", "Count (R1)": 2}]),
    )

    wrote = {}
    monkeypatch.setattr(
        barcodes_module,
        "write_frame",
        lambda frame, path, **kwargs: wrote.update({"frame": frame, "path": path, "kwargs": kwargs}),
    )

    runner = barcodes_module.sample_start_reads(
        tmp_path / "r1.fastq",
        None,
        tmp_path / "out.tsv",
    )
    out = runner()
    assert isinstance(out, pd.DataFrame)
    assert wrote["path"].name == "out.tsv"
    assert wrote["kwargs"]["index"] is False


def test_sample_start_reads_paired_with_barcodes_and_column_fallbacks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Cover paired-end classification with barcode-based column fallback logic."""
    monkeypatch.setattr(barcodes_module, "Runnable", DummyRunnable)
    monkeypatch.setattr(barcodes_module, "CallSpec", DummyCallSpec)
    monkeypatch.setattr(barcodes_module, "_input_fastq_name", lambda _p: "sample")

    def fake_count(path, **kwargs):
        return [("AA", 2)] if "r1" in path.name else [("TT", 2)]

    monkeypatch.setattr(barcodes_module, "_count_fastq_prefixes", fake_count)
    monkeypatch.setattr(
        barcodes_module,
        "_paired_end_frame",
        lambda name, r1, r2: pd.DataFrame(
            [{"input": name, "Sequence (R1)": "AA", "Count (R1)": 2, "Sequence (R2)": "TT", "Count (R2)": 2}]
        ),
    )
    monkeypatch.setattr(barcodes_module, "read_frame", lambda _p: pd.DataFrame([{"sample": "S", "b": "AA"}]))

    captured = {}

    def fake_classify(frame, barcodes, sample_col, start_cols, end_cols, *rest, **kwargs):
        captured["start"] = start_cols
        captured["end"] = end_cols
        return frame.assign(dummy=1)

    monkeypatch.setattr(barcodes_module, "_classify_paired_frame", fake_classify)
    monkeypatch.setattr(barcodes_module, "write_frame", lambda *args, **kwargs: None)

    runner = barcodes_module.sample_start_reads(
        tmp_path / "r1.fastq",
        tmp_path / "r2.fastq",
        tmp_path / "out.tsv",
        barcode_path=tmp_path / "barcodes.tsv",
        barcode_columns_r1=None,
        barcode_columns_r2=None,
    )
    _ = runner()
    assert captured["start"] == []
    assert captured["end"] == []


def test_flank_uniqueness_radius_delegates_to_details(monkeypatch: pytest.MonkeyPatch):
    """Ensure public wrapper returns first item from detailed helper result."""
    monkeypatch.setattr(
        barcodes_module,
        "_flank_uniqueness_radius_details",
        lambda *args, **kwargs: (7, ["x"]),
    )
    assert barcodes_module.flank_uniqueness_radius("AA", "AAAA") == 7


def test_barcode_check_uniqueness_radius_enriches_and_writes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Validate barcode uniqueness table annotation and final write call."""
    frame = pd.DataFrame({"start": ["AA"], "end": ["TT"]})
    monkeypatch.setattr(barcodes_module, "read_frame", lambda _p: frame.copy())
    monkeypatch.setattr(
        barcodes_module.SeqIO,
        "read",
        lambda *_args, **_kwargs: types.SimpleNamespace(seq="AATTGG"),
    )

    def fake_details(flank, genomic, max_hamming, direction):
        return (2, [f"hit:{flank}"])

    monkeypatch.setattr(barcodes_module, "_flank_uniqueness_radius_details", fake_details)
    wrote = {}
    monkeypatch.setattr(
        barcodes_module,
        "write_frame",
        lambda out_frame, out_path: wrote.update({"frame": out_frame, "path": out_path}),
    )

    out = barcodes_module.barcode_check_uniqueness_radius(
        tmp_path / "barcodes.tsv",
        tmp_path / "genome.fa",
        tmp_path / "out.tsv",
        3,
        ["start", "end"],
    )
    assert "Unique radius (start)" in out.columns
    assert "Non-exact hit examples (end)" in out.columns
    assert wrote["path"].name == "out.tsv"


def test_uniqueness_radius_returns_runnable_invoking_checker(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Cover uniqueness_radius factory and the callback forwarding arguments."""
    monkeypatch.setattr(barcodes_module, "Runnable", DummyRunnable)
    monkeypatch.setattr(barcodes_module, "CallSpec", DummyCallSpec)

    called = {}

    def fake_checker(*args):
        called["args"] = args
        return "done"

    monkeypatch.setattr(barcodes_module, "barcode_check_uniqueness_radius", fake_checker)
    runner = barcodes_module.uniqueness_radius(
        tmp_path / "barcodes.tsv",
        tmp_path / "genome.fa",
        tmp_path / "out.tsv",
        4,
        ["start"],
        direction="reverse",
    )
    assert runner() == "done"
    assert called["args"][3] == 4


def test_check_flank_code_returns_hits_and_uniqueness(monkeypatch: pytest.MonkeyPatch):
    """Verify hit collection and unique-flag behavior for flank scanning."""

    def dist(flank, window, max_distance):
        return 0 if flank == window else 2

    monkeypatch.setattr(barcodes_module, "hamming_early_break", dist)
    hits, unique = barcodes_module.check_flank_code("AA", "AATAA", max_hamming=0)
    assert 0 in hits and 3 in hits
    assert unique is False

    hits2, unique2 = barcodes_module.check_flank_code("AT", "GGATGG", max_hamming=0)
    assert len(hits2) == 1
    assert unique2 is True


def test_uniqueness_element_construction(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Exercise uniqueness element factory wiring using the real element model."""
    monkeypatch.setattr(barcodes_module, "uniqueness_radius", lambda *args, **kwargs: "runner")

    source_path = tmp_path / "barcodes.tsv"
    genomic_path = tmp_path / "genome.fa"
    source_path.write_text("start\nAAA\n", encoding="utf-8")
    genomic_path.write_text(">chr\nAAATTT\n", encoding="utf-8")

    source = elements_module.FileSource(source_path, tag=overlay_module.ElementTag(
        root="sample",
        level=0,
        stage=overlay_module.Stage.INPUT,
        method=overlay_module.Method.CHECK,
        state=overlay_module.State.RAW,
    ))
    source.pres = ()
    genomic = elements_module.FileSource(genomic_path, tag=overlay_module.ElementTag(
        root="genome",
        level=0,
        stage=overlay_module.Stage.INPUT,
        method=overlay_module.Method.CHECK,
        state=overlay_module.State.RAW,
    ))
    genomic.pres = ()

    out = barcodes_module.uniqueness(
        source,
        genomic,
        max_hamming=3,
        barcode_columns=["start"],
        direction="forward",
    )

    assert isinstance(out, elements_module.Element)
    assert out.run == "runner"
    assert out.pres == ()
    assert out.determinants == ("start", "3")
    assert out.artifacts.primary_name == "tsv"


def test_sampleadapters_element_construction_with_optional_barcodes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Validate sampleadapters wiring, determinants, and prereq aggregation."""
    monkeypatch.setattr(barcodes_module, "sample_start_reads", lambda *args, **kwargs: "runner2")

    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    r2.write_text("@r\nTGCA\n+\nIIII\n", encoding="utf-8")
    barcode_path = tmp_path / "barcodes.tsv"
    barcode_path.write_text("sample\tstart_barcode\nS1\tACGT\n", encoding="utf-8")

    source = elements_module.FastqSource(
        artifacts_module.FastqArtifact(r1, r2),
        tag=overlay_module.ElementTag(
            root="sample",
            level=0,
            stage=overlay_module.Stage.INPUT,
            method=overlay_module.Method.CHECK,
            state=overlay_module.State.RAW,
        ),
    )
    source.pres = (source,)

    barcode_src = elements_module.FileSource(
        barcode_path,
        tag=overlay_module.ElementTag(
            root="barcodes",
            level=1,
            stage=overlay_module.Stage.PREP,
            method=overlay_module.Method.CUSTOM,
            state=overlay_module.State.PREPROCESS,
        ),
    )
    barcode_src.pres = (barcode_src,)

    out = barcodes_module.sampleadapters(
        source,
        barcodes=barcode_src,
        start_length=12,
        sample_size=100,
        max_edit_distance=1,
        min_count=2,
        leven=True,
    )

    assert isinstance(out, elements_module.Element)
    assert out.run == "runner2"
    assert out.pres == (source, barcode_src)
    assert "12" in out.determinants and "100" in out.determinants


def test_sampleadapters_report_element_and_runner_callback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Ensure report element builds outputs and callback calls report writer."""
    called = {}

    def fake_report(table_path, report_path, summary_path):
        called["table"] = table_path
        called["report"] = report_path
        called["summary"] = summary_path
        return report_path

    monkeypatch.setattr(barcodes_module, "write_sampleadapters_report", fake_report)

    table_path = tmp_path / "in.tsv"
    table_path.write_text("sample\tstart\nS1\tAAA\n", encoding="utf-8")
    source = elements_module.Element(
        "source-ref",
        lambda: None,
        tag=overlay_module.ElementTag(
            root="sample",
            level=0,
            stage=overlay_module.Stage.INPUT,
            method=overlay_module.Method.CHECK,
            state=overlay_module.State.RAW,
        ),
        artifacts=artifacts_module.ArtifactSet(artifacts_module.TableArtifact(table_path)),
        pres=(),
    )

    out = barcodes_module.sampleadapters_report(source)

    assert isinstance(out, elements_module.Element)
    assert out.pres == (source,)
    assert called == {}

    result = out.runner()
    assert result == called["report"]
    assert called["table"] == table_path
