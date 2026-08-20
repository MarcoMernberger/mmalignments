from __future__ import annotations

import builtins
import gzip
import importlib
import io
import json
import subprocess
import sys
import types
from pathlib import Path

import pytest

# Mirror the repository's isolated import strategy so tests are resilient even
# when unrelated package imports are temporarily broken.
PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = PROJECT_ROOT / "src"
MMALIGNMENTS_ROOT = SRC_ROOT / "mmalignments"


def _ensure_namespace_package(name: str, path: Path) -> None:
    """Create a lightweight namespace package in sys.modules when absent."""
    if name in sys.modules:
        return
    module = types.ModuleType(name)
    module.__path__ = [str(path)]
    sys.modules[name] = module


_ensure_namespace_package("mmalignments", MMALIGNMENTS_ROOT)
_ensure_namespace_package("mmalignments.models", MMALIGNMENTS_ROOT / "models")
_ensure_namespace_package("mmalignments.services", MMALIGNMENTS_ROOT / "services")

io_module = importlib.import_module("mmalignments.services.io")


class DummyColumns:
    """Small helper that mimics pandas Index.tolist()."""

    def __init__(self, values):
        self._values = list(values)

    def tolist(self):
        return list(self._values)


class DummyFrameLike:
    """Tiny frame-like object with .columns for schema tests."""

    def __init__(self, columns):
        self.columns = DummyColumns(columns)


class RecordingFrame:
    """Record writer calls so write_frame behavior can be asserted."""

    def __init__(self):
        self.calls: list[tuple[str, Path, dict]] = []

    def to_csv(self, path: Path, **kwargs):
        self.calls.append(("csv", path, dict(kwargs)))

    def to_parquet(self, path: Path, **kwargs):
        self.calls.append(("parquet", path, dict(kwargs)))

    def to_excel(self, path: Path, **kwargs):
        self.calls.append(("excel", path, dict(kwargs)))


class ResolvablePath:
    """Mimic artifact-like objects that expose resolve()."""

    def __init__(self, path: Path):
        self.path = path

    def resolve(self) -> Path:
        return self.path.resolve()


def test_ensure_aggregates_results_from_all_paths(monkeypatch: pytest.MonkeyPatch):
    """Ensure ensure() calls ensure_path for every input and AND-aggregates results."""
    seen = []

    def fake_ensure_path(path):
        seen.append(path)
        return "ok" in str(path)

    monkeypatch.setattr(io_module, "ensure_path", fake_ensure_path)
    result = io_module.ensure("ok_a", "bad_b", "ok_c")

    # The reducer must keep processing all paths and return final conjunction.
    assert seen == ["ok_a", "bad_b", "ok_c"]
    assert result is False


def test_parents_forwards_parent_directories(monkeypatch: pytest.MonkeyPatch):
    """Validate parents() maps input files to their parent folders before ensure()."""
    captured = {}

    def fake_ensure(*paths):
        captured["paths"] = paths
        return True

    monkeypatch.setattr(io_module, "ensure", fake_ensure)

    out = io_module.parents(Path("/tmp/a/b.txt"), "/tmp/c/d.tsv")
    assert out is True
    assert captured["paths"] == (Path("/tmp/a"), Path("/tmp/c"))


def test_ensure_path_success_and_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Cover successful mkdir and exception-handling branch in ensure_path()."""
    target = tmp_path / "x" / "y"
    assert io_module.ensure_path(target) is True
    assert target.exists()

    def boom(self, parents=False, exist_ok=False):
        raise OSError("nope")

    monkeypatch.setattr(Path, "mkdir", boom)
    assert io_module.ensure_path(tmp_path / "z") is False


def test_open_fastq_uses_gzip_for_gz_and_open_for_plain(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Verify file opener dispatch based on extension for FASTQ inputs."""
    calls = []

    def fake_gzip_open(path, mode):
        calls.append(("gzip", Path(path), mode))
        return io.StringIO("gz")

    def fake_open(path, mode):
        calls.append(("open", Path(path), mode))
        return io.StringIO("plain")

    monkeypatch.setattr(gzip, "open", fake_gzip_open)
    monkeypatch.setattr(builtins, "open", fake_open)

    with io_module.open_fastq(tmp_path / "reads.fq.gz") as f1:
        assert f1.read() == "gz"
    with io_module.open_fastq(tmp_path / "reads.fq") as f2:
        assert f2.read() == "plain"

    assert calls == [
        ("gzip", tmp_path / "reads.fq.gz", "rt"),
        ("open", tmp_path / "reads.fq", "r"),
    ]


def test_write_json_persists_dictionary(tmp_path: Path):
    """Ensure write_json writes an indented JSON object to disk."""
    path = tmp_path / "result.json"
    payload = {"a": "1", "b": "2"}
    io_module.write_json(path, payload)

    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload


def test_write_fastq_check_results_renders_all_sections(tmp_path: Path):
    """Exercise report writer including checks, errors, and warnings blocks."""
    path = tmp_path / "check.txt"
    results = {
        "status": "FAIL",
        "checks": {"paired": True, "read_count": 10},
        "errors": ["missing R2"],
        "warnings": ["low complexity"],
    }
    io_module.write_fastq_check_results(path, "sampleA", "paired", results)
    text = path.read_text(encoding="utf-8")

    assert "FASTQ Input Check for sampleA" in text
    assert "Status: FAIL" in text
    assert "Pairing: paired" in text
    assert "paired: True" in text
    assert "- missing R2" in text
    assert "- low complexity" in text


def test_from_json_success_and_missing_file(tmp_path: Path):
    """Validate normal JSON load and FileNotFoundError path."""
    path = tmp_path / "params.json"
    path.write_text('{"x": 1, "y": "ok"}', encoding="utf-8")

    loaded = io_module.from_json(path)
    assert loaded["x"] == 1
    assert loaded["y"] == "ok"

    with pytest.raises(FileNotFoundError, match="JSON file not found"):
        io_module.from_json(tmp_path / "missing.json")


def test_load_param_json_success_and_error_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Cover valid object load, top-level-type validation, and wrapped missing-file error."""
    path = tmp_path / "spec.json"
    path.write_text('{"tool": {"param": "value"}}', encoding="utf-8")
    assert io_module.load_param_json(path) == {"tool": {"param": "value"}}

    monkeypatch.setattr(io_module, "from_json", lambda p: ["not", "object"])
    with pytest.raises(ValueError, match="Top-level JSON must be an object"):
        io_module.load_param_json(path)

    def raise_missing(_):
        raise FileNotFoundError("inner")

    monkeypatch.setattr(io_module, "from_json", raise_missing)
    with pytest.raises(FileNotFoundError, match="Param spec JSON not found"):
        io_module.load_param_json(path)


def test_open_target_none_path_append_and_file_object(tmp_path: Path):
    """Test open_target() behavior for None, path targets, append mode, and existing streams."""
    assert io_module.open_target(None, append=False) is None

    target = tmp_path / "out" / "log.txt"
    with io_module.open_target(target, append=False, encoding="utf-8") as f:
        f.write("first")
    with io_module.open_target(target, append=True, encoding="utf-8") as f:
        f.write("+second")
    assert target.read_text(encoding="utf-8") == "first+second"

    stream = io.StringIO()
    returned = io_module.open_target(stream, append=False)
    # Existing file-like objects should be returned untouched.
    assert returned is stream


def test_absolutize_converts_all_paths_to_absolute(tmp_path: Path):
    """Ensure absolutize() returns absolute Path values in input order."""
    rel = Path(".")
    abs_path = tmp_path.resolve()
    out = io_module.absolutize(rel, abs_path)
    assert isinstance(out, tuple)
    assert all(p.is_absolute() for p in out)


def test_paths_exists_and_paths_exists_raise(tmp_path: Path):
    """Cover closure factories that check path existence and raise on missing entries."""
    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    p1.write_text("a", encoding="utf-8")

    check_all = io_module.paths_exists(p1, ResolvablePath(p1), p2)
    assert check_all() is False

    p2.write_text("b", encoding="utf-8")
    assert check_all() is True

    check_raise = io_module.paths_exists_raise(p1, p2)
    check_raise()
    p2.unlink()
    with pytest.raises(FileNotFoundError, match="Required path does not exist"):
        check_raise()


def test_exists_factory_returns_callable_result(tmp_path: Path):
    """Verify exists() returns a callable that reflects filesystem state."""
    path = tmp_path / "flag.txt"
    checker = io_module.exists(path)
    assert checker() is False
    path.write_text("ok", encoding="utf-8")
    assert checker() is True


def test_write_fasta_and_write_fasta_from_list(tmp_path: Path):
    """Ensure both FASTA writers produce valid records and create parent directories."""
    p1 = tmp_path / "fa" / "one.fa"
    io_module.write_fasta(p1, {"s1": "AT", "s2": "CG"})
    text1 = p1.read_text(encoding="utf-8")
    assert ">s1\nAT\n" in text1
    assert ">s2\nCG\n" in text1

    p2 = tmp_path / "fa" / "two.fa"
    # zip() semantics: extra names are ignored if sequences are shorter.
    io_module.write_fasta_from_list(p2, ["n1", "n2"], ["AAA"])
    text2 = p2.read_text(encoding="utf-8")
    assert text2 == ">n1\nAAA\n"


def test_write_frames_delegates_to_write_frame(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Check write_frames() simply forwards each path to write_frame()."""
    seen = []

    def fake_write_frame(df, path, **kwargs):
        seen.append((df, path, dict(kwargs)))

    monkeypatch.setattr(io_module, "write_frame", fake_write_frame)
    frame = object()
    p1 = tmp_path / "a.tsv"
    p2 = tmp_path / "b.tsv"
    io_module.write_frames(frame, [p1, p2], index=False)
    assert seen == [(frame, p1, {"index": False}), (frame, p2, {"index": False})]


@pytest.mark.parametrize("ext", [".tsv", ".csv", ".txt"])
def test_write_frame_text_formats_use_to_csv(tmp_path: Path, ext: str):
    """Verify text-like formats route to to_csv with expected defaults."""
    frame = RecordingFrame()
    out = tmp_path / f"table{ext}"
    io_module.write_frame(frame, out)

    call = frame.calls[-1]
    assert call[0] == "csv"
    assert call[1] == out
    assert call[2]["sep"] == "\t"
    assert call[2]["index"] is True


def test_write_frame_parquet_and_excel_strip_sep(tmp_path: Path):
    """Ensure parquet/xls(x) writers remove text separator from params."""
    frame = RecordingFrame()

    pq_path = tmp_path / "table.parquet"
    io_module.write_frame(frame, pq_path, custom=1)
    kind1, _path1, kwargs1 = frame.calls[-1]
    assert kind1 == "parquet"
    assert "sep" not in kwargs1
    assert kwargs1["custom"] == 1

    xl_path = tmp_path / "table.xlsx"
    io_module.write_frame(frame, xl_path, custom=2)
    kind2, _path2, kwargs2 = frame.calls[-1]
    assert kind2 == "excel"
    assert "sep" not in kwargs2
    assert kwargs2["custom"] == 2


def test_write_frame_unsupported_extension_raises(tmp_path: Path):
    """Cover explicit ValueError for unsupported frame output format."""
    frame = RecordingFrame()
    with pytest.raises(ValueError, match="Unsupported file format"):
        io_module.write_frame(frame, tmp_path / "table.unknown")


def test_read_frame_routes_by_extension_and_default_usecols(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Validate read_frame() branch routing and unnamed-column filter default."""
    calls = []
    sentinel = object()

    def fake_read_csv(path, **kwargs):
        calls.append(("csv", Path(path), dict(kwargs)))
        return sentinel

    def fake_read_parquet(path, **kwargs):
        calls.append(("parquet", Path(path), dict(kwargs)))
        return sentinel

    def fake_read_excel(path, **kwargs):
        calls.append(("excel", Path(path), dict(kwargs)))
        return sentinel

    monkeypatch.setattr(io_module.pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(io_module.pd, "read_parquet", fake_read_parquet)
    monkeypatch.setattr(io_module.pd, "read_excel", fake_read_excel)

    assert io_module.read_frame(tmp_path / "a.tsv") is sentinel
    assert io_module.read_frame(tmp_path / "b.parquet", drop_unnamed_columns=False) is sentinel
    assert io_module.read_frame(tmp_path / "c.csv") is sentinel
    assert io_module.read_frame(tmp_path / "d.xlsx") is sentinel

    # For TSV/TXT/CSV branches, usecols should be defaulted unless provided.
    tsv_call = calls[0]
    assert tsv_call[0] == "csv"
    assert tsv_call[2]["sep"] == "\t"
    assert callable(tsv_call[2]["usecols"])
    assert tsv_call[2]["usecols"]("col") is True
    assert tsv_call[2]["usecols"]("Unnamed: 0") is False

    with pytest.raises(ValueError, match="Unsupported file format"):
        io_module.read_frame(tmp_path / "bad.bin")


def test_read_schema_routes_by_extension(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Cover read_schema() for tsv/parquet/csv/xlsx and unsupported extension."""
    monkeypatch.setattr(
        io_module.pd,
        "read_csv",
        lambda path, **kwargs: DummyFrameLike(["a", "b"]),
    )
    monkeypatch.setattr(
        io_module.pd,
        "read_excel",
        lambda path, **kwargs: DummyFrameLike(["x", "y"]),
    )
    monkeypatch.setattr(
        io_module.pq,
        "read_schema",
        lambda path: types.SimpleNamespace(names=["p", "q"]),
    )

    assert io_module.read_schema(tmp_path / "a.tsv") == ["a", "b"]
    assert io_module.read_schema(tmp_path / "b.parquet") == ["p", "q"]
    assert io_module.read_schema(tmp_path / "c.csv") == ["a", "b"]
    assert io_module.read_schema(tmp_path / "d.xlsx") == ["x", "y"]

    with pytest.raises(ValueError, match="Unsupported file format"):
        io_module.read_schema(tmp_path / "unknown.ext")


def test_concat_files_combines_inputs_in_order(tmp_path: Path):
    """Ensure concat_files preserves input order when writing combined output."""
    in1 = tmp_path / "in1.txt"
    in2 = tmp_path / "in2.txt"
    out = tmp_path / "out.txt"
    in1.write_text("a\n", encoding="utf-8")
    in2.write_text("b\n", encoding="utf-8")

    io_module.concat_files(out, in1, in2)
    assert out.read_text(encoding="utf-8") == "a\nb\n"


def test_get_paths_from_prefix_path_filters_by_stem_prefix(tmp_path: Path):
    """Verify prefix scanner returns only matching stems mapped to absolute paths."""
    base = tmp_path / "sample.tsv"
    (tmp_path / "sample.tsv").write_text("x", encoding="utf-8")
    (tmp_path / "sample.extra.txt").write_text("y", encoding="utf-8")
    (tmp_path / "other.tsv").write_text("z", encoding="utf-8")

    found = io_module.get_paths_from_prefix_path(base)
    assert "sample" in found
    assert "sample.extra" in found
    assert "other" not in found
    assert all(path.is_absolute() for path in found.values())


def test_concat_fastq_builds_cat_command_and_writes_binary_output(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Cover concat_fastq command construction and output-directory creation."""
    inp1 = tmp_path / "r1.fq"
    inp2 = tmp_path / "r2.fq"
    inp1.write_bytes(b"A")
    inp2.write_bytes(b"B")
    out = tmp_path / "nested" / "combined.fq"

    calls = {}

    def fake_run(cmd, stdout, check):
        calls["cmd"] = cmd
        calls["check"] = check
        # Simulate shell cat output to prove stdout is a writable binary stream.
        stdout.write(b"AB")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(io_module.subprocess, "run", fake_run)
    io_module.concat_fastq((inp1, inp2), out)

    assert calls["cmd"] == ["cat", str(inp1), str(inp2)]
    assert calls["check"] is True
    assert out.read_bytes() == b"AB"
