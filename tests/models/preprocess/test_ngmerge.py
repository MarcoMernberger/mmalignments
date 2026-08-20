from __future__ import annotations

import importlib
import sys
import types
import typing
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = PROJECT_ROOT / "src"
MMALIGNMENTS_ROOT = SRC_ROOT / "mmalignments"

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

# elements.py decorates non-Protocol classes with runtime_checkable.
typing.runtime_checkable = lambda cls: cls


def _ensure_namespace_package(name: str, path: Path) -> None:
    """Create a lightweight namespace package when missing."""
    if name in sys.modules:
        return
    module = types.ModuleType(name)
    module.__path__ = [str(path)]
    sys.modules[name] = module


_ensure_namespace_package("mmalignments", MMALIGNMENTS_ROOT)
_ensure_namespace_package("mmalignments.models", MMALIGNMENTS_ROOT / "models")
_ensure_namespace_package(
    "mmalignments.models.preprocess", MMALIGNMENTS_ROOT / "models" / "preprocess"
)
_ensure_namespace_package("mmalignments.services", MMALIGNMENTS_ROOT / "services")

artifacts_module = importlib.import_module("mmalignments.models.artifacts")
overlay_module = importlib.import_module("mmalignments.models.overlay")
elements_module = importlib.import_module("mmalignments.models.elements")

# ngmerge imports OutputSpec from artifacts, but this repository sometimes
# provides OutputSpec in overlay.
if not hasattr(artifacts_module, "OutputSpec"):
    artifacts_module.OutputSpec = overlay_module.OutputSpec

# Some source variants miss this helper; inject a safe fallback for imports.
if not hasattr(elements_module, "generate_element_key_name"):
    elements_module.generate_element_key_name = lambda *a, **k: ("k", "n")

ngmerge_module = importlib.import_module("mmalignments.models.preprocess.ngmerge")


def _fake_external_init(self, **kwargs):
    """Tiny External.__init__ stand-in to avoid heavyweight setup."""
    print(kwargs)
    self._name = kwargs.get("name")
    self._primary_binary = kwargs.get("primary_binary")
    self._version = kwargs.get("version")
    self._source = kwargs.get("source")
    self.parameters = kwargs.get("parameters")


def test_build_param_registry_contains_expected_keys():
    """Ensure inlined ParamRegistry includes thread spec and CLI flags."""
    registry = ngmerge_module._build_param_registry()
    specs = registry.default.specs

    # Validate a representative subset plus thread spec presence.
    assert "_thread_spec" in specs
    assert "r1" in specs
    assert "r2" in specs
    assert "o" in specs
    assert "n" in specs


def test_build_param_registry_as_mapping_returns_plain_dict():
    """Validate helper returns plain dict consumed by External.__init__."""
    out = ngmerge_module._build_param_registry_as_mapping()
    assert isinstance(out, dict)
    assert out == {}


def test_constructor_uses_mapping_fallback_and_overwrites_param_registry(
    monkeypatch: pytest.MonkeyPatch,
):
    """Cover NGmerge.__init__ branches for parameters None vs explicit."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    monkeypatch.setattr(ngmerge_module, "_build_param_registry", lambda: "REG")
    monkeypatch.setattr(
        ngmerge_module, "_build_param_registry_as_mapping", lambda: {"m": 1}
    )

    tool = ngmerge_module.NGmerge(parameters=None)
    assert tool.parameters == {"m": 1}
    assert tool.param_registry == "REG"

    tool2 = ngmerge_module.NGmerge(parameters={"x": "y"})
    assert tool2.parameters == {"x": "y"}
    assert tool2.param_registry == "REG"


def test_get_version_returns_cached_version_without_invoking_binary(
    monkeypatch: pytest.MonkeyPatch,
):
    """If _version is already set, get_version should return it immediately."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(version="1.2.3")
    called = {"n": 0}
    monkeypatch.setattr(tool, "ensure_binary", lambda: called.__setitem__("n", 1))

    assert tool.get_version("fallback") == "1.2.3"
    assert called["n"] == 0


def test_get_version_returns_fallback_when_binary_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    """Cover path where primary binary is unavailable or not ensured."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(primary_binary="NGmerge", version=None)

    monkeypatch.setattr(tool, "ensure_binary", lambda: False)
    assert tool.get_version("fallback") == "fallback"

    tool._primary_binary = ""
    assert tool.get_version("fallback2") == "fallback2"


def test_get_version_extracts_semver_and_handles_failures(
    monkeypatch: pytest.MonkeyPatch,
):
    """Validate semver parsing plus failure/exception fallbacks."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(primary_binary="NGmerge", version=None)
    monkeypatch.setattr(tool, "ensure_binary", lambda: True)

    ok = types.SimpleNamespace(returncode=0, stdout="NGmerge 2.4.1\n")
    monkeypatch.setattr(ngmerge_module.subprocess, "run", lambda *a, **k: ok)
    assert tool.get_version("fallback") == "2.4.1"

    bad = types.SimpleNamespace(returncode=1, stdout="")
    monkeypatch.setattr(ngmerge_module.subprocess, "run", lambda *a, **k: bad)
    assert tool.get_version("fallback") == "fallback"

    def boom(*args, **kwargs):
        raise RuntimeError("x")

    monkeypatch.setattr(ngmerge_module.subprocess, "run", boom)
    assert tool.get_version("fallback") == "fallback"


def test_default_output_dir_and_spec(monkeypatch: pytest.MonkeyPatch):
    """Check default path building and extension switching by compression."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)

    tool = ngmerge_module.NGmerge(version="0.1.0")
    monkeypatch.setattr(type(tool), "version_name", property(lambda self: "ngmerge0"))

    outdir = tool.default_output_dir("S01")
    assert outdir == Path("results") / "ngmerge0" / "S01"

    raw = tool.default_output_spec("S01", compression="Raw")
    gz = tool.default_output_spec("S01", compression="Gzip")
    assert raw.ext == "fq"
    assert gz.ext == "fq.gz"


def test_merge_raises_when_sample_missing_r1_or_r2(
    monkeypatch: pytest.MonkeyPatch,
):
    """Validate explicit input-shape guards in merge() for sample object."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(version="1.0.0")
    tag = overlay_module.ElementTag(
        root="S01",
        level=0,
        stage=overlay_module.Stage.INPUT,
        method=overlay_module.Method.CHECK,
        state=overlay_module.State.RAW,
    )

    with pytest.raises(ValueError, match="must have 'r1'"):
        tool.merge(types.SimpleNamespace(tag=tag, root="S01"))

    with pytest.raises(ValueError, match="must have 'r2'"):
        tool.merge(types.SimpleNamespace(tag=tag, root="S01", r1=Path("a.fq")))


def test_merge_constructs_element_for_paired_and_single_end(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Exercise merge() wiring using the real overlay and artifact data classes."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)

    tool = ngmerge_module.NGmerge(version="1.0.0")
    monkeypatch.setattr(type(tool), "version_name", property(lambda self: "v1"))
    monkeypatch.setattr(tool, "run_ngmerge", lambda **kwargs: "RUNNER")

    tag = overlay_module.ElementTag(
        root="S01",
        level=0,
        stage=overlay_module.Stage.INPUT,
        method=overlay_module.Method.CHECK,
        state=overlay_module.State.RAW,
    )
    sample = types.SimpleNamespace(
        tag=tag,
        root="S01",
        r1=tmp_path / "r1.fq",
        r2=tmp_path / "r2.fq",
        pres=(),
    )

    paired = tool.merge(sample, mode="adapter-removal")
    assert paired.tag.method == overlay_module.Method.NGMERGE
    assert paired.tag.state == overlay_module.State.MERGED
    assert paired.run == "RUNNER"
    assert paired.inputs == (sample.r1, sample.r2)
    assert paired.determinants == ("adapter-removal",)
    assert paired.artifacts.primary_name == "fastq"
    assert isinstance(paired.artifacts.primary, ngmerge_module.FastqArtifact)

    sample_single = types.SimpleNamespace(
        tag=tag,
        root="S02",
        r1=tmp_path / "single.fq",
        r2=None,
        pres=(),
    )
    single = tool.merge(sample_single, mode="stitch")
    assert single.inputs == (sample_single.r1,)


def test_run_ngmerge_builds_arguments_for_modes_and_compression(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    """Verify command argument assembly for stitch/adapter and gzip/raw."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(version="1.0.0")

    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    out = tmp_path / "out.fq"

    runnable = tool.run_ngmerge(
        r1,
        r2,
        out,
        mode="stitch",
        compression="Gzip",
    )
    args1 = runnable.command
    assert args1[1:7] == ["-1", str(r1), "-2", str(r2), "-o", str(out.resolve())]
    assert "-z" in args1
    assert "-a" not in args1
    assert args1[2] == str(r1)
    assert args1[4] == str(r2)
    assert args1[6] == str(out.resolve())

    runnable2 = tool.run_ngmerge(
        r1,
        r2,
        out,
        mode="adapter-removal",
        compression="Raw",
    )
    args2 = runnable2.command
    assert "-z" not in args2
    assert "-a" in args2
    assert "-d" in args2


def test_get_version_invokes_subprocess_with_expected_arguments(
    monkeypatch: pytest.MonkeyPatch,
):
    """Check subprocess.run call contract for version probing."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(primary_binary="NGmerge", version=None)
    monkeypatch.setattr(tool, "ensure_binary", lambda: True)

    captured = {}

    def fake_run(cmd, capture_output, text, timeout):
        captured["cmd"] = cmd
        captured["capture_output"] = capture_output
        captured["text"] = text
        captured["timeout"] = timeout
        return types.SimpleNamespace(returncode=0, stdout="NGmerge 9.9.9")

    monkeypatch.setattr(ngmerge_module.subprocess, "run", fake_run)
    assert tool.get_version("fallback") == "9.9.9"
    assert captured["cmd"] == ["NGmerge", "--version"]
    assert captured["capture_output"] is True
    assert captured["text"] is True
    assert captured["timeout"] == 10


def test_get_version_returns_fallback_when_stdout_lacks_semver(
    monkeypatch: pytest.MonkeyPatch,
):
    """If output has no semver token, fallback should be returned."""
    monkeypatch.setattr(ngmerge_module.External, "__init__", _fake_external_init)
    tool = ngmerge_module.NGmerge(primary_binary="NGmerge", version=None)
    monkeypatch.setattr(tool, "ensure_binary", lambda: True)
    monkeypatch.setattr(
        ngmerge_module.subprocess,
        "run",
        lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="version unknown"),
    )

    assert tool.get_version("fallback") == "fallback"
