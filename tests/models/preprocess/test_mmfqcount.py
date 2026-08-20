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

# Some source variants miss this helper; inject a safe fallback for imports.
if not hasattr(elements_module, "generate_element_key_name"):
    elements_module.generate_element_key_name = lambda *a, **k: ("k", "n")

# The compatibility shim should expose OutputSpec/FileSpec from the public module.
if not hasattr(artifacts_module, "OutputSpec"):
    artifacts_module.OutputSpec = overlay_module.OutputSpec
if not hasattr(artifacts_module, "FileSpec"):
    artifacts_module.FileSpec = overlay_module.FileSpec

mmfqcount_module = importlib.import_module("mmalignments.models.preprocess.mmfqcount")


def _fake_external_init(self, **kwargs):
    """Tiny External.__init__ stand-in to avoid heavyweight setup."""
    self._name = kwargs.get("name")
    self._primary_binary = kwargs.get("primary_binary")
    self._version = kwargs.get("version")
    self._source = kwargs.get("source")
    self.parameters = kwargs.get("parameters")


def test_build_param_registry_contains_expected_keys():
    """Ensure the inlined registry exposes the expected subcommand-specific flags."""
    registry = mmfqcount_module._build_param_registry()

    assert "trim_start" in registry.by_subcommand["count"].specs
    assert "trim_length" in registry.by_subcommand["count"].specs
    assert "seq_col" in registry.by_subcommand["match"].specs
    assert registry.by_subcommand["match"].specs["id_col"].default == "Name"


def test_build_param_registry_as_mapping_returns_plain_dict():
    """Validate the helper keeps the external init path simple."""
    out = mmfqcount_module._build_param_registry_as_mapping()
    assert isinstance(out, dict)
    assert out == {}


def test_constructor_uses_mapping_fallback_and_overwrites_param_registry(
    monkeypatch: pytest.MonkeyPatch,
):
    """Cover MmFqCount.__init__ branches for parameters None vs explicit."""
    monkeypatch.setattr(mmfqcount_module.External, "__init__", _fake_external_init)
    monkeypatch.setattr(mmfqcount_module, "_build_param_registry", lambda: "REG")
    monkeypatch.setattr(
        mmfqcount_module,
        "_build_param_registry_as_mapping",
        lambda: {"trim_start": "ACGT"},
    )

    tool = mmfqcount_module.MmFqCount(parameters=None)
    assert tool.parameters == {"trim_start": "ACGT"}
    assert tool.param_registry == "REG"

    tool2 = mmfqcount_module.MmFqCount(parameters={"x": "y"})
    assert tool2.parameters == {"x": "y"}
    assert tool2.param_registry == "REG"


def test_get_version_returns_cached_version_without_invoking_binary(
    monkeypatch: pytest.MonkeyPatch,
):
    """If _version is already set, get_version should return it immediately."""
    monkeypatch.setattr(mmfqcount_module.External, "__init__", _fake_external_init)
    tool = mmfqcount_module.MmFqCount(version="1.2.3")
    called = {"n": 0}
    monkeypatch.setattr(tool, "ensure_binary", lambda: called.__setitem__("n", 1))

    assert tool.get_version("fallback") == "1.2.3"
    assert called["n"] == 0


def test_get_version_returns_fallback_when_binary_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    """Cover the path where the binary is unavailable or not ensured."""
    monkeypatch.setattr(mmfqcount_module.External, "__init__", _fake_external_init)
    tool = mmfqcount_module.MmFqCount(primary_binary="mmfqcount", version=None)

    monkeypatch.setattr(tool, "ensure_binary", lambda: False)
    assert tool.get_version("fallback") == "fallback"

    tool._primary_binary = ""
    assert tool.get_version("fallback2") == "fallback2"


def test_get_version_extracts_semver_and_handles_failures(
    monkeypatch: pytest.MonkeyPatch,
):
    """Validate semver parsing plus failure/exception fallbacks."""
    monkeypatch.setattr(mmfqcount_module.External, "__init__", _fake_external_init)
    tool = mmfqcount_module.MmFqCount(primary_binary="mmfqcount", version=None)
    monkeypatch.setattr(tool, "ensure_binary", lambda: True)

    ok = types.SimpleNamespace(returncode=0, stdout="mmfqcount 2.4.1\n")
    monkeypatch.setattr(mmfqcount_module.subprocess, "run", lambda *a, **k: ok)
    assert tool.get_version("fallback") == "2.4.1"

    bad = types.SimpleNamespace(returncode=1, stdout="")
    monkeypatch.setattr(mmfqcount_module.subprocess, "run", lambda *a, **k: bad)
    assert tool.get_version("fallback") == "fallback"

    def boom(*args, **kwargs):
        raise RuntimeError("x")

    monkeypatch.setattr(mmfqcount_module.subprocess, "run", boom)
    assert tool.get_version("fallback") == "fallback"
