from __future__ import annotations

import hashlib
import json
from contextlib import contextmanager
from contextvars import ContextVar
from enum import Enum
from functools import cached_property, wraps
from pathlib import Path
from subprocess import CompletedProcess
from types import MappingProxyType
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Mapping,
    Optional,
    Protocol,
    TypeVar,
    cast,
)

from mmalignments.services.io import exists, parents

from .tags import ElementTag


def file_sig(p: Path) -> dict[str, Any]:
    if not p.exists():
        return {"path": str(p), "missing": True}
    st = p.stat()
    return {"path": str(p), "mtime_ns": st.st_mtime_ns, "size": st.st_size}


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def short_hash(signature: str, n: int = 8) -> str:
    return signature[:n]


# @dataclass(frozen=True)
# class ElementTag:
#     root: str  # sample
#     level: str  # No. in pipeline chain
#     omics: str  #
#     stage: str
#     method: str
#     state: str
#     ext: str | None = None  # e.g. "bam"

#     def suffix(self) -> str:
#         suffix = "[" + self.op
#         if self.tool and self.tool != "":
#             suffix += f",{self.tool}"
#         suffix += "]"
#         return suffix

#     def __repr__(self) -> str:
#         return f"OpTag(step='{self.step}', op='{self.op}', sig='{self.sig}', tool='
# {self.tool}', ext='{self.ext}')"


class ValidationPolicy(Enum):
    CHECK = "check"  # default behaviour
    FORCE_RUN = "force_run"
    FORCE_SKIP = "force_skip"


class Element:
    def __init__(
        self,
        key: str,
        tag: ElementTag,
        run: Callable[[], CompletedProcess],
        *,
        determinants: tuple[Any, ...] | None = None,
        inputs: tuple[Path, ...] | None = None,
        artifacts: Mapping[str, Any] | None = None,
        validator: Callable[[], tuple[bool, str]] | None = None,
        pres: tuple["Element", ...] | None = None,
        store_attributes: Mapping[str, Any] | None = None,
        name: str | None = None,
    ):
        self.key = key
        self.run = run
        self.tag = tag
        if tag is None:
            raise ValueError("A tag must be provided, was None.")
        self.artifacts = MappingProxyType(dict(artifacts or {}))
        self.pres = tuple(pres or [])
        self.validator = validator or (lambda: (True, "No validator, default to True"))
        self.determinants = tuple(str(det) for det in (determinants or ()))
        if self.determinants:
            self.key = f"{self.key}:" + self.determinants_as_str()
        self.inputs = tuple(sorted(inputs or [], key=lambda p: str(p)))
        self._signature: Optional[str] = None  # cache
        self._store_attributes = dict(store_attributes or {"multi_core": True})
        self._init_store_attributes()
        self._validation_policy = ValidationPolicy.CHECK
        self.name = (
            name or self.default_name
        )  # supplied name overrided tag-based default name
        self.build_provenance()
        self.root = self.tag.root if self.tag else self.name

    @property
    def validation_policy(self) -> ValidationPolicy:
        return self._validation_policy

    @validation_policy.setter
    def validation_policy(self, value: ValidationPolicy) -> None:
        self._validation_policy = value

    @property
    def output_files(self) -> list[Path]:
        return sorted(
            [v for v in self.artifacts.values() if isinstance(v, Path)], key=str
        )

    def _init_store_attributes(self):
        reserved = {
            "name",
            "key",
            "run",
            "validator",
            "determinants",
            "inputs",
            "artifacts",
            "pres",
        }  # do not overwrite by accident
        for attr, value in self._store_attributes.items():
            if attr in reserved:
                raise ValueError(
                    f"store_attributes key '{attr}' collides with reserved attribute"
                )
            setattr(self, attr, value)

    def __call__(self) -> Any:
        return self.run()

    def determinants_as_str(self) -> str:
        return ",".join(str(d) for d in self.determinants)

    @cached_property
    def default_name(self) -> str:
        return self.tag.default_name()

    @cached_property
    def default_output_file(self) -> str:
        return self.default_name + +(f".{self.tag.ext}" if self.tag.ext else "")

    def build_provenance(self) -> str:
        if not self.pres:
            prefix = ""
        elif len(self.pres) == 1:
            prefix = self.pres[0].provenance + "->"
        else:
            prefix = "(" + ",".join(p.provenance for p in self.pres) + ")->"
        return prefix + f"{self.name}"

    @cached_property
    def provenance(self) -> str:
        return self.build_provenance()

    def _compute_signature(self) -> str:
        sig_data = {
            "key": self.key,
            "determinants": tuple(self.determinants),
            "inputs": [file_sig(p) for p in self.inputs],
            "artifacts": self._artifact_sig(),
            "pre_sigs": sorted([pre.signature for pre in self.pres]),
        }
        return stable_hash(sig_data)

    def print_sig_data(self) -> None:
        sig_data = {
            "key": self.key,
            "determinants": tuple(self.determinants),
            "inputs": [file_sig(p) for p in self.inputs],
            "artifacts": self._artifact_sig(),
            "pre_sigs": sorted([pre.signature for pre in self.pres]),
        }
        print(json.dumps(sig_data, indent=2))

    # def invalidate(self) -> None:
    #     self._signature = None

    @cached_property
    def signature(self) -> str:
        return self._compute_signature()

    def outputs_ok(self) -> bool:
        return all(p.exists() and p.stat().st_size > 0 for p in self.output_files)

    def _artifact_sig(self) -> dict[str, Any]:
        sig: dict[str, Any] = {}
        for k, v in sorted(self.artifacts.items()):
            if isinstance(v, Path):
                sig[k] = str(v)
            elif isinstance(v, (str, int, float, bool)) or v is None:
                sig[k] = v
            else:
                raise TypeError(
                    f"artifact '{k}' has unsupported type for signature: {type(v)}"
                )
        return sig

    def skip(self, cached_signature: Optional[str] = None) -> tuple[bool, str]:
        if self.validation_policy == ValidationPolicy.FORCE_RUN:
            return False, "Validation policy forces run"

        if self.validation_policy == ValidationPolicy.FORCE_SKIP:
            return True, "Validation policy forces skip"

        if cached_signature is None:
            return False, "First run"

        if cached_signature != self.signature:
            return False, "Cached signature does not match"

        if not self.outputs_ok():
            return False, "Output files are not OK"
        if self.validator:
            check, reason = self.validator()
            if not check:
                return check, reason
        return True, "No relevant changes detected"

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"name='{self.name}', "
            f"key='{self.key}', "
            f"pres={len(self.pres)}, "
            f"inputs={len(self.inputs)}, "
            f"artifacts={list(self.artifacts.keys())}"
            f")"
        )

    def describe(self) -> str:
        artifactlist = [f"{k}: {v}" for k, v in self.artifacts.items()]
        return (
            f"{self.__class__.__name__}:\n"
            f"  key: {self.key}\n"
            f"  signature: {self.signature}\n"
            f"  determinants: {dict(self.determinants)}\n"
            f"  inputs: {[str(p) for p in self.inputs]}\n"
            f"  artifacts: {{ {', '.join(artifactlist)}}}\n"
            f"  pres: {[p.key for p in self.pres]}"
        )

    def __getattr__(self, name: str) -> Any:
        artifacts = self.__dict__.get("artifacts", {})
        if name in artifacts:
            return artifacts[name]
        raise AttributeError(f"{self.__class__.__name__} has no attribute '{name}'")

    def __hash__(self) -> int:
        return hash(self.key)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Element) and self.key == other.key


class MappedElement(Element):

    @property
    def bam(self) -> Path:
        v = self.artifacts["bam"]
        if not isinstance(v, Path):
            raise TypeError("artifact 'bam' must be a Path")
        return v

    # @property
    # def bai(self) -> Path:
    #     v = self.artifacts["bai"]
    #     if not isinstance(v, Path):
    #         raise TypeError("artifact 'bai' must be a Path")
    #     return v


class VcfElement(Element):
    """An Element that wraps a VCF/BCF file."""

    @property
    def vcf(self) -> Path:
        v = self.artifacts["vcf"]
        if not isinstance(v, Path):
            raise TypeError("artifact 'vcf' must be a Path")
        return v

    # This does not necessarily exist
    # @property
    # def tbi(self) -> Optional[Path]:
    #     """Path to the tabix index (.tbi) if present in artifacts."""
    #     v = self.artifacts.get("tbi")
    #     return Path(v) if v is not None else None
    #     return Path(v) if v is not None else None


F = TypeVar("F", bound=Callable[..., Element])


def _as_path(x: Any) -> Path | None:
    if isinstance(x, Path):
        return x
    if isinstance(x, str) and x.strip():
        return Path(x)
    return None


def _looks_like_filepath(p: Path) -> bool:
    s = str(p)
    # relativ oder absolut mit "/" oder Windows "\" -> likely a path
    if ("/" in s) or ("\\" in s):
        return True
    # oder hat eine Endung -> likely a file
    if p.suffix:
        return True
    return False


def element(
    fn: F | None = None,
    *,
    outputs: str | Iterable[str] | None = None,
    auto_outputs: bool = True,
) -> F | Callable[[F], F]:
    """
    Usable as:
      @element
      def builder(...): ...

    or:
      @element(outputs="bam")
      def builder(...): ...
    """

    def deco(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            e = func(*args, **kwargs)

            # intern (optional)
            reg = current_element_registry()
            if reg is not None:
                e = reg.intern(e)

            # mkdir parents for outputs (independent of registry)
            arts = getattr(e, "artifacts", None) or {}
            output_files = getattr(e, "output_files", [])

            if outputs is not None:
                keys = [outputs] if isinstance(outputs, str) else list(outputs)
                candidates = [arts.get(k) for k in keys]
            elif auto_outputs:
                candidates = list(
                    output_files
                )  # your Element.output_files uses artifacts Paths
            else:
                candidates = []

            out_paths: list[Path] = []
            for c in candidates:
                p = _as_path(c)
                if p is None:
                    continue
                if _looks_like_filepath(p):
                    out_paths.append(p)

            if out_paths:
                parents(*out_paths)

            return e

        return cast(F, wrapper)

    # called as @element
    if fn is not None:
        return deco(fn)

    # called as @element(...)
    return deco


@element
def FileElement(
    path: Path | str,
    name: str | None = None,
    attribute: str | None = None,
    tag: ElementTag | None = None,
) -> Element:
    """Helper to create an Element that represents a single file."""
    path = Path(path).absolute()
    stem = path.stem
    ext = path.suffix.lstrip(".")
    tag = tag or ElementTag(
        root=stem,
        level=0,
        omics=None,
        stage="file",
        method="",
        state="raw",
        ext=path.suffix,
    )
    attribute = attribute or ext  # e.g. "bam" for "sample.bam"
    runner = exists  # no-op runner, validation is done by skip logic
    runner.command = [f"check_exists({path})"]  # type: ignore[attr-defined]

    return Element(
        key=str(path),
        run=runner,
        tag=tag,
        name=name,  # name=name,
        inputs=[path],
        artifacts={attribute: path},
    )


class HasKey(Protocol):
    key: str


class HasSignature(Protocol):
    signature: str


T = TypeVar("T", bound=HasKey)

_current_registry: ContextVar["ElementRegistry[HasKey] | None"] = ContextVar(
    "current_element_registry", default=None
)


def current_element_registry() -> "ElementRegistry[HasKey] | None":
    return _current_registry.get()


@contextmanager
def element_build_context(registry: "ElementRegistry[HasKey]"):
    token = _current_registry.set(registry)
    try:
        yield
    finally:
        _current_registry.reset(token)


class ElementRegistry:
    """Intern registry: ensures one instance per key."""

    def __init__(self) -> None:
        self._by_key: Dict[str, HasKey] = {}

    def intern(self, e: HasKey) -> HasKey:
        existing = self._by_key.get(e.key)
        if existing is None:
            self._by_key[e.key] = e
            return e

        # optional sanity check, only if both have .signature
        ex_sig = getattr(existing, "signature", None)
        e_sig = getattr(e, "signature", None)
        if ex_sig is not None and e_sig is not None and ex_sig != e_sig:
            raise ValueError(f"Key collision with different signature: {e.key}")

        return existing
