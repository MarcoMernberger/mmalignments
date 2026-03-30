from __future__ import annotations

import hashlib
import json
import traceback
from enum import Enum
from functools import cached_property, wraps
from pathlib import Path
from subprocess import CompletedProcess
from types import MappingProxyType
from typing import Any, Callable, Iterable, Mapping, TypeVar, cast

from mmalignments.models.data import Pairing
from mmalignments.services.io import parents, paths_exists

from .registry import current_element_registry
from .tags import ElementTag, Method, PartialElementTag, Stage, State

# def file_sig(p: Path) -> dict[str, Any]:
#     if not p.exists():
#         return {"path": str(p), "missing": True}
#     st = p.stat()
#     return {"path": str(p), "mtime_ns": st.st_mtime_ns, "size": st.st_size}


def file_sig(p: Path, head_bytes: int = 65_536) -> dict[str, Any]:
    if not p.exists():
        return {"path": str(p), "missing": True}
    st = p.stat()
    # Content-hash of first 64kb plus file size for speed
    h = hashlib.sha256()
    with p.open("rb") as f:
        h.update(f.read(head_bytes))
    return {
        "path": str(p),
        "size": st.st_size,
        "head_sha256": h.hexdigest(),
    }


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def short_hash(signature: str, n: int = 8) -> str:
    return signature[:n]


class ValidationPolicy(Enum):
    CHECK = "check"  # default behaviour
    FORCE_RUN = "force_run"
    FORCE_SKIP = "force_skip"


class Element:
    def __init__(
        self,
        key: str,
        run: Callable[[], CompletedProcess | None | bool],
        tag: ElementTag,
        *,
        determinants: tuple[str, ...] | None = None,
        inputs: tuple[Path, ...] | None = None,
        artifacts: Mapping[str, Any] | None = None,
        validator: Callable[[], tuple[bool, str]] | None = None,
        pres: tuple["Element", ...] | None = None,
        empty_ok: bool = False,
        name: str | None = None,
    ):
        self.key = key
        self.run = run
        self.tag = tag
        self.threads = getattr(
            run, "threads", None
        )  # optional attribute for external runners
        if tag is None:
            raise ValueError("A tag must be provided, was None.")
        self.artifacts = MappingProxyType(dict(artifacts or {}))
        self.pres = tuple(pres or [])
        self.validator = validator
        self.determinants = tuple(str(det) for det in (determinants or ()))
        if self.determinants:
            self.key = f"{self.key}:" + self.determinants_as_str()
        self.inputs = tuple(sorted(inputs or [], key=lambda p: str(p)))
        self._signature: str | None = None  # cache
        # self._store_attributes = dict(store_attributes or {"multi_core": True})
        # self._init_store_attributes()
        self._validation_policy = ValidationPolicy.CHECK
        self.name = (
            name or self.default_name
        )  # supplied name overrided tag-based default name
        self.build_provenance()
        self.root = self.tag.root if self.tag else self.name
        self.empty_ok = empty_ok
        self.output_files  # trigger validation of output file paths
        self.validate_fields()
        self.creation_trace = "".join(traceback.format_stack()[:-1])

    def validate_fields(self) -> None:
        required_fields = ["key", "run", "tag"]
        for field in required_fields:
            if getattr(self, field) is None:
                raise ValueError(
                    f"Element '{self.name}' is missing required field: {field}."
                )
        if self.output_files is not None:
            for path in self.output_files:
                try:
                    assert isinstance(path, Path) or isinstance(path, str)
                except AssertionError:
                    raise AssertionError(
                        f"Output file '{path}' for {self.name} does not exist or is not a valid path."  # noqa: E501
                    )
        if self.inputs is not None:
            for path in self.inputs:
                try:
                    assert isinstance(path, Path) or (isinstance(path, str))
                except AssertionError:
                    raise AssertionError(
                        f"Input file '{path}' does not exist or is not a valid path."
                    )

    @property
    def validation_policy(self) -> ValidationPolicy:
        return self._validation_policy

    @validation_policy.setter
    def validation_policy(self, value: ValidationPolicy) -> None:
        self._validation_policy = value

    def force(self) -> Element:
        self.validation_policy = ValidationPolicy.FORCE_RUN
        return self

    @cached_property
    def output_files(self) -> Iterable[Path] | None:
        files = sorted(
            [v for v in self.artifacts.values() if isinstance(v, Path)], key=str
        )
        return files if files else None

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
        return self.tag.default_name

    @cached_property
    def default_output_file(self) -> str:
        return self.tag.default_output

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

    def sig_data(self) -> dict[str, Any]:
        sig_data = {
            "key": self.key,
            "determinants": tuple(self.determinants),
            "inputs": [file_sig(p) for p in self.inputs],
            "artifacts": self._artifact_sig(),
            "pre_sigs": sorted([pre.signature for pre in self.pres]),
        }
        return sig_data

    def _compute_signature(self) -> str:
        try:
            return stable_hash(self.sig_data())
        except Exception as e:
            print(e)
            raise RuntimeError(f"Failed to compute signature for {self.key!r}")

    def print_sig_data(self) -> None:
        sig_data = self.sig_data()
        print(json.dumps(sig_data, indent=2))

    # def invalidate(self) -> None:
    #     self._signature = None

    @cached_property
    def signature(self) -> str:
        return self._compute_signature()

    def outputs_ok(self) -> tuple[bool, str]:
        missing = []
        empty = []
        check = True
        for p in self.output_files or ():
            if not p.exists():
                missing.append(str(p))
                check = False
        if not check:
            return check, "Missing output files: " + ", ".join(missing)
        elif not self.empty_ok:
            for p in self.output_files or ():
                if p.stat().st_size == 0 and not self.empty_ok:
                    empty.append(str(p))
                    check = False
            if not check:
                return check, "Empty output files: " + ", ".join(empty)
        return True, "Output files are OK"

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

    def skip(self, cached_signature: str | None = None) -> tuple[bool, str]:
        # TODO turn this into Validation Policy
        if self.validation_policy == ValidationPolicy.FORCE_RUN:
            return False, "Validation policy forces run"

        elif self.validation_policy == ValidationPolicy.FORCE_SKIP:
            return True, "Validation policy forces skip"

        if cached_signature is None:
            return False, "First run"

        if cached_signature != self.signature:
            return False, "Cached signature does not match"

        check_output, reason = self.outputs_ok()
        if not check_output:
            return check_output, reason
        if self.validator:
            return self.validator()
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
        siglist = [f"{k}: {v}" for k, v in self.sig_data().items()]
        return (
            f"{self.__class__.__name__}:\n"
            f"  key:          {self.key}\n"
            f"  threads:      {self.threads}\n"
            f"  signature:    {self.signature}\n"
            f"  determinants: {', '.join(self.determinants)}\n"
            f"  inputs:       {[str(p) for p in self.inputs]}\n"
            f"  artifacts:    {{ {', '.join(artifactlist)}}}\n"
            f"  pres:         {[p.key for p in self.pres]}"
            f"  sigdata:      {{ {', '.join(siglist)}}}\n"
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


def get_candidates(
    arts: Mapping[str, Any],
    outputs: str | Iterable[str] | None = None,
    output_files: Iterable[Path] | None = None,
    auto_outputs: bool = True,
) -> list[Any]:
    if outputs is not None:
        keys = [outputs] if isinstance(outputs, str) else list(outputs)
        candidates = [arts.get(k) for k in keys]
    elif output_files and auto_outputs:
        candidates = list(
            output_files
        )  # your Element.output_files uses artifacts Paths
    else:
        candidates = []
    return candidates


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
            arts = e.artifacts
            output_files = e.output_files

            candidates = get_candidates(arts, outputs, output_files, auto_outputs)

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


class FilesElement(Element):

    def __init__(
        self,
        path: str | Path | Mapping[str, Path | str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        root: str | None = None,
        ext: str | None = None,
    ):

        paths = (
            path
            if isinstance(path, Mapping)
            else {Path(path).stem: Path(path).absolute()}
        )
        artifacts = {k: Path(v).absolute() for k, v in paths.items()}
        first = artifacts.values().__iter__().__next__()
        root = root or first.stem
        ext = ext or first.suffix.lstrip(".")
        tag = ElementTag(
            root=root,
            level=0,
            omics=None,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
            ext=ext,
        ).merge(tag)

        runner = paths_exists(
            *paths.values()
        )  # no-op runner, validation is done by skip logic
        runner.threads = 1
        runner.command = [f"check_exists({artifacts.values()})"]  # type: ignore[attr-defined]
        key = (
            f"{tag.default_name}_files_{'::'.join(str(p) for p in artifacts.values())}"
        )
        super().__init__(
            key=key,
            run=runner,
            tag=tag,
            validator=self.validate,
            inputs=tuple(artifacts.values()),
            artifacts=artifacts,
        )
        self.ext = tag.ext

    @cached_property
    def files(self) -> tuple[Path, ...]:
        return tuple(
            sorted([v for v in self.artifacts.values() if isinstance(v, Path)], key=str)
        )

    @cached_property
    def output_files(self) -> tuple[Path] | None:
        """
        overrides the Element output_files, the artifacts are not the output
        but the input files, so we return an empty tuple to avoid confusion
        """
        return None

    def validate(self) -> tuple[bool, str]:
        return True, "bypassed validation"
        md5sum = self.calc_md5sum()
        check = md5sum == self.md5sum
        return check, f"MD5 check {'passed' if check else 'failed'}"

    @cached_property
    def md5sum(self) -> str | None:
        return self.calc_md5sum()

    def calc_md5sum(self) -> str | None:
        md5 = ""
        for k, path in self.artifacts.items():
            if isinstance(path, Path) and path.exists():
                md5 += f"{k}:{hashlib.md5(path.read_bytes()).hexdigest()};"
        return md5


class FileElement(FilesElement):

    def __init__(
        self,
        filepath: Path | str,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        root: str | None = None,
    ):
        path = Path(filepath).absolute()
        ext = path.suffix.lstrip(".")
        by_suffix = {ext: path}
        super().__init__(by_suffix, root=root, tag=tag)


class Sample(FilesElement):

    def __init__(
        self,
        path: Path | str | Mapping[str, Path | str],
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        name: str | None = None,
    ):
        super().__init__(path, root=root, tag=tag)


class NextGenSampleElement(Sample):

    def __init__(
        self,
        path: Path | str | Mapping[str, Path | str],
        *,
        root: str,
        tag: PartialElementTag | ElementTag | None = None,
        read_group: str | None = None,
        reverse_reads: bool = False,
        cache_dir: Path | None = None,
        result_dir: Path | None = None,
        name: str | None = None,
    ):
        super().__init__(path, root=root, tag=tag)
        self.reverse_reads = reverse_reads
        self.read_group = read_group
        self.cache_dir = cache_dir or Path("cache") / "samples" / self.name
        self.result_dir = result_dir or Path("results") / "samples" / self.name
        self.pairing: Pairing = (
            Pairing.PAIRED if len(self.artifacts) > 1 else Pairing.SINGLE
        )

    @property
    def input_files(self) -> list[Path]:
        return sorted(self.artifacts.values(), key=str)


def sample_fastqs(
    sample: Sample | Element,
) -> tuple[Path, Path | None, str, str | None]:
    if isinstance(sample, Sample):
        r1 = Path(sample.input_files[0]).absolute()
        r2 = (
            Path(sample.input_files[1]).absolute()
            if len(sample.input_files) > 1
            else None
        )
        name = sample.name
        rg = getattr(sample, "read_group", None)
        return r1, r2, name, rg
    else:
        r1 = Path(sample.artifacts["fastq_r1"]).absolute()
        fastq_r2 = sample.artifacts.get("fastq_r2", None)
        r2 = Path(fastq_r2).absolute() if fastq_r2 else None
        name = sample.artifacts.get("sample_name", sample.name)
        rg = sample.artifacts.get("read_group")
        return r1, r2, name, rg


def generate_element_key_name(
    tag: ElementTag,
    tool_name: str,
    tool_version: str | None,
    subcommand: str | None = None,
    suffix: str | None = None,
    **params_str: Any,
) -> tuple[str, str]:
    """
    generate_element_key generates a unique key for an Element based on its tag 
    and other parameters. We want the keys to be somewhat informative for 
    debugging and provenance, but also unique and stable for caching.

    Parameters
    ----------
    tag : ElementTag
        The tag associated with the element.
    tool_version_name : str
        The name of the tool version.
    subcommand : str | None, optional
        The subcommand used, by default None
    suffix : str | None, optional
        The suffix to append to the key, by default None
    params : str | None, optional
        Additional parameters to include in the key, by default None

    Returns
    -------
    tuple[str, str]
        The generated element key.
    """
    parts = [tag.default_name, tool_name]
    if subcommand:
        parts.append(subcommand)
    short = "_".join(parts)
    if tool_version:
        parts.append(tool_version)
    if suffix:
        parts.append(suffix)
    for k, v in params_str.items():
        if v is not None:
            parts.append(f"{k}-{v}")
    key = "_".join(parts)
    ret = key, short
    return ret
