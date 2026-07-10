from __future__ import annotations

import hashlib
import json
import traceback
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property, wraps
from pathlib import Path
from subprocess import CompletedProcess
from typing import (
    Any,
    Callable,
    Iterable,
    Mapping,
    ParamSpec,
    TypeAlias,
    TypeVar,
    cast,
    overload,
)

import pandas as pd
from numpy import ndarray
from pandas import DataFrame, Series

from mmalignments.models.data import Pairing
from mmalignments.services.dependencies import (
    collect_code_dependency,
    file_sig,
    function_hash,
    stable_hash,
)
from mmalignments.services.io import parents, paths_exists

from .artifacts import Artifact, ArtifactSet, TableArtifact, TransientArtifact
from .registry import current_element_registry
from .tags import (
    ElementTag,
    Method,
    Omics,
    PartialElementTag,
    Stage,
    State,
)


class Runnable:
    def __init__(
        self,
        fn: Callable[[], Any],
        cmd: list[str] | None = None,
        display: str | None = None,
    ):
        self._fn = fn
        self.command_display = display
        self.command = cmd or [display]
        self.last_result: CompletedProcess | None = None
        self._fingerprint = self._compute_fingerprint()

    @property
    def fingerprint(self) -> str:
        return self._fingerprint

    def __call__(self) -> Any:
        result = self._fn()
        self.last_result = result
        return result

    def __name__(self) -> str:
        return getattr(self._fn, "__name__", repr(self._fn))

    def _compute_fingerprint(self) -> str:
        all_fns = collect_code_dependency(self._fn)

        safe_fns = (fn for fn in all_fns if fn is not None)

        hashes = [function_hash(fn) for fn in safe_fns]

        return hashlib.sha256("|".join(hashes).encode()).hexdigest()


@dataclass(frozen=True)
class CallSpec:
    path: tuple[str, ...]
    args: tuple[Any, ...] = field(default_factory=tuple)
    kwargs: dict[str, Any] = field(default_factory=dict)

    def render(self) -> str:
        callargs = ", ".join(repr(a) for a in self.args)
        callargs += ", ".join(f"{k}={repr(v)}" for k, v in self.kwargs.items())
        return f"{'.'.join(self.path)}({callargs})"


################################################################################
# Aliases
################################################################################

RunType: TypeAlias = Callable[[], CompletedProcess | None | bool | Any] | Runnable
FilesSourceType: TypeAlias = Path | str | Mapping[str, Path | str] | RunType

###############################################################################
# Elements
###############################################################################


class ValidationPolicy(Enum):
    CHECK = "check"  # default behaviour
    FORCE_RUN = "force_run"
    FORCE_SKIP = "force_skip"


class Element:
    def __init__(
        self,
        key: str,
        run: RunType,
        tag: ElementTag,
        *,
        determinants: tuple[str, ...] | None = None,
        inputs: tuple[Path, ...] | None = None,
        artifacts: Mapping[str, Any] | ArtifactSet | None = None,
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
        self.artifacts = ArtifactSet.from_any(artifacts)
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
            [v.resolve() for v in self.artifacts.values() if hasattr(v, "resolve")],
            key=str,
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
        if isinstance(self.run, Runnable):
            sig_data["run_hash"] = self.run.fingerprint

        return sig_data

    def _compute_signature(self) -> str:
        try:
            return stable_hash(self.sig_data())
        except Exception as e:
            print(e)
            print(self.sig_data())
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
            if isinstance(v, TransientArtifact):
                continue
            if isinstance(v, Path):
                sig[k] = str(v)
            elif isinstance(v, (str, int, float, bool)) or v is None:
                sig[k] = v
            elif isinstance(v, Artifact):
                sig[k] = v.signature()
            else:
                raise TypeError(
                    f"artifact '{k}' has unsupported type for signature: {type(v)}"
                )
        return sig

    def _explain_signature_diff(self, cached_sig_data: dict[str, Any] | None) -> str:
        """Compare current sig_data against a previously stored one and return
        a human-readable diff.  When *cached_sig_data* is ``None`` or the two
        dicts share no keys the method falls back to a simple mismatch message.
        """
        current = self.sig_data()
        if not cached_sig_data:
            return "Cached signature does not match (no cached sig_data available)"
        lines: list[str] = []
        all_keys = sorted(set(current) | set(cached_sig_data))
        for k in all_keys:
            old_v = cached_sig_data.get(k, "<missing>")
            new_v = current.get(k, "<missing>")
            if stable_hash(old_v) != stable_hash(new_v):
                lines.append(f"  {k}: {old_v!r} → {new_v!r}")
        if lines:
            return "Signature mismatch in:\n" + "\n".join(lines)
        # hashes of individual keys matched but overall hash differs (ordering/encoding)
        return "Cached signature does not match (unknown diff)"

    def skip(
        self,
        cached_signature: str | None = None,
        cached_sig_data: dict[str, Any] | None = None,
    ) -> tuple[bool, str]:
        # TODO turn this into Validation Policy
        if self.validation_policy == ValidationPolicy.FORCE_RUN:
            return False, "Validation policy forces run"

        elif self.validation_policy == ValidationPolicy.FORCE_SKIP:
            return True, "Validation policy forces skip"

        if cached_signature is None:
            return False, "First run"

        if cached_signature != self.signature:
            return False, self._explain_signature_diff(cached_sig_data)

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
        return (
            f"{self.__class__.__name__}:\n"
            f"  key:          {self.key}\n"
            f"  threads:      {self.threads}\n"
            f"  determinants: {', '.join(self.determinants)}\n"
            f"  inputs:       {[str(p) for p in self.inputs]}\n"
            f"  artifacts:    {{ {', '.join(artifactlist)}}}\n"
            f"  pres:         {[p.key for p in self.pres]}\n"
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

    @property
    def primary(self) -> Any:
        return self.artifacts.primary

    @property
    def file(self) -> Path | None:
        if "primary" in self.artifacts:
            return self.artifacts["primary"].resolve()
        if "parquet" in self.artifacts:
            return self.artifacts["parquet"]
        elif "tsv" in self.artifacts:
            return self.artifacts["tsv"]
        return next(iter(self.output_files), None) if self.output_files else None


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


P = ParamSpec("P")
R = TypeVar("R", bound=Element)


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


@overload
def element(fn: Callable[P, R]) -> Callable[P, R]: ...


@overload
def element(
    *,
    outputs: str | Iterable[str] | None = None,
    auto_outputs: bool = True,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def element(
    fn: Callable[P, R] | None = None,
    *,
    outputs: str | Iterable[str] | None = None,
    auto_outputs: bool = True,
) -> Callable[P, R] | Callable[[Callable[P, R]], Callable[P, R]]:
    """
        Usable as:
          @element
          def builder(...): ...
    Callable[[], CompletedProcess | None | bool | Any] | Runnable
        or:
          @element(outputs="bam")
          def builder(...): ...
    """

    def deco(func: Callable[P, R]) -> Callable[P, R]:
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            e = func(*args, **kwargs)

            # intern (optional)
            reg = current_element_registry()
            if reg is not None:
                e = cast(R, reg.intern(e))

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

        return wrapper

    # called as @element
    if fn is not None:
        return deco(fn)

    # called as @element(...)
    return deco


class FilesElement(Element):

    def __init__(
        self,
        path: Path | str | Mapping[str, Path | TableArtifact],
        *,
        # runner: Runnable | None = None,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        ext: str | None = None,
        is_prefix: bool = False,
        pres: tuple[Element, ...] | None = None,
    ):

        if is_prefix and isinstance(path, (str, Path)):
            artifacts = self.artifacts_from_prefix(path)
        else:
            artifacts = (
                path
                if isinstance(path, Mapping)
                else {Path(path).stem: Path(path).resolve()}
            )
        inputs = tuple(v.resolve() for v in artifacts.values() if hasattr(v, "resolve"))
        first = artifacts.values().__iter__().__next__().resolve()
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
        runner = self.exist
        key = f"{tag.default_name}::{'::'.join(str(p) for p in artifacts.values())}"

        super().__init__(
            key=key,
            run=runner,
            tag=tag,
            validator=self.validate,
            inputs=inputs,
            artifacts=artifacts,
            pres=pres,
        )
        self.ext = tag.ext

    @classmethod
    def artifacts_from_prefix(cls, prefix: str | Path) -> dict[str, Path]:
        artifacts = {}
        p = Path(prefix)
        file_dir = p.parent
        file_prefix = p.stem
        for p in file_dir.iterdir():
            if p.stem.startswith(file_prefix):
                artifacts[p.stem] = p.resolve()

        return artifacts

    @cached_property
    def files(self) -> tuple[Path, ...]:
        return tuple(
            sorted(
                [v.resolve() for v in self.artifacts.values() if hasattr(v, "resolve")],
                key=str,
            )
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

    def exist(self):
        return Runnable(
            paths_exists(*self.artifacts.values()),
            display=CallSpec(
                path=("io", "paths_exists"), args=tuple(self.artifacts.values())
            ).render(),
        )  # no-op runner, validation is done by skip logic


class FileElement(FilesElement):

    def __init__(
        self,
        source: Path | str | TableArtifact,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        root: str | None = None,
        pres: tuple[Element, ...] | None = None,
    ):
        if isinstance(source, str):
            source = Path(source)

        self.ext = source.resolve().suffix.lstrip(".")
        by_suffix = {self.ext: source}
        super().__init__(by_suffix, root=root, tag=tag, pres=pres)

    @property
    def file(self) -> Path:
        print("FileElement.file called, artifacts:", self.artifacts)
        if "primary" in self.artifacts:
            return self.artifacts["primary"].resolve()
        else:
            return self.artifacts[self.ext].resolve()

    @classmethod
    def resolve_path(cls, path: Path | str | TableArtifact) -> Path:
        if isinstance(path, str):
            path = Path(path)
        return path.resolve()


class Sample(FilesElement):
    #

    def __init__(
        self,
        path: Path | str | Mapping[str, Path | str],
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        is_prefix: bool = False,
        pres: tuple[Element, ...] | None = None,
    ):
        super().__init__(path, root=root, tag=tag, is_prefix=is_prefix, pres=pres)


class NextGenSampleElement(Sample):

    def __init__(
        self,
        path: Path | str | Mapping[str, Path | str],
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        read_group: str | None = None,
        reverse_reads: bool = False,
        cache_dir: Path | None = None,
        result_dir: Path | None = None,
        is_prefix: bool = False,
        pres: tuple[Element, ...] | None = None,
    ):
        if isinstance(path, (str, Path)):
            root = root or Path(path).stem
        else:
            first = next(iter(path.values()))
            root = root or Path(first).stem
        tag = ElementTag(
            root=root,
            level=0,
            omics=Omics.DNA,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
        ).merge(tag)
        super().__init__(path, root=root, tag=tag, is_prefix=is_prefix, pres=pres)
        self.reverse_reads = reverse_reads
        self.read_group = read_group
        self.cache_dir = cache_dir or Path("cache") / "samples" / self.name
        self.result_dir = result_dir or Path("results") / "samples" / self.name
        self.pairing: Pairing = (
            Pairing.PAIRED if len(self.artifacts) > 1 else Pairing.SINGLE
        )
        self.name = self.tag.default_name or root

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


@element
def register(element: Element):
    return element


# class TableElement(Element):

#     def __init__(
#         self,
#         key: str,
#         run: RunType,
#         tag: ElementTag,
#         *,
#         determinants: tuple | None = None,
#         inputs: tuple[Path, ...] | None = None,
#         artifacts: Mapping[str, Any] | None = None,
#         pres: tuple[Element, ...] | None = None,
#         name: str | None = None,
#     ):

#         super().__init__(
#             key=key,
#             run=run,
#             tag=tag,
#             determinants=determinants,
#             inputs=inputs,
#             artifacts=artifacts,
#             pres=pres,
#             name=name,
#         )

#     @property
#     def table(self) -> Path:
#         return self.artifacts.primary


# deprectated - this will in future be an TableArtifact
class TableElement(Element):
    """An Element that produces a TSV (for humans) and a Parquet file (for the pipeline).

    Parameters
    ----------
    key : str
        Unique cache key.
    run : Callable
        Zero-argument callable that executes the computation and writes both files.
    tag : ElementTag
        Tag used for naming / provenance.
    tsv : Path
        Output TSV path.
    parquet : Path
        Output Parquet path.
    column_roles : Mapping[str, str] | None
        Optional semantic role per column, e.g.::

            {
                "raw_A": "raw_counts",
                "raw_B": "raw_counts",
                "vst_A": "vst",
                "gene_desc": "annotation",
            }

        Well-known roles: ``"raw_counts"``, ``"cpm"``, ``"vst"``, ``"rlog"``,
        ``"log2fc"``, ``"fdr"``, ``"p"``, ``"stat"``, ``"annotation"``,
        ``"metadata"``.
    determinants, inputs, pres, name
        Forwarded to :class:`Element`.

    Artefacts
    ---------
    ``"tsv"``     → *tsv*
    ``"parquet"`` → *parquet*

    The pipeline should prefer ``element.parquet`` for chaining.
    Users / downstream tools consume ``element.tsv``.

    Examples
    --------
    raw = TableElement(..., column_roles={
    "KO_1": "raw_counts", "KO_2": "raw_counts",
    "WT_1": "raw_counts", "WT_2": "raw_counts",
    "gene_name": "annotation", "gene_desc": "annotation",
    })

    # VST — Annotationens get propagated automatically via propagate()
    vst = deseq2.vst(raw, sample_columns=["KO_1","KO_2","WT_1","WT_2"])
    # vst.column_roles == {"KO_1": "vst", ..., "gene_name": "annotation", ...}

    # Views
    vst.view("vst")          # only VST columns
    vst.view("annotation")   # only annotation columns
    vst.matrix("vst", include_annotations=True)  # VST + annotation columns together

    # Multiple roles at once
    vst.roles("vst", "annotation")

    # FFor PCA — only the matrix
    clustering.pca(vst.matrix("vst"))
    """

    # Roles that are considered *orthogonal* annotations and are propagated
    # automatically when a derived TableElement is created via
    # :meth:`propagate`.
    ANNOTATION_ROLES: frozenset[str] = frozenset({"annotation", "metadata"})

    def __init__(
        self,
        key: str,
        run: Runnable | Callable[[], CompletedProcess] | Path | str | None = None,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        root: str | None = None,
        tsv: Path | None = None,
        parquet: Path | None = None,
        column_roles: Mapping[str, str] | None = None,
        determinants: tuple | None = None,
        inputs: tuple[Path, ...] | None = None,
        pres: tuple[Element, ...] | None = None,
        artifacts: Mapping[str, Any] | None = None,
        name: str | None = None,
        index: str | None = None,
    ) -> None:
        # --- resolve file-based construction ---
        if isinstance(run, (Path, str)):
            # run is actually a file path
            p = Path(run).absolute()
            if p.suffix in (".parquet",):
                parquet = parquet or p
            elif p.suffix in (".tsv", ".csv", ".txt"):
                tsv = tsv or p
            else:
                raise ValueError(
                    f"Cannot infer file type from suffix {p.suffix!r}. "
                    "Pass tsv= or parquet= explicitly."
                )
            actual_run = paths_exists(p)
            actual_run.threads = 1
            actual_run.command = [f"check_exists({p})"]  # type: ignore[attr-defined]
            inputs = inputs or (p,)
            default_root = p.stem

        elif run is None:
            # tsv or parquet must be supplied explicitly
            if parquet is None and tsv is None:
                raise ValueError(
                    "Provide at least one of tsv= or parquet= when run=None."
                )
            paths_to_check = [p for p in (tsv, parquet) if p is not None]
            actual_run = paths_exists(*paths_to_check)
            actual_run.threads = 1
            actual_run.command = ["check_exists(...)"]  # type: ignore[attr-defined]
            inputs = inputs or tuple(paths_to_check)
            default_root = paths_to_check[0].stem

        else:
            actual_run = run
            default_root = root or (name or key)

        if tsv is None and parquet is None:
            raise ValueError("At least one of tsv= or parquet= must be provided.")

        # sentinel so we can tell callers "no TSV / no parquet"
        self._tsv = Path(tsv) if tsv is not None else None
        self._parquet = Path(parquet) if parquet is not None else None
        self.index_column = index
        self.column_roles: dict[str, str] = dict(column_roles or {})

        artifacts: dict[str, Any] = artifacts or {}
        if self._tsv:
            artifacts["tsv"] = self._tsv
        if self._parquet:
            artifacts["parquet"] = self._parquet
        root = root or default_root
        tag = ElementTag(
            root=root,
            level=0,
            omics=None,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
            ext=".tsv",
        ).merge(tag)

        super().__init__(
            key=key,
            run=actual_run,
            tag=tag,
            determinants=determinants,
            inputs=inputs,
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    @property
    def tsv(self) -> Path:
        if self._tsv is None:
            raise AttributeError(f"TableElement {self.name!r} has no TSV file.")
        return self._tsv

    @property
    def parquet(self) -> Path:
        if self._parquet is None:
            raise AttributeError(f"TableElement {self.name!r} has no Parquet file.")
        return self._parquet

    @property
    def index(self) -> pd.Index:
        """Return the DataFrame index (row labels) of this table."""
        return self.df.index

    @property
    def index_name(self) -> str:
        """Return the name of the DataFrame index column."""
        return self.df.index.name

    @property
    def file(self) -> Path:
        return self.parquet if self._parquet is not None else self.tsv

    # ------------------------------------------------------------------
    # column_roles helpers
    # ------------------------------------------------------------------

    def columns_for_role(self, role: str) -> list[str]:
        """Return all column names assigned to *role*."""
        return [col for col, r in self.column_roles.items() if r == role]

    def roles_present(self) -> set[str]:
        """Return the set of distinct roles registered for this table."""
        return set(self.column_roles.values())

    # ------------------------------------------------------------------
    # Views
    # ------------------------------------------------------------------

    def view(self, role: str) -> pd.DataFrame:
        """Return a :class:`~pandas.DataFrame` with only the columns that
        have the given *role*.

        The index of the full table is preserved.

        Parameters
        ----------
        role : str
            Semantic role to filter for, e.g. ``"vst"``, ``"annotation"``.

        Returns
        -------
        pd.DataFrame
            Subset of :attr:`df` with only the matching columns.

        Raises
        ------
        KeyError
            If no columns are registered for *role*.
        """
        cols = self.columns_for_role(role)
        if not cols:
            raise KeyError(
                f"No columns with role {role!r} in {self.name!r}.\n"
                f"Available roles: {sorted(self.roles_present())}"
            )
        return self.df[cols]

    def roles(
        self,
        *roles: str,
        extra_columns: Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return a DataFrame with columns matching any of the given *roles*.

        Parameters
        ----------
        *roles : str
            One or more role names to include.
        extra_columns : Iterable[str] | None
            Additional column names to include regardless of their role.

        Returns
        -------
        pd.DataFrame
            Combined subset; column order follows the original table.
        """
        wanted: set[str] = set()
        for role in roles:
            wanted.update(self.columns_for_role(role))
        if extra_columns:
            wanted.update(extra_columns)
        # preserve original column order
        ordered = [c for c in self.df.columns if c in wanted]
        return self.df[ordered]

    def propagate(
        self,
        roles: Iterable[str] | None = None,
    ) -> dict[str, str]:
        """Return a ``column_roles`` dict containing only the *annotation*
        (or custom) columns from this table — ready to be merged into a
        derived element's ``column_roles``.

        This lets transformations like VST or TMM carry annotation columns
        forward without having to enumerate them explicitly.

        Parameters
        ----------
        roles : Iterable[str] | None
            Roles to propagate.  Defaults to :attr:`ANNOTATION_ROLES`
            (``{"annotation", "metadata"}``).

        Returns
        -------
        dict[str, str]
            Mapping ``{column: role}`` for all matching columns.

        Examples
        --------
        >>> derived_roles = {"vst_A": "vst", "vst_B": "vst"}
        >>> derived_roles.update(raw_counts.propagate())
        >>> vst_element = TableElement(..., column_roles=derived_roles)
        """
        keep = frozenset(roles) if roles is not None else self.ANNOTATION_ROLES
        return {col: role for col, role in self.column_roles.items() if role in keep}

    # ------------------------------------------------------------------
    # Convenience: expression-matrix view
    # ------------------------------------------------------------------

    def matrix(
        self,
        kind: str,
    ) -> ndarray:
        """Return the expression matrix for a given normalisation *kind*.

        A thin convenience wrapper around :meth:`view`.

        Parameters
        ----------
        kind : str
            Normalisation kind / role, e.g. ``"vst"``, ``"rlog"``,
            ``"raw_counts"``, ``"cpm"``.

        Returns
        -------
        np.ndarray
        """
        return self.view(kind).to_numpy()

    # ------------------------------------------------------------------
    # Load from disk
    # ------------------------------------------------------------------

    @cached_property
    def df(self) -> DataFrame:
        """Load the result table (prefers Parquet, falls back to TSV)."""
        if self._parquet is not None and self._parquet.exists():
            df = pd.read_parquet(self._parquet)
        elif self._tsv is not None and self._tsv.exists():
            df = pd.read_csv(self._tsv, sep="\t")  # , index_col=self.index_column)
        else:
            raise FileNotFoundError(
                f"Neither Parquet nor TSV output found for {self.name!r}.\n"
                f"  parquet: {self._parquet}\n"
                f"  tsv:     {self._tsv}"
            )
        if self.index_column and self.index_column in df.columns:
            df["index"] = df[self.index_column].copy()
            df.set_index("index", inplace=True)
        return df


class AdataElement(Element):
    """An Element that wraps an :class:`anndata.AnnData` object persisted as
    an ``.h5ad`` file.

    Analogous to :class:`TableElement` for AnnData objects used in single-cell
    or bulk RNA-seq workflows.

    Parameters
    ----------
    key : str
        Unique cache key.
    run : Callable | Path | str | None
        Zero-argument callable that executes the computation and writes the
        ``.h5ad`` file.  If a *Path* or *str* pointing to an existing
        ``.h5ad`` file is passed, the element is created as an input-only
        node that validates the file exists.  Pass ``None`` together with an
        explicit *h5ad* path for the same effect.
    tag : PartialElementTag | ElementTag | None
        Tag used for naming / provenance.
    h5ad : Path | None
        Output ``.h5ad`` path.
    obs_roles : Mapping[str, str] | None
        Optional semantic role per ``obs`` column, e.g.::

            {
                "cell_type": "annotation",
                "sample":    "metadata",
                "UMAP_1":    "embedding",
                "UMAP_2":    "embedding",
            }

        Well-known roles: ``"annotation"``, ``"metadata"``, ``"embedding"``,
        ``"qc"``, ``"cluster"``.
    var_roles : Mapping[str, str] | None
        Optional semantic role per ``var`` column, e.g.::

            {"highly_variable": "feature_selection", "gene_name": "annotation"}
    determinants, inputs, pres, name
        Forwarded to :class:`Element`.

    Artefacts
    ---------
    ``"h5ad"`` → *h5ad*

    Examples
    --------
    >>> raw = AdataElement(
    ...     key="raw_counts",
    ...     run=compute_raw_counts,
    ...     tag=my_tag,
    ...     h5ad=Path("results/raw.h5ad"),
    ...     obs_roles={"cell_type": "annotation", "sample": "metadata"},
    ... )
    >>> raw.adata          # loads AnnData from disk
    >>> raw.obs_view("annotation")   # obs columns with role "annotation"
    >>> raw.obsm_key("embedding")    # first obsm key associated with role
    """

    ANNOTATION_ROLES: frozenset[str] = frozenset({"annotation", "metadata"})

    def __init__(
        self,
        key: str,
        run: Callable[[], Any] | Path | str | None = None,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        root: str | None = None,
        h5ad: Path | None = None,
        # obs_roles: Mapping[str, str] | None = None,
        # var_roles: Mapping[str, str] | None = None,
        determinants: tuple | None = None,
        inputs: tuple[Path, ...] | None = None,
        pres: tuple[Element, ...] | None = None,
        name: str | None = None,
    ) -> None:
        # --- resolve file-based construction ---
        if isinstance(run, (Path, str)):
            p = Path(run).absolute()
            if p.suffix not in (".h5ad",):
                raise ValueError(
                    f"Expected a .h5ad file, got suffix {p.suffix!r}. "
                    "Pass h5ad= explicitly."
                )
            h5ad = h5ad or p
            actual_run = paths_exists(p)
            actual_run.threads = 1
            actual_run.command = [f"check_exists({p})"]  # type: ignore[attr-defined]
            inputs = inputs or (p,)
            default_root = p.stem

        elif run is None:
            if h5ad is None:
                raise ValueError("Provide h5ad= when run=None.")
            actual_run = paths_exists(h5ad)
            actual_run.threads = 1
            actual_run.command = [f"check_exists({h5ad})"]  # type: ignore[attr-defined]
            inputs = inputs or (h5ad,)
            default_root = h5ad.stem

        else:
            if h5ad is None:
                raise ValueError("h5ad= must be provided when run is a callable.")
            actual_run = run
            default_root = root or (name or key)

        self._h5ad = Path(h5ad).absolute()
        # self.obs_roles: dict[str, str] = dict(obs_roles or {})
        # self.var_roles: dict[str, str] = dict(var_roles or {})

        artifacts: dict[str, Any] = {"h5ad": self._h5ad}

        root = root or default_root
        tag = ElementTag(
            root=root,
            level=0,
            omics=None,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
            ext=".h5ad",
        ).merge(tag)

        super().__init__(
            key=key,
            run=actual_run,
            tag=tag,
            determinants=determinants,
            inputs=inputs,
            artifacts=artifacts,
            pres=pres,
            name=name,
        )

    # ------------------------------------------------------------------
    # Core property
    # ------------------------------------------------------------------

    @property
    def h5ad(self) -> Path:
        return self._h5ad

    @cached_property
    def adata(self):
        """Load and return the :class:`anndata.AnnData` from disk."""
        try:
            import anndata as ad  # type: ignore[import]
        except ImportError as exc:
            raise ImportError(
                "anndata is required for AdataElement. "
                "Install it with: pip install anndata"
            ) from exc
        if not self._h5ad.exists():
            raise FileNotFoundError(
                f"h5ad file not found for {self.name!r}: {self._h5ad}"
            )
        return ad.read_h5ad(self._h5ad)

    def view(
        self,
        layer: str | None = None,
        obs_roles: dict[str, list[str]] = None,
        var_roles: dict[str, list[str]] = None,
    ) -> DataFrame:
        """Return a DataFrame of a given layer filtered by obs columns and/or var columns.

            Parameters
            ----------
            layer : str | None
                The layer to filter for, e.g. ``"raw"``, ``"cpm"``, ``"vst"``.
            obs_roles : dict[str, list[str]] | None
                Dictionary mapping obs roles to lists of column names.
            var_roles : dict[str, list[str]] | None
                Dictionary mapping var roles to lists of column names.

            Returns
            -------
            DataFrame
                A DataFrame containing the filtered data based on the specified layer and roles.

            Raises
            ------
        # --- Result ---
        if not diffs:
            result = True
            msg = f"[✓] No differences in sig_data for {key!r}"
        else:
            result = False
            msg = f"[✗] Sig_data differs for {key!r}:"
            for diff in diffs:
            KeyError
                If no ``obs`` or ``var`` columns are registered.
        """
        var_mask = Series(True, index=self.adata.var_names)
        if var_roles:
            for col, wanted in var_roles.items():
                var_mask &= self.adata.var[col].isin(wanted)
        obs_mask = Series(True, index=self.adata.obs_names)
        if obs_roles:
            for col, wanted in obs_roles.items():
                obs_mask &= self.adata.obs[col].isin(wanted)
        df_view = self.adata[layer][obs_mask, var_mask]
        return df_view

    # # ------------------------------------------------------------------
    # # obs helpers
    # # ------------------------------------------------------------------

    # def obs_columns_for_role(self, role: str) -> list[str]:
    #     """Return all ``obs`` column names assigned to *role*."""
    #     return [col for col, r in self.obs_roles.items() if r == role]

    # def obs_roles_present(self) -> set[str]:
    #     """Return the set of distinct roles registered for ``obs``."""
    #     return set(self.obs_roles.values())

    # def obs_view(self, role: str) -> DataFrame:
    #     """Return a DataFrame of ``obs`` columns that have the given *role*.

    #     Parameters
    #     ----------
    #     role : str
    #         Semantic role to filter for, e.g. ``"annotation"``, ``"qc"``.

    #     Returns
    #     -------
    #     DataFrame

    #     Raises
    #     ------
    #     KeyError
    #         If no ``obs`` columns are registered for *role*.
    #     """
    #     cols = self.obs_columns_for_role(role)
    #     if not cols:
    #         raise KeyError(
    #             f"No obs columns with role {role!r} in {self.name!r}.\n"
    #             f"Available roles: {sorted(self.obs_roles_present())}"
    #         )
    #     available = [c for c in cols if c in self.adata.obs.columns]
    #     return self.adata.obs[available]

    # # ------------------------------------------------------------------
    # # var helpers
    # # ------------------------------------------------------------------

    # def var_columns_for_role(self, role: str) -> list[str]:
    #     """Return all ``var`` column names assigned to *role*."""
    #     return [col for col, r in self.var_roles.items() if r == role]

    # def var_roles_present(self) -> set[str]:
    #     """Return the set of distinct roles registered for ``var``."""
    #     return set(self.var_roles.values())

    # def var_view(self, role: str) -> DataFrame:
    #     """Return a DataFrame of ``var`` columns that have the given *role*.

    #     Parameters
    #     ----------
    #     role : str
    #         Semantic role to filter for, e.g. ``"annotation"``, ``"feature_selection"``.

    #     Returns
    #     -------
    #     DataFrame

    #     Raises
    #     ------
    #     KeyError
    #         If no ``var`` columns are registered for *role*.
    #     """
    #     cols = self.var_columns_for_role(role)
    #     if not cols:
    #         raise KeyError(
    #             f"No var columns with role {role!r} in {self.name!r}.\n"
    #             f"Available roles: {sorted(self.var_roles_present())}"
    #         )
    #     available = [c for c in cols if c in self.adata.var.columns]
    #     return self.adata.var[available]

    # # ------------------------------------------------------------------
    # # obsm helpers
    # # ------------------------------------------------------------------

    # def obsm_keys_for_role(self, role: str) -> list[str]:
    #     """Return ``obsm`` keys whose name is registered under *role*.

    #     Convention: register obsm keys the same way as obs columns in
    #     ``obs_roles``, e.g. ``{"X_umap": "embedding", "X_pca": "embedding"}``.
    #     """
    #     return [k for k, r in self.obs_roles.items() if r == role and k in self.adata.obsm]

    # def matrix(self, obsm_key: str) -> ndarray:
    #     """Return a low-dimensional embedding / matrix stored in ``obsm``.

    #     Parameters
    #     ----------
    #     obsm_key : str
    #         Key in ``adata.obsm``, e.g. ``"X_umap"`` or ``"X_pca"``.

    #     Returns
    #     -------
    #     np.ndarray
    #     """
    #     return self.adata.obsm[obsm_key]

    # # ------------------------------------------------------------------
    # # Propagate annotation roles to derived elements
    # # ------------------------------------------------------------------

    # def propagate(
    #     self,
    #     obs_roles: Iterable[str] | None = None,
    #     var_roles: Iterable[str] | None = None,
    # ) -> tuple[dict[str, str], dict[str, str]]:
    #     """Return ``(obs_roles, var_roles)`` dicts containing only the
    #     *annotation* columns from this element — ready to be merged into a
    #     derived element's role mappings.

    #     Parameters
    #     ----------
    #     obs_roles : Iterable[str] | None
    #         ``obs`` roles to propagate.  Defaults to :attr:`ANNOTATION_ROLES`.
    #     var_roles : Iterable[str] | None
    #         ``var`` roles to propagate.  Defaults to :attr:`ANNOTATION_ROLES`.

    #     Returns
    #     -------
    #     tuple[dict[str, str], dict[str, str]]
    #         ``(obs_role_mapping, var_role_mapping)`` for annotation columns.

    #     Examples
    #     --------
    #     >>> obs_r, var_r = raw.propagate()
    #     >>> derived = AdataElement(..., obs_roles={**obs_r, "cluster": "cluster"})
    #     """
    #     keep_obs = frozenset(obs_roles) if obs_roles is not None else self.ANNOTATION_ROLES
    #     keep_var = frozenset(var_roles) if var_roles is not None else self.ANNOTATION_ROLES
    #     return (
    #         {col: role for col, role in self.obs_roles.items() if role in keep_obs},
    #         {col: role for col, role in self.var_roles.items() if role in keep_var},
    #     )

    # # ------------------------------------------------------------------
    # # Convenience: filter cells / genes by role
    # # ------------------------------------------------------------------

    # def subset_obs(self, role: str) -> Any:
    #     """Return an AnnData view filtered to cells where *role* columns are
    #     not null/NaN.

    #     Useful for e.g. restricting to annotated cells only.
    #     """
    #     cols = self.obs_columns_for_role(role)
    #     if not cols:
    #         raise KeyError(f"No obs columns with role {role!r}")
    #     mask = self.adata.obs[cols].notna().all(axis=1)
    #     return self.adata[mask]

    # def subset_var(self, role: str, value: Any = True) -> Any:
    #     """Return an AnnData view filtered to genes where a *role* column
    #     equals *value* (default ``True``).

    #     Useful for e.g. restricting to highly variable genes.
    #     """
    #     cols = self.var_columns_for_role(role)
    #     if not cols:
    #         raise KeyError(f"No var columns with role {role!r}")
    #     col = cols[0]
    #     mask = self.adata.var[col] == value
    #     return self.adata[:, mask]


def generate_element_key_name(
    tag: ElementTag,
    tool_name: str,
    tool_version: str | None = None,
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
    parts = [tag.default_name]
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


def explain_signature_diff(
    element,
    store: dict[str, Any],
) -> tuple[bool, str]:
    """
    Compares the sig_data of an Element with the stored entry.

    Returns True if identical, False if different.
    Provides a detailed explanation in case of differences.
    """
    key = element.key
    current = element.sig_data()

    if key not in store:
        raise ValueError(f"Key not found in store: {key!r}")

    stored = store[key]

    diffs: list[str] = []

    # --- key ---
    if current.get("key") != stored.get("key"):
        diffs.append(
            f"Keys differ:\n"
            f"    current: {current.get('key')!r}\n"
            f"    stored:  {stored.get('key')!r}"
        )

    # --- determinants ---
    cur_det = list(current.get("determinants", []))
    sto_det = list(stored.get("determinants", []))
    if cur_det != sto_det:
        diffs.append(
            f"Determinants differ:\n"
            f"    current: {cur_det}\n"
            f"    stored:  {sto_det}"
        )

    # --- inputs ---
    cur_inputs = current.get("inputs", [])
    sto_inputs = stored.get("inputs", [])
    if cur_inputs != sto_inputs:
        detail_lines = [
            f"Inputs differ: ({len(cur_inputs)} current vs {len(sto_inputs)} stored)"
        ]  # noqa: E501
        all_paths = {inp.get("path") for inp in cur_inputs + sto_inputs}
        for path in sorted(all_paths):
            cur_entry = next((i for i in cur_inputs if i.get("path") == path), None)
            sto_entry = next((i for i in sto_inputs if i.get("path") == path), None)
            if cur_entry != sto_entry:
                detail_lines.append(f"    path: {path!r}")
                detail_lines.append(f"      current: {cur_entry}")
                detail_lines.append(f"      stored:  {sto_entry}")
        diffs.append("\n".join(detail_lines))

    # --- artifacts ---
    cur_art = current.get("artifacts", {})
    sto_art = stored.get("artifacts", {})
    if cur_art != sto_art:
        all_artifact_keys = set(cur_art) | set(sto_art)
        detail_lines = ["Artifacts differ:"]
        for k in sorted(all_artifact_keys):
            if cur_art.get(k) != sto_art.get(k):
                detail_lines.append(f"    [{k!r}]")
                detail_lines.append(f"      current: {cur_art.get(k)!r}")
                detail_lines.append(f"      stored:  {sto_art.get(k)!r}")
        diffs.append("\n".join(detail_lines))

    # --- pre_sigs ---
    cur_pre = sorted(current.get("pre_sigs", []))
    sto_pre = sorted(stored.get("pre_sigs", []))
    if cur_pre != sto_pre:
        only_current = set(cur_pre) - set(sto_pre)
        only_stored = set(sto_pre) - set(cur_pre)
        detail_lines = ["Predecessor signatures differ:"]
        if only_current:
            detail_lines.append(f"    only in current: {sorted(only_current)}")
        if only_stored:
            detail_lines.append(f"    only in stored:  {sorted(only_stored)}")
        diffs.append("\n".join(detail_lines))

    # --- Result ---
    if not diffs:
        result = True
        msg = f"[✓] No differences in sig_data for {key!r}"
    else:
        result = False
        msg = f"[✗] Sig_data differs for {key!r}:"
        for diff in diffs:
            msg += f"\n{diff}"

    return result, msg
