import hashlib
import json
from functools import cached_property
from pathlib import Path
from typing import Any, Callable, Protocol, runtime_checkable


@runtime_checkable
class Signifiable(Protocol):

    @cached_property
    def signature(self) -> str: ...


class DynamicValue(Signifiable):
    def __init__(
        self, resolver: Callable[[], Any], dependency: list[Any] | None = None
    ):
        self.resolver = resolver
        self.dependency = dependency or []

    def resolve(self) -> Any:
        return self.resolver()

    @cached_property
    def signature(self):
        return combined_signature(function_hash(self.resolver), *self.dependency)


def depends(*deps):
    def decorator(fn):
        existing = list(getattr(fn, "__dependencies__", []))

        for dep in deps:
            if dep not in existing:
                existing.append(dep)

        fn.__dependencies__ = existing
        return fn

    return decorator


def function_hash(fn) -> str:
    if fn is None:
        return "none"
    try:
        code = fn.__code__
        payload = (
            code.co_code,
            # code.co_consts,   # this invalidates?
            code.co_names,
        )
        return hashlib.sha256(repr(payload).encode()).hexdigest()
    except Exception:
        return hashlib.sha256(fn.__name__.encode()).hexdigest()


def collect_code_dependency(fn, visited=None):
    if visited is None:
        visited = set()
    if fn is None or fn in visited:
        return []
    visited.add(fn)
    deps = getattr(fn, "__dependencies__", None) or []
    result = [fn]
    for d in deps:
        result.extend(collect_code_dependency(d, visited))
    return result


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


def file_signature(p: Path, head_bytes: int = 65_536) -> str:
    return str(file_sig(p, head_bytes))


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def short_hash(signature: str, n: int = 8) -> str:
    return signature[:n]


def combined_signature(*args) -> str:
    sigs = []
    for arg in args:
        if isinstance(arg, Signifiable):
            sigs.append(arg.signature)  # all Signifiable
        elif isinstance(arg, Path) and arg.is_file():  # all paths
            sigs.append(file_signature(arg))
        elif callable(arg):  # Callable, runnable and functions
            sigs.append(function_hash(arg))
        elif isinstance(arg, (str, int, float, bool)):
            sigs.append(str(arg))
        elif isinstance(arg, DynamicValue):  # all dynamic values
            sigs.append(arg.signature)
        else:
            sigs.append(stable_hash(arg))
    sig = stable_hash(sigs)
    return sig


def _compute_signature(self) -> str:
    try:
        return stable_hash(self.signature_data)
    except Exception as e:
        raise RuntimeError(
            f"Failed to compute signature for {self.key!r}, signature_data: {self.signature_data!r}"
        ) from e


def try_cast(v: str) -> Any:
    for cast in (int, float):
        try:
            return cast(v)
        except ValueError:
            pass
    return v
