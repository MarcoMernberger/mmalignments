import hashlib
import inspect
import json
from pathlib import Path
from typing import Any
from functools import wraps, lru_cache


def depends(*deps):
    def decorator(fn):
        fn.__dependencies__ = set(getattr(fn, "__dependencies__", set()))
        fn.__dependencies__.update(deps)
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


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def short_hash(signature: str, n: int = 8) -> str:
    return signature[:n]
