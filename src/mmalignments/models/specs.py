from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from enum import Enum
from subprocess import CompletedProcess
from typing import (
    Any,
    Callable,
)

from mmalignments.services.dependencies import (
    collect_code_dependency,
    function_hash,
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
    def signature(self) -> str:
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


class ValidationPolicy(Enum):
    CHECK = "check"  # default behaviour
    FORCE_RUN = "force_run"
    FORCE_SKIP = "force_skip"
