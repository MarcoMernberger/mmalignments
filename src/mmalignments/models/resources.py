from __future__ import annotations

import os
from contextvars import ContextVar
from dataclasses import dataclass, field


def _detect_cpus() -> int:
    for env_var in ("MMALIGNMENTS_THREADS", "SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"):
        val = os.environ.get(env_var)
        if val is not None:
            try:
                return max(1, int(val))
            except ValueError:
                pass
    cpus = os.cpu_count() or 1
    return max(1, cpus // 2)


# ContextVar: wird vom Executor gesetzt, von External gelesen
_current_resources: ContextVar[ResourceConfig | None] = ContextVar(
    "current_resources", default=None
)


def current_resources() -> ResourceConfig:
    """Return the ResourceConfig set by the active Executor build context.

    Falls back to auto-detection if called outside a build context.
    """
    r = _current_resources.get()
    return r if r is not None else ResourceConfig.detect()


@dataclass
class ResourceConfig:
    """Central resource configuration – single source of truth for the machine."""

    threads: int = field(default_factory=_detect_cpus)

    @classmethod
    def detect(cls) -> ResourceConfig:
        return cls(threads=_detect_cpus())

    @classmethod
    def from_env(cls, default_threads: int = 4) -> ResourceConfig:
        val = os.environ.get("MMALIGNMENTS_THREADS")
        if val is not None:
            try:
                return cls(threads=max(1, int(val)))
            except ValueError:
                pass
        return cls(threads=default_threads)

    def for_tool(self, fraction: float = 1.0) -> int:
        """Return thread count for a specific tool (e.g. 0.5 = half the CPUs)."""
        return max(1, int(self.threads * fraction))
