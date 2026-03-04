# mmalignments/models/registry.py
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from typing import Dict, Protocol, TypeVar


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
