from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from typing import Dict, Iterator, TYPE_CHECKING

if TYPE_CHECKING:
    from .elements import Element


_current_registry: ContextVar["ElementRegistry | None"] = ContextVar(
    "current_element_registry", default=None
)


def current_element_registry() -> "ElementRegistry | None":
    return _current_registry.get()


@contextmanager
def element_build_context(registry: "ElementRegistry"):
    token = _current_registry.set(registry)
    try:
        yield
    finally:
        _current_registry.reset(token)


class ElementRegistry:
    """Intern registry: ensures one instance per key and tracks all Elements
    created inside an :func:`element_build_context`."""

    def __init__(self) -> None:
        self._by_key: Dict[str, Element] = {}

    def intern(self, e: Element) -> Element:
        """Register *e* and return the canonical instance for its key.

        If an Element with the same key already exists, the existing instance
        is returned (and a signature check is performed).  Otherwise *e* is
        stored and returned as-is.

        Parameters
        ----------
        e : Element
            The Element to intern.

        Raises
        ------
        ValueError
            If a key collision with a different signature is detected.
        """
        existing = self._by_key.get(e.key)
        if existing is None:
            self._by_key[e.key] = e
            return e

        ex_sig = getattr(existing, "signature", None)
        e_sig = getattr(e, "signature", None)
        if ex_sig is not None and e_sig is not None and ex_sig != e_sig:
            raise ValueError(f"Key collision with different signature: {e.key}")

        return existing

    def keys(self) -> Iterator[str]:
        """Iterate over all registered Element keys."""
        return iter(self._by_key.keys())

    def get(self, key: str) -> Element | None:
        """Return the Element registered under *key*, or ``None``."""
        return self._by_key.get(key)

    def print(self) -> None:
        """Print all registered Elements to stdout."""
        print("ElementRegistry:")
        for key, element in sorted(self._by_key.items()):
            sig = getattr(element, "signature", None)
            if sig:
                print(f"  {key}  [{sig}]")
            else:
                print(f"  {key}")
