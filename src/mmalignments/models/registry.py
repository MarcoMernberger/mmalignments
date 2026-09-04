from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterator

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
        """Register *e* and return it if the key is unique.

        If an Element with the same key already exists, registration is only
        accepted when the exact same object is passed again. Registering a
        different Element object under an existing key raises immediately.

        Parameters
        ----------
        e : Element
            The Element to intern.

        Raises
        ------
        ValueError
            If a different Element object is registered with an existing key.
        """
        existing = self._by_key.get(e.key)
        if existing is None:
            self._by_key[e.key] = e
            return e

        if existing is e:
            return e

        raise ValueError(
            "Duplicate key collision in ElementRegistry.intern:\n"
            f"  key: {e.key!r}\n"
            f"  existing element: {existing!r}\n"
            f"  new element:      {e!r}\n"
            "Element keys must be unique per distinct Element object."
        )

    def keys(self) -> Iterator[str]:
        """Iterate over all registered Element keys."""
        return iter(self._by_key.keys())

    def get(self, key: str) -> Element | None:
        """Return the Element registered under *key*, or ``None``."""
        return self._by_key.get(key)

    def print(self) -> None:
        """Print all registered Elements to stdout."""
        print("ElementRegistry:")
        for key in sorted(self._by_key):
            print(f"  {key}")

    def write_registry(self, outfile: Path) -> None:
        with open(outfile, "w") as outp:
            for key, element in sorted(self._by_key.items()):
                outp.write(f"{key}:\n{element.describe()}\n")
