from __future__ import annotations

from collections.abc import Callable, Hashable
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from threading import RLock
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
    """Track Elements created inside an :func:`element_build_context`.

    Each element key belongs to one object. Shared resources also have a
    caller-supplied identity for looking up their Element before construction.
    """

    def __init__(self) -> None:
        """Initialize the registry for a pipeline."""
        self._by_key: Dict[str, Element] = {}
        self._shared: dict[Hashable, Element] = {}
        self._lock = RLock()

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
        with self._lock:
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

    def get_or_create(
        self, resource: Hashable, factory: Callable[[], Element]
    ) -> Element:
        """Get a shared prerequisite or build and register it on first use.

        Parameters
        ----------
        resource : Hashable
            Namespaced identity covering every setting that affects the
            resource. Requests with the same identity share a single Element.
        factory : Callable[[], Element]
            Builds the Element on a cache miss.

        Returns
        -------
        Element
            The registered shared prerequisite.

        Raises
        ------
        ValueError
            If the new Element's key belongs to a different Element.
        """
        with self._lock:
            existing = self._shared.get(resource)
            if existing is not None:
                return existing
            element = self.intern(factory())
            self._shared[resource] = element
            return element

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
