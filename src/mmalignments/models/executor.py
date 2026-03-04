from __future__ import annotations

import inspect
import json
import logging
import os
import shlex
from contextlib import contextmanager
from logging import Logger
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterable, Mapping

from mmalignments.models.tasks import Element
from mmalignments.models.registry import ElementRegistry, element_build_context
from mmalignments.services.io import from_json, parents


class Executor:
    """
    Executor for Element DAGs with:
      - record()/play() auto-collection from caller scope
      - deterministic toposort
      - dry-run planning
      - per-node signature cache persistence (robust against later crashes)
      - optional continue_on_error
      - Variant A: Element interning via ElementRegistry + contextvars
    """

    def __init__(
        self,
        logger: Logger | None = None,
        cache_path: Path | None = None,
        verbose_level: int = logging.INFO,
        registry: ElementRegistry | None = None,
    ):
        self.log_dir = cache_path or Path(".runlog")
        self.cache_dir = cache_path or Path("cache/run")
        self.signature_store_path = self.cache_dir / "signatures.json"
        self.logger = logger or logging.getLogger("pipeline.run")
        self.dot_path = self.cache_dir / "dag.dot"
        self.verbose_level = verbose_level

        # NEW: registry for Variant A (interning)
        self.registry = registry or ElementRegistry()

        self._baseline_ids: set[int] | None = None

    # -------------------------------------------------------------------------
    # Build context (Variant A)
    # -------------------------------------------------------------------------

    @contextmanager
    def build(self):
        """
        Activate ElementRegistry for interning while building Elements.

        Usage:
            with executor.build():
                a = FileElement(...)
                b = tool.step(a, ...)
        """
        with element_build_context(self.registry):
            yield

    # -------------------------------------------------------------------------
    # Cache
    # -------------------------------------------------------------------------

    def load_cache(self, path: Path) -> dict[str, str]:
        if not path.exists():
            return {}
        data = from_json(path)
        return {str(k): str(v) for k, v in (data or {}).items()}

    def save_cache(self, path: Path, cache: dict[str, str]) -> None:
        parents(path)
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            prefix=path.name + ".",
            dir=str(path.parent),
            delete=False,
        ) as tmp:
            json.dump(cache, tmp, indent=2, sort_keys=True)
            tmp.write("\n")
            tmp.flush()
            os.fsync(tmp.fileno())
            os.replace(tmp.name, path)

    # -------------------------------------------------------------------------
    # record/play auto-collection from caller scope
    # -------------------------------------------------------------------------

    def record(self) -> None:
        """
        Snapshot current Elements in caller scope,
        so play() can collect only new ones.
        """
        frame = inspect.currentframe()
        if frame is None or frame.f_back is None:
            self._baseline_ids = set()
            return

        caller = frame.f_back
        baseline = self._collect_elements_from_mapping(
            {**caller.f_globals, **caller.f_locals}
        )
        self._baseline_ids = {id(e) for e in baseline}

    def play(
        self,
        *,
        dry_run: bool = False,
        verbose: bool = False,
        dot_path: str | Path | None = None,
        print_plan: bool = True,
        continue_on_error: bool = False,
        raise_on_error: bool = True,
    ) -> None:
        """
        Collect new Elements from caller scope (relative to record()) and run.
        """
        dot_path = dot_path or self.dot_path

        frame = inspect.currentframe()
        if frame is None or frame.f_back is None:
            raise RuntimeError("Cannot inspect caller frame for auto-collection")
        caller = frame.f_back

        all_found = self._collect_elements_from_mapping(
            {**caller.f_globals, **caller.f_locals}
        )

        if self._baseline_ids is None:
            targets = all_found
        else:
            targets = [e for e in all_found if id(e) not in self._baseline_ids]

        if not targets:
            self.msg("No Elements found to run.", verbose=verbose)
            return

        self.run_targets(
            targets,
            dry_run=dry_run,
            verbose=verbose,
            dot_path=dot_path,
            print_plan=print_plan,
            continue_on_error=continue_on_error,
            raise_on_error=raise_on_error,
        )

    def _collect_elements_from_mapping(self, ns: Mapping[str, Any]) -> list[Element]:
        """Recursively find Elements in a namespace mapping (locals/globals)."""
        out: dict[str, Element] = {}

        def visit(x: Any) -> None:
            if x is None:
                return
            if isinstance(x, Element):
                out[x.key] = x
                return
            if isinstance(x, dict):
                for v in x.values():
                    visit(v)
                return
            if isinstance(x, (list, tuple, set)):
                for v in x:
                    visit(v)
                return

        for v in ns.values():
            visit(v)
        return list(out.values())

    # -------------------------------------------------------------------------
    # Logging helper
    # -------------------------------------------------------------------------

    def msg(self, message: str, verbose: bool = False, level: str = "INFO") -> None:
        if level == "INFO":
            self.logger.info(message)
        elif level == "WARNING":
            self.logger.warning(message)
        elif level == "ERROR":
            self.logger.error(message)
        else:
            self.logger.debug(message)
        if verbose:
            print(message)

    # -------------------------------------------------------------------------
    # Core execution
    # -------------------------------------------------------------------------

    def run_targets(
        self,
        targets: Iterable[Element],
        *,
        dry_run: bool = False,
        verbose: bool = False,
        dot_path: str | Path | None = None,
        print_plan: bool = True,
        continue_on_error: bool = False,
        raise_on_error: bool = True,
    ) -> None:
        cache = self.load_cache(self.signature_store_path)

        nodes = self.collect(targets)
        order = self.toposort(nodes)

        self.write_dot(nodes, dot_path, verbose)

        failures: list[tuple[str, str, Exception]] = []

        for node in order:
            cached_sig: str | None = cache.get(node.key)
            skip, reason = node.skip(cached_signature=cached_sig)
            will_run = not skip

            if print_plan:
                self.print_plan(will_run, node, reason, verbose)

            if dry_run or not will_run:
                continue

            try:
                node()
                cache[node.key] = node.signature
                self.save_cache(self.signature_store_path, cache)
            except Exception as e:
                failures.append((node.name, node.key, e))
                self.msg(
                    f"Error executing {node.name} ({node.key}): {e}",
                    verbose=verbose,
                    level="ERROR",
                )
                self.msg(
                    f"  provenance: {node.provenance}", verbose=True, level="ERROR"
                )
                if not continue_on_error:
                    self.save_cache(self.signature_store_path, cache)
                    raise

        if not dry_run:
            self.save_cache(self.signature_store_path, cache)

        if failures and raise_on_error:
            lines = [
                f"{len(failures)} element(s) failed:",
                *[
                    f"- {name}: {key} ({type(err).__name__}: {err})"
                    for name, key, err in failures
                ],
            ]
            raise RuntimeError("\n".join(lines))

    def print_plan(self, will_run, node, reason, verbose) -> None:
        status = "RUN " if will_run else "SKIP"
        print(node)
        self.msg(f"[{status}] {node.name}", verbose=verbose)
        self.msg(f"  key: {node.key}", verbose=verbose)
        self.msg(f"  reason: {reason}", verbose=verbose)
        self.msg(f"  prov: {node.provenance}", verbose=verbose)  #   remove
        self.msg(f"  output: {node.default_output_file}", verbose=verbose)  # remove
        self.msg(f"  artifacts: {node.artifacts}", verbose=verbose)  # remove

        cmd = getattr(node.run, "command", None)
        display = getattr(node.run, "command_display", None)
        if will_run:
            if cmd:
                self.msg(f"  cmd:  {shlex.join(cmd)}", verbose=verbose)
            else:
                self.msg("  cmd:  <no command available>", verbose=verbose)
            if display and display != cmd:
                self.msg(f"  display:  {display}", verbose=verbose)
            if verbose:
                self.msg(f"  prov: {node.provenance}", verbose=verbose)
        self.msg("", verbose=verbose)

    # -------------------------------------------------------------------------
    # DAG helpers
    # -------------------------------------------------------------------------

    def collect(self, targets: Iterable[Element]) -> list[Element]:
        """
        Collect closure of prerequisites. With Variant A, duplicates by key
        should never appear as different objects (interning guarantees that).
        """
        by_key: dict[str, Element] = {}
        stack = list(targets)

        while stack:
            e = stack.pop()
            existing = by_key.get(e.key)

            if existing is None:
                by_key[e.key] = e
                stack.extend(e.pres)
                continue

            # With interning, this should be the same object:
            if existing is e:
                continue

            # If this happens, Elements were built outside executor.build() or without _intern()
            raise ValueError(
                f"Duplicate key with different Element objects: {e.key}\n"
                "This indicates missing Element interning. Build Elements inside "
                "`with executor.build():` and ensure factories call _intern()."
            )

        return list(by_key.values())

    def toposort(self, nodes: Iterable[Element]) -> list[Element]:
        nodes = list(nodes)

        indeg: dict[str, int] = {n.key: 0 for n in nodes}
        by_key: dict[str, Element] = {n.key: n for n in nodes}
        succ: dict[str, list[str]] = {n.key: [] for n in nodes}

        for n in nodes:
            for pre in n.pres:
                if pre.key not in by_key:
                    continue
                succ[pre.key].append(n.key)
                indeg[n.key] += 1

        q: list[str] = [k for k, d in indeg.items() if d == 0]
        q.sort()

        out: list[Element] = []
        while q:
            k = q.pop(0)
            out.append(by_key[k])
            for s in succ[k]:
                indeg[s] -= 1
                if indeg[s] == 0:
                    q.append(s)
                    q.sort()

        if len(out) != len(nodes):
            cyc = [k for k, d in indeg.items() if d > 0]
            raise RuntimeError(f"Dependency cycle detected among: {cyc}")

        return out

    # -------------------------------------------------------------------------
    # DOT graph generation
    # -------------------------------------------------------------------------

    def to_dot(self, targets: Iterable[Element]) -> str:
        nodes = self.collect(targets)
        by_key: dict[str, Element] = {n.key: n for n in nodes}

        edges: set[tuple[str, str]] = set()
        for n in nodes:
            for pre in n.pres:
                if pre.key in by_key:
                    edges.add((pre.key, n.key))

        def label(e: Element) -> str:
            return getattr(e, "name", e.key)

        lines = [
            "digraph DAG {",
            "  rankdir=LR;",
            "  node [shape=box];",
        ]

        for k in sorted(by_key):
            ll = label(by_key[k]).replace('"', '\\"')
            lines.append(f'  "{k}" [label="{ll}"];')

        for a, b in sorted(edges):
            lines.append(f'  "{a}" -> "{b}";')

        lines.append("}")
        return "\n".join(lines)

    def write_dot(
        self, targets: Iterable[Element], dot_path: str | Path, verbose: bool = False
    ) -> Path:
        dot_path = Path(dot_path)
        dot_path.parent.mkdir(parents=True, exist_ok=True)
        dot_path.write_text(self.to_dot(targets), encoding="utf-8")
        self.msg(f"Wrote DAG dot file to {dot_path}", verbose=verbose)
        return dot_path
