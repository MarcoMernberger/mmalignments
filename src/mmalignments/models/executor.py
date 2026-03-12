from __future__ import annotations

import json
import logging
import os
import shlex
import signal
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, Future
from contextlib import contextmanager
from logging import Logger
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable
from subprocess import CalledProcessError
from mmalignments.models.elements import Element
from mmalignments.models.registry import ElementRegistry, element_build_context
from mmalignments.models.resources import ResourceConfig, _current_resources
from mmalignments.models.status import (  # type: ignore[import]
    NodeState,
    ProgressReporter,
)
from mmalignments.services.io import from_json, parents


class Executor:
    """
    Executor for Element DAGs with:
      - build() context manager that activates the ElementRegistry
      - record() snapshots the registry state so play() only runs new Elements
      - deterministic toposort
      - dry-run planning
      - per-node signature cache persistence (robust against later crashes)
      - optional continue_on_error
    """

    def __init__(
        self,
        logger: Logger | None = None,
        cache_path: Path | None = None,
        verbose_level: int = logging.INFO,
        registry: ElementRegistry | None = None,
        resources: ResourceConfig | None = None,
    ):
        self.resources = resources or ResourceConfig.detect()
        self.log_dir = cache_path or Path(".runlog")
        self.cache_dir = cache_path or Path("cache/run")
        self.signature_store_path = self.cache_dir / "signatures.json"
        self.logger = logger or logging.getLogger("pipeline.run")
        self.dot_path = self.cache_dir / "dag.dot"
        self.verbose_level = verbose_level
        self.registry = registry or ElementRegistry()
        self.stop_event = threading.Event()
        self.state = "init"
        signal.signal(signal.SIGTERM, self._handle_interrupt)
        # Keys present in the registry at the time record() was called.
        # play() will only run elements whose keys were added after record().
        self._baseline_keys: set[str] | None = None
        self.active_processes: set[subprocess.Popen] = set()

    def _handle_interrupt(self, signum, frame):
        self.logger.warning("Interrupt received — stopping pipeline...")
        self.stop_event.set()
        for proc in list(self.active_processes):  # terminate running subprocesses
            proc.terminate()

    def _keyboard_listener(self):
        for line in sys.stdin:
            if line.strip().lower() == "q":
                self.logger.warning("User requested stop.")
                self.stop_event.set()
                for proc in list(self.active_processes):
                    proc.terminate()
                break

    ###########################################################################
    # Build context  ← injects BOTH registry and resources into contextvars
    ###########################################################################

    @contextmanager
    def build(self):
        """Activate the ElementRegistry and ResourceConfig for this pipeline.

        Every ``Element`` created inside this block is automatically interned.
        Every ``External.apply_threads()`` call automatically uses
        ``self.resources`` without any explicit passing.

        Usage
        -----
        >>> executor = Executor(resources=ResourceConfig(threads=16))
        >>> with executor.build():
        ...     mapped = bwa.alignsort(sample, genome)
        ...     marked = gatk.mark(mapped)
        """
        resources_token = _current_resources.set(self.resources)
        try:
            with element_build_context(self.registry):
                yield
        finally:
            _current_resources.reset(resources_token)

    ###########################################################################
    # record / play
    ###########################################################################

    def record(self) -> None:
        """Snapshot the set of currently registered Element keys.

        After calling this, :meth:`play` will only execute Elements that were
        registered *after* this snapshot was taken.  Call it right after
        ``executor.build()`` has been entered and before defining any pipeline
        Elements.
        """
        self._baseline_keys = set(self.registry.keys())

    def play(
        self,
        *,
        verbose: bool = False,
        log_run_only: bool = False,
        progress: bool = True,
        parallel: bool = False,
        continue_on_error: bool = False,
        final_raise_on_error: bool = True,
        dry_run: bool = False,
        max_workers: int | None = None,
        dot_path: str | Path | None = None,
    ) -> None:
        """Run all Elements registered since :meth:`record` was called."""
        dot_path = dot_path or self.dot_path
        all_keys = set(self.registry.keys())
        new_keys = (
            all_keys
            if self._baseline_keys is None
            else all_keys - self._baseline_keys  # noqa: E501
        )
        targets = [e for k in new_keys if (e := self.registry.get(k)) is not None]
        if not targets:
            self.log(None, "No Elements found to run.", level="INFO")
            return

        self.run_targets(
            targets,
            dry_run=dry_run,
            verbose=verbose,
            dot_path=dot_path,
            continue_on_error=continue_on_error,
            final_raise_on_error=final_raise_on_error,
            log_run_only=log_run_only,
            progress=progress,
            parallel=parallel,
            max_workers=max_workers or self.resources.threads,
        )

    ###########################################################################
    # Cache
    ###########################################################################

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

    ###########################################################################
    # Logging helper
    ###########################################################################

    def msg(self, message: str, level: str = "INFO") -> None:
        """Log messages to the logger"""
        if level == "INFO":
            self.logger.info(message)
        elif level == "WARNING":
            self.logger.warning(message)
        elif level == "ERROR":
            self.logger.error(message)
        else:
            self.logger.debug(message)
        print(message)

    def log(
        self,
        reporter: ProgressReporter | None,
        message: str,
        level: str = "INFO",
    ) -> None:
        """If we have a reporter for live progress, we cannot log directly"""
        if reporter:
            # Route through the live display only — suppress StreamHandler out
            reporter.push_log(f"[{level}] {message}")
        else:
            self.msg(message, level=level)

    def compose_element_message(self, will_run, node, reason, verbose) -> str:
        """for dry runs or runs without ProgressReporter"""
        if will_run:
            status = "RUN "
            t = f"  [{node.threads}t]" if hasattr(node, "threads") else ""
        else:
            status = "SKIP"
            t = ""
        message = f"{status} {node.name}  ({reason}){t}"
        if verbose:
            message += f"\n  key: {node.key}"
            message += f"\n  actual output: {node.output_files}"
            message += f"\n  artifacts: {node.artifacts}"
            message += f"\n  tag: {node.tag}"
            message += f"\n  pre: {node.pres}"
            cmd = getattr(node.run, "command", None)
            display = getattr(node.run, "command_display", None)
            if will_run:
                command = shlex.join(cmd) if cmd else "<no command available>"
                message += f"\n  cmd:  {command}"
                if display and display != cmd:
                    message += f"\n  display:  {display}"
                    message += f"\n  prov: {node.provenance}"
                message += f"\n  default output: {node.default_output_file}"
        return message

    ###########################################################################
    # Core execution
    ###########################################################################

    def run_targets(
        self,
        targets: list[Element],
        *,
        dry_run: bool = False,
        verbose: bool = False,
        dot_path: str | Path | None = None,
        continue_on_error: bool = True,
        final_raise_on_error: bool = True,
        log_run_only: bool = False,
        progress: bool = True,
        parallel: bool = False,  # ← neu
        max_workers: int | None = None,
    ) -> None:
        self.state = "running"
        if self.stop_event.is_set():
            self.log(None, "Pipeline already aborted before start.", level="WARNING")
            return
        cache = self.load_cache(self.signature_store_path)
        nodes = self.collect(targets)
        order = self.toposort(nodes)
        self.check_duplicate_outputs(nodes)
        self.write_dot(nodes, dot_path or self.dot_path, verbose)

        reporter: ProgressReporter | None = (
            ProgressReporter() if (progress and not dry_run) else None
        )
        if reporter:
            reporter.register(order)
            reporter.start_live()

        failures: list[tuple[str, str, Exception]] = []

        try:
            if sys.stdin.isatty():
                threading.Thread(
                    target=self._keyboard_listener,
                    daemon=True,
                ).start()
            if parallel and not dry_run:
                self.log(None, "Starting parallel execution...", level="INFO")
                self._run_parallel(
                    order,
                    cache=cache,
                    reporter=reporter,
                    continue_on_error=continue_on_error,
                    failures=failures,
                    max_workers=max_workers,
                    log_run_only=log_run_only,
                    verbose=verbose,
                )
            else:
                self.log(None, "Starting sequential execution...", level="INFO")
                self._run_sequential(
                    order,
                    cache=cache,
                    reporter=reporter,
                    continue_on_error=continue_on_error,
                    failures=failures,
                    dry_run=dry_run,
                    log_run_only=log_run_only,
                    verbose=verbose,
                )
        finally:
            if reporter:
                # Mark any nodes that were never reached (e.g. after error)
                for node in order:
                    np = reporter._nodes.get(node.key)
                    if np is not None and np.state == NodeState.PENDING:
                        reporter.mark_skip(node.key, "not reached")
                reporter.stop_live()

        if not dry_run:
            self.save_cache(self.signature_store_path, cache)

        if failures and final_raise_on_error:
            lines = [
                f"{len(failures)} element(s) failed:",
                *[
                    f"- {name}: {key} ({type(err).__name__}: {err})"
                    for name, key, err in failures
                ],
            ]
            raise RuntimeError("\n".join(lines))

    ###########################################################################
    # Sequential execution
    ###########################################################################

    def _run_sequential(
        self,
        order: list[Element],
        *,
        cache: dict[str, str],
        reporter: ProgressReporter | None,
        continue_on_error: bool,
        failures: list[tuple[str, str, Exception]],
        dry_run: bool,
        log_run_only: bool,
        verbose: bool,
    ) -> None:
        for node in order:
            if self.stop_event.is_set():
                self.log(reporter, "Pipeline aborted by user.", level="WARNING")
                if reporter:
                    reporter.stop_live()
                break
            cached_sig = cache.get(node.key)
            skip, reason = node.skip(cached_signature=cached_sig)
            will_run = not skip
            message = self.compose_element_message(will_run, node, reason, verbose)
            if will_run or (not log_run_only):
                self.log(reporter, message, level="INFO")

            if dry_run:
                continue

            if skip:
                if reporter:
                    reporter.mark_skip(node.key, reason)
                continue

            if reporter:
                reporter.mark_start(node.key)
            try:
                node()
                # update cache immediately after each node finishes, so we don't lose progress if the pipeline crashes later  # noqa: E501
                cache[node.key] = node.signature
                self.save_cache(self.signature_store_path, cache)
                if reporter:
                    reporter.mark_done(node.key)
            except CalledProcessError as e:
                failures.append((node.name, node.key, e))
                if reporter:
                    reporter.mark_failed(node.key)
                self.log(
                    reporter,
                    f"Error: {node.name}: failed{e}",
                    level="ERROR",
                )
                if not continue_on_error:
                    self.save_cache(self.signature_store_path, cache)
                    raise

    ###########################################################################
    # Parallel execution  (DAG-aware: runs independent nodes concurrently)
    ###########################################################################

    def _run_parallel(
        self,
        order: list[Element],
        *,
        cache: dict[str, str],
        reporter: ProgressReporter | None,
        continue_on_error: bool,
        failures: list[tuple[str, str, Exception]],
        max_workers: int | None,
        log_run_only: bool,
        verbose: bool,
    ) -> None:

        # Build predecessor-count map from toposorted order
        by_key: dict[str, Element] = {n.key: n for n in order}
        # How many prerequisites are not yet done
        pending_deps: dict[str, int] = {n.key: 0 for n in order}
        successors: dict[str, list[str]] = {n.key: [] for n in order}

        for node in order:
            for pre in node.pres:
                if pre.key in by_key:
                    pending_deps[node.key] += 1
                    successors[pre.key].append(node.key)

        cache_lock = threading.Lock()
        done_event = threading.Event()
        active: set[str] = set()
        completed: set[str] = set()
        lock = threading.Lock()
        abort = self.stop_event

        # nodes whose deps are all done
        ready: list[str] = [k for k, d in pending_deps.items() if d == 0]

        # get available workers
        workers = max_workers or max(1, (self.resources.threads // 4))

        # let the ThreadPoolExecutionor handle the threads
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures: dict[Future, str] = {}

            def submit_ready() -> None:
                if abort.is_set():
                    return
                for key in list(ready):
                    if abort.is_set():
                        return
                    if key not in active and key not in completed:
                        ready.remove(key)
                        active.add(key)
                        f = pool.submit(_run_one, key)
                        futures[f] = key

            def _run_one(key: str) -> None:
                if abort.is_set():
                    return
                node = by_key[key]
                cached_sig = cache.get(node.key)
                skip, reason = node.skip(cached_signature=cached_sig)
                message = self.compose_element_message(not skip, node, reason, verbose)
                if not skip or (not log_run_only):
                    self.log(reporter, message)
                if skip:
                    if reporter:
                        reporter.mark_skip(node.key, reason)
                else:
                    if reporter:
                        reporter.mark_start(node.key)
                    try:
                        node()
                        with cache_lock:
                            cache[node.key] = node.signature
                            self.save_cache(self.signature_store_path, cache)
                        if reporter:
                            reporter.mark_done(node.key)
                    except Exception as e:
                        failures.append((node.name, node.key, e))
                        if reporter:
                            reporter.mark_failed(node.key)
                        self.log(reporter, f"Error: {node.name}: {e}", level="ERROR")
                        if not continue_on_error:
                            abort.set()
                        return

                # unblock successors
                with lock:
                    if abort.is_set():
                        done_event.set()
                        return
                    completed.add(key)
                    active.discard(key)
                    for s in successors[key]:
                        pending_deps[s] -= 1
                        if pending_deps[s] == 0:
                            ready.append(s)
                    submit_ready()
                    if len(completed) == len(order):
                        done_event.set()

            with lock:
                submit_ready()

            while not done_event.is_set():
                if abort.is_set():
                    break
                done_event.wait(0.5)

            if abort.is_set():
                pool.shutdown(wait=False, cancel_futures=True)

    ###########################################################################
    # DAG helpers
    ###########################################################################

    def check_duplicate_outputs(self, nodes: Iterable[Element]) -> None:
        """Raise an error if two different Elements share an output file path.

        Parameters
        ----------
        nodes : Iterable[Element]
            All Elements in the execution graph (as returned by
            :meth:`collect`).

        Raises
        ------
        ValueError
            If two Elements with different keys share at least one output path.
        """
        seen: dict[Path, str] = {}
        conflicts: list[str] = []
        for node in nodes:
            for path in node.output_files or ():
                owner = seen.get(path)
                if owner is None:
                    seen[path] = node.key
                elif owner != node.key:
                    conflicts.append(
                        f"  {path}\n"
                        f"    claimed by: {owner}\n"
                        f"    and also by: {node.key}"
                    )
        if conflicts:
            raise ValueError(
                f"{len(conflicts)} output path conflict(s) detected:\n"
                + "\n".join(conflicts)
            )

    def collect(self, targets: Iterable[Element]) -> list[Element]:
        """Collect the full transitive closure of prerequisites.

        Parameters
        ----------
        targets : Iterable[Element]
            Root Elements to start collection from.

        Raises
        ------
        ValueError
            If two different Element objects share the same key (indicates
            Elements were built outside ``executor.build()``).
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

            if existing is e:
                continue

            raise ValueError(
                f"Duplicate key with different Element objects: {e.key}\n"
                "This indicates missing Element interning. Build Elements "
                "inside with executor.build() and make sure _intern() is called"
            )

        return list(by_key.values())

    def toposort(self, nodes: Iterable[Element]) -> list[Element]:
        """
        Return Elements in topological order (prerequisites before
        dependents).

        Parameters
        ----------
        nodes : Iterable[Element]
            All collected Elements.

        Raises
        ------
        RuntimeError
            If a dependency cycle is detected.
        """
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

    ###########################################################################
    # DOT graph generation
    ###########################################################################

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

    def write_dot(self, targets: Iterable[Element], dot_path: str | Path) -> Path:
        dot_path = Path(dot_path)
        dot_path.parent.mkdir(parents=True, exist_ok=True)
        dot_path.write_text(self.to_dot(targets), encoding="utf-8")
        self.log(None, f"Wrote DAG dot file to {dot_path}")
        return dot_path

    ###########################################################################
    # Capture
    ###########################################################################

    def capture_store_state(self, outstore: Path) -> None:
        """Save the current state of the workspace to a file."""
        all_keys = set(self.registry.keys())
        new_keys = (
            all_keys
            if self._baseline_keys is None
            else all_keys - self._baseline_keys  # noqa: E501
        )
        targets = [e for k in new_keys if (e := self.registry.get(k)) is not None]
        cache = self.load_cache(self.signature_store_path)
        nodes = self.collect(targets)
        self.check_duplicate_outputs(nodes)
        order = self.toposort(nodes)
        for node in order:
            cached_sig = cache.get(node.key)
            skip, reason = node.skip(cached_signature=cached_sig)
            if not skip and reason == "first run":
                if node.outputs_ok():
                    print(node)
                    raise ValueError(
                        f"Node {node.key} appears to be a first run but its outputs already exist. This may indicate a problem with the node's skip logic or signature computation."
                    )
            if skip:
                cache[node.key] = node.signature
                self.save_cache(outstore, cache)
