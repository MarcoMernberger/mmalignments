from __future__ import annotations

import json
import logging
import os
import shlex
import signal
import sys
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from logging import Logger
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterable, Literal

from rich.console import Console  # type: ignore[import]

from mmalignments.models.elements import Element, explain_signature_diff
from mmalignments.models.externals import External
from mmalignments.models.registry import ElementRegistry, element_build_context
from mmalignments.models.resources import ResourceConfig, _current_resources
from mmalignments.models.status import (
    NodeState,  # type: ignore[import]
    ProgressReporter,
)
from mmalignments.services.errors import PipelineError
from mmalignments.services.io import ensure, from_json, parents
from mmalignments.services.logging import initlog
from mmalignments.services.time import now_as_str


@dataclass
class NodeState:
    key: str
    state: Literal["DONE", "SKIPPED", "FAILED"] = "SKIPPED"
    reason: str = ""
    start: datetime | None = None
    stop: datetime | None = None
    elapsed: float | None = None

    def mark_start(self):
        self.start = datetime.now()

    def mark_done(self):
        self.stop = datetime.now()
        self.elapsed = (self.stop - self.start).total_seconds() if self.start else None
        self.state = "DONE"

    def mark_failed(self, reason: str):
        self.stop = datetime.now()
        self.elapsed = (self.stop - self.start).total_seconds() if self.start else None
        self.state = "FAILED"
        self.reason = reason

    def mark_skipped(self, reason: str):
        self.state = "SKIPPED"
        self.reason = reason
        self.start = self.stop = datetime.now()
        self.elapsed = 0


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
        self.main_dir = cache_path or Path("cache/.run")
        self.cache_dir = self.main_dir / "store"
        self.resources = resources or ResourceConfig.detect()
        self.log_dir = self.main_dir / "logs"
        self.elements_file = self.log_dir / f"order_{now_as_str()}.txt"
        ensure(self.log_dir, self.cache_dir)
        self.logger = logger or initlog(console=True, log_dir=self.log_dir)
        self.signature_store_path = self.cache_dir / "signatures.json"
        self.signature_data_path = self.cache_dir / "signature_data.json"
        # self.logger = logger or logging.getLogger("pipeline.run")
        self.dot_path = self.cache_dir / "dag.dot"
        self.verbose_level = verbose_level
        self.registry = registry or ElementRegistry()
        self.stop_event = threading.Event()
        self.state = "init"
        if threading.current_thread() is threading.main_thread():
            signal.signal(signal.SIGTERM, self._handle_interrupt)
            # signal.signal(signal.SIGINT, self._handle_interrupt)
        # Keys present in the registry at the time record() was called.
        # play() will only run elements whose keys were added after record().
        self._baseline_keys: set[str] | None = None

    def _handle_interrupt(self, signum, frame):
        self.logger.warning("Interrupt received — stopping pipeline...")
        self.stop_event.set()

    def _keyboard_listener(self):
        for line in sys.stdin:
            if line.strip().lower() == "q":
                self.logger.warning("User requested stop.")
                self.stop_event.set()
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
            External.add_main_logger(self.logger)
            with element_build_context(self.registry):
                yield
        finally:
            External.remove_main_logger()
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
        log_signatures: bool = False,
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
        if log_run_only:
            self.log(
                None,
                "Logging only Elements that need to run",
                level="INFO",
            )
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
            log_signatures=log_signatures,
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

    def load_sigstore(self, path: Path) -> dict[str, Any]:
        if not path.exists():
            return {}
        data = from_json(path)
        return data

    def save_sigstore(self, path: Path, cache: dict[str, Any]) -> None:
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
        # print(message)

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

    def compose_element_message(self, skip, node, reason, verbose) -> str:
        """for dry runs or runs without ProgressReporter"""
        if not skip:
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
            if not skip:
                command = shlex.join(cmd) if cmd else "<no command available>"
                message += f"\n  cmd:  {command}"
                if display and display != command:
                    message += f"\n  display:  {display}"
                    message += f"\n  prov: {node.provenance}"
                message += f"\n  default output: {node.default_output_file}"
        return message

    # def initlog(self, console: bool = False) -> Logger:
    #     timestamp = self.get_timestamp()
    #     return self.setup_run_logger(timestamp, console=console)

    # def setup_run_logger(self, timestamp: str, console: bool = False) -> Logger:
    #     run_log = self.log_dir / f"run_{timestamp}.log"
    #     logger = self.logger  # logging.getLogger("pipeline.run")
    #     logger.setLevel(logging.INFO)
    #     logger.propagate = False

    #     formatter = logging.Formatter(
    #         "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    #     )

    #     if not any(
    #         isinstance(h, logging.FileHandler) and h.baseFilename == str(run_log)
    #         for h in logger.handlers
    #     ):
    #         fh = logging.FileHandler(run_log)
    #         fh.setLevel(logging.INFO)
    #         fh.setFormatter(formatter)
    #         logger.addHandler(fh)

    #     if console and not any(
    #         isinstance(h, logging.StreamHandler)
    #         and not isinstance(h, logging.FileHandler)
    #         for h in logger.handlers
    #     ):
    #         ch = logging.StreamHandler()
    #         ch.setLevel(logging.INFO)
    #         ch.setFormatter(formatter)
    #         logger.addHandler(ch)

    #     return logger

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
        log_signatures: bool = False,
    ) -> None:
        self.state = "running"
        if self.stop_event.is_set():
            self.log(None, "Pipeline already aborted before start.", level="WARNING")
            return
        cache = self.load_cache(self.signature_store_path)
        nodes = self.collect(targets)
        order = self.toposort(nodes)
        self.registry.write_registry(self.elements_file)
        self.check_duplicate_outputs(nodes)
        self.write_dot(nodes, dot_path or self.dot_path)
        reporter: ProgressReporter | None = (
            ProgressReporter(logger=self.logger) if (progress and not dry_run) else None
        )
        if reporter:
            reporter.register(order)
            reporter.start_live()

        failures: list[tuple[str, str, PipelineError]] = []

        try:
            if sys.stdin.isatty():
                threading.Thread(
                    target=self._keyboard_listener,
                    daemon=True,
                ).start()
            if parallel and not dry_run:
                self.log(None, "Starting parallel execution...", level="INFO")
                failures, reporter = self._run_parallel(
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
                failures, reporter = self._run_sequential(
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

                # if sys.stdin.isatty():
                #     while True:
                #         Console().print("\nPipeline finished. Choose an option:")
                #         Console().print("  [1] Show all error panels")
                #         Console().print("  [2] Show overview of failed nodes")
                #         Console().print("  [3] Exit to command line")
                #         choice = input("Enter 1/2/3: ").strip()
                #         if choice == "1":
                #             for panel in reporter.error_panels():
                #                 Console().print(panel)
                #         elif choice == "2":
                #             failed = [
                #                 n
                #                 for n in reporter._nodes.values()
                #                 if n.state == NodeState.FAILED
                #             ]
                #             if failed:
                #                 Console().print(
                #                     f"[bold red]Failed Nodes ({len(failed)}):[/bold red]"  # noqa: E501
                #                 )
                #                 for n in failed:
                #                     Console().print(f"- {n.name}")
                #             else:
                #                 Console().print("[green]No failed nodes[/green]")
                #         elif choice == "3":
                #             break
                #         else:
                #             Console().print(
                #                 "[yellow]Invalid choice, try again[/yellow]"
                #             )

                reporter.stop_live()

        if not dry_run:
            self.save_cache(self.signature_store_path, cache)

        if failures:
            if reporter:
                for error_panel in reporter.error_panels():
                    Console().print(error_panel)
            for name, key, err in failures:
                self.log(
                    None,
                    f"Element failed: {name} (key: {key}): {type(err).__name__}: {err.info.log_text}",  # noqa: E501
                    level="ERROR",
                )
            if final_raise_on_error:
                # lines = [
                #     f"{len(failures)} element(s) failed:",
                #     *[
                #         f"- {name}: {key} ({type(err).__name__}: {err})"
                #         for name, key, err in failures
                #     ],
                # ]
                raise RuntimeError("Pipeline failed with errors in elements: ")
        if log_signatures:
            self.log_sigs(order)
        # self.prune()

    def log_sigs(self, order) -> None:
        cache = self.load_cache(self.signature_store_path)
        for node in order:
            if node is not None:
                sig = cache.get(node.key, None)
                self.log(None, f"SIG {node.name}: {sig}")
        self.registry.write_registry(Path("test.reagistry.txt"))

    def prune(self) -> None:
        """Save the current state of the workspace to a file."""
        all_keys = set(self.registry.keys())
        cache = self.load_cache(self.signature_store_path)
        clean_cache = {}
        pruned_from_cache = {}
        for key in cache:
            cached_sig = cache.get(key, None)
            if key in all_keys:
                clean_cache[key] = cached_sig
            else:
                pruned_from_cache[key] = cached_sig
        outpruned = self.signature_store_path.with_suffix(".stale.json")
        self.save_cache(outpruned, pruned_from_cache)
        self.save_cache(
            self.signature_store_path.with_suffix(".clean.json"), clean_cache
        )

    ###########################################################################
    # Sequential execution
    ###########################################################################

    def check_and_stop(self, reporter) -> bool:
        if self.stop_event.is_set():
            self.log(reporter, "Pipeline aborted by user.", level="WARNING")
            if reporter:
                reporter.stop_live()
                return True
        return False

    def _evaluate_node(self, node, cache) -> tuple[bool, str]:
        cached_sig = cache.get(node.key)
        skip, reason = node.skip(cached_signature=cached_sig)
        return skip, reason

    def _log_node(self, node, reporter, skip, reason, log_run_only, verbose) -> str:
        message = self.compose_element_message(skip, node, reason, verbose)
        if (not skip) or (not log_run_only):
            self.log(reporter, message, level="INFO")
        return message

    def update_cache(self, cache, node) -> None:
        cache[node.key] = node.signature
        self.save_cache(self.signature_store_path, cache)

    def _run_node(self, node, cache, reporter, continue_on_error, failures):
        try:
            node()
            # update cache immediately after each node finishes, so we don't lose progress if the pipeline crashes later  # noqa: E501
            self.update_cache(cache, node)
            if reporter:
                reporter.mark_done(node.key)
        except PipelineError as e:
            e.info.add_creation_trace(node)
            # Console().print(e.info.panel)
            self.log(reporter, f"{node.name} failed!", level="ERROR")
            if reporter:
                reporter.push_error_panel(e.info.panel)
                reporter.mark_failed(node.key)
            else:
                Console().print(e.info.panel)
            failures.append((node.name, node.key, e))
            if not continue_on_error:
                self.save_cache(self.signature_store_path, cache)
                raise

        except Exception as e:
            err = PipelineError(e, phase="EXECUTION")
            self.log(
                reporter,
                f"Error: {node.name}: failed with {type(e).__name__}: {e}",
                level="ERROR",
            )
            if reporter:
                reporter.push_error_panel(err)
                reporter.mark_failed(node.key)
            else:
                Console().print(err.info.panel)
            failures.append((node.name, node.key, err))
            if not continue_on_error:
                self.save_cache(self.signature_store_path, cache)
                raise
        return cache, reporter, failures

    def write_element(self, node) -> None:
        with self.elements_file.open("a", encoding="utf-8") as f:
            f.write(f"{node.key}\t{node.name}\n{node.describe()}\n")
            f.flush()
            os.fsync(f.fileno())

    def _run_sequential(
        self,
        order: list[Element],
        *,
        cache: dict[str, str],
        reporter: ProgressReporter | None,
        continue_on_error: bool,
        failures: list[tuple[str, str, PipelineError]],
        dry_run: bool,
        log_run_only: bool,
        verbose: bool,
    ) -> tuple[list[tuple[str, str, PipelineError]], ProgressReporter | None]:
        sigstore = self.load_sigstore(self.signature_data_path)
        # Keys of nodes that failed or were blocked by an upstream failure.
        # Used to propagate UPSTREAM_FAILED to their dependents.
        blocked_keys: set[str] = set()
        reasons_to_run = {}
        need_to_run: list[Element] = []
        skipped: set[str] = set()
        for node in order:
            sigdata = dict(node.sig_data())  #
            sigdata["name"] = node.tag.default_name
            sigstore[node.key] = sigdata
            self.save_sigstore(self.signature_data_path, sigstore)
            # check if the node must run or skipped
            skip, reason = self._evaluate_node(node, cache)
            if skip:
                skipped.add(node.key)
                # log skipped node
                self._log_node(node, reporter, skip, reason, log_run_only, verbose)
                if reporter:
                    reporter.mark_skip(node.key, reason)
                continue
            else:
                reasons_to_run[node.key] = reason
                need_to_run.append(node)
            # check if we need to stop
            if self.check_and_stop(reporter):
                break
        for node in need_to_run:
            # check if a direct prerequisite failed or was blocked
            upstream_failed = any(pre.key in blocked_keys for pre in node.pres)
            if upstream_failed and not dry_run:
                blocked_keys.add(node.key)
                self._log_node(
                    node, reporter, True, "upstream failed", log_run_only, verbose
                )
                if reporter:
                    reporter.mark_upstream_failed(node.key, "upstream failed")
                continue
            # log running node
            self._log_node(
                node, reporter, False, reasons_to_run[node.key], log_run_only, verbose
            )
            if dry_run:
                continue

            if reporter:
                reporter.mark_start(node.key)
            cache, reporter, failures = self._run_node(
                node, cache, reporter, continue_on_error, failures
            )
            # If the node just failed, add it to blocked_keys so its
            # dependents are skipped with UPSTREAM_FAILED.
            if failures and failures[-1][1] == node.key:
                blocked_keys.add(node.key)
        return failures, reporter

    ###########################################################################
    # Parallel execution  (DAG-aware: runs independent nodes concurrently)
    ###########################################################################

    def _init_pending_and_successors(
        self, order
    ) -> tuple[dict[str, Element], dict[str, int], dict[str, list[str]]]:
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
        return by_key, pending_deps, successors

    def _run_parallel(
        self,
        order: list[Element],
        *,
        cache: dict[str, str],
        reporter: ProgressReporter | None,
        continue_on_error: bool,
        failures: list[tuple[str, str, PipelineError]],
        max_workers: int | None,
        log_run_only: bool,
        verbose: bool,
    ) -> tuple[list[tuple[str, str, PipelineError]], ProgressReporter | None]:

        by_key, pending_deps, successors = self._init_pending_and_successors(order)
        cache_lock = threading.Lock()
        done_event = threading.Event()
        active: set[str] = set()
        completed: set[str] = set()
        # Keys of nodes that failed or were blocked — their successors get
        # UPSTREAM_FAILED instead of being submitted.
        blocked_keys: set[str] = set()
        lock = threading.Lock()
        abort = self.stop_event

        # nodes whose deps are all done
        ready: list[str] = [k for k, d in pending_deps.items() if d == 0]

        # get available workers
        workers = max_workers or max(1, (self.resources.threads // 4))

        # let the ThreadPoolExecutionor handle the threads
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures: dict[Future, str] = {}

            def _mark_successors_blocked(key: str) -> None:
                """Recursively mark all transitive successors as upstream-failed.
                Must be called with *lock* held."""
                for s in successors.get(key, []):
                    if s in blocked_keys or s in completed:
                        continue
                    blocked_keys.add(s)
                    completed.add(s)  # treat as "done" so the done_event fires
                    self.log(
                        reporter,
                        f"SKIP {by_key[s].name}  (upstream failed)",
                        level="WARNING",
                    )
                    if reporter:
                        reporter.mark_upstream_failed(s, "upstream failed")
                    _mark_successors_blocked(s)

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
                skip, reason = self._evaluate_node(node, cache)
                self._log_node(node, reporter, skip, reason, log_run_only, verbose)
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
                        if not isinstance(e, PipelineError):
                            e = PipelineError(e, phase="EXECUTION")
                        failures.append((node.name, node.key, e))

                        if reporter:
                            reporter.mark_failed(node.key)
                        self.log(reporter, f"Error: {node.name}: {e}", level="ERROR")
                        # Mark all transitive successors as upstream-failed.
                        with lock:
                            blocked_keys.add(key)
                            _mark_successors_blocked(key)
                            if len(completed) + len(active) - 1 >= len(order):
                                done_event.set()
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
        return failures, reporter

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

    def ordered_nodes_for_run(self) -> list[Element]:
        all_keys = set(self.registry.keys())
        new_keys = (
            all_keys
            if self._baseline_keys is None
            else all_keys - self._baseline_keys  # noqa: E501
        )
        targets = [e for k in new_keys if (e := self.registry.get(k)) is not None]
        nodes = self.collect(targets)
        self.check_duplicate_outputs(nodes)
        order = self.toposort(nodes)
        return order

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

    def capture_store_state(
        self, outstore: Path, accept_existing: bool = False
    ) -> None:
        """Save the current state of the workspace to a file."""
        # all_keys = set(self.registry.keys())
        # new_keys = (
        #     all_keys
        #     if self._baseline_keys is None
        #     else all_keys - self._baseline_keys  # noqa: E501
        # )
        # targets = [e for k in new_keys if (e := self.registry.get(k)) is not None]
        # cache = self.load_cache(self.signature_store_path)
        # nodes = self.collect(targets)
        # self.check_duplicate_outputs(nodes)
        # order = self.toposort(nodes)
        sig_data_store = outstore.with_suffix(".data.json")
        order = self.ordered_nodes_for_run()
        cache = self.load_cache(self.signature_store_path)
        data_cache = self.load_sigstore(self.signature_data_path)
        for node in order:
            cached_sig = cache.get(node.key)
            skip, reason = node.skip(cached_signature=cached_sig)
            sigdata = dict(node.sig_data())  #
            sigdata["name"] = node.tag.default_name
            data_cache[node.key] = sigdata
            self.save_sigstore(sig_data_store, data_cache)

            if not skip and reason == "First run":
                if node.outputs_ok():
                    if accept_existing:
                        cache[node.key] = node.signature
                        self.save_cache(outstore, cache)
                        self.log(
                            None,
                            f"Node {node.key} is being captured for the first time, but its outputs already exist and are valid. Accepting existing outputs as-is.",  # noqa: E501
                            level="WARNING",
                        )
                    else:
                        raise RuntimeError(
                            f"Node {node.key} is being captured for the first time, but its outputs already exist and are valid. This likely means the node was run before without the cache, or the cache was cleared."  # noqa: E501
                        )
            elif not skip and reason == "Cached signature does not match":
                _, check_result = explain_signature_diff(node, data_cache)
                self.log(
                    None,
                    f"Node {node.key} cache signature is different:\n{check_result}",  # noqa: E501
                    level="WARNING",
                )

                if node.outputs_ok():
                    if accept_existing:
                        cache[node.key] = node.signature
                        self.save_cache(outstore, cache)
                        self.log(
                            None,
                            f"Node {node.key} cache signature differs, but its outputs already exist and are valid. Accepting existing outputs as-is.",  # noqa: E501
                            level="WARNING",
                        )
                    else:
                        raise RuntimeError(
                            f"Node {node.key} cache signature differs, but its outputs already exist and are valid. This likely means the node was run before without the cache, or the cache was cleared."  # noqa: E501
                        )

            elif skip:
                cache[node.key] = node.signature
                self.save_cache(outstore, cache)

    def merge_signature_stores(
        self,
        store_paths: list[Path],
        out_path: Path,
        conflict: Literal["override", "raise"] = "raise",
    ) -> None:
        """Merge multiple signature stores into one, checking for conflicts."""
        merged: dict[str, str] = {}
        for path in store_paths:
            cache = self.load_cache(path)
            data_cache = self.load_sigstore(path.with_suffix(".data.json"))
            for key, sig in cache.items():
                if key in merged and merged[key] != sig:
                    if conflict == "override":
                        self.log(
                            None,
                            f"Conflict for key {key}: {merged[key]} vs {sig} (from {path}). Overriding with new value.",  # noqa: E501
                            level="WARNING",
                        )
                    else:
                        raise RuntimeError(
                            f"Conflict for key {key}: {merged[key]} vs {sig} (from {path})"  # noqa: E501
                        )
                merged[key] = sig
        self.save_cache(out_path, merged)
        self.save_sigstore(out_path.with_suffix(".data.json"), data_cache)

    def update_store(
        self, store_paths: list[Path], conflict: Literal["override", "raise"] = "raise"
    ) -> None:
        """Update a signature store with new key-signature pairs."""
        out_path = self.signature_store_path
        store_paths = [out_path] + store_paths
        self.merge_signature_stores(store_paths, out_path, conflict)

    ############################################################################
    # Reporting
    ############################################################################
    # def write_report(self, order: list[Element], path: Path):
    #     with path.open("w", encoding="utf-8") as f:
    #         f.write("FILE LINEAGE\n")
    #         f.write("=" * 80 + "\n")

    #         for file, producer in produced_by.items():
    #             f.write(f"{file}  <--  {producer}\n")

    #         f.write("PIPELINE REPORT\n")
    #         f.write("=" * 80 + "\n\n")

    #         for i, node in enumerate(order, 1):
    #             f.write(f"[{i}] {node.name}\n")
    #             f.write(f"  Key: {node.key}\n")

    #             # Inputs
    #             f.write("  Inputs:\n")
    #             for inp in node.inputs:
    #                 f.write(f"    - {inp}\n")

    #             # Outputs
    #             f.write("  Outputs:\n")
    #             for name, out in node.artifacts.items():
    #                 f.write(f"    - {name}: {out}\n")

    #             # Dependencies
    #             f.write("  Depends on:\n")
    #             for pre in node.pres:
    #                 f.write(f"    - {pre.name}\n")

    #             f.write("\n")

    # def write_execution_order(self, order, path):
    # with open(path, "w") as f:
    #     for i, node in enumerate(order, 1):
    #         f.write(f"{i:03d}  {node.name}\n")

    # def render_dag(self, dot_path: Path):
    #     subprocess.run(["dot", "-Tpng", str(dot_path), "-o", str(dot_path.with_suffix(".png"))])

    # def build_lineage(self, nodes):
    #     produced_by = {}

    #     for node in nodes:
    #         for out in node.output_files:
    #             produced_by[str(out)] = node.name

    #     return produced_by
