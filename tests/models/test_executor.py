from __future__ import annotations

import importlib
import io
import json
import logging
import sys
import types
import typing
from pathlib import Path

import pytest

# Keep module imports isolated from package-level side effects in this repository.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
MMALIGNMENTS_ROOT = SRC_ROOT / "mmalignments"


def _ensure_namespace_package(name: str, path: Path) -> None:
    if name in sys.modules:
        return
    module = types.ModuleType(name)
    module.__path__ = [str(path)]
    sys.modules[name] = module


_ensure_namespace_package("mmalignments", MMALIGNMENTS_ROOT)
_ensure_namespace_package("mmalignments.models", MMALIGNMENTS_ROOT / "models")
_ensure_namespace_package("mmalignments.services", MMALIGNMENTS_ROOT / "services")

# elements.py applies @runtime_checkable to non-Protocol classes; neutralize it.
typing.runtime_checkable = lambda cls: cls

elements_module = importlib.import_module("mmalignments.models.elements")
if not hasattr(elements_module, "generate_element_key_name"):
    setattr(elements_module, "generate_element_key_name", lambda *args, **kwargs: "key")

executor_module = importlib.import_module("mmalignments.models.executor")
resources_module = importlib.import_module("mmalignments.models.resources")

Executor = executor_module.Executor
NodeState = executor_module.NodeState


class DummyTag:
    def __init__(self, default_name: str):
        self.default_name = default_name


class DummyRun:
    def __init__(self, command: list[str] | None = None, display: str | None = None):
        self.command = command or ["echo", "ok"]
        self.command_display = display

    def __call__(self):
        return True


class DummyInfo:
    def __init__(self):
        self.panel = "panel"
        self.log_text = "log_text"
        self.traces: list[str] = []

    def add_creation_trace(self, node):
        self.traces.append(node.key)


class DummyPipelineError(Exception):
    def __init__(self, cause: Exception | str = "x", phase: str | None = None):
        super().__init__(str(cause))
        self.info = DummyInfo()


class DummyReporter:
    def __init__(self):
        self.logs: list[str] = []
        self.started = False
        self.stopped = False
        self._nodes: dict[str, types.SimpleNamespace] = {}
        self.upstream_failed: list[tuple[str, str]] = []
        self.failed: list[str] = []
        self.skipped: list[tuple[str, str]] = []
        self.done: list[str] = []
        self.started_nodes: list[str] = []
        self.error_panels_list: list[str] = []

    def register(self, nodes):
        for n in nodes:
            self._nodes[n.key] = types.SimpleNamespace(state=NodeState.PENDING)

    def start_live(self):
        self.started = True

    def stop_live(self):
        self.stopped = True

    def push_log(self, message: str):
        self.logs.append(message)

    def mark_skip(self, key: str, reason: str = ""):
        self.skipped.append((key, reason))
        self._nodes[key].state = NodeState.SKIPPED

    def mark_start(self, key: str):
        self.started_nodes.append(key)
        self._nodes[key].state = NodeState.RUNNING

    def mark_done(self, key: str):
        self.done.append(key)
        self._nodes[key].state = NodeState.DONE

    def mark_failed(self, key: str):
        self.failed.append(key)
        self._nodes[key].state = NodeState.FAILED

    def mark_upstream_failed(self, key: str, reason: str = "upstream failed"):
        self.upstream_failed.append((key, reason))
        self._nodes[key].state = NodeState.UPSTREAM_FAILED

    def push_error_panel(self, panel):
        self.error_panels_list.append(panel)

    def error_panels(self):
        return list(self.error_panels_list)


class DummyNode:
    def __init__(
        self,
        key: str,
        *,
        pres: tuple[DummyNode, ...] = (),
        skip_result: tuple[bool, str] = (True, "cached"),
        output_files: tuple[Path, ...] = (),
        outputs_ok_result: bool = True,
        signature: str | None = None,
        sig_data: dict | None = None,
        run: DummyRun | None = None,
        call_error: Exception | None = None,
    ):
        self.key = key
        self.name = f"node_{key}"
        self.threads = 2
        self.pres = tuple(pres)
        self.output_files = tuple(output_files)
        self.files = tuple(output_files)
        self.artifacts = {"primary": "artifact"}
        self.tag = DummyTag(default_name=f"tag_{key}")
        self.provenance = f"prov_{key}"
        self.default_output_file = f"default_{key}.txt"
        self.run = run or DummyRun()
        self.skip_result = skip_result
        self._sig_data = sig_data or {
            "key": key,
            "determinants": [],
            "inputs": [],
            "artifacts": {},
            "pre_sigs": [],
        }
        self.signature = signature or f"sig_{key}"
        self._outputs_ok_result = outputs_ok_result
        self._call_error = call_error
        self.called = 0
        self.creation_trace = f"trace_{key}"

    def skip(self, cached_signature=None, cached_sig_data=None):
        return self.skip_result

    def sig_data(self):
        return dict(self._sig_data)

    def outputs_ok(self):
        return self._outputs_ok_result

    def describe(self):
        return f"describe_{self.key}"

    def __call__(self):
        self.called += 1
        if self._call_error is not None:
            raise self._call_error
        return True


def _logger() -> logging.Logger:
    logger = logging.getLogger("test.executor")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)
    return logger


def _executor(tmp_path: Path) -> Executor:
    # Every test uses an isolated cache/log directory under pytest temp path.
    return Executor(
        logger=_logger(),
        cache_path=tmp_path / "run_cache",
    )


def test_init_creates_dirs_and_defaults(tmp_path: Path):
    """Validate constructor bootstraps directories and baseline state."""
    ex = _executor(tmp_path)
    assert ex.main_dir.exists()
    assert ex.cache_dir.exists()
    assert ex.log_dir.exists()
    assert ex.tool_log_dir.exists()
    assert ex.state == "init"
    assert ex._baseline_keys is None


def test_cleanup_old_run_logs_keeps_latest_and_ignores_unlink_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Ensure log pruning keeps newest files and never crashes on unlink failures."""
    ex = _executor(tmp_path)
    for i in range(8):
        (ex.log_dir / f"run_{i}.log").write_text("x", encoding="utf-8")
        (ex.log_dir / f"order_{i}.txt").write_text("y", encoding="utf-8")

    target = ex.log_dir / "run_0.log"
    original_unlink = Path.unlink

    def flaky_unlink(self):
        if self == target:
            raise OSError("boom")
        return original_unlink(self)

    monkeypatch.setattr(Path, "unlink", flaky_unlink)
    ex._cleanup_old_run_logs(keep=5)
    assert len(list(ex.log_dir.glob("run_*.log"))) <= 6
    assert len(list(ex.log_dir.glob("order_*.txt"))) <= 5


def test_interrupt_and_keyboard_listener_set_stop_event(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Check both signal and keyboard-triggered stop paths."""
    ex = _executor(tmp_path)
    warnings: list[str] = []
    monkeypatch.setattr(ex.logger, "warning", lambda m: warnings.append(m))
    ex._handle_interrupt(15, None)
    assert ex.stop_event.is_set()
    assert any("Interrupt received" in w for w in warnings)

    ex2 = _executor(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO("q\n"))
    warnings2: list[str] = []
    monkeypatch.setattr(ex2.logger, "warning", lambda m: warnings2.append(m))
    ex2._keyboard_listener()
    assert ex2.stop_event.is_set()
    assert any("User requested stop" in w for w in warnings2)


def test_build_context_sets_contextvars_and_external_logger_hooks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Verify build() sets resources/tool-log context and calls External logger hooks."""
    ex = _executor(tmp_path)
    calls: list[str] = []
    monkeypatch.setattr(
        executor_module.External, "add_main_logger", lambda logger: calls.append("add")
    )
    monkeypatch.setattr(
        executor_module.External, "remove_main_logger", lambda: calls.append("remove")
    )

    with ex.build():
        current = resources_module._current_resources.get()
        assert current is ex.resources
        assert resources_module._current_tool_log_dir.get() == ex.tool_log_dir

    assert calls == ["add", "remove"]
    assert resources_module._current_tool_log_dir.get() is None


def test_record_and_play_selects_new_targets_and_handles_empty_registry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Ensure record/play baseline filtering works and empty runs are logged."""
    ex = _executor(tmp_path)
    messages: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": messages.append(message)
    )

    ex.record()
    ex.play()
    assert any("No Elements found to run" in m for m in messages)

    node_old = DummyNode("old")
    node_new = DummyNode("new")
    ex.registry.intern(node_old)
    ex.record()
    ex.registry.intern(node_new)

    captured = {}

    def fake_run_targets(targets, **kwargs):
        captured["targets"] = targets
        captured["kwargs"] = kwargs

    monkeypatch.setattr(ex, "run_targets", fake_run_targets)
    ex.play(log_run_only=True, progress=False)
    assert captured["targets"] == [node_new]
    assert captured["kwargs"]["log_run_only"] is True


def test_load_and_save_cache_and_sigstore_roundtrip(tmp_path: Path):
    """Cover JSON cache/sigstore load-save code paths including missing-file branch."""
    ex = _executor(tmp_path)
    missing = tmp_path / "nope.json"
    assert ex.load_cache(missing) == {}
    assert ex.load_sigstore(missing) == {}

    cache_path = tmp_path / "c.json"
    ex.save_cache(cache_path, {"a": "1", "b": "2"})
    assert ex.load_cache(cache_path) == {"a": "1", "b": "2"}

    sigstore_path = tmp_path / "s.json"
    ex.save_sigstore(sigstore_path, {"x": {"k": 1}})
    assert ex.load_sigstore(sigstore_path) == {"x": {"k": 1}}


def test_msg_and_log_paths_for_levels_and_reporter(tmp_path: Path):
    """Check log routing for INFO/WARNING/ERROR/DEBUG and reporter redirection."""
    ex = _executor(tmp_path)
    seen: list[tuple[str, str]] = []
    ex.logger = types.SimpleNamespace(
        info=lambda m: seen.append(("info", m)),
        warning=lambda m: seen.append(("warning", m)),
        error=lambda m: seen.append(("error", m)),
        debug=lambda m: seen.append(("debug", m)),
    )

    ex.msg("i", level="INFO")
    ex.msg("w", level="WARNING")
    ex.msg("e", level="ERROR")
    ex.msg("d", level="OTHER")
    assert seen == [("info", "i"), ("warning", "w"), ("error", "e"), ("debug", "d")]

    reporter = DummyReporter()
    ex.log(reporter, "hello", level="INFO")
    assert reporter.logs == ["[INFO] hello"]


def test_compose_element_message_skip_and_verbose_details(tmp_path: Path):
    """Validate run/skip message formatting and verbose branch details."""
    ex = _executor(tmp_path)
    run = DummyRun(command=["echo", "run"], display="echo run")
    node = DummyNode("a", run=run)

    # Non-skip includes RUN and thread suffix.
    plain = ex.compose_element_message(False, node, "why", verbose=False)
    assert plain.startswith("RUN ")
    assert "[2t]" in plain

    # Verbose mode emits key, outputs, artifacts, tag, pres, and command lines.
    verbose = ex.compose_element_message(False, node, "why", verbose=True)
    assert "key:" in verbose
    assert "actual output:" in verbose
    assert "artifacts:" in verbose
    assert "cmd:" in verbose

    # Skip mode emits SKIP and no thread suffix.
    skipped = ex.compose_element_message(True, node, "cached", verbose=True)
    assert skipped.startswith("SKIP")


def test_check_and_stop_behavior(tmp_path: Path):
    """Ensure check_and_stop reports abort and handles reporter stop path."""
    ex = _executor(tmp_path)
    reporter = DummyReporter()
    ex.stop_event.set()
    assert ex.check_and_stop(reporter) is True
    assert reporter.stopped is True

    ex2 = _executor(tmp_path)
    assert ex2.check_and_stop(None) is False


def test_evaluate_log_update_and_write_element_helpers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover helper methods used by sequential/parallel execution."""
    ex = _executor(tmp_path)
    node = DummyNode("n")

    # _evaluate_node delegates to node.skip using cache key lookup.
    skip, reason = ex._evaluate_node(node, {"n": "sig"}, {"x": 1})
    assert (skip, reason) == node.skip_result

    # _log_node returns compose output and conditionally emits via log().
    logs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": logs.append(message)
    )
    msg = ex._log_node(node, None, True, "cached", log_run_only=False, verbose=False)
    assert "SKIP" in msg
    assert logs

    # update_cache persists both in-memory map and file.
    cache = {}
    ex.update_cache(cache, node)
    assert cache[node.key] == node.signature
    on_disk = json.loads(ex.signature_store_path.read_text(encoding="utf-8"))
    assert on_disk[node.key] == node.signature

    # write_element appends entry and fsyncs.
    ex.write_element(node)
    text = ex.elements_file.read_text(encoding="utf-8")
    assert node.key in text
    assert node.describe() in text


def test_run_node_success_pipelineerror_and_generic_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Exercise _run_node for success, PipelineError branch, and generic exception wrapping."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    reporter = DummyReporter()
    cache = {}
    failures: list[tuple[str, str, Exception]] = []

    ok_node = DummyNode("ok")
    reporter.register([ok_node])
    cache, reporter, failures = ex._run_node(ok_node, cache, reporter, True, failures)
    assert ok_node.key in cache
    assert reporter.done == [ok_node.key]
    assert failures == []

    pe_node = DummyNode("pe", call_error=DummyPipelineError("boom"))
    reporter.register([pe_node])
    cache, reporter, failures = ex._run_node(pe_node, cache, reporter, True, failures)
    assert failures[-1][0] == pe_node.name
    assert reporter.failed[-1] == pe_node.key

    ge_node = DummyNode("ge", call_error=RuntimeError("oops"))
    reporter.register([ge_node])
    cache, reporter, failures = ex._run_node(ge_node, cache, reporter, True, failures)
    assert failures[-1][0] == ge_node.name
    assert reporter.failed[-1] == ge_node.key


def test_run_node_raises_when_continue_on_error_false(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Verify _run_node re-raises errors when continue_on_error is disabled."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)
    with pytest.raises(DummyPipelineError):
        ex._run_node(
            DummyNode("x", call_error=DummyPipelineError("fail")),
            {},
            None,
            False,
            [],
        )


def test_init_pending_and_successors_counts_dependencies(tmp_path: Path):
    """Confirm dependency count and successor map initialization for DAG scheduling."""
    ex = _executor(tmp_path)
    a = DummyNode("a")
    b = DummyNode("b", pres=(a,))
    c = DummyNode("c", pres=(a, b))

    by_key, pending, succ = ex._init_pending_and_successors([a, b, c])
    assert set(by_key) == {"a", "b", "c"}
    assert pending["a"] == 0
    assert pending["b"] == 1
    assert pending["c"] == 2
    assert succ["a"] == ["b", "c"]


def test_run_sequential_skip_run_and_upstream_blocking(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover sequential planner for skipped nodes, runnable nodes, and upstream-failed propagation."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    # Build a 2-node chain where first fails and second should be blocked.
    n1 = DummyNode(
        "n1", skip_result=(False, "First run"), call_error=RuntimeError("boom")
    )
    n2 = DummyNode("n2", pres=(n1,), skip_result=(False, "First run"))
    order = [n1, n2]

    reporter = DummyReporter()
    reporter.register(order)
    failures: list[tuple[str, str, Exception]] = []

    failures, rep = ex._run_sequential(
        order,
        cache={},
        reporter=reporter,
        continue_on_error=True,
        failures=failures,
        dry_run=False,
        log_run_only=False,
        verbose=False,
    )
    assert failures
    assert any(k == "n2" for k, _reason in rep.upstream_failed)


def test_run_parallel_executes_ready_nodes_and_records_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Exercise DAG-aware parallel scheduler on a successful multi-node run."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    # Two runnable nodes without failures should complete and be marked done.
    a = DummyNode("a", skip_result=(False, "First run"))
    b = DummyNode("b", skip_result=(False, "First run"))
    order = [a, b]
    reporter = DummyReporter()
    reporter.register(order)

    failures, rep = ex._run_parallel(
        order,
        cache={},
        reporter=reporter,
        continue_on_error=True,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []
    assert set(rep.done) == {"a", "b"}


def test_check_duplicate_outputs_detects_conflicts(tmp_path: Path):
    """Verify duplicate output detection for conflicting and non-conflicting nodes."""
    ex = _executor(tmp_path)
    p = tmp_path / "x.txt"
    n1 = DummyNode("1", output_files=(p,))
    n2 = DummyNode("2", output_files=(p,))
    with pytest.raises(ValueError, match="output path conflict"):
        ex.check_duplicate_outputs([n1, n2])

    ex.check_duplicate_outputs(
        [
            DummyNode("3", output_files=(tmp_path / "a",)),
            DummyNode("4", output_files=(tmp_path / "b",)),
        ]
    )


def test_collect_handles_identity_and_duplicate_key_conflict(tmp_path: Path):
    """Test transitive prerequisite collection and duplicate-key conflict guard."""
    ex = _executor(tmp_path)
    a = DummyNode("a")
    b = DummyNode("b", pres=(a,))
    nodes = ex.collect([b])
    assert {n.key for n in nodes} == {"a", "b"}

    # Same key but distinct object is rejected.
    x1 = DummyNode("dup")
    x2 = DummyNode("dup")
    y = DummyNode("y", pres=(x1, x2))
    with pytest.raises(
        ValueError, match="Duplicate key with different Element objects"
    ):
        ex.collect([y])


def test_toposort_orders_nodes_and_detects_cycle(tmp_path: Path):
    """Validate deterministic topological order and cycle detection branch."""
    ex = _executor(tmp_path)
    a = DummyNode("a")
    b = DummyNode("b", pres=(a,))
    c = DummyNode("c", pres=(b,))
    order = ex.toposort([c, b, a])
    assert [n.key for n in order] == ["a", "b", "c"]

    # Create a cycle explicitly: a <- b <- a.
    a.pres = (b,)
    b.pres = (a,)
    with pytest.raises(RuntimeError, match="Dependency cycle"):
        ex.toposort([a, b])


def test_ordered_nodes_for_run_uses_baseline_and_dag_helpers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Ensure ordered_nodes_for_run uses baseline filtering and helper pipeline."""
    ex = _executor(tmp_path)
    old = DummyNode("old")
    new = DummyNode("new")
    ex.registry.intern(old)
    ex.record()
    ex.registry.intern(new)

    monkeypatch.setattr(ex, "collect", lambda targets: targets)
    monkeypatch.setattr(ex, "check_duplicate_outputs", lambda nodes: None)
    monkeypatch.setattr(ex, "toposort", lambda nodes: list(nodes))

    out = ex.ordered_nodes_for_run()
    assert out == [new]


def test_to_dot_and_write_dot_emit_graph(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover DOT graph serialization and file writing path."""
    ex = _executor(tmp_path)
    a = DummyNode("a")
    b = DummyNode("b", pres=(a,))

    dot = ex.to_dot([b])
    assert "digraph DAG" in dot
    assert '"a" -> "b"' in dot

    log_msgs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": log_msgs.append(message)
    )
    outpath = ex.write_dot([b], tmp_path / "dag" / "graph.dot")
    assert outpath.exists()
    assert "Wrote DAG dot file" in log_msgs[-1]


def test_run_targets_dispatches_to_sequential_and_parallel_and_raises_on_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Test run_targets orchestration for sequential/parallel paths, failure reporting, and final raise."""
    ex = _executor(tmp_path)
    nodes = [DummyNode("a")]

    monkeypatch.setattr(ex, "collect", lambda targets: nodes)
    monkeypatch.setattr(ex, "toposort", lambda n: list(n))
    monkeypatch.setattr(ex, "check_duplicate_outputs", lambda n: None)
    monkeypatch.setattr(ex, "write_dot", lambda n, p: p)
    monkeypatch.setattr(ex.registry, "write_registry", lambda path: None)

    # Sequential path with no failures should complete and save cache.
    monkeypatch.setattr(ex, "_run_sequential", lambda *args, **kwargs: ([], None))
    ex.run_targets(nodes, progress=False, dry_run=False, parallel=False)

    # Parallel path with failures and final_raise_on_error should raise RuntimeError.
    dummy_err = DummyPipelineError("e")
    monkeypatch.setattr(
        ex, "_run_parallel", lambda *args, **kwargs: ([("n", "k", dummy_err)], None)
    )
    with pytest.raises(RuntimeError, match="Pipeline failed"):
        ex.run_targets(
            nodes,
            progress=False,
            dry_run=False,
            parallel=True,
            final_raise_on_error=True,
        )


def test_log_sigs_and_prune_cache_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Verify signature logging and cache pruning outputs."""
    ex = _executor(tmp_path)
    n = DummyNode("n")
    ex.save_cache(ex.signature_store_path, {"n": "sig_n", "stale": "old"})
    ex.registry.intern(n)

    msgs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": msgs.append(message)
    )
    monkeypatch.setattr(
        ex.registry, "write_registry", lambda path: msgs.append(f"wrote:{path}")
    )
    ex.log_sigs([n])
    assert any("SIG" in m for m in msgs)

    ex.prune()
    assert ex.signature_store_path.with_suffix(".stale.json").exists()
    assert ex.signature_store_path.with_suffix(".clean.json").exists()


def test_capture_store_state_branches_for_first_run_mismatch_and_skip(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover capture_store_state branches: first-run accept/reject, mismatch handling, and skipped updates."""
    ex = _executor(tmp_path)

    # Prepare nodes that hit different capture branches.
    n_first = DummyNode(
        "first", skip_result=(False, "First run"), outputs_ok_result=True
    )
    n_diff = DummyNode(
        "diff",
        skip_result=(False, "Cached signature does not match"),
        outputs_ok_result=True,
    )
    n_skip = DummyNode("skip", skip_result=(True, "cached"))

    monkeypatch.setattr(ex, "ordered_nodes_for_run", lambda: [n_first, n_diff, n_skip])
    monkeypatch.setattr(
        executor_module,
        "explain_signature_diff",
        lambda node, store: (False, "details"),
    )

    outstore = tmp_path / "sigstore.json"
    ex.save_cache(ex.signature_store_path, {})
    ex.save_sigstore(ex.signature_data_path, {})

    # Rejecting existing outputs on first-run must raise.
    with pytest.raises(RuntimeError, match="being captured for the first time"):
        ex.capture_store_state(outstore, accept_existing=False)

    # Accepting existing outputs should update caches without raising.
    ex.capture_store_state(outstore, accept_existing=True)
    data = json.loads(outstore.read_text(encoding="utf-8"))
    assert data["first"] == n_first.signature
    assert data["diff"] == n_diff.signature
    assert data["skip"] == n_skip.signature


def test_merge_signature_stores_conflict_override_and_raise(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Check merge conflict behavior for override and strict-raise modes."""
    ex = _executor(tmp_path)
    s1 = tmp_path / "s1.json"
    s2 = tmp_path / "s2.json"

    ex.save_cache(s1, {"k": "v1"})
    ex.save_sigstore(s1.with_suffix(".data.json"), {"k": {"a": 1}})
    ex.save_cache(s2, {"k": "v2", "x": "z"})
    ex.save_sigstore(s2.with_suffix(".data.json"), {"k": {"a": 2}, "x": {"a": 3}})

    out = tmp_path / "out.json"
    warnings: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": warnings.append(message)
    )
    ex.merge_signature_stores([s1, s2], out, conflict="override")
    merged = json.loads(out.read_text(encoding="utf-8"))
    assert merged["k"] == "v2"
    assert "x" in merged
    assert any("Conflict for key" in w for w in warnings)

    with pytest.raises(RuntimeError, match="Conflict for key"):
        ex.merge_signature_stores([s1, s2], out, conflict="raise")


def test_update_store_prepends_default_store(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Ensure update_store merges default store plus provided stores in order."""
    ex = _executor(tmp_path)
    seen = {}

    def fake_merge(paths, out_path, conflict):
        seen["paths"] = paths
        seen["out"] = out_path
        seen["conflict"] = conflict

    monkeypatch.setattr(ex, "merge_signature_stores", fake_merge)
    extra = tmp_path / "extra.json"
    ex.update_store([extra], conflict="override")
    assert seen["paths"][0] == ex.signature_store_path
    assert seen["paths"][1] == extra
    assert seen["out"] == ex.signature_store_path
    assert seen["conflict"] == "override"


def test_constructor_non_main_thread_skips_signal_registration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover constructor branch where current thread is not the main thread."""

    class DummyThread:
        pass

    monkeypatch.setattr(
        executor_module.threading,
        "current_thread",
        lambda: DummyThread(),
    )
    called = {"count": 0}
    monkeypatch.setattr(
        executor_module.signal,
        "signal",
        lambda *args, **kwargs: called.__setitem__("count", called["count"] + 1),
    )
    _ = _executor(tmp_path)
    assert called["count"] == 0


def test_keyboard_listener_non_q_and_empty_input_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Exercise listener branch for non-q lines and immediate EOF."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO("x\n"))
    ex._keyboard_listener()
    assert ex.stop_event.is_set() is False

    ex2 = _executor(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))
    ex2._keyboard_listener()
    assert ex2.stop_event.is_set() is False


def test_compose_element_message_display_and_provenance_verbose_branch(tmp_path: Path):
    """Hit verbose sub-branch where command_display differs from command string."""
    ex = _executor(tmp_path)
    run = DummyRun(command=["echo", "run"], display="python run.py --display")
    node = DummyNode("branch", run=run)
    msg = ex.compose_element_message(False, node, "reason", verbose=True)
    assert "display:" in msg
    assert "prov:" in msg


def test_run_targets_abort_and_reporter_lifecycle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover run_targets early-abort and progress reporter lifecycle paths."""
    ex = _executor(tmp_path)
    ex.stop_event.set()
    logs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": logs.append(message)
    )
    ex.run_targets([DummyNode("a")], progress=False)
    assert any("already aborted before start" in m for m in logs)

    ex2 = _executor(tmp_path)
    node = DummyNode("n")
    monkeypatch.setattr(ex2, "collect", lambda targets: [node])
    monkeypatch.setattr(ex2, "toposort", lambda nodes: nodes)
    monkeypatch.setattr(ex2, "check_duplicate_outputs", lambda nodes: None)
    monkeypatch.setattr(ex2, "write_dot", lambda nodes, path: path)
    monkeypatch.setattr(ex2.registry, "write_registry", lambda path: None)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(
        executor_module.threading,
        "Thread",
        lambda *args, **kwargs: types.SimpleNamespace(start=lambda: None),
    )

    reporter = DummyReporter()
    monkeypatch.setattr(executor_module, "ProgressReporter", lambda logger: reporter)
    monkeypatch.setattr(ex2, "_run_sequential", lambda *args, **kwargs: ([], reporter))
    ex2.run_targets([node], progress=True, dry_run=False, parallel=False)
    assert reporter.started is True
    assert reporter.stopped is True


def test_run_targets_failure_reporting_with_reporter_and_log_signatures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover reporter error panel printing and log_signatures hook in run_targets."""
    ex = _executor(tmp_path)
    node = DummyNode("a")
    monkeypatch.setattr(ex, "collect", lambda targets: [node])
    monkeypatch.setattr(ex, "toposort", lambda nodes: nodes)
    monkeypatch.setattr(ex, "check_duplicate_outputs", lambda nodes: None)
    monkeypatch.setattr(ex, "write_dot", lambda nodes, path: path)
    monkeypatch.setattr(ex.registry, "write_registry", lambda path: None)

    reporter = DummyReporter()
    reporter.register([node])
    reporter.error_panels_list = ["panel1"]
    err = DummyPipelineError("boom")
    monkeypatch.setattr(executor_module, "ProgressReporter", lambda logger: reporter)
    monkeypatch.setattr(
        ex,
        "_run_sequential",
        lambda *args, **kwargs: ([(node.name, node.key, err)], reporter),
    )
    printed: list[str] = []
    monkeypatch.setattr(
        executor_module.Console,
        "print",
        lambda self, x=None: printed.append(str(x)),
    )
    calls = {"sigs": 0}
    monkeypatch.setattr(
        ex,
        "log_sigs",
        lambda order: calls.__setitem__("sigs", calls["sigs"] + 1),
    )

    ex.run_targets(
        [node],
        progress=True,
        dry_run=False,
        parallel=False,
        final_raise_on_error=False,
        log_signatures=True,
    )
    assert printed
    assert calls["sigs"] == 1


def test_check_and_stop_stop_event_without_reporter(tmp_path: Path):
    """Cover stop-event branch in check_and_stop when reporter is None."""
    ex = _executor(tmp_path)
    ex.stop_event.set()
    assert ex.check_and_stop(None) is False


def test_log_node_skip_suppressed_when_log_run_only_true(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover _log_node branch where skipped nodes are hidden under log_run_only mode."""
    ex = _executor(tmp_path)
    node = DummyNode("s")
    logs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": logs.append(message)
    )
    ex._log_node(node, None, True, "cached", log_run_only=True, verbose=False)
    assert logs == []


def test_run_node_console_paths_without_reporter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover _run_node exception branches that print directly via Console when reporter is absent."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)
    printed: list[str] = []
    monkeypatch.setattr(
        executor_module.Console,
        "print",
        lambda self, panel=None: printed.append(str(panel)),
    )

    _, _, failures = ex._run_node(
        DummyNode("pe", call_error=DummyPipelineError("x")),
        {},
        None,
        True,
        [],
    )
    assert failures and printed

    printed.clear()
    _, _, failures2 = ex._run_node(
        DummyNode("ge", call_error=RuntimeError("x")),
        {},
        None,
        True,
        [],
    )
    assert failures2 and printed

    with pytest.raises(RuntimeError):
        ex._run_node(
            DummyNode("ge2", call_error=RuntimeError("x")),
            {},
            None,
            False,
            [],
        )


def test_run_sequential_additional_branches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover sequential branches: stop break, dry-run continue, and upstream fail propagation."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    skip_node = DummyNode("skip", skip_result=(True, "cached"))
    run_node = DummyNode("run", skip_result=(False, "First run"))

    calls = {"n": 0}

    def stop_once(reporter):
        calls["n"] += 1
        return calls["n"] == 1

    monkeypatch.setattr(ex, "check_and_stop", stop_once)
    failures, _ = ex._run_sequential(
        [skip_node, run_node],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        dry_run=False,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []

    monkeypatch.setattr(ex, "check_and_stop", lambda reporter: False)
    run_dry = DummyNode("run_dry", skip_result=(False, "First run"))
    ex._run_sequential(
        [run_dry],
        cache={},
        reporter=DummyReporter(),
        continue_on_error=True,
        failures=[],
        dry_run=True,
        log_run_only=False,
        verbose=False,
    )
    assert run_dry.called == 0

    fail = DummyNode(
        "fail",
        skip_result=(False, "First run"),
        call_error=RuntimeError("boom"),
    )
    dep = DummyNode("dep", skip_result=(False, "First run"), pres=(fail,))
    rep = DummyReporter()
    rep.register([fail, dep])
    failures3, rep3 = ex._run_sequential(
        [fail, dep],
        cache={},
        reporter=rep,
        continue_on_error=True,
        failures=[],
        dry_run=False,
        log_run_only=False,
        verbose=False,
    )
    assert failures3
    assert any(k == "dep" for k, _reason in rep3.upstream_failed)


def test_init_pending_successors_ignores_uncollected_prereq(tmp_path: Path):
    """Cover pending/successor branch where prerequisite key is outside by_key."""
    ex = _executor(tmp_path)
    external = DummyNode("ext")
    local = DummyNode("loc", pres=(external,))
    _, pending, succ = ex._init_pending_and_successors([local])
    assert pending["loc"] == 0
    assert succ["loc"] == []


def test_run_parallel_skip_abort_and_dependency_release(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover parallel skip branch, global abort branch, and dependency-release scheduling."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    sk = DummyNode("sk", skip_result=(True, "cached"))
    rep = DummyReporter()
    rep.register([sk])
    failures, rep2 = ex._run_parallel(
        [sk],
        cache={},
        reporter=rep,
        continue_on_error=True,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []
    assert any(k == "sk" for k, _r in rep2.skipped)

    ex2 = _executor(tmp_path)
    ex2.stop_event.set()
    failures2, _ = ex2._run_parallel(
        [DummyNode("a", skip_result=(False, "First run"))],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures2 == []

    ex3 = _executor(tmp_path)
    a = DummyNode("a", skip_result=(False, "First run"))
    b = DummyNode("b", skip_result=(False, "First run"), pres=(a,))
    rep3 = DummyReporter()
    rep3.register([a, b])
    failures3, rep3 = ex3._run_parallel(
        [a, b],
        cache={},
        reporter=rep3,
        continue_on_error=True,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert failures3 == []
    assert set(rep3.done) == {"a", "b"}


def test_run_parallel_failure_with_continue_on_error_false(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Exercise parallel failure path while forcing abort to avoid waiting on done_event."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)
    a = DummyNode(
        "a",
        skip_result=(False, "First run"),
        call_error=RuntimeError("boom"),
    )
    b = DummyNode("b", pres=(a,), skip_result=(False, "First run"))
    rep = DummyReporter()
    rep.register([a, b])
    failures, rep = ex._run_parallel(
        [a, b],
        cache={},
        reporter=rep,
        continue_on_error=False,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert failures
    assert rep.failed


def test_check_duplicate_outputs_same_owner_no_conflict(tmp_path: Path):
    """Cover duplicate-path branch where same owner key is repeated and should not conflict."""
    ex = _executor(tmp_path)
    p = tmp_path / "same.txt"
    node = DummyNode("same", output_files=(p, p))
    ex.check_duplicate_outputs([node])


def test_collect_same_object_duplicate_identity_branch(tmp_path: Path):
    """Cover collect branch where same object appears twice and is ignored by identity check."""
    ex = _executor(tmp_path)
    shared = DummyNode("s")
    root = DummyNode("r", pres=(shared, shared))
    nodes = ex.collect([root])
    assert {n.key for n in nodes} == {"r", "s"}


def test_toposort_external_pre_and_multidegree(tmp_path: Path):
    """Cover toposort external-pre skip and indegree decrement-not-zero path."""
    ex = _executor(tmp_path)
    ext = DummyNode("ext")
    a = DummyNode("a")
    b = DummyNode("b")
    c = DummyNode("c", pres=(ext, a, b))
    order = ex.toposort([a, b, c])
    assert order[-1].key == "c"


def test_to_dot_skips_uncollected_pre_edges(tmp_path: Path):
    """Cover to_dot branch where prerequisites absent from by_key do not emit edges."""
    ex = _executor(tmp_path)
    ext = DummyNode("ext")
    local = DummyNode("loc", pres=(ext,))
    # Force collect() to return only the local node so the pre-key check can fail.
    ex.collect = lambda targets: [local]  # type: ignore[method-assign]
    dot = ex.to_dot([local])
    assert '"ext" -> "loc"' not in dot


def test_run_targets_dry_run_skips_cache_save_and_pending_false_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover run_targets dry-run branch and pending-check false path in finalization."""
    ex = _executor(tmp_path)
    node = DummyNode("n")
    monkeypatch.setattr(ex, "collect", lambda targets: [node])
    monkeypatch.setattr(ex, "toposort", lambda nodes: nodes)
    monkeypatch.setattr(ex, "check_duplicate_outputs", lambda nodes: None)
    monkeypatch.setattr(ex, "write_dot", lambda nodes, path: path)
    monkeypatch.setattr(ex.registry, "write_registry", lambda path: None)

    reporter = DummyReporter()
    reporter.register([node])
    # Mark as already done so "pending" condition is false in finally-loop check.
    reporter._nodes[node.key].state = NodeState.DONE
    monkeypatch.setattr(executor_module, "ProgressReporter", lambda logger: reporter)
    monkeypatch.setattr(ex, "_run_sequential", lambda *args, **kwargs: ([], reporter))

    saved = {"count": 0}
    monkeypatch.setattr(
        ex,
        "save_cache",
        lambda path, cache: saved.__setitem__("count", saved["count"] + 1),
    )
    ex.run_targets([node], progress=True, dry_run=True, parallel=False)
    assert saved["count"] == 0


def test_log_sigs_handles_none_nodes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Cover log_sigs branch where an order entry is None and should be ignored."""
    ex = _executor(tmp_path)
    logs: list[str] = []
    monkeypatch.setattr(
        ex, "log", lambda reporter, message, level="INFO": logs.append(message)
    )
    monkeypatch.setattr(ex.registry, "write_registry", lambda path: None)
    ex.log_sigs([None])
    assert logs == []


def test_run_sequential_skip_with_reporter_and_upstream_without_reporter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover reporter mark_skip and upstream_failed branch when reporter is None."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    sk = DummyNode("sk", skip_result=(True, "cached"))
    rep = DummyReporter()
    rep.register([sk])
    ex._run_sequential(
        [sk],
        cache={},
        reporter=rep,
        continue_on_error=True,
        failures=[],
        dry_run=False,
        log_run_only=False,
        verbose=False,
    )
    assert rep.skipped and rep.skipped[0][0] == "sk"

    fail = DummyNode(
        "fail", skip_result=(False, "First run"), call_error=RuntimeError("boom")
    )
    dep = DummyNode("dep", skip_result=(False, "First run"), pres=(fail,))
    failures, _ = ex._run_sequential(
        [fail, dep],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        dry_run=False,
        log_run_only=False,
        verbose=False,
    )
    assert failures


def test_run_parallel_additional_failure_and_abort_branches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover additional _run_parallel branches for abort checks and PipelineError passthrough."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    # Skip node with no reporter covers skip path branch that bypasses mark_skip.
    sk = DummyNode("sk", skip_result=(True, "cached"))
    failures, _ = ex._run_parallel(
        [sk],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []

    # Single-node failure with continue_on_error=False must set abort and return.
    pe = DummyNode(
        "pe", skip_result=(False, "First run"), call_error=DummyPipelineError("x")
    )
    failures2, _ = ex._run_parallel(
        [pe],
        cache={},
        reporter=None,
        continue_on_error=False,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures2

    # Abort during node execution should trigger abort check in the unblock-successors section.
    ex2 = _executor(tmp_path)

    class AbortNode(DummyNode):
        def __call__(self):
            self.called += 1
            ex2.stop_event.set()
            return True

    ab = AbortNode("ab", skip_result=(False, "First run"))
    failures3, _ = ex2._run_parallel(
        [ab],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures3 == []


def test_run_parallel_pending_deps_not_zero_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover branch where successor dependency counter decrements but is still non-zero."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)
    a = DummyNode("a", skip_result=(False, "First run"))
    b = DummyNode("b", skip_result=(False, "First run"))
    c = DummyNode("c", skip_result=(False, "First run"), pres=(a, b))
    rep = DummyReporter()
    rep.register([a, b, c])
    failures, rep = ex._run_parallel(
        [a, b, c],
        cache={},
        reporter=rep,
        continue_on_error=True,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []


def test_run_parallel_submit_ready_internal_abort_and_early_worker_return(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover submit loop abort re-check and worker early-return when abort is already set."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    a = DummyNode("a", skip_result=(False, "First run"))
    b = DummyNode("b", skip_result=(False, "First run"))

    monkeypatch.setattr(
        ex,
        "_init_pending_and_successors",
        lambda order: (
            {"a": a, "b": b},
            {"a": 0, "b": 0},
            {"a": [], "b": []},
        ),
    )

    class FakePool:
        def __init__(self, *args, **kwargs):
            self.calls = 0

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, key):
            self.calls += 1
            if self.calls == 1:
                # Ensure the worker sees abort already set (line 756 path),
                # then the submit loop's next-iteration abort re-check fires.
                ex.stop_event.set()
                fn(key)
            return object()

        def shutdown(self, wait=False, cancel_futures=True):
            return None

    monkeypatch.setattr(executor_module, "ThreadPoolExecutor", FakePool)

    failures, _ = ex._run_parallel(
        [a, b],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []


def test_run_parallel_blocked_successor_already_seen_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover recursive block propagation when a shared successor is already blocked/completed."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    class LenZeroList(list):
        def __len__(self):
            # Forces the failure-path done_event guard to trigger, avoiding hangs
            # while still executing the blocking recursion in worker threads.
            return 0

    child = DummyNode("c", skip_result=(False, "First run"))
    left = DummyNode(
        "a", skip_result=(False, "First run"), call_error=RuntimeError("boom")
    )
    right = DummyNode(
        "b", skip_result=(False, "First run"), call_error=RuntimeError("boom")
    )
    child.pres = (left, right)
    order = LenZeroList([left, right, child])

    failures, _ = ex._run_parallel(
        order,
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=2,
        log_run_only=False,
        verbose=False,
    )
    assert len(failures) >= 1


def test_run_parallel_duplicate_ready_key_hits_active_completed_guard(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover submit_ready branch where a key is encountered again while active/completed."""
    ex = _executor(tmp_path)
    monkeypatch.setattr(executor_module, "PipelineError", DummyPipelineError)

    node = DummyNode("dup", skip_result=(True, "cached"))

    class DuplicatePending(dict):
        def items(self):
            # Return the same key twice so submit_ready re-checks active/completed.
            return [("dup", 0), ("dup", 0)]

    pending = DuplicatePending({"dup": 0})
    monkeypatch.setattr(
        ex,
        "_init_pending_and_successors",
        lambda order: ({"dup": node}, pending, {"dup": []}),
    )

    failures, _ = ex._run_parallel(
        [node],
        cache={},
        reporter=None,
        continue_on_error=True,
        failures=[],
        max_workers=1,
        log_run_only=False,
        verbose=False,
    )
    assert failures == []


def test_capture_store_state_mismatch_outputs_not_ok_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover mismatch branch where outputs_ok is false and therefore no accept/raise action occurs."""
    ex = _executor(tmp_path)
    mm = DummyNode(
        "mm",
        skip_result=(False, "Cached signature does not match"),
        outputs_ok_result=False,
    )
    monkeypatch.setattr(ex, "ordered_nodes_for_run", lambda: [mm])
    monkeypatch.setattr(
        executor_module, "explain_signature_diff", lambda node, store: (False, "x")
    )
    outstore = tmp_path / "out_mm.json"
    ex.save_cache(ex.signature_store_path, {})
    ex.save_sigstore(ex.signature_data_path, {})
    ex.capture_store_state(outstore, accept_existing=False)
    # No cache write for mm happens in this path; output store remains absent or empty.
    assert (
        not outstore.exists()
        or json.loads(outstore.read_text(encoding="utf-8") or "{}") == {}
    )


def test_capture_store_state_additional_reason_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Cover capture_store_state branches for outputs_ok false and mismatch-raise behavior."""
    ex = _executor(tmp_path)
    n_first_no_outputs = DummyNode(
        "fno",
        skip_result=(False, "First run"),
        outputs_ok_result=False,
    )
    n_mismatch_raise = DummyNode(
        "mm",
        skip_result=(False, "Cached signature does not match"),
        outputs_ok_result=True,
    )
    n_other = DummyNode(
        "other", skip_result=(False, "Other reason"), outputs_ok_result=True
    )
    monkeypatch.setattr(
        ex,
        "ordered_nodes_for_run",
        lambda: [n_first_no_outputs, n_mismatch_raise, n_other],
    )
    monkeypatch.setattr(
        executor_module,
        "explain_signature_diff",
        lambda node, store: (False, "x"),
    )
    outstore = tmp_path / "out2.json"
    ex.save_cache(ex.signature_store_path, {})
    ex.save_sigstore(ex.signature_data_path, {})

    with pytest.raises(RuntimeError, match="cache signature differs"):
        ex.capture_store_state(outstore, accept_existing=False)

    ex.capture_store_state(outstore, accept_existing=True)
    data = json.loads(outstore.read_text(encoding="utf-8"))
    assert data["mm"] == n_mismatch_raise.signature
