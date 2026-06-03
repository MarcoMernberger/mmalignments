from __future__ import annotations

import sys
import threading
import time
from collections import deque
from dataclasses import dataclass, field
from enum import Enum
from logging import Logger
from typing import IO

from rich.console import Console  # type: ignore[import]
from rich.layout import Layout  # type: ignore[import]
from rich.live import Live  # type: ignore[import]
from rich.panel import Panel  # type: ignore[import]
from rich.style import Style  # type: ignore[import]
from rich.table import Table  # type: ignore[import]
from rich.text import Text  # type: ignore[import]

from mmalignments.services.logging import enable_console


class NodeState(Enum):
    RUNNING = "running"
    SKIPPED = "skipped"
    PENDING = "pending"
    FAILED = "failed"
    UPSTREAM_FAILED = "upstream_failed"
    DONE = "done"


STATE_PRIORITY = {
    NodeState.RUNNING: 0,
    NodeState.FAILED: 1,
    NodeState.UPSTREAM_FAILED: 2,
    NodeState.DONE: 3,
    NodeState.PENDING: 4,
    NodeState.SKIPPED: 5,
}


_STYLE = {
    NodeState.PENDING: Style(color="cyan"),
    NodeState.SKIPPED: Style(color="bright_black"),
    NodeState.RUNNING: Style(color="yellow", bold=True),
    NodeState.DONE: Style(color="green"),
    NodeState.FAILED: Style(color="red"),
    NodeState.UPSTREAM_FAILED: Style(color="dark_orange"),
}


_ICON = {
    NodeState.PENDING: "○",
    NodeState.SKIPPED: "⊘",
    NodeState.RUNNING: "◎",
    NodeState.DONE: "✓",
    NodeState.FAILED: "✗",
    NodeState.UPSTREAM_FAILED: "↯",
}


@dataclass
class NodeProgress:
    name: str
    key: str
    threads: int
    state: NodeState = NodeState.PENDING
    reason: str = ""
    elapsed: float = 0.0
    _start: float = field(default=0.0, repr=False, compare=False)

    def start(self) -> None:
        self._start = time.monotonic()
        self.state = NodeState.RUNNING

    def finish(self, failed: bool = False) -> None:
        self.elapsed = time.monotonic() - self._start
        self.state = NodeState.FAILED if failed else NodeState.DONE

    def skip(self, reason: str = "") -> None:
        self.state = NodeState.SKIPPED
        self.reason = reason

    def upstream_failed(self, reason: str = "upstream failed") -> None:
        self.state = NodeState.UPSTREAM_FAILED
        self.reason = reason

    def fmt_elapsed(self) -> str:
        if self.state == NodeState.RUNNING:
            secs = time.monotonic() - self._start
        elif self.elapsed > 0:
            secs = self.elapsed
        else:
            return ""
        m, s = divmod(int(secs), 60)
        h, m = divmod(m, 60)
        if h:
            return f"[{h}h {m:02d}m {s:02d}s]"
        if m:
            return f"[{m:02d}m {s:02d}s]"
        return f"[{s}s]"


_SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


def _spinner_char() -> str:
    """Return current spinner frame based on wall-clock time. Thread-safe."""
    return _SPINNER[int(time.monotonic() * 10) % len(_SPINNER)]


class _LiveLayout:
    """Thin wrapper so rich re-calls _build_layout on every auto-refresh."""

    def __init__(self, reporter: ProgressReporter) -> None:
        self._reporter = reporter
        self.final = False

    def __rich_console__(self, console, options):
        yield self._reporter._build_layout(self.final)


class ProgressReporter:
    """Split-pane live progress display using rich.

    Left pane  → node list with status icons + elapsed time
    Right pane → intercepted log messages (push_log)

    Parameters
    ----------
    stream :
        Target stream (default ``sys.stderr``).
    refresh_interval :
        Live refresh rate in seconds.
    max_log_lines :
        How many log lines to keep visible on the right.
    """

    def __init__(
        self,
        logger: Logger,
        stream: IO | None = None,
        refresh_per_second: int = 10,
        max_lines: int = 55,
    ) -> None:
        self._console = Console(file=stream or sys.stderr, highlight=False)
        self._refresh_per_second = refresh_per_second
        self.max_visible = max_lines
        self._logger = logger
        self._nodes: dict[str, NodeProgress] = {}
        self._order: list[str] = []
        self._log_buffer: deque[str] = deque(maxlen=self.max_visible)
        self._final_errors: list[Panel] = []
        self._lock = threading.Lock()
        self._live: Live | None = None
        self._renderable: _LiveLayout | None = None
        self.log_width = self._console.size.width# // 2 - 2

    ###########################################################################
    # Node Status handling
    ###########################################################################

    def error_panels(self) -> list[Panel]:
        with self._lock:
            return list(self._final_errors)

    def register(self, nodes: list) -> None:
        """Register nodes. Call before start_live()."""
        with self._lock:
            for node in nodes:
                if node.key in self._nodes:
                    continue
                threads = getattr(node, "threads", 1)
                self._nodes[node.key] = NodeProgress(
                    name=node.name, key=node.key, threads=threads
                )
                self._order.append(node.key)

    def push_log(self, message: str) -> None:
        """Add a log line to the right-hand pane."""
        with self._lock:
            self._log_buffer.append(message.rstrip())

    def push_error_panel(self, renderable: Panel):
        with self._lock:
            self._final_errors.append(renderable)

    def mark_skip(self, key: str, reason: str = "") -> None:
        with self._lock:
            if key in self._nodes:
                self._nodes[key].skip(reason)

    def mark_start(self, key: str) -> None:
        with self._lock:
            if key in self._nodes:
                self._nodes[key].start()

    def mark_done(self, key: str) -> None:
        with self._lock:
            if key in self._nodes:
                self._nodes[key].finish(failed=False)

    def mark_failed(self, key: str) -> None:
        with self._lock:
            if key in self._nodes:
                self._nodes[key].finish(failed=True)

    def mark_upstream_failed(self, key: str, reason: str = "upstream failed") -> None:
        with self._lock:
            if key in self._nodes:
                self._nodes[key].upstream_failed(reason)

    def start_live(self) -> None:
        self._renderable = _LiveLayout(self)
        self.layout = Layout()
        enable_console(
            self._logger, False
        )  # disable console echo when using live display
        self._live = Live(
            self._renderable,
            console=self._console,
            refresh_per_second=self._refresh_per_second,
            transient=False,
        )
        if self._live:
            self._live.start()

    def stop_live(self) -> None:
        enable_console(
            self._logger, True
        )  # re-enable console echo when stopping live display
        if self._live:
            if self._renderable:
                self._renderable.final = True
            self._live.refresh()  # one final render with final=True
            self._live.stop()
            self._live = None

    def summary(self) -> Text:
        from collections import Counter

        with self._lock:
            counts: Counter[NodeState] = Counter(
                np.state for np in self._nodes.values()
            )
            total = len(self._nodes)
        t = Text(f"Total {total}")
        for state, label in [
            (NodeState.DONE, "done"),
            (NodeState.SKIPPED, "skipped"),
            (NodeState.FAILED, "failed"),
            (NodeState.UPSTREAM_FAILED, "upstream failed"),
            (NodeState.RUNNING, "running"),
            (NodeState.PENDING, "pending"),
        ]:
            n = counts[state]
            t.append(" │ ")
            t.append(f"{_ICON[state]} {n} {label}", style=_STYLE[state])
        return t

    ###########################################################################
    # Panel Construction
    ###########################################################################

    def _build_node_panel(self, final: bool = False) -> Panel:
        """Build the left-hand node table."""
        spin = _spinner_char()  # time-based; no background thread needed

        # Snapshot all display values inside the lock to avoid races.
        with self._lock:
            nodes = [
                (k, np.state, np.name, np.threads, np.fmt_elapsed())
                for k in self._order
                if (np := self._nodes.get(k)) is not None
            ]

        running = [n for n in nodes if n[1] == NodeState.RUNNING]
        others = [n for n in nodes if n[1] != NodeState.RUNNING]
        remaining_slots = self.max_visible - len(running)
        if remaining_slots > 0:
            others.sort(key=lambda x: STATE_PRIORITY[x[1]])
            snapshot = running + others[:remaining_slots]
        else:
            snapshot = running[: self.max_visible]
        table = Table.grid(padding=(0, 1))
        table.add_column(width=2)  # icon
        table.add_column(min_width=self.log_width - 19)  # name
        table.add_column(width=5)  # threads
        table.add_column(width=12)  # elapsed

        running_count = 0
        for _, state, name, threads, elapsed in snapshot:
            if state == NodeState.RUNNING:
                running_count += 1
                icon_char = spin if not final else _ICON[state]
            else:
                icon_char = _ICON[state]
            table.add_row(
                Text(icon_char, style=_STYLE[state]),
                Text(name[: self.log_width - 19], style=_STYLE[state]),
                Text(
                    f"({threads}t)" if threads and threads > 1 else "",
                    style=Style(color="bright_black"),
                ),
                Text(elapsed, style=Style(color="bright_black")),
            )

        title = (
            Text(f"Pipeline running {spin}", style="bold yellow")
            if running_count and not final
            else Text("Pipeline", style="bold")
        )
        return Panel(
            table,
            title=title,
            subtitle=self.summary(),
            border_style="bright_black",
        )

    def _build_log_panel(self) -> Panel:
        """Build the right-hand log panel."""
        with self._lock:
            lines = list(self._log_buffer)

        log_text = Text()
        for i, line in enumerate(lines):
            if i:
                log_text.append("\n")
            log_text.append(
                line[: self.log_width], style=Style(color="bright_white", dim=True)
            )

        return Panel(
            log_text,
            title=Text("Logs", style="bold"),
            border_style="bright_black",
        )

    def _build_layout(self, final: bool = False) -> Layout:
        self.layout.split_row(
            Layout(self._build_node_panel(final=final), name="nodes", ratio=1),
            Layout(self._build_log_panel(), name="logs", ratio=1),
        )
        return self.layout
