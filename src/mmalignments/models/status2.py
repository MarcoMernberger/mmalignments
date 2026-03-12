# from __future__ import annotations

# import sys
# import threading
# import time
# from collections import deque
# from dataclasses import dataclass, field
# from enum import Enum
# from typing import IO


# class NodeState(Enum):
#     PENDING = "pending"
#     SKIPPED = "skipped"
#     RUNNING = "running"
#     DONE = "done"
#     FAILED = "failed"


# # ANSI colour for the whole line
# _LINE_COLOUR = {
#     NodeState.PENDING: "\033[90m",  # grey
#     NodeState.SKIPPED: "\033[36m",  # cyan
#     NodeState.RUNNING: "\033[33m",  # yellow / bold handled separately
#     NodeState.DONE: "\033[32m",  # green
#     NodeState.FAILED: "\033[31m",  # red
# }
# _RESET = "\033[0m"
# _BOLD = "\033[1m"
# _DIM = "\033[2m"

# _ICON = {
#     NodeState.PENDING: "○",
#     NodeState.SKIPPED: "⊘",
#     NodeState.RUNNING: "◎",
#     NodeState.DONE: "✓",
#     NodeState.FAILED: "✗",
# }

# # Summary colour per state
# _SUMMARY_COLOUR = {
#     NodeState.DONE: "\033[32m",
#     NodeState.SKIPPED: "\033[36m",
#     NodeState.FAILED: "\033[31m",
#     NodeState.RUNNING: "\033[33m",
#     NodeState.PENDING: "\033[90m",
# }


# def _use_colour(stream: IO) -> bool:
#     return hasattr(stream, "isatty") and stream.isatty()


# @dataclass
# class NodeProgress:
#     name: str
#     key: str
#     threads: int
#     state: NodeState = NodeState.PENDING
#     reason: str = ""
#     elapsed: float = 0.0
#     _start: float = field(default=0.0, repr=False, compare=False)

#     def start(self) -> None:
#         self._start = time.monotonic()
#         self.state = NodeState.RUNNING

#     def finish(self, failed: bool = False) -> None:
#         self.elapsed = time.monotonic() - self._start
#         self.state = NodeState.FAILED if failed else NodeState.DONE

#     def skip(self, reason: str = "") -> None:
#         self.state = NodeState.SKIPPED
#         self.reason = reason

#     def fmt_elapsed(self) -> str:
#         """Format elapsed time, shown for RUNNING/DONE/FAILED."""
#         if self.state == NodeState.RUNNING:
#             # live elapsed
#             secs = time.monotonic() - self._start
#         elif self.elapsed > 0:
#             secs = self.elapsed
#         else:
#             return ""
#         m, s = divmod(int(secs), 60)
#         h, m = divmod(m, 60)
#         if h:
#             return f" [{h}h {m:02d}m {s:02d}s]"
#         if m:
#             return f" [{m:02d}m {s:02d}s]"
#         return f" [{s}s]"

#     def render(self, width: int = 52, colour: bool = True) -> str:
#         c = _LINE_COLOUR.get(self.state, "") if colour else ""
#         icon = _ICON[self.state]
#         label = self.name[:width].ljust(width)
#         t = f" ({self.threads}t)" if (self.threads and self.threads > 1) else "     "
#         elapsed = self.fmt_elapsed()
#         reset = _RESET if colour else ""

#         if self.state == NodeState.RUNNING and colour:
#             # bold + colour for running
#             return f"  {_BOLD}{c}{icon} {label}{t}{elapsed}{reset}"
#         return f"  {c}{icon} {label}{t}{elapsed}{reset}"


# class ProgressReporter:
#     """Live progress display that intercepts logging so nothing else writes
#     to the terminal while the spinner is active.

#     Parameters
#     ----------
#     stream :
#         Target stream (default ``sys.stderr``).
#     refresh_interval :
#         Spinner refresh rate in seconds.
#     plain :
#         Force plain (no ANSI) output.
#     max_log_lines :
#         How many intercepted log lines to show above the node list.
#     """

#     _SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"

#     def __init__(
#         self,
#         stream: IO | None = None,
#         refresh_interval: float = 0.3,
#         plain: bool = False,
#         max_log_lines: int = 6,
#     ) -> None:
#         self._stream = stream or sys.stderr
#         self._interval = refresh_interval
#         self._colour = _use_colour(self._stream) and not plain
#         self._plain = plain or not _use_colour(self._stream)

#         self._nodes: dict[str, NodeProgress] = {}
#         self._order: list[str] = []
#         self._log_buffer: deque[str] = deque(maxlen=max_log_lines)

#         self._lock = threading.Lock()
#         self._spinner = 0
#         self._stop = threading.Event()
#         self._thread: threading.Thread | None = None
#         self._printed = 0  # ANSI lines currently on screen
#         self._last_render: float = 0.0

#     # ------------------------------------------------------------------
#     # Public API
#     # ------------------------------------------------------------------

#     def register(self, nodes: list) -> None:
#         """Register nodes. Call before start_live()."""
#         with self._lock:
#             for node in nodes:
#                 if node.key in self._nodes:
#                     continue
#                 threads = getattr(node, "threads", 1)
#                 np = NodeProgress(name=node.name, key=node.key, threads=threads)
#                 self._nodes[node.key] = np
#                 self._order.append(node.key)

#     def push_log(self, message: str) -> None:
#         """Intercept a log/print message – shown inside the progress display."""
#         with self._lock:
#             # strip any trailing newline for clean rendering
#             self._log_buffer.append(message.rstrip())

#     def mark_skip(self, key: str, reason: str = "") -> None:
#         with self._lock:
#             if key in self._nodes:
#                 self._nodes[key].skip(reason)
#             self._refresh_locked()

#     def mark_start(self, key: str) -> None:
#         with self._lock:
#             if key in self._nodes:
#                 self._nodes[key].start()
#             self._refresh_locked()

#     def mark_done(self, key: str) -> None:
#         with self._lock:
#             if key in self._nodes:
#                 self._nodes[key].finish(failed=False)
#             self._refresh_locked()

#     def mark_failed(self, key: str) -> None:
#         with self._lock:
#             if key in self._nodes:
#                 self._nodes[key].finish(failed=True)
#             self._refresh_locked()

#     def start_live(self) -> None:
#         if self._plain:
#             return
#         self._stop.clear()
#         self._thread = threading.Thread(target=self._loop, daemon=True)
#         self._thread.start()

#     def stop_live(self) -> None:
#         self._stop.set()
#         if self._thread:
#             self._thread.join(timeout=2)
#         with self._lock:
#             self._erase()
#             self._render_full(final=True)

#     def summary(self, colour: bool | None = None) -> str:
#         c = self._colour if colour is None else colour
#         counts = dict.fromkeys(NodeState, 0)
#         for np in self._nodes.values():
#             counts[np.state] += 1
#         total = len(self._nodes)

#         parts = [f"Total {total}"]
#         for state, label in [
#             (NodeState.DONE, "done"),
#             (NodeState.SKIPPED, "skipped"),
#             (NodeState.FAILED, "failed"),
#             (NodeState.RUNNING, "running"),
#             (NodeState.PENDING, "pending"),
#         ]:
#             n = counts[state]
#             col = _SUMMARY_COLOUR[state] if c else ""
#             rst = _RESET if c else ""
#             ico = _ICON[state]
#             parts.append(f"{col}{ico} {n} {label}{rst}")

#         return " | ".join(parts)

#     # ------------------------------------------------------------------
#     # Internal
#     # ------------------------------------------------------------------

#     def _loop(self) -> None:
#         while not self._stop.is_set():
#             self._spinner = (self._spinner + 1) % len(self._SPINNER)
#             with self._lock:
#                 # Skip render if a mark_* already refreshed very recently
#                 if time.monotonic() - self._last_render >= self._interval * 0.5:
#                     self._erase()
#                     self._render_full()
#             time.sleep(self._interval)

#     def _refresh_locked(self) -> None:
#         """Immediate refresh. Caller must hold self._lock."""
#         if self._plain:
#             return
#         self._erase()
#         self._render_full()

#     def _refresh(self) -> None:
#         """Immediate refresh (called from mark_* methods)."""
#         if self._plain:
#             return
#         with self._lock:
#             self._refresh_locked()

#     def _erase(self) -> None:
#         """Erase all lines written by the last render. Caller must hold lock."""
#         if self._printed == 0:
#             return
#         self._stream.write(f"\033[{self._printed}A\033[J")
#         self._stream.flush()
#         self._printed = 0

#     def _render_full(self, final: bool = False) -> None:
#         """Render the complete progress display. Caller must hold lock."""
#         c = self._colour
#         lines: list[str] = []

#         # ── header ──────────────────────────────────────────────────────────
#         running_count = sum(
#             1 for np in self._nodes.values() if np.state == NodeState.RUNNING
#         )
#         if running_count and not final:
#             spin = self._SPINNER[self._spinner]
#             header = (
#                 f"  {_BOLD}Pipeline running {spin}{_RESET}"
#                 if c
#                 else "  Pipeline running"
#             )
#         else:
#             header = f"  {_BOLD}Pipeline{_RESET}" if c else "  Pipeline"
#         lines.append(header)
#         lines.append("")

#         # ── intercepted log lines ────────────────────────────────────────────
#         if self._log_buffer:
#             dim = _DIM if c else ""
#             rst = _RESET if c else ""
#             for msg in self._log_buffer:
#                 lines.append(f"  {dim}{msg}{rst}")
#             lines.append("")

#         # ── node list ────────────────────────────────────────────────────────
#         for key in self._order:
#             np = self._nodes[key]
#             line = np.render(colour=c)
#             # replace ◎ with live spinner char for running nodes
#             if np.state == NodeState.RUNNING and not final:
#                 spin = self._SPINNER[self._spinner]
#                 line = line.replace(_ICON[NodeState.RUNNING], spin, 1)
#             lines.append(line)

#         lines.append("")
#         lines.append("  " + self.summary())
#         lines.append("")

#         text = "\n".join(lines) + "\n"
#         self._stream.write(text)
#         self._stream.flush()
#         self._printed = len(lines) + 1
#         self._last_render = time.monotonic()
