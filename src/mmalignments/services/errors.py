from __future__ import annotations

import shlex
import subprocess
import sys
import traceback
from datetime import datetime
from pathlib import Path
from subprocess import CalledProcessError
from typing import TYPE_CHECKING, Iterable, Optional

from rich.console import Group  # type: ignore[import]
from rich.panel import Panel  # type: ignore[import]
from rich.syntax import Syntax  # type: ignore[import]
from rich.text import Text  # type: ignore[import]
from rich.traceback import Traceback  # type: ignore[import]

from mmalignments.utils import constants  # type: ignore[import]

if TYPE_CHECKING:
    from mmalignments.models.elements import Element  # type: ignore[import]


def _isatty(stream) -> bool:
    try:
        return stream.isatty()
    except Exception:
        return False


def _ansi(s: str, code: str, enabled: bool) -> str:
    if not enabled:
        return s
    return f"\x1b[{code}m{s}\x1b[0m"


def _truncate(s: str | None, max_chars: int) -> str:
    s = s or ""
    if max_chars <= 0 or len(s) <= max_chars:
        return s
    head = s[:max_chars]
    return head + f"\n… <truncated, {len(s) - max_chars} chars omitted>"


def _build_cmd_list(e: Exception) -> list[str]:
    c = getattr(e, "cmd", None)
    if isinstance(c, (list, tuple)):
        cmd_list = list(c)
    elif (
        isinstance(e.args, (list, tuple))
        and e.args
        and isinstance(e.args[0], (list, tuple))
    ):
        cmd_list = list(e.args[0])  # rarely
    else:
        cmd_list = [str(getattr(e, "cmd", e.args))]
    return cmd_list


def _resolve_err_out(
    stderr: str | None, stdout: str | None, e: Exception
) -> tuple[str, str]:
    err = ""
    out = ""
    if stderr is not None:
        err = stderr
    elif isinstance(e, CalledProcessError):
        err = e.stderr or ""
    if stdout is not None:
        out = stdout
    elif isinstance(e, CalledProcessError):
        out = e.stdout or ""
    return err, out


def _resolve_returncode(e: Exception) -> int | None:
    return getattr(e, "returncode", None)


def format_called_process_error(
    *,
    e: Exception,
    phase: str = "RUN",
    cmd: Iterable[str] | None = None,
    stdout: str | None = None,
    stderr: str | None = None,
    show_stack: bool = False,
    trace: str | None = None,
    color: bool | None = None,
    width: int = 120,
    max_stream_chars: int = 1000,
) -> str:
    """
    Returns (print_text, log_text).
    - print_text: pretty, optionally colored, possibly truncated stdout/stderr
    - log_text: no color, full stdout/stderr (or provided), not truncated unless
      you choose to truncate it by passing max_stream_chars, which defaults to
      1000 chars per stream but can be set to 0 or negative for no truncation.
    """
    # Determine command
    cmd_list = list(cmd) if cmd is not None else None
    if cmd_list is None:
        # CalledProcessError typically has .cmd; fallback to .args
        if isinstance(e, CalledProcessError):
            cmd_list = _build_cmd_list(e)
        else:
            cmd_list = [str(getattr(e, "cmd", e.args))]

    cmd_str = shlex.join(cmd_list)

    # Streams: prefer explicit, else use e.stderr/e.stdout/e.output
    err, out = _resolve_err_out(stderr, stdout, e)

    # Decide color
    if color is None:
        color = _isatty(sys.stderr) or _isatty(sys.stdout)
    # textwrap.fill(cmd_str, width=width, break_long_words=False, break_on_hyphens=False)  # noqa: E501
    cmd_wrapped = cmd_str

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = f"[{ts}] External command failed\nPhase: {phase}"
    rc = _resolve_returncode(e)
    rc_line = f"Return code: {rc}"

    # Colored pieces for terminal
    header_c = _ansi(header, "1;31", color)  # bold red
    cmd_label_c = _ansi("Command:", "1;36", color)  # bold cyan
    cmd_c = _ansi(cmd_wrapped, "1;33", color)  # bold yellow
    rc_c = _ansi(rc_line, "1;31", color)

    sep = "=" * min(width, 120)

    print_parts = [
        sep,
        header_c,
        rc_c,
        f"{cmd_label_c} {cmd_c}",
    ]

    if show_stack:
        print_parts.append(_ansi("Stack (caller):", "1;35", color))
        if trace:
            print_parts.append(trace.strip())

    # Truncate for terminal readability
    out_short = _truncate(out, max_stream_chars)
    err_short = _truncate(err, max_stream_chars)

    if err_short.strip():
        print_parts.append(_ansi("STDERR:", "1;31", color))
        print_parts.append(err_short.rstrip())
    if out_short.strip():
        print_parts.append(_ansi("STDOUT:", "1;34", color))
        print_parts.append(out_short.rstrip())

    print_parts.append(sep)
    print_text = "\n".join(print_parts) + "\n"

    # Log text (no ANSI), usually keep full streams
    log_parts = [
        sep,
        header,
        rc_line,
        f"Command:\n{cmd_str}",
    ]
    if err.strip():
        log_parts.append("STDERR:")
        log_parts.append(err.rstrip())
    if out.strip():
        log_parts.append("STDOUT:")
        log_parts.append(out.rstrip())
    log_parts.append(sep)

    return print_text


def handle_called_process_error(
    *,
    e: subprocess.CalledProcessError,
    cmd: Optional[Iterable[str]] = None,
    stdout: Optional[str] = None,
    stderr: Optional[str] = None,
    creation_trace: str | None = None,
    logfile: Path | None = None,
    trace: str | None = None,
    show_stack: bool = True,
    color: bool | None = None,
) -> tuple[str, str]:
    """
    Build a rich error message (with optional colors + stacktrace) and return
    it and a non-colored log text.
    Does NOT exit or print.
    """

    print_text, log_text = format_called_process_error(
        e=e,
        cmd=cmd,
        stdout=stdout,
        stderr=stderr,
        show_stack=False,  # we handle stack manually
        trace=trace,
        color=color,
    )
    # -------- STACKTRACE BLOCK --------
    stack_block = ""
    if show_stack:
        stack_sep = "===========STACKTRACE==========="
        stack_content = trace or traceback.format_exc()

        if color is None:
            color = _isatty(sys.stderr) or _isatty(sys.stdout)

        stack_sep_c = _ansi(stack_sep, "1;35", color)  # magenta
        stack_content_c = _ansi(stack_content.strip(), "0;37", color)

        stack_block = (
            f"{stack_sep_c}\n" f"{stack_content_c}\n" f"{stack_sep_c}\n\n"
        )  # noqa: E501
    creation_block = ""
    if creation_trace:
        creation_block = (
            f"===========ELEMENT CREATION===========\n"
            f"{creation_trace}\n"
            f"===========ELEMENT CREATION===========\n\n"
        )
    # -------- OPTIONAL LOGFILE HINT --------
    logfile_block = ""
    if logfile:
        logfile_block = f"\nFull log file: {logfile}\n"

    # -------- FINAL TEXT --------
    final_text = stack_block + creation_block + print_text + logfile_block

    return final_text, log_text


class ErrorInfo:
    def __init__(
        self,
        *,
        e: Exception,
        cmd: Iterable[str] | None = None,
        stdout: str | None = None,
        stderr: str | None = None,
        trace: str | None = None,
        logfile: Path | None = None,
        creation_trace: str | None = None,
        phase: str = "RUN",
        color: bool = True,
        show_stack: bool = True,
    ):
        """
        Initialize an ErrorInfo object. This object encapsulates all relevant
        information about an error that occurred during pipeline execution, including
        the original exception, command, outputs, stack trace, and more. It also
        provides methods to format this information for logging and display purposes.

        Parameters
        ----------
        e : Exception
            The original exception that was caught.
        cmd : Iterable[str] | None, optional
            The command that was being executed when the error occurred, by default None
        stdout : str | None, optional
            The standard output captured from the command, by default None
        stderr : str | None, optional
            The standard error captured from the command, by default None
        trace : str | None, optional
            The traceback of the exception, by default None
        logfile : Path | None, optional
            The path to the log file, by default None
        creation_trace : str | None, optional
            The trace of the element creation, by default None
        phase : str, optional
            The phase of the error, by default "RUN"
        color : bool, optional
            Whether to use color in the output, by default True
        show_stack : bool, optional
            Whether to show the stack trace, by default True
        """
        self.e = e
        self.phase = phase
        self.type = type(e)
        self._creation_trace = creation_trace
        self.cmd = list(cmd) if cmd is not None else _build_cmd_list(e)
        err, out = _resolve_err_out(stderr, stdout, e)
        self.stderr = err
        self.stdout = out
        self.trace = trace
        self.logfile = logfile
        self._creation_trace = None
        self.color = color
        self.show_stack = show_stack
        self.linewidth = 120
        self.sepchar = constants.SEPCHAR

    def as_dict(self):
        """
        Convert the ErrorInfo object to a dictionary representation.

        Returns
        -------
        dict
            A dictionary containing the error information.
        """
        return {
            "cmd": self.cmd,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "trace": self.trace,
        }

    @property
    def pretty_ansi_text(self):
        """
        Get the pretty ANSI text representation of the error.

        Returns
        -------
        str
            The pretty ANSI text representation of the error.
        """
        if not hasattr(self, "_pretty_ansi_text"):
            self._pretty_ansi_text = self.build_pretty_ansi_text()
        return self._pretty_ansi_text

    @property
    def log_text(self) -> str:
        """Get the log text representation of the error."""
        if not hasattr(self, "_log_text"):
            self._log_text = self.build_log_text()
        return self._log_text

    @property
    def panel(self):
        """Get the error panel representation of the error."""
        if not hasattr(self, "_error_panel"):
            self._error_panel = self.build_error_panel()
        return self._error_panel

    @property
    def creation_trace(self) -> str | None:
        """Get the creation trace of the error."""
        return self._creation_trace

    def add_creation_trace(self, node: Element):
        """Add the creation trace from an element to the error."""
        creation_trace = getattr(node, "creation_trace", None)
        self._creation_trace = creation_trace

    def build_error_panel(self):
        elements = []

        # Header
        elements.append(Text("External command failed", style="bold red"))
        elements.append(Text(f"Phase: {self.phase}", style="bold yellow"))

        # Stacktrace
        if self.show_stack:
            if self.trace:
                elements.append(
                    Panel(
                        self.trace, title="Stacktrace", border_style="magenta"
                    )  # noqa: E501
                )
            else:
                elements.append(
                    Traceback.from_exception(
                        type(self.e), self.e, self.e.__traceback__
                    )  # noqa: E501
                )

        # Command
        elements.append(
            Panel(Syntax(shlex.join(self.cmd), "bash"), title="Command")
        )  # noqa: E501

        # Creation Trace
        if self.creation_trace:
            elements.append(
                Panel(
                    Syntax(self.creation_trace, "python"),
                    title="Element Creation",
                    border_style="yellow",
                )
            )
        # STDERR
        elements.append(
            Panel(self.stderr, title="STDERR", border_style="red")
        )  # noqa: E501

        # STDOUT
        elements.append(
            Panel(self.stdout, title="STDOUT", border_style="blue")
        )  # noqa: E501

        return Panel(
            Group(*elements), title="Pipeline Error", border_style="red"
        )  # noqa: E501

    def build_log_text(self) -> str:
        parts = []
        section_end = self.sepchar * self.linewidth + "\n"

        parts.append("\n")
        section = "External command failed".center(self.linewidth, self.sepchar)
        parts.append(section)
        parts.append(f"Phase: {self.phase}")
        parts.append(f"Return code: {_resolve_returncode(self.e)}\n")
        parts.append(section_end)

        if self.trace:
            section = "STACKTRACE".center(self.linewidth, self.sepchar)
            parts.append("\n" + section)
            parts.append(self.trace.strip())
            parts.append(section_end)

        if self.creation_trace:
            section = "ELEMENT CREATION".center(self.linewidth, self.sepchar)
            parts.append("\n" + section)
            parts.append(self.creation_trace.strip())
            parts.append(section_end)

        section = "PROCESS".center(self.linewidth, self.sepchar)
        parts.append("\n" + section)
        if self.cmd:
            parts.append("\n ---Command---:")
            parts.append(shlex.join(self.cmd))

        if self.stderr:
            parts.append("\n ---STDERR---:")
            parts.append(self.stderr.rstrip())

        if self.stdout:
            parts.append("\n ---STDOUT---:")
            parts.append(self.stdout.rstrip())
        parts.append(section_end)

        section = "LOGFILE".center(self.linewidth, self.sepchar)
        if self.logfile:
            parts.append("\n" + section)
            parts.append(f"\nFull log file: {self.logfile}")
            parts.append(section_end)

        return "\n".join(parts) + "\n"

    def build_pretty_ansi_text(self) -> str:
        print_text = format_called_process_error(
            e=self.e,
            cmd=self.cmd,
            stdout=self.stdout,
            stderr=self.stderr,
            show_stack=False,  # we handle stack manually
            trace=self.trace,
            color=self.color,
        )
        # -------- STACKTRACE BLOCK --------
        stack_block = ""
        if self.show_stack:
            stack_sep = "STACKTRACE".center(self.linewidth, self.sepchar)
            stack_content = self.trace or traceback.format_exc()

            color = self.color
            if self.color is None:
                color = _isatty(sys.stderr) or _isatty(sys.stdout)

            stack_sep_c = _ansi(stack_sep, "1;35", color)  # magenta
            stack_content_c = _ansi(stack_content.strip(), "0;37", color)

            stack_block = (
                f"{stack_sep_c}\n" f"{stack_content_c}\n" f"{stack_sep_c}\n\n"
            )  # noqa: E501
        creation_block = ""
        if self.creation_trace:
            creation_block = (
                f"{'ELEMENT CREATION'.center(self.linewidth, self.sepchar)}\n"
                f"{self.creation_trace}\n"
                f"{'ELEMENT CREATION'.center(self.linewidth, self.sepchar)}\n\n"
            )
        # -------- OPTIONAL LOGFILE HINT --------
        logfile_block = ""
        if self.logfile:
            logfile_block = f"\nFull log file: {self.logfile}\n"

        # -------- FINAL TEXT --------
        final_text = stack_block + creation_block + print_text + logfile_block
        return final_text


class PipelineError(Exception):
    def __init__(
        self,
        e: Exception,
        cmd=None,
        stdout=None,
        stderr=None,
        logfile=None,
        phase="RUN",
    ):
        """
        A PipelineError wraps any exception that occurs during pipeline execution,
        especially subprocess.CalledProcessError from external commands, and enriches it
        with additional context like the command, stdout/stderr, stack trace, and phase
        of execution.

        Parameters
        ----------
        e : Exception
            The original exception that was caught.
        cmd : _type_, optional
            The command that was being executed when the exception occurred, by default
            None
        stdout : _type_, optional
            The standard output captured from the command, by default None
        stderr : _type_, optional
            The standard error captured from the command, by default None
        logfile : _type_, optional
            The path to the log file, by default None
        phase : str, optional
            The phase of execution when the error occurred, by default "RUN"
        """
        self.info = PipelineError.build_error_info_from_exception(
            e=e,
            cmd=cmd,
            stdout=stdout,
            stderr=stderr,
            logfile=logfile,
            phase=phase,
        )
        self.exception = e
        self.exception_type: type | None = (
            type(self.exception) if self.exception else None
        )
        super().__init__(self.info.log_text)

    @staticmethod
    def build_error_info_from_exception(
        *,
        e: Exception,
        cmd=None,
        stdout=None,
        stderr=None,
        logfile=None,
        phase="RUN",
    ) -> ErrorInfo:
        """
        Build error information from an exception.

        Parameters
        ----------
        e : Exception
            The original exception that was caught.
        cmd : _type_, optional
            The command that was being executed when the exception occurred, by default
            None
        stdout : _type_, optional
            The standard output captured from the command, by default None
        stderr : _type_, optional
            The standard error captured from the command, by default None
        logfile : _type_, optional
            The path to the log file, by default None
        phase : str, optional
            The phase of execution when the error occurred, by default "RUN"

        Returns
        -------
        ErrorInfo
            An ErrorInfo object containing detailed information about the error.
        """
        stack = "".join(traceback.format_stack())
        exc = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        trace = stack + "\n--- Exception ---\n" + exc

        return ErrorInfo(
            e=e,
            cmd=cmd,
            stdout=stdout,
            stderr=stderr,
            trace=trace,
            logfile=logfile,
            phase=phase,
        )


# def build_error_info(
#     *,
#     e: subprocess.CalledProcessError,
#     cmd: Optional[Iterable[str]] = None,
#     stdout: Optional[str] = None,
#     stderr: Optional[str] = None,
#     logfile: Path | None = None,
#     trace: str | None = None,
#     show_stack: bool = True,
#     color: bool = True,
# ) -> ErrorInfo:
#     return ErrorInfo(
#         e=e,
#         cmd=cmd,
#         stdout=stdout,
#         stderr=stderr,
#         trace=trace,
#         logfile=logfile,
#         color=color,
#         show_stack=show_stack,
#     )


# def log_text(self) -> str:
#     cmd_list = list(self.cmd) if self.cmd else getattr(self.e, "cmd", [])
#     cmd_str = shlex.join(cmd_list)

#     out = self.stdout or getattr(self.e, "stdout", "") or ""
#     err = self.stderr or getattr(self.e, "stderr", "") or ""

#     ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     header = f"[{ts}] External command failed"
#     rc_line = f"Return code: {self.e.returncode}"

#     log_parts = [
#         "=" * 120,
#         header,
#         rc_line,
#         f"Command:\n{cmd_str}",
#     ]
#     if err.strip():
#         log_parts.append("STDERR:")
#         log_parts.append(err.rstrip())
#     if out.strip():
#         log_parts.append("STDOUT:")
#         log_parts.append(out.rstrip())
#     log_parts.append("=" * 120)
#     return "\n".join(log_parts) + "\n"


# def format_called_process_error_rich(
#     *,
#     e,
#     cmd=None,
#     stdout=None,
#     stderr=None,
#     trace=None,
#     logfile=None,
# ):
#     cmd_list = list(cmd) if cmd else getattr(e, "cmd", [])
#     cmd_str = shlex.join(cmd_list)

#     out = stdout or getattr(e, "stdout", "") or ""
#     err = stderr or getattr(e, "stderr", "") or ""

#     # --- HEADER ---
#     header = Text()
#     header.append("External command failed\n", style="bold red")
#     header.append(f"Return code: {e.returncode}\n", style="red")

#     # --- COMMAND ---
#     cmd_block = Syntax(
#         cmd_str,
#         "bash",
#         theme="monokai",
#         word_wrap=True,
#     )

#     # --- STDERR / STDOUT ---
#     blocks = []

#     if err.strip():
#         blocks.append(
#             Panel(
#                 err.strip(),
#                 title="STDERR",
#                 border_style="red",
#             )
#         )

#     if out.strip():
#         blocks.append(
#             Panel(
#                 out.strip(),
#                 title="STDOUT",
#                 border_style="blue",
#             )
#         )

#     # --- STACKTRACE ---
#     tb = None
#     if trace:
#         tb = Panel(
#             Text(trace.strip(), style="dim"),
#             title="STACKTRACE",
#             border_style="magenta",
#         )
#     else:
#         tb = Traceback.from_exception(type(e), e, e.__traceback__)

#     # --- LOGFILE ---
#     logfile_text = None
#     if logfile:
#         logfile_text = Text(f"\nLog file: {logfile}", style="dim")

#     # --- FINAL GROUP ---
#     return Panel(
#         Group(
#             header,
#             Panel(cmd_block, title="Command", border_style="cyan"),
#             tb,
#             *blocks,
#             logfile_text if logfile_text else Text(""),
#         ),
#         border_style="red",
#         title="Pipeline Error",
#     )


# def generate_rich_error_message_for_node(node, rich_msg):
#     # ---------Creation TRACE BLOCK -----------------
#     creation_trace = getattr(node, "creation_trace", None)
#     if creation_trace:
#         creation_block = (
#             "===========ELEMENT CREATION===========\n"
#             + creation_trace.strip()
#             + "\n===========ELEMENT CREATION===========\n\n"
#         )


#    # -------- LOGGING (no ANSI) --------#
#    if logger:
#        logger.error(log_text)
