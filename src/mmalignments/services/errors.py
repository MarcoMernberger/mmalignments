import shlex
import subprocess
import sys
import textwrap
import traceback
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional


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


def format_called_process_error(
    *,
    e: subprocess.CalledProcessError,
    cmd: Optional[Iterable[str]] = None,
    stdout: Optional[str] = None,
    stderr: Optional[str] = None,
    show_stack: bool = False,
    trace: str | None = None,
    color: bool | None = None,
    width: int = 120,
    max_stream_chars: int = 1000,
) -> tuple[str, str]:
    """
    Returns (print_text, log_text).
    - print_text: pretty, optionally colored, possibly truncated stdout/stderr
    - log_text: no color, full stdout/stderr (or provided), not truncated unless you
      choose
    """
    # Determine command
    cmd_list = list(cmd) if cmd is not None else None
    if cmd_list is None:
        # CalledProcessError typically has .cmd; fallback to .args
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

    cmd_str = shlex.join(cmd_list)

    # Streams: prefer explicit, else use e.stderr/e.stdout/e.output
    err = stderr if stderr is not None else (e.stderr or "")
    out = (
        stdout if stdout is not None else (getattr(e, "stdout", None) or e.output or "")
    )

    # Decide color
    if color is None:
        color = _isatty(sys.stderr) or _isatty(sys.stdout)

    cmd_wrapped = cmd_str  # textwrap.fill(cmd_str, width=width, break_long_words=False, break_on_hyphens=False)  # noqa: E501

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = f"[{ts}] External command failed"
    rc_line = f"Return code: {e.returncode}"

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
        # Note: stack printed separately in handler (we keep formatting text-only here)

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
    log_text = "\n".join(log_parts) + "\n"

    return print_text, log_text


def handle_called_process_error(
    *,
    e: subprocess.CalledProcessError,
    cmd: Optional[Iterable[str]] = None,
    stdout: Optional[str] = None,
    stderr: Optional[str] = None,
    logger=None,
    logfile: Path | None = None,
    trace: str | None = None,
    exit_code: int = 1,
    show_stack: bool = True,
    color: bool | None = None,
) -> None:
    """
    Print + log a pretty error message and exit without Python traceback.
    """
    print_text, log_text = format_called_process_error(
        e=e,
        cmd=cmd,
        stdout=stdout,
        stderr=stderr,
        show_stack=False,
        trace=trace,
        color=color,
    )

    if show_stack:
        sys.stderr.write("\n")
        sys.stderr.write(trace or traceback.format_exc())
        sys.stderr.write("\n")

    if logfile:
        log_line = f"\nFull log file: {logfile}\n"
        print_text += log_line
        log_text += log_line

    # Print to stderr (better for errors)
    sys.stderr.write(print_text)

    if logger:
        logger.error(log_text)

    raise SystemExit(exit_code)
