import logging
import os
import subprocess
from datetime import datetime
from logging import FileHandler, Logger
from pathlib import Path
from typing import IO, Any, Iterator

from .errors import ErrorInfo
from .io import ensure, open_target
from .time import now_as_str, str_to_timestamp, timestamp_to_str


def initlog(
    console: bool = False, log_dir: Path | str | None = None, level=logging.INFO
) -> Logger:
    timestamp = now_as_str()
    return setup_run_logger(timestamp, console=console, log_dir=log_dir, level=level)


def setup_run_logger(
    timestamp: str,
    console: bool = False,
    log_dir: Path | str | None = None,
    level=logging.INFO,
) -> Logger:
    log_dir = log_dir or Path(".logs")
    ensure(log_dir)
    run_log = Path(log_dir) / f"run_{timestamp}.log"
    logger = logging.getLogger("pipeline")
    logger.setLevel(level)
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    if not any(
        isinstance(h, logging.FileHandler) and h.baseFilename == str(run_log)
        for h in logger.handlers
    ):
        fh = logging.FileHandler(run_log)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    if console and not any(
        isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        for h in logger.handlers
    ):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger


def enable_console(logger: logging.Logger, enabled: bool) -> None:
    for h in logger.handlers:
        if isinstance(h, logging.StreamHandler) and not isinstance(
            h, logging.FileHandler
        ):
            if enabled:
                h.setLevel(logging.INFO)
            else:
                h.setLevel(logging.CRITICAL + 1)


class ExternalLogger:
    def __init__(
        self,
        log_dir: Path,
        name: str,
        timestamp: datetime,
        log_level=logging.INFO,
    ):
        """
        __init__ _summary_

        _extended_summary_

        Parameters
        ----------
        log_dir : Path
            Directory where log files will be stored.
        name : str
            Name prefix for log files (e.g. tool name).
        timestamp : datetime
            Timestamp to include in log filenames for uniqueness.
        log_level : int, optional
            Logging level for the logger, by default logging.INFO
        """
        self.log_dir = log_dir
        self.name = name
        self.timestamp = timestamp
        self._loglevel = log_level
        self.logger = logging.getLogger(f"external.{self.name}.{id(self)}")
        self.logger.propagate = False
        self.logger.handlers.clear()
        self.combined_log_path: Path | None = None
        self.stdout_log_path: Path | None = None
        self.stderr_log_path: Path | None = None
        self.file_handler: FileHandler | None = None
        self.stdout_stream = None
        self.stderr_stream = None
        self.success = False
        self.setup()

    def current_suffix(self) -> str:
        return f"{timestamp_to_str(self.timestamp)}.log"

    def extract_timestamp(self, filepath: Path) -> datetime | None:
        try:
            return str_to_timestamp(str(filepath).split(".")[-2])
        except Exception:
            return None

    def _get_log_files_for_run(self) -> tuple[Path, Path, Path]:
        """
        Generate log file paths for the current run based on the name and
        timestamp.

        Returns
        -------
        tuple[Path, Path, Path]
            Paths for the combined, stdout, and stderr log files.
        """
        pid = os.getpid()
        timestamp = timestamp_to_str(self.timestamp)
        suffix = f".{pid}.{timestamp}.log"
        combined_log_path = self.log_dir / f"{self.name}{suffix}"
        stdout_log_path = self.log_dir / f"{self.name}.stdout{suffix}"
        stderr_log_path = self.log_dir / f"{self.name}.stderr{suffix}"
        return combined_log_path, stdout_log_path, stderr_log_path

    def setup(self) -> None:
        """Set up log files and file handler for a run."""
        ensure(self.log_dir)
        (combined_log_path, stdout_log_path, stderr_log_path) = (
            self._get_log_files_for_run()
        )
        self._cleanup_old_logs()
        # Create file handler for logger messages
        with open(combined_log_path, "w", encoding="utf-8") as f:
            f.write("\n" + +"\n")
            f.write("LOG\n")
            f.write("=" * 80 + "\n\n")
        file_handler = FileHandler(combined_log_path, mode="a", encoding="utf-8")
        file_handler.setLevel(self._loglevel)
        fmt = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        file_handler.setFormatter(fmt)
        self.logger.setLevel(self._loglevel)
        self.logger.addHandler(file_handler)
        # return the generated base (name + timestamp) so callers can
        # identify and manage related per-run files
        self.combined_log_path = combined_log_path
        self.stdout_log_path = stdout_log_path
        self.stderr_log_path = stderr_log_path
        self.file_handler = file_handler

    def is_old_log(self, path: Path) -> bool:
        if path.is_file() and path.stem.startswith(self.name) and path.suffix == ".log":
            ts = self.extract_timestamp(path)
            if ts is None:
                self.logger.debug(
                    "Could not parse log timestamp '%s', skipping", path.stem
                )
                return False
            return ts < self.timestamp
        return False

    def _collect_old_log_files(self) -> Iterator[Path]:
        """Collect all log files related to the current run."""
        for filepath in self.log_dir.iterdir():
            if self.is_old_log(filepath):
                yield filepath

    def _cleanup_old_logs(
        self,
    ) -> None:
        """
        Remove older per-run logs for this tool in *log_dir*.

        Files that match the pattern ``{name}_*.log`` are examined. If their
        timestamp (encoded in the filename) is older than the one in
        ``current_base`` they are removed along with their corresponding
        ``.stdout.log`` and ``.stderr.log`` files.

        Parameters
        ----------
        current_prefix : str
            The prefix used in log filenames for the current run
            (e.g. tool name).
        """
        for path in self._collect_old_log_files():
            try:
                path.unlink()
            except Exception as e:
                self.log_logging_error(
                    f"Failed to remove old log {path} at {self.timestamp}\nException was: {e}"
                )

    def resolve_stream(
        self, from_cfg: Path | None | IO, default_path: Path | None
    ) -> Path | IO | None:
        if from_cfg is not None:
            return from_cfg
        return default_path

    def prepare_output_streams(
        self,
        output: Path | None = None,
        pipe_output: bool = False,
        append: bool = False,
        stderr_from_cfg: Path | None | IO = None,
        stdout_from_cfg: Path | None | IO = None,
    ) -> None:
        self.stdout_stream = None
        self.stderr_stream = None
        # decide stdout
        # Open per-stream files if not capturing output

        if pipe_output:
            self.stdout_stream = subprocess.PIPE
            self.stderr_stream = subprocess.PIPE
        else:
            out = self.resolve_stream(stdout_from_cfg, self.stdout_log_path)
            err = self.resolve_stream(stderr_from_cfg, self.stderr_log_path)
            self.stdout_stream = open_target(out, append=append)
            self.stderr_stream = open_target(err, append=append)
        if output:
            # If output file is specified, redirect to it always
            self.stdout_stream = open_target(output, append=append)

    def _finalize_streams(self, stdout_fh: Any, stderr_fh: Any) -> None:
        for fh in (stdout_fh, stderr_fh):
            try:
                if fh is None or fh is subprocess.PIPE:
                    continue
                close = getattr(fh, "close", None)
                if callable(close):
                    close()
            except Exception:
                self.log_logging_error("Failed to close stdout file object")

    def get_streams(self):
        return self.stdout_stream, self.stderr_stream

    def write_streams(self, cp):
        try:
            if self.stdout_log_path and cp.stdout:
                self.stdout_log_path.write_text(cp.stdout)

            if self.stderr_log_path and cp.stderr:
                self.stderr_log_path.write_text(cp.stderr)
        except Exception as e:
            self.log_logging_error(f"Failed to write streams: {e}")

    def finalize(self) -> None:
        """Finalize logging by combining log files and adding success marker.

        Parameters
        ----------
        combined_log_path : Path
            Path to the combined log file.
        stdout_log_path : Path
            Path to the stdout log file.
        stderr_log_path : Path
            Path to the stderr log file.
        file_handler : FileHandler | None
            File handler to remove from logger.
        success : bool
            Whether the subprocess completed successfully.
        """
        try:
            # Remove and close the file handler first
            if self.file_handler is not None:
                self.logger.removeHandler(self.file_handler)
                self.file_handler.close()

            # Append subprocess streams to combined log
            if self.combined_log_path:
                with open(self.combined_log_path, mode="a", encoding="utf-8") as outf:
                    outf.write("\n" + "=" * 80 + "\n")
                    outf.write("SUBPROCESS OUTPUT\n")
                    outf.write("=" * 80 + "\n\n")

                    # add stderr to main log file
                    if self.stderr_log_path and self.stderr_log_path.exists():
                        outf.write("--- STDERR ---\n")
                        with open(self.stderr_log_path, "r", encoding="utf-8") as f:
                            outf.write(f.read())
                        outf.write("\n")

                    # add stdout to main log file
                    if self.stdout_log_path and self.stdout_log_path.exists():
                        outf.write("--- STDOUT ---\n")
                        with open(self.stdout_log_path, "r", encoding="utf-8") as f:
                            outf.write(f.read())
                        outf.write("\n")

                    # Add success/failure marker
                    outf.write("\n" + "=" * 80 + "\n")
                    if self.success:
                        outf.write("STATUS: SUCCESS\n")
                        outf.write("The command completed successfully.\n")
                    else:
                        outf.write("STATUS: FAILED\n")
                        outf.write("The command failed or was interrupted.\n")
                    outf.write("=" * 80 + "\n")

        except Exception as e:
            self.log_logging_error(f"Failed while finalizing run log files: {e}")
            # Clean up temporary per-stream log files
        try:
            if self.stderr_log_path and self.stderr_log_path.exists():
                self.stderr_log_path.unlink()
            if self.stdout_log_path and self.stdout_log_path.exists():
                self.stdout_log_path.unlink()
        except Exception as e:
            self.log_logging_error(
                f"Failed to remove temporary per-stream log files: {e}"
            )
        self._finalize_streams(self.stdout_stream, self.stderr_stream)

    def log_logging_error(self, msg: str) -> None:
        safe_log_error(self.logger, msg)

    def log_error(self, error_info: ErrorInfo) -> None:
        safe_log_error(self.logger, error_info.log_text)

    def mark_success(self) -> None:
        self.success = True


def safe_log_error(logger, msg: str):
    try:
        if logger:
            logger.error(msg)
    except Exception:
        print(f"[LOGGER FAILED]\n{msg}\n")
