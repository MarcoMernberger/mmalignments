"""
Contains a class that serves as an interface for external tools. It provides
a method to run the tool with specified parameters and handles the execution
and error checking.

Each External tool or algorithm needs to have at least a name, a version, a
source (e.g. an url, github or anything), a primary binary, a dictionary of
parameters as a property, and a method to build the command for execution.

Some tools will need to call other tools in a sequential manner.
if so, the run method allows for specifying a callback function that will be
called before and/or after the execution of the main command. This allows for
chaining multiple tools together in a flexible way.

Mostly that involves running a command in a subprocess, capturing stdout and
stderr, and checking the return code to determine if the execution was
successful. It also allows for specifying a callback function to be called
after successful execution, and it manages the creation of output files and
directories as needed.

Specific tools will be subclassing the External class and implementing the
required properties and methods to define how to run that specific tool and
how to build the command for it.
"""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from functools import cached_property, wraps
from inspect import signature
from pathlib import Path
from subprocess import CompletedProcess
from typing import IO, Any, Callable, Iterable, Mapping, TypeAlias

from mmalignments.models.parameters import (
    Params,
    ParamSet,
    ToolThreadSpec,
    initialize_param_registry,
)
from mmalignments.models.resources import (
    ResourceConfig,  # type: ignore[import]
    current_resources,
)
from mmalignments.services.errors import (
    PipelineError,
    build_error_info_from_exception,
)
from mmalignments.services.io import parents
from mmalignments.services.logging import ExternalLogger

from .elements import ElementTag, generate_element_key_name

logger = logging.getLogger(__name__)


###############################################################################
# Forwarding handler
###############################################################################


class ForwardingHandler(logging.Handler):
    """Forward WARNING and above records to a target logger (e.g. the Executor
    main logger) without re-entering the source logger.

    Attach this to the External module logger so that warnings/errors that
    happen inside any External tool are also visible in the pipeline's central
    log / console output, while INFO-level chatter stays in the per-run file
    only.

    Usage::

        # in Executor.build() or wherever the main logger is available:
        External.add_main_logger(captain)

        # tear down when the pipeline finishes:
        External.remove_main_logger()
    """

    def __init__(
        self, target: logging.Logger, level: int = logging.WARNING
    ) -> None:  # noqa: E501
        super().__init__(level)
        self._target = target

    def emit(self, record: logging.LogRecord) -> None:
        # Avoid infinite loops if target and source share a handler chain
        if record.name == self._target.name:
            return
        try:
            self._target.handle(record)
        except Exception:
            self.handleError(record)


###############################################################################
# Config
###############################################################################


@dataclass
class ExternalRunConfig:
    """Configuration for running an External tool."""

    cwd: Path | None = None
    env: dict[str, str] | None = None
    pipe_output: bool = True
    check: bool = True
    timeout: float | None = None
    threads: int = field(
        default_factory=lambda: ResourceConfig.detect().threads
    )  # noqa: E501
    multi: bool = True
    stdout: Path | None | IO = None
    stderr: Path | None | IO = None
    append: bool = False
    log_dir: Path | None = None

    def __post_init__(self) -> None:
        if self.stdout and self.stderr and self.stdout == self.stderr:
            raise ValueError(
                "stdout and stderr cannot be the same path or file object"
            )  # noqa: E501
        if self.log_dir and not isinstance(self.log_dir, Path):
            raise ValueError(f"log_dir must be a directory: {type(self.log_dir)}")
        if self.threads < 1:
            raise ValueError("threads must be a positive integer")


###############################################################################
# External Wrapper
###############################################################################


class Runnable:
    def __init__(
        self, fn: Callable[[], CompletedProcess | None], cmd: list[str], display: str
    ):
        self._fn = fn
        self.command = cmd
        self.command_display = display
        self.last_result: CompletedProcess | None = None

    def __call__(self) -> CompletedProcess | None:
        result = self._fn()
        self.last_result = result
        return result

    def __name__(self) -> str:
        return self._fn.__name__


###############################################################################
# External Wrapper
###############################################################################


class External:
    """Base class for external tools.

    Subclasses should provide at least the ``primary_binary`` property and
    may override ``build_cmd`` to change how the command line is assembled.

    The class offers a thin, well-tested wrapper around subprocess.run that
    builds a command from a parameters dictionary, runs the binary, captures
    stdout/stderr and (optionally) calls user-provided callbacks before/after
    execution.
    """

    def __init__(
        self,
        name: str,
        primary_binary: str | None = None,
        folder: Path | None = None,
        version: str | None = None,
        source: str | None = None,
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
        loglevel: int = logging.INFO,
    ) -> None:
        """Create a new External tool wrapper.

        Parameters
        ----------
        name : str
            Logical name of the tool (e.g. "bwa-mem2").
        primary_binary : str | None
            Executable name or path used to invoke the tool (may be None
            for abstract subclasses).
        folder: Path | None
            Optional path to a folder where the tool should write its output
        version : str | None
            Optional version string (used as a cache/fallback for
            get_version()).
        source : str | None
            Optional human-readable source/URL for the tool (documentation
            or homepage).
        parameters : Mapping[str, ParamSet] | ParamSet | str | Path | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
            If a file path or string is provided, it will be loaded from JSON
            and converted to ParamSet.
        loglevel : int
            Logging level for the tool's logger.
        """
        self._name = name
        self._primary_binary = primary_binary
        self._version = version
        self._folder = folder
        if version is None:
            self._version = self.get_version(version)
        self._source = source
        self._loglevel = loglevel
        if parameters is None:
            parameters = (
                Path(__file__).parent / f"{self.name}.json"
            )  # default path  TODO: put it somewhere else

        self.__init_parameters(parameters)

    def __init_parameters(
        self, parameters: Mapping[str, ParamSet] | ParamSet | str | Path
    ) -> None:
        """
        Initialize the parameters for the External tool.

        Override in Subclasses.

        Parameters
        ----------
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Known parameters for the External tool.
        """

        self.param_registry = initialize_param_registry(self.name, parameters)

    ###########################################################################
    # Properties
    ###########################################################################

    @property
    def name(self) -> str:
        return self._name

    @property
    def primary_binary(self) -> str | None:
        return self._primary_binary

    @property
    def version(self) -> str | None:
        return self._version

    @property
    def folder(self) -> Path | None:
        return self._folder

    @version.setter
    def version(self, value: str) -> None:
        self._version = self.get_version(value)

    @property
    def version_name(self) -> str:
        """Convenience property for versioned tool names, e.g. "bwa-mem2_2.2.1"."""  # noqa: E501
        if self._version:
            return f"{self._name}_{self._version}"
        return self._name

    @property
    def source(self) -> str | None:
        return self._source

    @cached_property
    def parameters(self) -> Mapping[str, ParamSet]:
        """
        Return the ParamSet for this tool.
        """
        return {
            "default": self.param_registry.default,
            **self.param_registry.by_subcommand,
        }

    @property
    def ts_format(self) -> str:
        """
        Return the format string used for timestamps in log filenames.

        Returns
        -------
        str
            The format string for timestamps.
        """
        return "%Y-%m-%d-%H-%M-%S"

    def __repr__(self) -> str:  # pragma: no cover - trivial
        return f"<External name={self.name} binary={self.primary_binary}>"

    ###########################################################################
    # Main-logger forwarding (class-level, shared across all External instances)
    ###########################################################################

    _forwarding_handler: ForwardingHandler | None = None

    @classmethod
    def add_main_logger(
        cls,
        main_logger: logging.Logger,
        level: int = logging.WARNING,
    ) -> None:
        """Attach a :class:`ForwardingHandler` to the External module logger.

        After this call every WARNING / ERROR / CRITICAL record emitted by any
        ``External`` subclass (or the module-level ``logger``) is forwarded to
        *main_logger* — typically the ``Executor`` captain logger — so that
        problems are visible in the central pipeline log and on the console
        without duplicating INFO noise.

        Parameters
        ----------
        main_logger:
            The target logger to forward records to (e.g. ``captain``).
        level:
            Minimum level to forward; defaults to ``logging.WARNING``.
        """
        # Remove any previously registered forwarding handler first
        cls.remove_main_logger()
        handler = ForwardingHandler(main_logger, level=level)
        handler.setLevel(level)
        logger.addHandler(handler)
        cls._forwarding_handler = handler

    @classmethod
    def remove_main_logger(cls) -> None:
        """Detach the forwarding handler installed by :meth:`add_main_logger`."""
        if cls._forwarding_handler is not None:
            logger.removeHandler(cls._forwarding_handler)
            cls._forwarding_handler = None

    ###########################################################################
    # Helpers
    ###########################################################################

    def abs(self, path: Path | str) -> Path:
        """Absolutize a path."""
        return Path(path).absolute()

    def strabs(self, path: Path | str) -> str:
        """Absolutize a path."""
        return str(self.abs(path))

    def ensure_binary(self) -> bool:
        """Check whether the configured primary binary is available.

        Returns
        -------
        bool
            True when the executable can be found on PATH; False otherwise.
        """
        if not self.primary_binary:
            return False
        return shutil.which(self.primary_binary) is not None

    def get_version(self, fallback: str | None = None) -> str | None:
        """Return a cached or discovered version string for the tool.

        The method first returns the value provided at construction time
        (``self._version``). If that is not set it will attempt to run the
        binary with common version flags and return the first non-empty line
        from stdout/stderr. If the binary is missing or no information can be
        determined, ``fallback`` is returned.

        Parameters
        ----------
        fallback : str | None
                Value to return when no version string can be determined.

        Returns
        -------
        str | None
                Detected version string or the provided fallback.
        """
        if self._version:
            return self._version

        if not self.primary_binary or not self.ensure_binary():
            return fallback

        for flag in ("--version", "-v", "-V", "version"):
            try:
                cp = subprocess.run(
                    [self.primary_binary, flag],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False,
                )
                out = (cp.stdout or "").strip()
                if out:
                    # return the first non-empty line
                    return out.splitlines()[0]
            except Exception:
                continue

        return fallback

    def extract_timestamp_part(
        self, current_base: str, cur_prefix: str
    ) -> datetime | None:
        """
        Extract the timestamp part from a log filename base.

        Parameters
        ----------
        current_base : str
            The base name of the log file.
        cur_prefix : str
            The expected prefix of the log file (e.g. "{toolname}_").

        Returns
        -------
        str | None
            The extracted timestamp part, or None if it cannot be parsed.
        """
        cur_ts = None
        try:
            cur_ts_with_pid = current_base[len(cur_prefix) + 1 :]
            cur_ts_str = cur_ts_with_pid.split("_")[0]
            cur_ts = self._str_to_timestamp(cur_ts_str)
        except Exception:
            # if we cannot parse current timestamp return None
            pass
        return cur_ts

    def _timestamp_to_str(self, timestamp: datetime) -> str:
        """
        Convert a datetime object to a string format used in log filenames.

        Parameters
        ----------
        timestamp : datetime
            The timestamp to convert.

        Returns
        -------
        str
            The formatted timestamp string.
        """
        return timestamp.strftime(self.ts_format)

    def _str_to_timestamp(self, cur_ts_str: str) -> datetime:
        """
        Convert a timestamp string from log filenames back to a datetime object.

        Parameters
        ----------
        cur_ts_str : str
            The timestamp string to convert.

        Returns
        -------
        datetime
            The corresponding datetime object.
        """
        ret = datetime.strptime(cur_ts_str, self.ts_format)
        return ret

    def build_element_name(
        self,
        tag: ElementTag,
        subcommand: str | None = None,
        suffix: str | None = None,
        **param_str: Any,
    ) -> tuple[str, str]:
        """
        build_element_name returns a tuple of (element key, element short_key)
        based on the provided parameters.

        Parameters
        ----------
        tag : ElementTag
            the tag for the element containing meta information about the
            identity of the element.
        subcommand : str | None, optional
            The subcommand associated with the element, by default None
        suffix : str | None, optional
            An optional suffix for the element name, by default None

        Returns
        -------
        tuple[str, str]
            A tuple containing the element key and short key.
        """
        return generate_element_key_name(
            tag, self.name, self.version, subcommand, suffix, **param_str
        )

    ###########################################################################
    # Multi-Threading
    ###########################################################################

    def apply_threads(
        self,
        params: Params | None,
        cfg: ExternalRunConfig | None,
        resources: ResourceConfig | None = None,
        subroutine: str | None = None,
    ) -> tuple[Params, ExternalRunConfig]:
        """Inject the resolved thread count into *params* and *cfg*.

        This is the **single place** where machine resources × tool capability
        → concrete thread values. Call it at the top of every high-level
        method before building the runner.

        Parameters
        ----------
        params : Params | None
            Caller-supplied params. Thread flag is added unless already set.
        cfg : ExternalRunConfig | None
            Caller-supplied run config. ``cfg.available_threads`` is set.
        resources : ResourceConfig | None
            Machine resource config. Falls back to ``ResourceConfig.detect()``.
        subroutine : str | None
            Optional subroutine name for per-subroutine specs
            (resolved via ``_thread_spec_for``).

        Returns
        -------
        tuple[Params, ExternalRunConfig]
            Updated (params, cfg) with thread values injected.
        """
        resources = resources or ResourceConfig.detect()
        thread_spec = self._thread_spec_for(subroutine)
        cfg = cfg or ExternalRunConfig()
        params = params or Params()

        if thread_spec is not None:
            # get the specific parameter set for the subroutine or default
            paramset = self.get_paramset(subroutine)
            n_threads = thread_spec.resolve(resources)

            # Always tell the scheduler how many CPUs this job needs
            cfg = ExternalRunConfig(
                **{
                    **cfg.__dict__,
                    "threads": n_threads,
                }
            )
            # Only inject if the tool has a flag AND the caller hasn't set it already   #   noqa: E501
            if thread_spec.multi and thread_spec.flag is not None:
                spec = paramset.get_spec(thread_spec.flag)
                if spec and spec.flag and spec.flag not in params:
                    params = params.override(**{spec.name: n_threads})

        return params, cfg

    def _thread_spec_for(self, subroutine: str | None) -> ToolThreadSpec | None:
        """Return the ToolThreadSpec for *subroutine* (or the default).

        Override in subclasses to provide per-subroutine specs:

        .. code-block:: python

            _thread_specs = {
                "align": ToolThreadSpec(param_flag="-t"),
                "sort":  ToolThreadSpec(param_flag="-@", fraction=0.5),
            }

            def _thread_spec_for(self, subroutine):
                return self._thread_specs.get(subroutine, self.thread_spec)
        """
        thread_spec = self.get_paramset(subroutine).get_spec(
            "_thread_spec"
        )  # ensure subroutine exists and has a spec
        if isinstance(thread_spec, ToolThreadSpec):
            return thread_spec
        else:
            return None

    ###########################################################################
    # Logging
    ###########################################################################

    def _get_log_dir(
        self, cfg: ExternalRunConfig, outpaths: list[Path] | None = None
    ) -> Path:
        """Determine the log directory based on the configuration.

        Parameters
        ----------
        cfg : ExternalRunConfig
            Configuration for the run, which may specify log_dir or cwd.

        Returns
        -------
        Path
            The directory where logs should be stored.
        """
        current = Path.cwd().resolve()
        fallback = current / ".logs"
        if self.folder:
            return self.folder / ".logs"
        elif outpaths:
            for p in outpaths:
                pp = Path(p)
                if pp.is_dir():
                    log_dir = (pp / ".logs").resolve()
                else:
                    log_dir = (pp.parent / ".logs").resolve()
                if current in log_dir.parents or current == log_dir:
                    return log_dir
        else:
            if cfg.cwd is not None:
                log_dir = (Path(cfg.cwd) / ".logs").resolve()
                if current in log_dir.parents or current == log_dir:
                    return log_dir
        return fallback

    # def get_timestamp_with_pid(self, timestamp: datetime) -> str:
    #     """Generate a timestamp string for log file naming."""
    #     timestamp_pid = f"{self._timestamp_to_str(timestamp)}_{os.getpid()}"
    #     return timestamp_pid

    # def _delete_associated_logs(
    #     self, log_dir: Path, path: Path, timestamp: datetime
    # ) -> None:
    #     stem = path.stem
    #     combined, stdout, stderr = self._get_log_files_for_run(log_dir, stem, timestamp)
    #     try:
    #         combined.unlink()
    #     except Exception as e:
    #         logger.info(
    #             f"Failed to remove old combined log {path} with {timestamp}\nException was: {e}"
    #         )
    #     try:
    #         if stdout.exists():
    #             stdout.unlink()
    #     except Exception as e:
    #         logger.info(
    #             f"Failed to remove old stdout log {stdout} with {timestamp}\nException was: {e}"
    #         )
    #     try:
    #         if stderr.exists():
    #             stderr.unlink()
    #     except Exception as e:
    #         logger.info(
    #             f"Failed to remove old stderr log {stderr} with {timestamp}\nException was: {e}"
    #         )

    # def _get_log_files_for_run(
    #     self, log_dir: Path, base: str, timestamp: datetime
    # ) -> tuple[Path, Path, Path]:
    #     combined_log_path = (
    #         log_dir / f"{base}_{self._timestamp_to_str(timestamp)}_{os.getpid()}.log"
    #     )
    #     stdout_log_path = combined_log_path.with_suffix(".stdout.log")
    #     stderr_log_path = combined_log_path.with_suffix(".stderr.log")
    #     return combined_log_path, stdout_log_path, stderr_log_path

    # def _setup_run_logging(
    #     self,
    #     log_dir: Path,
    #     filebase: str,
    #     timestamp: datetime,
    # ) -> tuple[Path, Path, Path, FileHandler | None]:
    #     """Set up log files and file handler for a run.

    #     Parameters
    #     ----------
    #     cwd : Path | None
    #         Working directory for logs; if None, uses current directory.

    #     Returns
    #     -------
    #     tuple[Path, Path, Path, FileHandler | None, str]
    #         (combined_log_path, stdout_log_path, stderr_log_path, file_handler, filebase)
    #     """
    #     ensure(log_dir)
    #     combined_log_path, stdout_log_path, stderr_log_path = (
    #         self._get_log_files_for_run(log_dir, filebase, timestamp)
    #     )
    #     # Create file handler for logger messages
    #     with open(combined_log_path, "w", encoding="utf-8") as f:
    #         f.write("\n" + "=" * 80 + "\n")
    #         f.write("LOG\n")
    #         f.write("=" * 80 + "\n\n")
    #     file_handler = FileHandler(combined_log_path, mode="a", encoding="utf-8")
    #     file_handler.setLevel(self._loglevel)
    #     fmt = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    #     file_handler.setFormatter(fmt)
    #     logger.setLevel(self._loglevel)
    #     logger.addHandler(file_handler)
    #     # return the generated base (name + timestamp) so callers can
    #     # identify and manage related per-run files
    #     return combined_log_path, stdout_log_path, stderr_log_path, file_handler

    # def _finalize_run_logging(
    #     self,
    #     combined_log_path: Path,
    #     stdout_log_path: Path,
    #     stderr_log_path: Path,
    #     file_handler: FileHandler | None,
    #     success: bool,
    # ) -> None:
    #     """Finalize logging by combining log files and adding success marker.

    #     Parameters
    #     ----------
    #     combined_log_path : Path
    #         Path to the combined log file.
    #     stdout_log_path : Path
    #         Path to the stdout log file.
    #     stderr_log_path : Path
    #         Path to the stderr log file.
    #     file_handler : FileHandler | None
    #         File handler to remove from logger.
    #     success : bool
    #         Whether the subprocess completed successfully.
    #     """
    #     try:
    #         # Remove and close the file handler first
    #         if file_handler is not None:
    #             logger.removeHandler(file_handler)
    #             file_handler.close()

    #         # Append subprocess streams to combined log
    #         with open(combined_log_path, "a", encoding="utf-8") as outf:
    #             outf.write("\n" + "=" * 80 + "\n")
    #             outf.write("SUBPROCESS OUTPUT\n")
    #             outf.write("=" * 80 + "\n\n")

    #             if stderr_log_path.exists():
    #                 outf.write("--- STDERR ---\n")
    #                 with open(stderr_log_path, "r", encoding="utf-8") as f:
    #                     outf.write(f.read())
    #                 outf.write("\n")

    #             if stdout_log_path.exists():
    #                 outf.write("--- STDOUT ---\n")
    #                 with open(stdout_log_path, "r", encoding="utf-8") as f:
    #                     outf.write(f.read())
    #                 outf.write("\n")

    #             # Add success/failure marker
    #             outf.write("\n" + "=" * 80 + "\n")
    #             if success:
    #                 outf.write("STATUS: SUCCESS\n")
    #                 outf.write("The command completed successfully.\n")
    #             else:
    #                 outf.write("STATUS: FAILED\n")
    #                 outf.write("The command failed or was interrupted.\n")
    #             outf.write("=" * 80 + "\n")

    #         # Clean up temporary per-stream log files
    #         try:
    #             if stderr_log_path.exists():
    #                 stderr_log_path.unlink()
    #             if stdout_log_path.exists():
    #                 stdout_log_path.unlink()
    #         except Exception:
    #             logger.exception("Failed to remove temporary per-stream log files")

    #     except Exception:
    #         logger.exception("Failed while finalizing run log files")

    # def _cleanup_old_logs(
    #     self,
    #     log_dir: Path,
    #     current_prefix: str,
    #     current_timestamp: datetime,
    # ) -> None:
    #     """
    #     Remove older per-run logs for this tool in *log_dir*.

    #     Files that match the pattern ``{name}_*.log`` are examined. If their
    #     timestamp (encoded in the filename) is older than the one in
    #     ``current_base`` they are removed along with their corresponding
    #     ``.stdout.log`` and ``.stderr.log`` files.

    #     Parameters
    #     ----------
    #     log_dir : Path
    #         Directory where log files are stored.
    #     current_prefix : str
    #         The prefix used in log filenames for the current run (e.g. tool name).
    #     current_timestamp : datetime
    #         The timestamp of the current run, used to identify older logs.
    #     """

    #     # extract timestamp part from current_base (strip pid suffix)
    #     current_base = f"{current_prefix}_{self._timestamp_to_str(current_timestamp)}"
    #     for p in Path(log_dir).glob(f"{current_prefix}_*.log"):
    #         stem = p.stem  # base without suffix
    #         if stem == current_base:
    #             continue
    #         ts = self.extract_timestamp_part(stem, current_prefix)
    #         if ts is None:
    #             logger.debug("Could not parse log timestamp '%s', skipping", stem)
    #             continue
    #         if ts < current_timestamp:
    #             # remove combined and related per-stream logs
    #             self._delete_associated_logs(log_dir=log_dir, path=p, timestamp=ts)
    #         # except Exception:
    #         #     logger.exception(f"Failed to process log file {p} for cleanup")
    #         #     continue
    #         #     raise

    # def _prepare_output_streams(
    #     self,
    #     cfg: ExternalRunConfig,
    #     output: Path | None = None,
    #     stdout_log_path: Path | None = None,
    #     stderr_log_path: Path | None = None,
    # ) -> tuple[Any, Any]:
    #     stdout_stream = None
    #     stderr_stream = None
    #     # decide stdout
    #     # Open per-stream files if not capturing output
    #     if cfg.capture_output:
    #         if cfg.stdout is None:
    #             stdout_stream = subprocess.PIPE
    #         else:
    #             stdout_stream = open_target(cfg.stdout, append=cfg.append)
    #         if cfg.stderr is None:
    #             stderr_stream = subprocess.PIPE
    #         else:
    #             stderr_stream = open_target(cfg.stderr, append=cfg.append)
    #     else:
    #         if stdout_log_path:
    #             stdout_stream = open(stdout_log_path, "w", encoding="utf-8")
    #         else:
    #             stdout_stream = None
    #         if stderr_log_path:
    #             stderr_stream = open(stderr_log_path, "w", encoding="utf-8")
    #         else:
    #             stderr_stream = None
    #     if output:
    #         # If output file is specified, redirect to it always (append if specified)
    #         stdout_stream = open_target(output, append=cfg.append)
    #     return stdout_stream, stderr_stream

    # def _finalize_streams(self, stdout_fh: Any, stderr_fh: Any) -> None:
    #     for fh in (stdout_fh, stderr_fh):
    #         try:
    #             if fh is None or fh is subprocess.PIPE:
    #                 continue
    #             close = getattr(fh, "close", None)
    #             if callable(close):
    #                 close()
    #         except Exception:
    #             logger.exception("Failed to close stdout file object")

    ####################################################################################
    # Command-building helpers
    ####################################################################################

    def subcommand(self, arguments: list[str] | None) -> str | None:
        if not arguments:
            return None
        # only if there are registered subcommands
        if not getattr(self.param_registry, "by_subcommand", {}):
            return None
        cand = arguments[0]
        subcommand = cand if cand in self.param_registry.by_subcommand else None
        return subcommand

    def get_paramset(self, subroutine: str | None = None) -> ParamSet:
        """
        Return the ParamSet for this tool or a specific subroutine.

        Parameters
        ----------
        subroutine : str | None, optional
            Name of the subroutine, by default None (returns the main ParamSet).
            If the tool has multiple subroutines with different parameters,
            this can be used to retrieve the specific ParamSet for that
            subroutine.

        Returns
        -------
        ParamSet
            The ParamSet for the specified subroutine or the main ParamSet if
            no subroutine is specified.
        """
        return self.param_registry.for_subcommand(subroutine)

    def signature_determinants(
        self, params: Params | None, subroutine: str | None = None
    ) -> tuple[str, ...]:
        """
        Return a list of command-line tokens that represent the run signature
        based on the provided parameters. This can be used for checking
        re-runs.

        Parameters
        ----------
        override_params : Params
                Parameters upplied to instance run.
        subroutine : str | None
                Optional tool subroutine name to check the provided parameters
                (e.g. "align" or "qc").

        Returns
        -------
        list[str]
                List of command-line tokens representing the tool's signature.
        """
        if not params:
            return ()
        paramset = self.param_registry.for_subcommand(subroutine)
        # Include parameters that affect output in a deterministic order
        signature_determinants = paramset.signature_determinants(params)
        return tuple([str(d) for d in signature_determinants])

    def to_cli(self, params: Params, subroutine: str | None = None) -> list[str]:
        """
        Return a list of command-line tokens that represent the tool's CLI
        arguments based on the provided parameters. This is used for building
        CLI command.

        Parameters
        ----------
        override_params : Params
                Parameters upplied to instance run.
        subroutine : str | None
                Optional tool subroutine name to check the provided parameters
                (e.g. "align" or "qc").

        Returns
        -------
        list[str]
                List of command-line tokens representing the tool's signature.
        """
        paramset = self.param_registry.for_subcommand(subroutine)
        cli_arguments = paramset.to_cli(params)
        return cli_arguments

    def build_cmd(
        self,
        arguments: Iterable[str] | None = None,
        params: Params | None = None,
    ) -> list[str]:
        """Build the command-line to execute.

        The command is constructed from the configured ``primary_binary`` and
        a list of positional arguments plus appropriate additional parameters.
        Additional parameters are wrapped as Params object and converted to CLI
        arguments via the ``to_cli`` method.

        Parameters
        ----------
        arguments : list[str] | None
                Positional arguments appended to the command after flags.
        params : Params | None
                Parameter values that override the instance defaults for this
                invocation.

        Returns
        -------
        List[str]
                Full command-line as a list suitable for ``subprocess.run``.
        """
        if not self.primary_binary:
            raise ValueError("primary_binary must be set to build a command")

        args = list(arguments or [])
        cmd: list[str] = [self.primary_binary] + args
        # Keep order in params provided
        if params:
            subcommand = self.subcommand(args)
            cli_args = self.to_cli(params, subroutine=subcommand)
            cmd.extend(cli_args)
        return cmd

    ###########################################################################
    # Runner
    ###########################################################################
    # External.run expects a single callable(post(cp)). Wrap the list
    # of zero-argument post-callables in a single function that will be
    # called with the CompletedProcess after the subprocess finished.

    def runnable(
        self,
        *,
        arguments: Iterable[str] | None = None,
        subcommand: str | None = None,
        output: Path | None = None,
        pre: Callable[[], CompletedProcess | None] | Runnable | None = None,
        post: Callable[[], CompletedProcess | None] | Runnable | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:
        """Execute the external tool and optionally run callbacks.

        This is a thin wrapper around ``subprocess.run`` that builds the
        command via :meth:`build_cmd`, captures stdout/stderr by default and
        raises ``subprocess.CalledProcessError`` when ``check`` is True and the
        exit code is non-zero.

        Parameters
        ----------
        arguments : Iterable[str] | None
            Positional arguments to the command call appended after the
            primary binary.
        subcommand : str | None
            Optional subcommand name for selecting parameter sets and thread
            specs.
        output : Path | None
            Optional path to a file where the command's output should be
            written.
            This is intended for commands that write to stdout and allows
            redirecting that output to a file. If the command writes to files
            directly, this can be None.
        pre : Callable[[], CompletedProcess | None] | Runnable | None
            Optional callback executed before the subprocess is started.
        post : Callable[[], CompletedProcess | None] | Runnable | None
            Optional callback executed after successful subprocess completion.
        params : Params | None
            Additional parameter overrides passed to `build_cmd` for this run.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Runnable
            A zero-argument callable which, when invoked, runs the configured
            external command and returns the CompletedProcess. The
            callable stores the most recent result in the attribute
            ``last_result``.

        Raises
        ------
        subprocess.CalledProcessError
            When ``check`` is True and the subprocess returns a non-zero code.
        """
        cfg = cfg or ExternalRunConfig()
        params = params or Params()
        cmd = self.build_cmd(arguments, params)

        # Build a zero-argument callable that executes the configured run.
        def _runner() -> CompletedProcess:
            """Execute this External invocation and return CompletedProcess.

            The closure captures the arguments passed to :meth:`run` so the
            runner can be invoked later without parameters.
            """
            # Set up logging infrastructure (always create log files)
            name = (
                f"{self.name}_{subcommand}" if subcommand else self.name
            )  # noqa: E501
            now = datetime.now()
            log_dir = self._get_log_dir(cfg)
            call_logger = ExternalLogger(log_dir, name, now)
            call_logger.prepare_output_streams(
                output=output,
                pipe_output=cfg.pipe_output,
                append=cfg.append,
                stderr_from_cfg=cfg.stderr,
                stdout_from_cfg=cfg.stdout,
            )

            def try_add(
                fn: Callable[[], CompletedProcess | None] | Runnable | None,
                phase: str,
                cmd: list[str] | None = None,
                call_logger: ExternalLogger | None = None,
                name: str | None = None,
            ) -> Any:
                if fn:
                    name = getattr(fn, "__name__", None)
                    if not name and isinstance(fn, Runnable):
                        name = fn.command_display

                    phase = f"{phase} ({name})" if name else phase
                    try:
                        return fn()
                    except PipelineError:
                        raise
                    except Exception as e:
                        error_info = build_error_info_from_exception(
                            e=e,
                            cmd=fn.command if isinstance(fn, Runnable) else cmd,
                            logfile=(
                                call_logger.combined_log_path if call_logger else None
                            ),
                            phase=phase,
                        )
                        if call_logger:
                            call_logger.log_error(error_info)
                        raise PipelineError(error_info)

            def try_main(
                run_process: Callable[[], CompletedProcess] | Runnable,
                call_logger: ExternalLogger,
                name: str | None = None,
            ) -> Any:
                cp = run_process()
                if cfg.check and cp:
                    try:
                        cp.check_returncode()
                        call_logger.write_streams(cp)
                        return cp
                    except subprocess.CalledProcessError as e:
                        print("RAISING ERROR NOW")
                        # raise
                        error_info = build_error_info_from_exception(
                            e=e,
                            cmd=cmd,
                            stdout=e.output,
                            stderr=e.stderr,
                            logfile=call_logger.combined_log_path,
                            phase=f"RUN ({name})",
                        )
                        # this is the externals logger, not the main logger
                        call_logger.log_error(error_info)
                        raise PipelineError(error_info)

            try:
                # Execute pre-callback
                try_add(pre, phase="PRE", call_logger=call_logger)
                # if pre:
                #     try:
                #         pre()
                #     except PipelineError:
                #         raise
                #     except Exception as e:
                #         error_info = build_error_info_from_exception(
                #             e=e,
                #             cmd=pre.command if isinstance(pre, Runnable) else None,
                #             logfile=(
                #             phase=f"PRE ({pre.__name__})",
                #         )
                #         call_logger.log_error(error_info)
                #         raise PipelineError(error_info)

                # Execute the main subprocess command
                def run_process() -> CompletedProcess:
                    return subprocess.run(
                        cmd,
                        cwd=str(cfg.cwd) if cfg.cwd is not None else None,
                        env=cfg.env,
                        stdout=call_logger.stdout_stream,
                        stderr=call_logger.stderr_stream,
                        text=True,
                        timeout=cfg.timeout,
                        check=False,  # we'll handle this manually to log errors
                    )

                cp = try_main(run_process, call_logger=call_logger, name=name)
                # cp = subprocess.run(
                #     cmd,
                #     cwd=str(cfg.cwd) if cfg.cwd is not None else None,
                #     env=cfg.env,
                #     stdout=call_logger.stdout_stream,
                #     stderr=call_logger.stderr_stream,
                #     text=True,
                #     timeout=cfg.timeout,
                #     check=False,
                # )
                # cp = try_call(
                #     run_process,
                #     phase=f"RUN ({name})",
                #     cmd=cmd,
                #     call_logger=call_logger,
                # )
                # Write captured output to per-stream files
                # call_logger.write_streams(cp)
                # if cfg.check:
                #     try:
                #         cp.check_returncode()
                #     except subprocess.CalledProcessError as e:
                #         print("RAISING ERROR NOW")
                #         # raise
                #         error_info = build_error_info_from_exception(
                #             e=e,
                #             cmd=cmd,
                #             stdout=e.output,
                #             stderr=e.stderr,
                #             logfile=call_logger.combined_log_path,
                #             phase="RUN",
                #         )
                #         # this is the externals logger, not the main logger
                #         call_logger.log_error(error_info)
                #         raise PipelineError(error_info)

                # Execute post-callback
                try_add(post, phase="POST", call_logger=call_logger)
                # if post:
                #     try:
                #         post()
                #     except PipelineError:
                #         raise
                #     except Exception as e:
                #         error_info = build_error_info_from_exception(
                #             e=e,
                #             cmd=post.command if isinstance(post, Runnable) else None,
                #             logfile=call_logger.combined_log_path,
                #             phase=f"POST ({post.__name__})",
                #         )
                #         call_logger.log_error(error_info)
                #         raise PipelineError(error_info)
                # Mark as successful
                call_logger.mark_success()
                return cp

            finally:
                # Close stream file objects if they were opened
                call_logger.finalize()

        # execute immediately (preserve previous behavior) and return the
        # callable so caller may re-run later by calling the returned function.
        _display = shlex.join(cmd)
        if output:
            _display += f" > {output}"
        return Runnable(_runner, cmd, _display)


###############################################################################
# Decorator
###############################################################################

SubroutineIn: TypeAlias = tuple[
    list[str],  # arguments
    str | None,  # subcommand
    list[str | Path],  # in paths
    list[str | Path],  # out paths
    str | Path | None,  # output if piped
    Callable[[], CompletedProcess | None] | None,  # pre
    Callable[[], CompletedProcess | None] | None,  # post
]


def subroutine(
    fn: Callable[
        ...,
        SubroutineIn,
    ],
) -> Callable[..., Callable[[], CompletedProcess]]:
    """
    Decorator: wrapped function returns (arguments, outputs).
    wrapper returns a runner callable (Callable[[], CompletedProcess])
    from self.runnable and sets runner.command.
    """
    sig = signature(fn)

    @wraps(fn)
    def wrapper(self, *args, **kwargs) -> Callable[[], CompletedProcess]:
        # get the bound parameters from signature
        bound = sig.bind(self, *args, **kwargs)
        bound.apply_defaults()

        # get cfg and params supplied and make sure they are not None
        params = bound.arguments.get("params", None)
        params = params or Params()
        cfg = bound.arguments.get("cfg", None)
        cfg = cfg or ExternalRunConfig()
        # what yre our ressources?
        resources = current_resources()
        # ensure output dirs
        arguments, subcommand, inpaths, outpaths, pipeoutput, pre, post = fn(
            *bound.args, **bound.kwargs
        )
        outpaths = [Path(p) for p in outpaths if p is not None]
        paths = inpaths + outpaths
        if pipeoutput:
            paths = paths + [pipeoutput]
        if not subcommand:
            subcommand = self.subcommand(arguments)

        # make sure the folders exist
        parents(*paths)

        params, cfg = self.apply_threads(
            params, cfg=cfg, resources=resources, subroutine=subcommand
        )
        arguments = [str(arg) for arg in arguments]
        runner = self.runnable(
            arguments=arguments,
            subcommand=subcommand,
            output=pipeoutput,
            params=params,
            pre=pre,
            post=post,
            cfg=cfg,
        )
        runner.threads = cfg.threads
        return runner

    return wrapper
