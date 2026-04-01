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
from typing import IO, Any, Callable, Iterable, Mapping, Sequence, TypeAlias

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

    def _initalize_logger(
        self, cfg: ExternalRunConfig, name: str, now: datetime, output: Path | None
    ) -> ExternalLogger:
        """
        Initialize the ExternalLogger for the current run based on the provided
        configuration.

        Parameters
        ----------
        cfg : ExternalRunConfig
            Configuration for the run, which may specify log_dir or cwd.
        name : str
            Name of the logger.
        now : datetime
            Current timestamp for log file naming.
        output : Path | None
            An optional output path if the command call needs to be piped into an output
            file. Then the logger sets the output stream to this file.

        Returns
        -------
        ExternalLogger
            The initialized ExternalLogger for the current run.
        """
        log_dir = self._get_log_dir(cfg)
        call_logger = ExternalLogger(log_dir, name, now)
        call_logger.prepare_output_streams(
            output=output,
            pipe_output=cfg.pipe_output,
            append=cfg.append,
            stderr_from_cfg=cfg.stderr,
            stdout_from_cfg=cfg.stdout,
        )
        return call_logger

    ####################################################################################
    # Command-building helpers
    ####################################################################################

    def subcommand(self, arguments: list[str] | None) -> str | None:
        """
        A best-effor helper function that tries to resolve the subcommand of the main
        external function call. For this to work, the subcommand must have been supplied
        in the tool parameter json or have been inferred otherwise.
        If this does not work, it assumes there is no subcommand.

        Parameters
        ----------
        arguments : list[str] | None
            Command line arguments to be called.

        Returns
        -------
        str | None
            The specific subcommand if it can be resolved, else None.
        """
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
    # Run Helper
    ###########################################################################

    def _run_callback(
        self,
        fn: Callable[[], CompletedProcess | None] | Runnable | None,
        phase: str,
        cmd: list[str] | None = None,
        call_logger: ExternalLogger | None = None,
        name: str | None = None,
    ) -> Any:
        """
        A wrapper that executes pre and post callbacks that need to run prior or after
        the main command. In addition, it handles logging and error reporting.

        Parameters
        ----------
        fn : Callable[[], CompletedProcess  |  None] | Runnable | None
            the callback function to execute, which can be a simple callable or a
            Runnable instance.
        phase : str
            The phase of execution, used for logging and error reporting.
        cmd : list[str] | None, optional
            The command being executed, by default None
        call_logger : ExternalLogger | None, optional
            The logger instance used for logging errors, by default None
        name : str | None, optional
            The name of the callback function, by default None

        Returns
        -------
        Any
            The result of the callback function execution, if any.

        Raises
        ------
        PipelineError
            Pipeline-specific Error with more debug information.
        """
        if fn:
            name = getattr(fn, "__name__", None)
            if not name and isinstance(fn, Runnable):
                name = fn.command_display

            phase = f"{phase} ({name})" if name else phase
            try:
                return fn()
            except PipelineError:
                raise
            except KeyboardInterrupt:
                if call_logger:
                    call_logger.logger.warning("Execution aborted by user")
                raise
            except Exception as e:
                error = PipelineError(
                    e=e,
                    cmd=fn.command if isinstance(fn, Runnable) else cmd,
                    logfile=(call_logger.combined_log_path if call_logger else None),
                    phase=phase,
                )
                if call_logger:
                    call_logger.log_error(error.info)
                raise error

    def run_process(
        self, cmd: list[str], cfg: ExternalRunConfig, call_logger: ExternalLogger
    ) -> CompletedProcess | None:
        # the function that executes the subprocess run
        try:
            cp = subprocess.run(
                cmd,
                cwd=str(cfg.cwd) if cfg.cwd is not None else None,
                env=cfg.env,
                stdout=call_logger.stdout_stream,
                stderr=call_logger.stderr_stream,
                text=True,
                timeout=cfg.timeout,
                check=False,  # we'll handle this manually to log errors
            )
        except KeyboardInterrupt:
            call_logger.logger.warning("Execution aborted by user")
            raise
        call_logger.write_streams(cp)
        if cfg.check and cp:
            try:
                cp.check_returncode()  # may raise CalledProcessError
            except subprocess.CalledProcessError as e:
                # make sure we have the output stream for debugging
                e.stdout = cp.stdout
                e.stderr = cp.stderr
                raise
        return cp

    ###########################################################################
    # Runner
    ###########################################################################

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

        This is a wrapper around ``subprocess.run`` that builds the
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
        def _runner() -> CompletedProcess | None:
            """Execute this External invocation and return CompletedProcess.

            The closure captures the arguments passed to :meth:`run` so the
            runner can be invoked later without parameters.
            """
            # make a name for error report
            name = f"{self.name}_{subcommand}" if subcommand else self.name
            call_logger = None
            try:
                # timestamp
                now = datetime.now()

                # initialize the logger
                call_logger = self._initalize_logger(cfg, name, now, output)

                # Execute pre-callback if there is one
                self._run_callback(pre, phase="PRE", call_logger=call_logger)

                # Execute the main subprocess command

                # execute main call
                cp = self.run_process(cmd, cfg, call_logger)

                # Execute post-callback if there is one
                self._run_callback(post, phase="POST", call_logger=call_logger)

                # mark success if we got here without exceptions
                call_logger.mark_success()

                # return the CompletedProcess from the main subprocess call
                return cp

            except PipelineError:
                # catch and raise PipelineErrors from _run_callback
                raise

            # except subprocess.CalledProcessError as e:
            #     # catch CalledProcessError from main call
            #     error = PipelineError(
            #         e=e,
            #         cmd=cmd,
            #         stdout=e.stdout,
            #         stderr=e.stderr,
            #         logfile=call_logger.combined_log_path,
            #         phase=phase,
            #     )

            #     call_logger.log_error(error.info)
            #     raise error

            except Exception as e:
                # turn any Exception (including CalledProcessError) into a PipelineError
                phase = "RUN" if not name else f"RUN ({name})"
                error = PipelineError(
                    e=e,
                    cmd=cmd,
                    stdout=getattr(e, "stdout", None)
                    or getattr(call_logger, "stdout_stream", None),
                    stderr=getattr(e, "stderr", None)
                    or getattr(call_logger, "stderr_stream", None),
                    logfile=getattr(call_logger, "combined_log_path", None),
                    phase=phase,
                )
                # log if possible
                if call_logger:
                    call_logger.log_error(error.info)
                raise error

            finally:
                # Close stream file objects if they were opened
                # Errors in logging are currently just logged and repressed, so finalize
                # is safe. if that changes, implement try catch here.
                if call_logger:
                    call_logger.finalize()

        # add a display command to see the actual command call later
        _display = f"{shlex.join(cmd)} > {output}" if output else shlex.join(cmd)

        return Runnable(_runner, cmd, _display)


###############################################################################
# Decorator
###############################################################################

SubroutineIn: TypeAlias = tuple[
    list[str],  # arguments
    str | None,  # subcommand
    Sequence[str | Path],  # in paths
    Sequence[str | Path],  # out paths
    str | Path | None,  # output if piped
    Callable[[], CompletedProcess | None | Any] | None,  # pre
    Callable[[], CompletedProcess | None | Any] | None,  # post
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
        paths = list(inpaths) + list(outpaths)
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
