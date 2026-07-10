"""Base class for tools that delegate computation to an R script via Rscript."""

from __future__ import annotations

import json
import logging
from datetime import datetime
from functools import wraps
from inspect import signature
from pathlib import Path
from typing import Any, Callable, Mapping, TypeAlias

import pandas as pd  # type: ignore[import]
import rpy2.robjects as ro  # type: ignore[import]
from rpy2.rinterface import NULL  # type: ignore[import]
from rpy2.robjects import pandas2ri  # type: ignore[import]
from rpy2.robjects.vectors import (
    DataFrame as RDataFrame,  # type: ignore[import]
)
from rpy2.robjects.vectors import (  # type: ignore[import]
    FloatVector,
    ListVector,  # type: ignore[import]
    StrVector,
)

from mmalignments.models.externals import (
    External,
    ExternalRunConfig,
    Runnable,
    subroutine,
)
from mmalignments.models.parameters import Params, ParamSet
from mmalignments.models.resources import (
    current_resources,
)
from mmalignments.services.errors import (
    PipelineError,
)
from mmalignments.services.io import parents

PROJECT_ROOT = Path(__file__).resolve().parent.parent
R_SCRIPT_DIR = PROJECT_ROOT / "r"


logger = logging.getLogger(__name__)


def _convert_value(v: Any) -> Any:
    try:
        if isinstance(v, pd.DataFrame):
            with (ro.default_converter + pandas2ri.converter).context():
                return pandas2ri.py2rpy(v)
    except ImportError:
        pass
    """Recursively convert a Python value to an rpy2-compatible R type."""
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, dict):
        return ro.ListVector({dk: _convert_value(dv) for dk, dv in v.items()})
    if isinstance(v, list) and v and all(isinstance(i, str) for i in v):
        return StrVector(v)
    if (
        isinstance(v, list)
        and v
        and all(isinstance(i, (int, float)) for i in v)  # noqa: E501
    ):
        return FloatVector(v)
    if isinstance(v, list) and not v:
        return StrVector([])
    if v is None:
        return NULL
    return v


def listvector_to_dict(x: ListVector) -> dict[str, "pd.DataFrame"]:
    return {name: pandas2ri.rpy2py(x.rx2(name)) for name in x.names}


def _r_to_pandas(r_result: Any) -> Any:
    with (ro.default_converter + pandas2ri.converter).context():

        # CASE 1: actual R data.frame
        if isinstance(r_result, RDataFrame):
            return ro.conversion.get_conversion().rpy2py(r_result)

        # CASE 2: named list of data.frames
        if isinstance(r_result, ListVector):
            return {
                name: ro.conversion.get_conversion().rpy2py(r_result.rx2(name))
                for name in r_result.names
            }

        return ro.conversion.get_conversion().rpy2py(r_result)


RSubroutineIn: TypeAlias = tuple[
    str,  # r_function
    dict[str, Any] | Callable[[], dict[str, Any]],  # payload or context
    list[Path] | None,  # inpaths
    list[Path] | None,  # outpaths
    Any | None,  # pipeoutput
    Runnable | Callable | None,  # pre
    Runnable | Callable | None,  # post
]

logger = logging.getLogger(__name__)


class RScriptInternal(External):
    """External subclass that calls R functions directly via rpy2 (in-process).

    The contract:
    - A ``.R`` source file is loaded once per call via ``source()``.
    - Individual R functions are retrieved from the global environment and
      called with the converted *payload* dict as keyword arguments.
    - Output post-processing (writing TSV/Parquet, merging annotations) is
      done in Python via the *post* callback which receives the R result.

    Subclasses implement ``@element`` high-level methods that call
    ``@rsubroutine`` low-level methods returning
    ``(r_function, payload, inpaths, outpaths, pipeoutput, pre, post)``.
    """

    def __init__(
        self,
        name: str,
        r_source_file: str | Path,
        version: str | None = None,
        source: str | None = None,
        parameters: (
            Mapping[str, ParamSet] | ParamSet | str | Path | None
        ) = None,  # noqa: E501
        loglevel: int = logging.INFO,
    ) -> None:
        self.rsource = Path(r_source_file).resolve()
        if not self.rsource.exists():
            raise FileNotFoundError(f"R source file not found: {self.rsource}")
        super().__init__(
            name=name,
            primary_binary=None,
            version=version,
            source=source,
            parameters=parameters or {},
            loglevel=loglevel,
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        try:
            ver = str(ro.r("R.version.string")[0])
            return ver
        except Exception:
            return fallback

    def _load_r_function(self, r_function: str):
        """Source the R file and return the named function object."""
        ro.r["source"](str(self.rsource))
        fn = ro.globalenv[r_function]
        if fn is None:
            raise NameError(
                f"R function {r_function!r} not found in {self.rsource}"
            )  # noqa: E501
        return fn

    def runnable(
        self,
        *,
        r_function: str,
        context: dict[str, Any] | Callable[[], dict[str, Any]],
        output: Path | None = None,
        pre: Runnable | Callable | None = None,
        post: Callable | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:
        """Build a :class:`Runnable` that calls *r_function* with *payload*.

        Parameters
        ----------
        r_function : str
            Name of the R function to call (must exist in *r_source_file*).
        context : dict[str, Any] | Callable[[], dict[str, Any]]
            Either a dict of keyword arguments to pass to the R function, or a
            zero-argument callable that returns a dict with a payload-dict and
            an optional context-dict. The context-dict is passed to the *post*
            callable.
        output : Path | None
            Primary output path (used by the logger).
        pre : Runnable | Callable | None
            Optional zero-argument callable to run before the R call.
        post : Callable | None
            Optional single-argument callable ``post(result) -> result`` called
            with the pandas DataFrame returned by the R function.
        params, cfg
            Standard pipeline parameters.
        """
        cfg = cfg or ExternalRunConfig()

        def _runner():
            _name = f"{self.rsource.stem}.{r_function}"
            call_logger = None
            try:
                now = datetime.now()
                call_logger = self._initalize_logger(cfg, _name, now, output)

                self._run_callback(pre, phase="PRE", call_logger=call_logger)

                r_fn = self._load_r_function(r_function)
                if callable(context):
                    _context = context()
                    _payload = _context["payload"]
                    post_context = _context.get("context", None)
                else:
                    _payload, post_context = context, None
                _payload = adjust_r_types(_payload)
                r_result = r_fn(**_payload)
                result = _r_to_pandas(r_result)

                # post receives (and returns) the pandas result
                if post is not None:
                    result = post(result, post_context)

                call_logger.mark_success()
                return result

            except PipelineError:
                raise

            except Exception as e:
                _phase = f"RUN ({_name})"
                error = PipelineError(
                    e=e,
                    cmd=[r_function],
                    stdout=getattr(e, "stdout", None)
                    or getattr(call_logger, "stdout_stream", "None"),
                    stderr=getattr(e, "stderr", None)
                    or getattr(call_logger, "stderr_stream", "None"),
                    logfile=getattr(call_logger, "combined_log_path", None),
                    phase=_phase,
                )
                if call_logger:
                    call_logger.log_error(error.info)
                raise error

            finally:
                if call_logger:
                    call_logger.finalize()

        return Runnable(
            _runner, [r_function], f"R: {self.rsource.stem}::{r_function}"
        )  # noqa: E501


def adjust_r_types(payload: dict[str, Any]) -> dict[str, Any]:
    """Convert Python values in *payload* to rpy2-compatible R types."""
    return {k: _convert_value(v) for k, v in payload.items()}


def rsubroutine(
    fn: Callable[
        ...,
        RSubroutineIn,
    ],
) -> Callable[..., Callable[[], Any]]:
    """Decorator for :class:`RScriptInternal` methods.

    The decorated method must return a 7-tuple::

        (r_function, payload, inpaths, outpaths, pipeoutput, pre, post)

    where *post* is a optional ``Callable[[DataFrame, dict | None], DataFrame]``
    that receives the R result (already converted to pandas) and an optional
    context dictionary, and handles file-writing / annotation merging.
    """
    sig = signature(fn)

    @wraps(fn)
    def wrapper(self, *args, **kwargs) -> Callable[[], Any]:
        bound = sig.bind(self, *args, **kwargs)
        bound.apply_defaults()

        params = bound.arguments.get("params", None) or Params()
        cfg = bound.arguments.get("cfg", None) or ExternalRunConfig()
        resources = current_resources()

        r_function, context, inpaths, outpaths, _, pre, post = fn(
            *bound.args, **bound.kwargs
        )
        outpaths = [Path(p) for p in (outpaths or []) if p is not None]
        inpaths = list(inpaths or [])
        paths = inpaths + outpaths

        parents(*paths)
        params, cfg = self.apply_threads(
            params, cfg=cfg, resources=resources, subroutine=r_function
        )
        runner = self.runnable(
            r_function=r_function,
            context=context,
            output=outpaths[0] if outpaths else None,
            pre=pre,
            post=post,
            cfg=cfg,
        )
        runner.threads = cfg.threads
        return runner

    return wrapper


class RScriptExternal(External):
    """External subclass that calls ``Rscript <script.R> <params.json>``.

    The contract is simple:
    - The R script receives one positional argument: the path to a JSON file.
    - That JSON file contains all inputs, outputs, and algorithm parameters.
    - The R script writes TSV and Parquet files to the paths specified in JSON.

    Subclasses only need to provide the path to the bundled R script via
    :meth:`r_script` and implement the ``@element`` high-level methods.
    """

    # Path to the directory that contains the bundled R scripts (override in
    # subclass or let the default point to the package's ``r/`` folder).

    def __init__(
        self,
        name: str,
        r_script_dir: Path | None = None,
        version: str | None = None,
        source: str | None = None,
        parameters: (
            Mapping[str, ParamSet] | ParamSet | str | Path | None
        ) = None,  # noqa: E501
    ) -> None:
        self.r_script_dir = r_script_dir or (Path(__file__).parent / "r")
        super().__init__(
            name=name,
            primary_binary="Rscript",
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        import subprocess

        try:
            cp = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                check=False,
            )
            # Rscript --version prints to stderr
            out = (cp.stderr or cp.stdout or "").strip()
            if out:
                return out.splitlines()[0]
        except Exception:
            pass
        return fallback

    def r_script(self, script_name: str) -> Path:
        """Return the absolute path to a bundled R script.

        Parameters
        ----------
        script_name : str
            Filename, e.g. ``"deseq2_unpaired.R"``.
        """
        path = (self.r_script_dir / script_name).resolve()
        if not path.exists():
            raise FileNotFoundError(
                f"R script not found: {path}\n"
                f"Expected it in r_script_dir={self.r_script_dir}"
            )
        return path

    def write_params_json(
        self, params_path: Path, payload: dict[str, Any]
    ) -> None:  # noqa: E501
        """Serialise *payload* to *params_path* as UTF-8 JSON.

        All Path values are converted to strings so they survive JSON
        round-trips.
        """
        parents(params_path)

        def _default(obj: Any) -> Any:
            if isinstance(obj, Path):
                return str(obj)
            raise TypeError(
                f"Object of type {type(obj)} is not JSON serializable"
            )  # noqa: E501

        params_path.write_text(
            json.dumps(payload, indent=2, default=_default), encoding="utf-8"
        )
        logger.debug("Wrote params JSON to %s", params_path)

    @subroutine
    def run_r_script(
        self,
        script_name: str,
        params_path: Path,
        output_files: list[Path],
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple:
        """Low-level subroutine: call ``Rscript <script.R> <params.json>``.

        Parameters
        ----------
        script_name : str
            R script filename within the bundled ``r/`` directory.
        params_path : Path
            Path to the JSON file that was already written by the caller.
        output_files : list[Path]
            Paths to files the R script will produce (used for output dir
            creation and caching).
        params : Params | None
            Additional CLI parameters (usually unused for R scripts).
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        5-tuple consumed by the ``@subroutine`` decorator.
        """
        script_path = self.r_script(script_name)
        arguments = [str(script_path), str(params_path)]
        # subroutine decorator expects a 7-tuple:
        # (arguments, subcommand, inpaths, outpaths, pipeoutput, pre, post)
        return (
            arguments,
            None,
            [params_path],
            list(output_files),
            None,
            None,
            None,
        )  # noqa: E501
