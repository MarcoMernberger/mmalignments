"""Base class for tools that delegate computation to an R script via Rscript."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Mapping

from mmalignments.models.externals import External, ExternalRunConfig, subroutine
from mmalignments.models.parameters import Params, ParamSet
from mmalignments.services.io import parents

logger = logging.getLogger(__name__)


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
    _r_script_dir: Path | None = None

    def __init__(
        self,
        name: str,
        r_script_dir: Path | None = None,
        version: str | None = None,
        source: str | None = None,
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        self._r_script_dir = r_script_dir or (Path(__file__).parent / "r")
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
        path = (self._r_script_dir / script_name).resolve()
        if not path.exists():
            raise FileNotFoundError(
                f"R script not found: {path}\n"
                f"Expected it in r_script_dir={self._r_script_dir}"
            )
        return path

    def write_params_json(self, params_path: Path, payload: dict[str, Any]) -> None:
        """Serialise *payload* to *params_path* as UTF-8 JSON.

        All Path values are converted to strings so they survive JSON round-trips.
        """
        parents(params_path)

        def _default(obj: Any) -> Any:
            if isinstance(obj, Path):
                return str(obj)
            raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

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
        return arguments, output_files, None, None, None
