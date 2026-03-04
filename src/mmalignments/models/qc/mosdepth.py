"""Mosdepth wrapper for depth calculation.

Provides a Mosdepth class that wraps the `mosdepth` executable and exposes a
`run_mosdepth` method which returns a zero-argument callable (consistent with
other External-based wrappers in this project).

Usage example
-------------
from mmalignments.models.qc.mosdepth import Mosdepth
m = Mosdepth()
runner = m.run_mosdepth(
    prefix="sample", bam_file="sample.bam", threads=4, by="targets.bed"
)
runner()
"""

from __future__ import annotations

import gzip
import json
import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.tasks import Element, MappedElement, element

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet

logger = logging.getLogger(__name__)


class Mosdepth(External):
    """Wrapper for mosdepth depth calculation tool.

    The `run_mosdepth` method builds the command-line according to provided
    options and returns a zero-argument callable that executes the command.
    """

    def __init__(
        self,
        name: str = "mosdepth",
        primary_binary: str = "mosdepth",
        version: str | None = None,
        source: str = "https://github.com/brentp/mosdepth",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """
        Initialize Mosdepth wrapper.

        Parameters
        ----------
        name : str, optional
            Tool name (default: "mosdepth").
        primary_binary : str, optional
            Path to Mosdepth executable (default: "mosdepth").
        version : Optional[str], optional
            Version string override (default: None).
        source : str, optional
            URL/source for the tool (default: "https://github.com/brentp/mosdepth").
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        parameters_file = Path(__file__).parent / "mosdepth.json"
        parameters = parameters or parameters_file

        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            out = (cp.stdout or cp.stderr or "").strip()
            if out:
                # mosdepth prints something like: 'mosdepth 0.3.1'
                for token in out.split():
                    if token[0].isdigit():
                        return token
        except Exception:
            pass
        return fallback

    @subroutine
    def depth_coverage(
        self,
        input_bam: Path | str,
        output_prefix: Path | str,
        *,
        targets: Path | str | int | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Calculate coverage depth on target regions using mosdepth.

        Creates a zero-argument callable that runs mosdepth to calculate coverage depth
        statistics on specified target regions (typically from a BED file).

        Parameters
        ----------
        input_bam : Path | str
            Input BAM file.
        output_prefix : Path | str
            Output prefix for mosdepth files (e.g., "qc/sample.targets").
        targets : Path | str | int | None
            BED file defining target regions, or an integer for fixed-size windows.
        params : Params | None
            Additional mosdepth parameters as a dictionary (e.g., {"-n": True} for no-per-base).
        cfg : ExternalRunConfig | None
            Configuration for the external run, including threads and working directory.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes mosdepth.

        Examples
        --------
        >>> mosdepth = Mosdepth()
        >>> coverage = mosdepth.coverage(
        ...     input_bam="sample.bam",
        ...     output_prefix="qc/sample.targets",
        ...     by="targets.padded.bed",
        ...     cfg=ExternalRunConfig(threads=8)
        ... )
        >>> coverage()
        """
        prefix = str(Path(output_prefix).absolute())
        arguments: list[str] = [prefix]

        if isinstance(targets, (int, str, Path)):
            params = Params(by=str(targets), **params.to_dict())

        # Finally add the BAM/CRAM input
        arguments.append(self.strabs(input_bam))
        return arguments, [prefix], None, None, None

    @element
    def coverage(
        self,
        mapped: MappedElement,
        targets: Element | Path | str,
        *,
        output_file_prefix: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> "Element":
        """Calculate coverage depth from a MappedElement.

        Creates an Element that calculates coverage depth statistics on
        specified target regions using mosdepth.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        targets : Element | Path | str
            Element or BED file defining target regions.
        output_file_prefix : Path | str | None
            Output prefix for mosdepth files. If not provided, defaults to
            "qc/mosdepth/{bam_stem}.targets" in the BAM's parent directory.
        params : Params | None
            Additional mosdepth parameters as a dictionary (e.g., {"-n": True}
            for no-per-base).
        cfg : ExternalRunConfig | None
            Configuration for the external run, including threads and working
            directory.

        Returns
        -------
        Element
                An Element that executes mosdepth when run.

        Examples
        --------
        >>> mosdepth = Mosdepth()
        >>> depth_elem = mosdepth.depth(mapped, targets="targets.bed", threads=8)
        >>> depth_elem.run()
        """

        input_bam = mapped.bam
        qc_dir = input_bam.parent / "qc" / "mosdepth"
        output_prefix = Path(
            output_file_prefix or qc_dir / f"{input_bam.stem}.targets"
        ).absolute()
        output_summary = f"{output_prefix}.mosdepth.summary.txt"
        output_dist = f"{output_prefix}.mosdepth.global.dist.txt"
        artifacts = {
            "summary": Path(output_summary),
            "global_dist": Path(output_dist),
        }

        # mosdepth creates several output files with the prefix
        no_per_base = params.get("-n") or params.get("--no-per-base") or False
        if not no_per_base:
            artifacts["per_base"] = Path(f"{output_prefix}.per-base.bed.gz")
        if "-q" in params or "--quantize" in params:
            artifacts["quantized"] = Path(f"{output_prefix}.quantized.bed.gz")
        if "-b" in params or "--by" in params:
            artifacts["regions"] = Path(f"{output_prefix}.regions.bed.gz")
        if "-T" in params or "--thresholds" in params:
            artifacts["thresholds"] = Path(f"{output_prefix}.thresholds.bed.gz")

        # Extract BED path from Element if needed
        targets_path = targets.bed if isinstance(targets, Element) else targets

        # Build prerequisites list
        pres = [mapped]
        if isinstance(targets, Element):
            pres.append(targets)

        runner = self.depth_coverage(
            input_bam=input_bam,
            output_prefix=output_prefix,
            targets=targets_path,
            params=params,
            cfg=cfg,
        )

        determinants = self.signature_determinants(params)
        return Element(
            name=f"{input_bam.stem}_coverage_{self.name}",
            key=f"{input_bam.stem}_coverage_{self.version_name}",
            run=runner,
            determinants=determinants,
            inputs=[input_bam, targets_path],
            artifacts=artifacts,
            pres=pres,
            store_attributes={"multi_core": True},
        )

    @staticmethod
    def callable_mb_from_mosdepth_per_base(
        self,
        input_per_base_bed_gz: str | Path,
        min_dp: int = 10,
        output_json: str | Path | None = None,
    ) -> dict:
        """
        Calculate callable bases from a mosdepth per-base BED file.

        Parameters
        ----------
        input_per_base_bed_gz : str | Path
            Path to the mosdepth per-base BED file.
        min_dp : int, optional
            Minimum depth to consider a base callable, by default 10
        output_json : str | Path | None, optional
            Path to the output JSON file, by default None

        Returns
        -------
        dict
            Dictionary containing the callable bases and callable megabases.
        """

        def __calc_and_write():
            per_base_bed_gz = Path(input_per_base_bed_gz).absolute()
            callable_bases = 0

            with gzip.open(per_base_bed_gz, "rt") as f:
                for line in f:
                    if not line or line[0] == "#":
                        continue
                    chrom, start, end, depth = line.rstrip("\n").split("\t")[:4]
                    if int(depth) >= min_dp:
                        callable_bases += int(end) - int(start)
            ret = {
                "method": self.name,
                "min_dp": min_dp,
                "callable_bases": callable_bases,
                "callable_mb": callable_bases / 1_000_000,
            }
            open_json = Path(output_json).absolute() if output_json else None
            if open_json:
                with open(open_json, "w") as f:
                    json.dump(ret, f, indent=2)
            return ret

        return __calc_and_write

    @element
    def callable_mb(
        self,
        mosdepth_element: Element,
        *,
        output_json: str | Path = None,
        min_dp: int = 10,
    ) -> Element:

        if "per_base" not in mosdepth_element.artifacts:
            raise ValueError(
                "mosdepth_element needs per_base artifact. Run mosdepth with no_per_base=False."  # noqa: E501
            )

        per_base = mosdepth_element.artifacts["per_base"]
        out = Path(
            output_json
            or (
                Path(per_base)
                .with_suffix("")
                .with_suffix("")  # strip .gz and .bed
                .with_suffix(".callable.json")
            )
        ).absolute()
        runner = self.callable_mb_from_mosdepth_per_base(
            per_base, min_dp=min_dp, output_json=out
        )
        determinants = ([f"min_dp={min_dp}"],)
        return Element(
            name=f"{Path(per_base).stem}_callable_mb_dp{min_dp}",
            key=f"{mosdepth_element.key}_callable_mb_dp{min_dp}",
            run=runner,
            determinants=determinants,
            inputs=[per_base],
            artifacts={"callable": out},
            pres=[mosdepth_element],
        )

    @staticmethod
    def extract_callable_mb(summary_file: Path) -> tuple[int, float]:
        callable_bases = 0

        with open(summary_file) as f:
            for line in f:
                if line.startswith("CALLABLE"):
                    callable_bases = int(line.split()[1])
                    break

        callable_mb = callable_bases / 1_000_000
        return callable_bases, callable_mb
