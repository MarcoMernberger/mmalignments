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

from mmalignments.models.elements import Element, MappedElement, element
from mmalignments.models.tags import (
    ElementTag,
    Method,
    Stage,
    State,
    merge_tag,
    PartialElementTag,
    from_prior,
)

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
        version : str | None
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

    ####################################################################################
    # Elements and subroutines
    ####################################################################################

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
            Additional mosdepth parameters (e.g., Params(no_per_base=True)).
        cfg : ExternalRunConfig | None
            Configuration for the external run, including threads and working directory.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes mosdepth.

        Examples
        --------
        >>> mosdepth = Mosdepth()
        >>> coverage = mosdepth.depth_coverage(
        ...     input_bam="sample.bam",
        ...     output_prefix="qc/sample.targets",
        ...     targets="targets.padded.bed",
        ...     cfg=ExternalRunConfig(threads=8)
        ... )
        >>> coverage()
        """
        prefix = str(Path(output_prefix).absolute())
        arguments: list[str] = [prefix]

        # Finally add the BAM/CRAM input

        arguments.append(self.strabs(input_bam))
        return arguments, [prefix], None, None, None

    @element
    def coverage(
        self,
        mapped: MappedElement,
        targets: Element | Path | str,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        fileprefix: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Calculate coverage depth from a MappedElement.

        Creates an Element that calculates coverage depth statistics on
        specified target regions using mosdepth.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        targets : Element | Path | str
            Element or BED file defining target regions.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default
            naming. If not provided, a default tag will be generated based on
            the input Element's name.
        outdir : Path | str | None
            Directory for the sorted BAM file. If not provided, defaults to
            the same directory as the input BAM file.
        fileprefix : Path | str | None
            Filename override. If not provided, defaults to ``tag.default_output``.
        params : Params | None
            Additional mosdepth parameters.
        cfg : ExternalRunConfig | None
            Configuration for running the mosdepth command (e.g., working
            directory, threads).

        Returns
        -------
        Element
                An Element that executes mosdepth when run.

        Examples
        --------
        >>> mosdepth = Mosdepth()
        >>> from mmalignments.models.externals import ExternalRunConfig
        >>> depth_elem = mosdepth.coverage(
        ...     mapped, targets="targets.bed", cfg=ExternalRunConfig(threads=8)
        ... )
        >>> depth_elem.run()
        """
        default_tag = from_prior(
            mapped.tag,
            tag,
            stage=Stage.QC,
            method=Method.MOSDEPTH,
            state=State.REPORT,
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag

        input_bam = mapped.bam
        qc_dir = Path(outdir or input_bam.parent / "qc" / "mosdepth")
        outprefix = qc_dir / (fileprefix or tag.default_output)

        # these are the output files with fixed suffixes mosdepth will produce
        output_summary = Path(f"{outprefix}.mosdepth.summary.txt")
        output_dist = Path(f"{outprefix}.mosdepth.global.dist.txt")
        artifacts = {
            "summary": output_summary,
            "global": output_dist,
        }
        # mosdepth creates several output files with the prefix
        no_per_base = params.get("n") or params.get("no_per_base") or False
        if not no_per_base:
            artifacts["per_base"] = Path(f"{outprefix}.per-base.bed.gz")
        if "q" in params or "quantize" in params:
            artifacts["quantized"] = Path(f"{outprefix}.quantized.bed.gz")
        if "b" in params or "by" in params:
            artifacts["regions"] = Path(f"{outprefix}.regions.bed.gz")
        if "T" in params or "thresholds" in params:
            artifacts["thresholds"] = Path(f"{outprefix}.thresholds.bed.gz")

        # Extract BED path from Element if needed
        targets_path = targets.bed if isinstance(targets, Element) else targets

        # Build prerequisites list
        pres = [mapped]
        if isinstance(targets, Element):
            pres.append(targets)

        runner = self.depth_coverage(
            input_bam=input_bam,
            output_prefix=outprefix,
            targets=targets_path,
            params=params,
            cfg=cfg,
        )
        key = f"{tag.default_name}_coverage_{self.version_name}"
        determinants = self.signature_determinants(params)
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=[input_bam, targets_path],
            artifacts=artifacts,
            pres=pres,
        )

    @staticmethod
    def callable_mb_from_mosdepth_per_base(
        self,
        input_per_base_bed_gz: Path | str,
        min_dp: int = 10,
        output_json: Path | str | None = None,
    ) -> Callable[[], dict]:
        """
        Calculate callable bases from a mosdepth per-base BED file.

        Parameters
        ----------
        input_per_base_bed_gz : Path | str
            Path to the mosdepth per-base BED file.
        min_dp : int
            Minimum depth to consider a base callable (default: 10).
        output_json : Path | str | None
            Optional path to the output JSON file.

        Returns
        -------
        Callable[[], dict]
            Zero-argument callable that returns a dictionary containing
            the callable bases and callable megabases when executed.
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
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:

        params = Params or Params(min_dp=10)
        if "per_base" not in mosdepth_element.artifacts:
            raise ValueError(
                "mosdepth_element needs per_base artifact. Run mosdepth with no_per_base=False."  # noqa: E501
            )
        default_tag = ElementTag(
            root=mosdepth_element.tag.root,
            level=mosdepth_element.tag.level + 1,
            stage=Stage.QC,
            method=Method.MOSDEPTH,
            state=State.STAT,
            omics=mosdepth_element.tag.omics,
            ext="callabele.json",
            param=f"min_dp={params.min_dp}",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag
        per_base = mosdepth_element.artifacts["per_base"]
        outdir = Path(outdir) or per_base.parent
        filename = filename or tag.default_output
        out_json = outdir / filename
        runner = self.callable_mb_from_mosdepth_per_base(
            per_base, output_json=out_json, params=params, cfg=cfg
        )
        determinants = ([f"min_dp={params.min_dp}"],)
        key = (
            f"{tag.default_name}_callable_mb_{mosdepth_element.key}_dp_{params.min_dp}"
        )
        return Element(
            key=key,
            run=runner,
            determinants=determinants,
            inputs=[per_base],
            artifacts={"callable": out_json},
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
