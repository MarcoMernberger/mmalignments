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
from typing import Callable, Mapping, cast

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import Element, MappedElement, element, FileSource
from mmalignments.models.overlay import (
    CfgType,
    OutSpec,
    OutType,
    ParType,
    TagType,
)
from mmalignments.models.tags import (
    Method,
    Stage,
    State,
)

from ..externals import External, SubroutineIn, subroutine
from ..parameters import ParamSet

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
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
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
        parameters : Mapping[str, ParamSet] | ParamSet | str | Path | None
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
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
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
        par : ParType | None
            Additional mosdepth parameters (e.g., ParType(no_per_base=True)).
        cfg : CfgType | None
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
        subcommand = None
        return arguments, subcommand, [input_bam], [prefix], None, None, None

    @element
    def coverage(
        self,
        mapped: MappedElement,
        *,
        targets: Element | FileSource | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Calculate coverage depth from a MappedElement.

        Creates an Element that calculates coverage depth statistics on
        specified target regions using mosdepth.

        Parameters
        ----------
        mapped : MappedElement
            MappedElement containing the BAM file.
        targets : Element | FileSource | None
            Element with a ``bed`` artifact defining target regions, or a FileSource pointing to a BED file.
        tag : TagType | None
            Tag override for the output Element.
        out : OutType | None
            Output specification override. ``folder`` and ``stem`` define the
            mosdepth output prefix (``folder/stem``); defaults to
            ``<bam dir>/qc/mosdepth`` and the tag's default name.
        par : ParType | None
            Additional mosdepth parameters.
        cfg : CfgType | None
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
        par = par or ParType()
        resolved_tag = mapped.tag.bump(
            stage=Stage.QC,
            method=Method.MOSDEPTH,
            state=State.REPORT,
        ).resolve(tag)
        if targets:
            targets_bed = targets.artifacts.primary.resolve()
            par = par.override(by=self.strabs(targets_bed))
        input_bam = mapped.bam
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(
                folder=input_bam.parent / "qc" / "mosdepth",
                ext="mosdepth.summary.txt",
            ),
            out,
        )
        outprefix = Path(resolved_out.folder) / resolved_out.stem

        # these are the output files with fixed suffixes mosdepth will produce
        output_summary = Path(f"{outprefix}.mosdepth.summary.txt")
        output_dist = Path(f"{outprefix}.mosdepth.global.dist.txt")
        extras: dict[str, Path] = {"global": output_dist}
        # mosdepth creates several output files with the prefix
        if not par.get("no_per_base"):
            extras["per_base"] = Path(f"{outprefix}.per-base.bed.gz")
        if "q" in par or "quantize" in par:
            extras["quantized"] = Path(f"{outprefix}.quantized.bed.gz")
        if "b" in par or "by" in par:
            extras["regions"] = Path(f"{outprefix}.regions.bed.gz")
        if "T" in par or "thresholds" in par:
            extras["thresholds"] = Path(f"{outprefix}.thresholds.bed.gz")
        artifacts = ArtifactSet(output_summary, primary_name="summary", **extras)

        # Extract BED path from Element if needed
        inputs = [input_bam, targets.bed] if targets else [input_bam]
        # Build prerequisites list
        pres = (mapped, targets) if targets else (mapped,)

        runner = self.depth_coverage(
            input_bam=input_bam,
            output_prefix=outprefix,
            targets=targets.bed if targets else None,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(resolved_tag, self.version_name, "coverage")
        determinants = self.signature_determinants(par)
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=tuple(inputs),
            artifacts=artifacts,
            pres=pres,
        )

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
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:

        par = par or ParType(min_dp=10)
        min_dp = cast(int, par.get("min_dp", 10))
        if "per_base" not in mosdepth_element.artifacts:
            raise ValueError(
                "mosdepth_element needs per_base artifact. Run mosdepth with no_per_base=False."  # noqa: E501
            )
        resolved_tag = mosdepth_element.tag.bump(
            stage=Stage.QC,
            method=Method.MOSDEPTH,
            state=State.STAT,
            flag=f"min_dp{min_dp}",
        ).resolve(tag)
        per_base = mosdepth_element.artifacts["per_base"]
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=per_base.parent, ext="callable.json"), out
        )
        out_json = resolved_out.file
        runner = self.callable_mb_from_mosdepth_per_base(
            per_base, output_json=out_json, min_dp=min_dp
        )
        determinants = (f"min_dp={min_dp}",)
        key = Element.generate_key(
            resolved_tag, self.version_name, "callable_mb", min_dp=min_dp
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(per_base,),
            artifacts=ArtifactSet(out_json, primary_name="callable"),
            pres=(mosdepth_element,),
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
