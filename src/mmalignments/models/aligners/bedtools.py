"""Module contains a bedtools interface for BED file processing."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Mapping, Union

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import Element, element
from mmalignments.models.overlay import (
    CfgType,
    OutSpec,
    OutType,
    ParType,
    PartialElementTag,
    SpecValues,
    TagType,
    pick_spec,
)
from mmalignments.models.tags import Method, State

from ..externals import External, SubroutineIn, subroutine
from ..parameters import Params, ParamSet

logger = logging.getLogger(__name__)


class Bedtools(External):
    """Bedtools interface for BED file processing.

    Provides methods for sorting, extending (slop), and merging BED files
    using bedtools. Returns zero-argument callables that can be chained as
    post-processing steps or Elements for pipeline integration.

    Examples
    --------
    Sort and merge a BED file::

        bt = Bedtools()
        sort_runner = bt.sort_bed("targets.bed", "targets.sorted.bed")
        merge_runner = bt.merge_bed("targets.sorted.bed", "targets.merged.bed")

    Use convenience method for sort+merge:

        bt = Bedtools()
        merged_element = bt.mergesort(
            input_element,
            genome="mm10.genome",
            output_bed="targets.merged.bed"
        )
    """

    def __init__(
        self,
        name: str = "bedtools",
        primary_binary: str = "bedtools",
        version: str | None = None,
        source: str = "https://github.com/arq5x/bedtools2",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        """Initialize Bedtools wrapper.

        Parameters
        ----------
        name : str
            Tool name (default: "bedtools").
        primary_binary : str
            Binary executable name (default: "bedtools").
        version : Optional[str]
            Version string override.
        source : str
            URL/source for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        parameters_file = Path(__file__).parent / f"{Path(__file__).stem}.json"
        parameters = parameters or parameters_file

        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """
        Get bedtools version string.

        Parameters
        ----------
        fallback : str | None, optional
                Value to return if version cannot be determined (default: None).

        Returns
        -------
        str | None
                Version string (e.g., "2.30.0") or fallback if not found.
        """
        if self._version:
            return self._version

        if not self.primary_binary or not self.ensure_binary():
            return fallback

        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            output = cp.stdout.strip()
            if output:
                # bedtools output: "bedtools v2.30.0"
                parts = output.split()
                for part in parts:
                    if part.startswith("v"):
                        return part[1:]  # Strip 'v' prefix
                    elif part[0].isdigit() and "." in part:
                        return part
        except subprocess.CalledProcessError:
            return fallback

        return fallback

    ###########################################################################
    # Sort (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def sort(
        self,
        bed_element: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        Sort a BED file using bedtools sort.

        Creates an Element that sorts the BED file from the input element.

        Parameters
        ----------
        bed_element : Element
            Element containing the BED file to be sorted (expects 'bed' in artifacts).
        tag : TagType | None
            Partial or full Element tag override. Defaults are derived from the
            input Element's tag.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BED file and extension ``bed``.
        par : ParType | None
            Parameters for the sort command. Overrides instance defaults.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Element
            New Element representing the sorted BED file, with the
            sort command as its run callable and the input element
            as a pre-requisite.

        Examples
        --------
        >>> bt = Bedtools()
        >>> sorted_elem = bt.sort(input_element)
        >>> sorted_elem.run()  # Execute the sort
        """
        input_bed = bed_element.bed
        resolved_tag = bed_element.tag.bump(
            method=Method.BEDTOOLS, state=State.SORT
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bed.parent, ext="bed"), out
        )
        output_bed = resolved_out.file

        runner = self.sort_bed(
            input_bed=input_bed,
            output_bed=output_bed,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(resolved_tag, self.version_name, "sort")
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=self.signature_determinants(par, subroutine="sort"),
            inputs=(input_bed,),
            artifacts=ArtifactSet(output_bed, primary_name="bed"),
            pres=(bed_element,),
        )

    @subroutine
    def sort_bed(
        self,
        input_bed: Path | str,
        output_bed: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """
        Sort a BED file using bedtools sort.

        Creates a zero-argument callable that sorts the input BED file and
        writes the sorted output to the specified path.

        Parameters
        ----------
        input_bed : Union[Path, str]
            Path to the unsorted BED file.
        output_bed : Union[Path, str]
            Path for the sorted output BED file.
        par : ParType | None
            Parameters for sorting call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Callable[[], CompletedProcess]
                Zero-argument callable that performs the sort when invoked.

        Examples
        --------
        >>> bt = Bedtools()
        >>> sort_runner = bt.sort_bed("targets.bed", "targets.sorted.bed")
        >>> sort_runner()  # Execute the sort
        """
        # Build command: bedtools sort -i input.bed > output.bed
        arguments = ["sort", "-i", str(input_bed)]
        subcommand = "sort"
        return arguments, subcommand, [input_bed], [output_bed], output_bed, None, None

    ###########################################################################
    # Slop (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def slop(
        self,
        bed_element: Element,
        genomesizes: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
        slop_bp: int | None = None,
    ) -> Element:
        """
        Extend BED intervals using bedtools slop.

        Creates an Element that extends intervals in the BED file from the
        input element.

        Parameters
        ----------
        bed_element : Element
            Element containing the BED file to extend (expects 'bed' in artifacts).
        genomesizes : Element
            Element containing the chromosome sizes information.
        tag : TagType | None
            Partial or full Element tag override. Defaults are derived from the
            input Element's tag.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BED file and extension ``bed``.
        par : ParType | None
            Additional parameter overrides.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).
        slop_bp : int | None
            Number of base pairs to extend on both sides (default: 100).

        Returns
        -------
        Element
            New Element representing the extended BED file, with the
            slop command as its run callable and the input element
            as a pre-requisite.

        Examples
        --------
        >>> bt = Bedtools()
        >>> slop_elem = bt.slop(sorted_element, "mm10.genome", slop_bp=100)
        >>> slop_elem.run()  # Execute the slop
        """
        input_bed = bed_element.bed
        genome_file = genomesizes.sizes
        resolved_tag = bed_element.tag.bump(
            method=Method.BEDTOOLS, state=State.PADDED
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bed.parent, ext="bed"), out
        )
        output_bed = resolved_out.file

        if slop_bp:
            if slop_bp > 0:
                par = (par or Params()).override(b=slop_bp)
            if slop_bp < 0:
                par = (par or Params()).override(l=slop_bp)

        runner = self.slop_bed(
            input_bed=input_bed,
            genome_file=genome_file,
            output_bed=output_bed,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(
            resolved_tag, self.version_name, "slop", slop_bp=slop_bp
        )
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=self.signature_determinants(par, subroutine="slop"),
            inputs=(input_bed, genome_file),
            artifacts=ArtifactSet(output_bed, primary_name="bed"),
            pres=(bed_element,),
        )

    @subroutine
    def slop_bed(
        self,
        input_bed: Path | str,
        genome_file: Path | str,
        output_bed: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Extend BED intervals using bedtools slop.

        Creates a zero-argument callable that extends intervals in the input
        BED file by the specified amount and writes the output.

        Parameters
        ----------
        input_bed : Path | str
                Path to the input BED file.
        genome_file : Path | str
                Path to the genome file (chr\tlength format).
        output_bed : Path | str
                Path for the output BED file with extended intervals.
        par : ParType | None
                Additional parameter overrides (e.g., {"l": 50, "r": 150} for asymmetric
                extension).
        cfg : CfgType | None
                Configuration for running the command.

        Returns
        -------
        Callable[[], None]
                Zero-argument callable that performs the slop when invoked.

        Examples
        --------
        >>> bt = Bedtools()
        >>> slop_runner = bt.slop_bed(
        ...     "targets.bed", "mm10.genome", "targets.pad100.bed", par={"b": 100}
        ... )
        >>> slop_runner()  # Execute the slop
        """
        # Build command: bedtools slop -i input.bed -g genome -b 100 > output.bed
        arguments = [
            "slop",
            "-i",
            self.strabs(input_bed),
            "-g",
            self.strabs(genome_file),
        ]
        if par:
            for k, v in par.items():
                arguments.extend([f"-{k}", str(v)])

        subcommand = "slop"
        return (
            arguments,
            subcommand,
            [input_bed, genome_file],
            [self.abs(output_bed)],
            self.abs(output_bed),
            None,
            None,
        )  # noqa: E501

    ###########################################################################
    # Merge (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def merge(
        self,
        bed_element: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Merge overlapping intervals using bedtools merge.

        Creates an Element that merges overlapping intervals in the BED file
        from the input element (which should be sorted).

        Parameters
        ----------
        bed_element : Element
            Element containing the sorted BED file to merge
            (expects 'bed' in artifacts).
        tag : TagType | None
            Partial or full Element tag override. Defaults are derived from the
            input Element's tag.
        out : OutType | None
            Output specification override. Defaults to the directory of the
            input BED file and extension ``bed``.
        par : ParType | None
            Parameter for the merge command. Overrides instance defaults.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Element
            New Element representing the merged BED file, with the
            merge command as its run callable and the input element
            as a pre-requisite.

        Examples
        --------
        >>> bt = Bedtools()
        >>> merged_elem = bt.merge(sorted_element)
        >>> merged_elem.run()  # Execute the merge
        """
        input_bed = bed_element.bed
        resolved_tag = bed_element.tag.bump(
            method=Method.BEDTOOLS, state=State.MERGED
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_bed.parent, ext="bed"), out
        )
        output_bed = resolved_out.file

        runner = self.merge_bed(
            input_bed=input_bed,
            output_bed=output_bed,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(resolved_tag, self.version_name, "merge")
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=self.signature_determinants(par, subroutine="merge"),
            inputs=(input_bed,),
            artifacts=ArtifactSet(output_bed, primary_name="bed"),
            pres=(bed_element,),
        )

    @subroutine
    def merge_bed(
        self,
        input_bed: Union[Path, str],
        output_bed: Union[Path, str],
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Merge overlapping intervals using bedtools merge.

        Creates a zero-argument callable that merges overlapping intervals
        in the input BED file (which must be sorted) and writes the output.

        Parameters
        ----------
        input_bed : Union[Path, str]
            Path to the sorted input BED file.
        output_bed : Union[Path, str]
            Path for the merged output BED file.
        par : ParType | None
            Parameter overrides (e.g., {"d": 100} to merge intervals within 100bp).
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Callable[[], CompletedProcess]
                Zero-argument callable that performs the merge when invoked.

        Examples
        --------
        >>> bt = Bedtools()
        >>> merge_runner = bt.merge_bed("targets.sorted.bed", "targets.merged.bed")
        >>> merge_runner()  # Execute the merge
        """
        # Build command: bedtools merge -i input.bed > output.bed
        arguments = ["merge", "-i", str(input_bed)]
        subcommand = "merge"
        return arguments, subcommand, [input_bed], [output_bed], output_bed, None, None

    ###########################################################################
    # Conevenience methods for common workflows
    ###########################################################################

    def mergesort(
        self,
        bed_element: Element,
        *,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> Element:
        """Convenience method: merge and sort a BED file in one step.

        This combines bedtools merge and sort operations.

        Parameters
        ----------
        bed_element : Element
            Element containing the BED file to sort and merge.
        tag : TagType | SpecValues[TagType] | None
            Tag override, either for both steps or per step using
            ``SpecValues``. For per-step overrides, declare
            ``valid_steps=("merge", "sort")`` and put values in the
            ``values`` mapping.
        out : OutType | SpecValues[OutType] | None
            Output specification override, either for both steps or per step
            using ``SpecValues``. Use per-step values if a non-default ``out``
            is given, since both steps would otherwise write to the same file.
        par : ParType | SpecValues[ParType] | None
            Parameters, either for both steps or per step using ``SpecValues``.
        cfg : CfgType | SpecValues[CfgType] | None
            Run configuration, either for both steps or per step using
            ``SpecValues``.


        Returns
        -------
        Element
            New Element representing the sorted and merged BED file.

        Examples
        --------
        >>> bt = Bedtools()
        >>> merged_elem = bt.mergesort(input_element)
        >>> merged_elem.run()  # Execute sort+merge
        """
        merge_element = self.merge(  # type: ignore[reportCallIssue]
            bed_element,
            tag=pick_spec(tag, "merge"),
            out=pick_spec(out, "merge"),
            par=pick_spec(par, "merge"),
            cfg=pick_spec(cfg, "merge"),
        )
        sort_element = self.sort(  # type: ignore[reportCallIssue]
            merge_element,
            tag=pick_spec(tag, "sort"),
            out=pick_spec(out, "sort"),
            par=pick_spec(par, "sort"),
            cfg=pick_spec(cfg, "sort"),
        )
        return sort_element

    def pad(
        self,
        targets: Element,
        genomesizes: Element,
        *,
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
        slop_bp: int = 100,
    ) -> Element:
        """
        Convenience method to pad BED intervals by a specified number of base
        pairs.

        Parameters
        ----------
        targets : Element
            Element containing the targets BED file to pad.
        genomesizes : Element
            Element containing the genome sizes file.
        tag : TagType | SpecValues[TagType] | None
            Tag override, either for all steps or per step using
            ``SpecValues``. For per-step overrides, declare
            ``valid_steps=("sort", "slop", "merge")`` and put values in the
            ``values`` mapping.
        out : OutType | SpecValues[OutType] | None
            Output specification override, either for all steps or per step
            using ``SpecValues``. The final sort step always uses defaults.
        par : ParType | SpecValues[ParType] | None
            Parameters, either for all steps or per step using ``SpecValues``.
        cfg : CfgType | SpecValues[CfgType] | None
            Run configuration, either for all steps or per step using
            ``SpecValues``.
        slop_bp : int, optional
            Number of base pairs to pad each interval, by default 100

        Returns
        -------
        Element
            New Element representing the padded BED file.
        """
        tag_sort = pick_spec(tag, "sort")
        sorted_element = self.sort(  # type: ignore[reportCallIssue]
            bed_element=targets,
            tag=tag_sort,
            out=pick_spec(out, "sort"),
            par=pick_spec(par, "sort"),
            cfg=pick_spec(cfg, "sort"),
        )
        slopped = self.slop(  # type: ignore[reportCallIssue]
            bed_element=sorted_element,
            genomesizes=genomesizes,
            tag=pick_spec(tag, "slop"),
            out=pick_spec(out, "slop"),
            par=pick_spec(par, "slop"),
            cfg=pick_spec(cfg, "slop"),
            slop_bp=slop_bp,
        )
        padded = PartialElementTag(state=State.PADDED)
        tag_sort = padded if tag_sort is None else tag_sort.resolve(padded)
        merged_element = self.mergesort(
            bed_element=slopped,
            tag=SpecValues(
                None,
                valid_steps=("merge", "sort"),
                values={"sort": tag_sort, "merge": pick_spec(tag, "merge")},
            ),
            out=SpecValues(
                None,
                valid_steps=("merge", "sort"),
                values={"merge": pick_spec(out, "merge")},
            ),
            par=par,
            cfg=cfg,
        )
        return merged_element
