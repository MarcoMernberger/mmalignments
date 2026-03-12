"""Module contains a bedtools interface for BED file processing."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping, Union

from mmalignments.models.elements import Element, element
from mmalignments.models.tags import ElementTag, Method, State, merge_tag
from mmalignments.models.tags import Tag

from ..externals import External, ExternalRunConfig, subroutine
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
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
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
        parameters_file = Path(__file__).parent / "bedtools.json"
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
    # Helpers
    ###########################################################################

    def build_outfile(
        self,
        outdir: Path | str | None,
        filename: Path | str | None,
        tag: ElementTag,
        input_bed: Path,
    ) -> Path:
        """Determine the output path for a BED file based on provided parameters.

        Parameters
        ----------
        outdir : Path | str | None
            Directory for the output BED file. If None, defaults to the input BED file's
            directory.
        filename : Path | str | None
            Filename for the output BED file. If None, defaults to a name based on the
            tag.
        tag : ElementTag
            Tag containing metadata for the output file.
        input_bed : Path
            Path to the input BED file.

        Returns
        -------
        Path
            Path to the output BED file.
        """
        output_dir = outdir or input_bed.parent
        output_filename = filename or tag.default_output
        return Path(output_dir) / output_filename

    def get_tag_params_cfg(
        self,
        tags: Mapping[str, ElementTag | None] | None,
        parameters: Mapping[str, Params] | None,
        cfgs: Mapping[str, ExternalRunConfig | None] | None,
        subroutine_name: str,
    ) -> tuple[ElementTag | None, Params, ExternalRunConfig | None]:
        """Helper to get tag, params, and cfg for a given subroutine.

        This centralizes the logic for retrieving the appropriate tag, parameters,
        and configuration for a specific subroutine (e.g., "sort", "merge") from
        the provided mappings. It checks for the presence of the subroutine key
        in each mapping and returns the corresponding values or defaults.

        Parameters
        ----------
        tags : Mapping[str, ElementTag | None] | None
            Optional mapping of subroutine names to ElementTags.
        parameters : Mapping[str, Params] | None
            Optional mapping of subroutine names to Params.
        cfgs : Mapping[str, ExternalRunConfig | None] | None
            Optional mapping of subroutine names to ExternalRunConfig.
        subroutine_name : str
            Name of the subroutine to retrieve values for (e.g., "sort").

        Returns
        -------
        tuple[ElementTag | None, Params, ExternalRunConfig | None]
            A tuple containing the ElementTag (or None), Params (defaulting to empty),
            and ExternalRunConfig (or None) for the specified subroutine.
        """
        tag = tags[subroutine_name] if tags and subroutine_name in tags else None
        params = (
            parameters[subroutine_name]
            if parameters and subroutine_name in parameters
            else Params()
        )
        cfg = cfgs[subroutine_name] if cfgs and subroutine_name in cfgs else None
        return tag, params, cfg

    ###########################################################################
    # Sort (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def sort(
        self,
        bed_element: Element,
        *,
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """
        Sort a BED file using bedtools sort.

        Creates an Element that sorts the BED file from the input element.

        Parameters
        ----------
        bed_element : Element
            Element containing the BED file to be sorted (expects 'bed' in artifacts).
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default naming.
            If not provided, a default tag will be generated based on the input
            Element's tag.
        outdir : Path | str | None
            Directory for the sorted output BED file. If not provided, defaults to
            the same directory as the input BED file.
        filename : Path | str | None
            Filename for the sorted output BED file. If not provided, defaults to
            the default filename based on the tag. Overrides default filename from tag.
        params : Params | None
            Parameters for the sort command. Overrides instance defaults.
        cfg : ExternalRunConfig | None
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
        root = bed_element.tag.root or bed_element.bed.stem
        default_tag = ElementTag(
            root=root,
            level=bed_element.tag.level + 1,
            stage=bed_element.tag.stage,
            method=Method.BEDTOOLS,
            state=State.SORT,
            omics=bed_element.tag.omics,
            ext="bed",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag

        output_bed = self.build_outfile(outdir, filename, tag, input_bed)
        params = params or Params()
        cfg = cfg or ExternalRunConfig()
        determinants = self.signature_determinants(params, subroutine="sort")
        # the callable to create the file
        runner = self.sort_bed(
            input_bed=input_bed,
            output_bed=output_bed,
            params=params,
            cfg=cfg,
        )
        # build element
        element = Element(
            key=f"{tag.default_name}_sorted_{self.version_name}",
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=[input_bed],
            artifacts={"bed": output_bed},
            pres=[bed_element],
        )
        return element

    @subroutine
    def sort_bed(
        self,
        input_bed: Path | str,
        output_bed: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
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
        params : Params | None
            Parameters for sorting call.
        cfg : ExternalRunConfig | None
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
        return arguments, [output_bed], output_bed, None, None

    ###########################################################################
    # Slop (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def slop(
        self,
        bed_element: Element,
        genomesizes: Element,
        *,
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
        slop_bp: int | None = None,
    ) -> Element:
        """
        Extend BED intervals using bedtools slop.

        Creates an Element that extends intervals in the BED file from the
        input element.

        Parameters
        ----------
        input_element : Element
            Element containing the BED file to extend (expects 'bed' in artifacts).
        genomesizes : Element
            Element containing the chromosome sizes information.
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default naming.
            If not provided, a default tag will be generated based on the input
            Element's tag.
        outdir : Path | str | None
            Directory for the sorted output BED file. If not provided, defaults to
            the same directory as the input BED file.
        filename : Path | str | None
            Filename for the sorted output BED file. If not provided, defaults to
            the default filename based on the tag. Overrides default filename from tag.
        params : Params | None
            Additional parameter overrides.
        cfg : ExternalRunConfig | None
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
        default_tag = ElementTag(
            root=bed_element.tag.root or bed_element.bed.stem,
            level=bed_element.tag.level + 1,
            stage=bed_element.tag.stage,
            method=Method.BEDTOOLS,
            state=State.PADDED,
            omics=bed_element.tag.omics,
            ext="bed",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag
        output_bed = self.build_outfile(outdir, filename, tag, input_bed)

        print(f"Slop parameters ", self.get_paramset("slop"))
        if slop_bp:
            if slop_bp > 0:
                params = Params(b=slop_bp, **(params.to_dict() if params else {}))
            if slop_bp < 0:
                params = Params(l=slop_bp, **(params.to_dict() if params else {}))

        runner = self.slop_bed(
            input_bed=input_bed,
            genome_file=genome_file,
            output_bed=output_bed,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params, subroutine="slop")
        element = Element(
            key=f"{tag.default_name}_slop_{slop_bp}_{self.version_name}",
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=[input_bed, genome_file],
            artifacts={"bed": output_bed},
            pres=[bed_element],
        )
        return element

    @subroutine
    def slop_bed(
        self,
        input_bed: Path | str,
        genome_file: Path | str,
        output_bed: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
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
        params : Params | None
                Additional parameter overrides (e.g., {"l": 50, "r": 150} for asymmetric
                extension).
        cfg : ExternalRunConfig | None
                Configuration for running the command.

        Returns
        -------
        Callable[[], None]
                Zero-argument callable that performs the slop when invoked.

        Examples
        --------
        >>> bt = Bedtools()
        >>> slop_runner = bt.slop_bed(
        ...     "targets.bed", "mm10.genome", "targets.pad100.bed", slop_bp=100
        ... )
        >>> slop_runner()  # Execute the slop
        """
        paths = [input_bed, genome_file, output_bed]

        # Build command: bedtools slop -i input.bed -g genome -b 100 > output.bed
        arguments = [
            "slop",
            "-i",
            self.strabs(input_bed),
            "-g",
            self.strabs(genome_file),
        ]

        return arguments, paths, self.abs(output_bed), None, None

    ###########################################################################
    # Merge (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def merge(
        self,
        bed_element: Element,
        *,
        tag: Tag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Merge overlapping intervals using bedtools merge.

        Creates an Element that merges overlapping intervals in the BED file
        from the input element (which should be sorted).

        Parameters
        ----------
        bed_element : Element
            Element containing the sorted BED file to merge
            (expects 'bed' in artifacts).
        tag : Tag | ElementTag | None
            Partial or full Element tag for the output Element, used for default naming.
            If not provided, a default tag will be generated based on the input
            Element's tag.
        outdir : Path | str | None
            Directory for the sorted output BED file. If not provided, defaults to
            the same directory as the input BED file.
        filename : Path | str | None
            Filename for the sorted output BED file. If not provided, defaults to
            the default filename based on the tag. Overrides default filename from tag.
        params : Params | None
            Parameter for the merge command. Overrides instance defaults.
        cfg : ExternalRunConfig | None
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
        default_tag = ElementTag(
            root=bed_element.tag.root or bed_element.bed.stem,
            level=bed_element.tag.level + 1,
            stage=bed_element.tag.stage,
            method=Method.BEDTOOLS,
            state=State.MERGED,
            omics=bed_element.tag.omics,
            ext="bed",
        )
        tag = merge_tag(default_tag, tag) if tag is not None else default_tag
        output_bed = self.build_outfile(outdir, filename, tag, input_bed)
        params = params or Params()
        cfg = cfg or ExternalRunConfig()
        determinants = self.signature_determinants(params, subroutine="merge")
        # build callable
        runner = self.merge_bed(
            input_bed=input_bed,
            output_bed=output_bed,
            params=params,
            cfg=cfg,
        )
        # nuild element
        element = Element(
            key=f"{tag.default_name}_merged_{self.version_name}",
            run=runner,
            tag=tag,
            determinants=determinants,
            inputs=[input_bed],
            artifacts={"bed": output_bed},
            pres=[bed_element],
            name=None,  # name=f"{input_element.name}.merged",
        )
        return element

    @subroutine
    def merge_bed(
        self,
        input_bed: Union[Path, str],
        output_bed: Union[Path, str],
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], CompletedProcess]:
        """Merge overlapping intervals using bedtools merge.

        Creates a zero-argument callable that merges overlapping intervals
        in the input BED file (which must be sorted) and writes the output.

        Parameters
        ----------
        input_bed : Union[Path, str]
            Path to the sorted input BED file.
        output_bed : Union[Path, str]
            Path for the merged output BED file.
        params : Params | None
            Parameter overrides (e.g., {"d": 100} to merge intervals within 100bp).
        cfg : ExternalRunConfig | None
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
        return arguments, [output_bed], output_bed, None, None

    ###########################################################################
    # Conevenience methods for common workflows
    ###########################################################################

    def mergesort(
        self,
        bed_element: Element,
        *,
        tags: Mapping[str, ElementTag | None] | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        parameters: Mapping[str, Params] | None = None,
        cfgs: Mapping[str, ExternalRunConfig | None] | None = None,
    ) -> Element:
        """Convenience method: sort and merge a BED file in one step.

        This combines bedtools sort and merge operations
        (equivalent to: bedtools sort | bedtools merge).

        Parameters
        ----------
        bed_element : Element
            Element containing the BED file to sort and merge.
        tags : Mapping[str, ElementTag | None] | None
            Mandatory tags for the output Element for default naming. If not provided,
            default tags will be generated based on each input Element's tag.
        outdir : Path | str | None
            Directory for the sorted output BED file. If not provided, defaults to
            the same directory as the input BED file.
        filename : Path | str | None
            Filename for the sorted output BED file. If not provided, defaults to
            the default filename based on the tag. Overrides default filename from tag.
        parameters : Mapping[str, Params] | None
            Parameter for sort and merge commands. Each one gets its own Params
            object in the mapping (e.g., {"sort": Params(...), "merge": Params(...)}).
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Element
            New Element representing the sorted and merged BED file.

        Examples
        --------
        >>> bt = Bedtools()
        >>> merged_elem = bt.mergesort(input_element, "targets.merged.bed")
        >>> merged_elem.run()  # Execute sort+merge
        """
        tag, params, cfg = self.get_tag_params_cfg(tags, parameters, cfgs, "merge")
        merge_element = self.merge(
            bed_element,
            tag=tag,
            outdir=outdir,
            filename=filename,
            params=params,
            cfg=cfg,
        )
        tag, params, cfg = self.get_tag_params_cfg(tags, parameters, cfgs, "sort")
        sort_element = self.sort(merge_element, tag=tag, params=params, cfg=cfg)
        return sort_element

    def pad(
        self,
        targets: Element,
        genomesizes: Element,
        *,
        tags: Mapping[str, ElementTag | None] | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        parameters: Mapping[str, Params] | None = None,
        cfgs: Mapping[str, ExternalRunConfig | None] | None = None,
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
        tags : Mapping[str, ElementTag | None] | None
            Optional tags for the sort, slop, and merge steps. If not provided,
            default tags will be generated based on the input Element's tag.
        outdir : Path | str | None
            Directory for the output BED file. If not provided, defaults to the same
            directory as the input BED file.
        filename : Path | str | None
            Filename for the output BED file. If not provided, defaults to the default
            filename based on the tag. Overrides default filename from tag.
        parameters : Mapping[str, Params] | None, optional
            Parameters for sort, slop, and merge commands. Each one gets its
            own Params object in the mapping (e.g.,
            {"sort": Params(...), "slop": Params(...), "merge": Params(...)}
            ), by default None.
        cfgs : Mapping[str, ExternalRunConfig | None] | None, optional
            Configuration for running the command, by default None
        slop_bp : int, optional
            Number of base pairs to pad each interval, by default 100

        Returns
        -------
        Element
            New Element representing the padded BED file.
        """
        tag_sort, params, cfg = self.get_tag_params_cfg(tags, parameters, cfgs, "sort")
        sorted_element = self.sort(
            bed_element=targets,
            tag=tag_sort,
            outdir=outdir,
            params=params,
            cfg=cfg,
        )
        tag, params, cfg = self.get_tag_params_cfg(tags, parameters, cfgs, "slop")
        slopped = self.slop(
            bed_element=sorted_element,
            genomesizes=genomesizes,
            tag=tag,
            outdir=outdir,
            params=params,
            cfg=cfg,
            slop_bp=slop_bp,
        )
        tag_merge, _, _ = self.get_tag_params_cfg(tags, parameters, cfgs, "merge")
        tag_sort = Tag(tag_sort, state=State.PADDED)
        merged_element = self.mergesort(
            bed_element=slopped,
            outdir=outdir,
            filename=filename,
            tags={"sort": tag_sort, "merge": tag_merge},
            parameters=parameters,
            cfgs=cfgs,
        )
        return merged_element
