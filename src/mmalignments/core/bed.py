"""Some heler functionality for bed files."""

from pathlib import Path
from subprocess import CompletedProcess
from typing import Callable, Mapping

from mmalignments.models.elements import Element, element
from mmalignments.models.tags import (
    ElementTag,
    Method,
    State,
    merge_tag,
    PartialElementTag,
    from_prior,
)
from mmalignments.services.io import absolutize, parents


def ensemblmap_bed(
    input_bed: Path | str,
    output_bed: Path | str,
    *,
    mapping: Mapping[str, str] | None = None,
) -> Callable[[], None]:
    """
    Fix a BED file by renaming chromosomes using a mapping or removing the
    'chr' prefix.

    Parameters
    ----------
    input_bed : Path | str
        Path to the input BED file.
    output_bed : Path | str
        Path to the output BED file.
    mapping : dict[str, str] | None, optional
        Dictionary mapping old chromosome names to new names. If None, 'chr'
        prefix will be removed.
    chromosome_sizes : Path | str | None, optional
        Path to a file containing chromosome sizes. If provided, can be used
        for additional validation or processing.

    Returns
    -------
    Callable[[], CompletedProcess]
        A zero-argument callable that performs the BED file fixing when invoked.
    """
    input_bed, output_bed = absolutize(input_bed, output_bed)

    def fix_name(chr: str) -> str:
        if mapping and chr in mapping:
            return mapping[chr]
        else:
            if "_random" in chr:
                parts = chr.split("_")
                return f"{parts[1]}.1"
            elif chr.startswith("chr"):
                return chr[3:]
            else:
                raise ValueError(
                    f"Chromosome name '{chr}' not found in mapping and does not start with 'chr'."  # noqa: E501
                )

    def __fix():
        parents(output_bed)
        with input_bed.open("r") as inp, output_bed.open("w") as out:
            for line in inp:
                if line.startswith(("#", "track", "browser")):
                    out.write(line)
                    continue

                parts = line.rstrip("\n").split("\t")

                if parts:
                    chrom = parts[0]
                    parts[0] = fix_name(chrom)

                out.write("\t".join(parts) + "\n")

    return __fix


@element
def ensemblmap(
    bed_element: Element,
    *,
    tag: PartialElementTag | ElementTag | None = None,
    outdir: Path | str | None = None,
    filename: Path | str | None = None,
    mapping: Mapping[str, str] | None = None,
) -> Element:
    """
    Fix a BED file by renaming chromosomes using a mapping or removing the
    'chr' prefix.

    Parameters
    ----------
    bed_element : Element
        Element containing the input BED file.
    output_bed : Path | str
        Path to the output BED file.
    tag : Tag | ElementTag | None
        Partial or full Element tag for the output Element, used for default naming. If
        not provided, a default tag will be generated based on the input Element's tag.
    outdir : Optional[Union[Path, str]]
        Directory for the sorted output BED file. If not provided, defaults to
        the same directory as the input BED file.
    filename : Optional[Union[Path, str]]
        Filename for the sorted output BED file. If not provided, defaults to
        the default filename based on the tag. Overrides default filename from tag.
    mapping : dict[str, str] | None, optional
        Dictionary mapping old chromosome names to new names. If None, 'chr'
        prefix will be removed.

    Returns
    -------
    Element
        New Element representing the fixed BED file.
    """
    input_bed = bed_element.bed
    tag = from_prior(
        bed_element.tag,
        tag,
        method=Method.CONTIGMAP,
        state=State.HARMONIZED,
        ext="bed",
    )
    # else: full ElementTag passed – use as-is
    outdir = Path(outdir or input_bed.parent)
    filename = filename or tag.default_output
    output_bed = outdir / filename

    runner = ensemblmap_bed(
        input_bed=input_bed,
        output_bed=output_bed,
        mapping=mapping,
    )
    runner.command = ["ensemblmap_bed"]
    runner.threads = 1
    determinants = [
        "".join([f"{k}={v}" for k, v in (mapping or {}).items()]) if mapping else []
    ]
    key = f"{tag.default_name}_mapped_chromosomes_to_ensembl"
    element = Element(
        key=key,
        run=runner,
        tag=tag,
        determinants=determinants,
        inputs=[input_bed],
        artifacts={"bed": output_bed},
        pres=[bed_element],
    )
    return element
