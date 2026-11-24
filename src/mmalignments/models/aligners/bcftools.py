"""
Module contains a bcftools interface for VCF file processing and filtering.
"""

from __future__ import annotations

import gzip
import json
import logging
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Callable, Iterable, Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.data import HardFilterThresholds
from mmalignments.models.elements import (
    Element,
    VcfElement,
    element,
)
from mmalignments.models.tags import (
    Method,
    Stage,
    State,
)
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
from mmalignments.services.errors import PipelineError
from mmalignments.services.io import ensure, from_json
from mmalignments.models.genome import (
    EnsemblGenome,  # type: ignore[import]
)

from ..externals import External, SubroutineIn, subroutine
from ..parameters import ParamSet

logger = logging.getLogger(__name__)


_TOTAL_FAD_HEADER = (
    '##INFO=<ID=TOTAL_FAD,Number=1,Type=Integer,Description="Total informative '
    'fragment depth across all alleles at locus prior to multiallelic splitting">'
)
_FRAG_VAF_HEADER = (
    '##FORMAT=<ID=FRAG_VAF,Number=A,Type=Float,Description="Fragment-based variant '
    "allele fraction per alt allele (alt_FAD / sum(FAD)) computed prior to "
    'multiallelic splitting">'
)


def write_vcf(
    header_lines: list[str], records: Iterable[str], output_vcf: Path
) -> Path:
    """Write VCF header lines and records, using gzip for a .gz output."""
    opener = gzip.open if output_vcf.suffix == ".gz" else open
    with opener(output_vcf, "wt", encoding="utf-8") as handle:
        for lines in (header_lines, records):
            for line in lines:
                handle.write(line.rstrip("\r\n") + "\n")
    return output_vcf


def _annotate_fragment_record(line: str) -> str:
    """Add site-wide depth and per-sample fractions from complete FAD vectors."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 10:
        return line
    format_keys = fields[8].split(":")
    if "FAD" not in format_keys:
        return line

    fad_index = format_keys.index("FAD")
    if "FRAG_VAF" not in format_keys:
        format_keys.append("FRAG_VAF")
    vaf_index = format_keys.index("FRAG_VAF")
    allele_count = 1 + len(fields[4].split(","))
    totals = []
    for index in range(9, len(fields)):
        values = fields[index].split(":")
        values.extend(["."] * (len(format_keys) - len(values)))
        fad = values[fad_index].split(",")
        values[vaf_index] = "."
        if len(fad) == allele_count and all(value.isdecimal() for value in fad):
            depths = list(map(int, fad))
            total = sum(depths)
            totals.append(total)
            values[vaf_index] = ",".join(
                str(round(depth / total, 6) if total else 0.0) for depth in depths[1:]
            )
        fields[index] = ":".join(values)

    if not totals:
        return line
    info = [
        value
        for value in fields[7].split(";")
        if value not in {"", "."} and value.split("=", 1)[0] != "TOTAL_FAD"
    ]
    fields[7] = ";".join([*info, f"TOTAL_FAD={sum(totals)}"])
    fields[8] = ":".join(format_keys)
    return "\t".join(fields)


def _annotate_fragment_vafs(input_vcf: Path, annotated_vcf: Path) -> Path:
    """Stream a VCF with pre-split TOTAL_FAD and FRAG_VAF annotations."""
    opener = gzip.open if input_vcf.suffix == ".gz" else open
    header_lines = []
    with opener(input_vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                for prefix, definition in (
                    ("##INFO=<ID=TOTAL_FAD,", _TOTAL_FAD_HEADER),
                    ("##FORMAT=<ID=FRAG_VAF,", _FRAG_VAF_HEADER),
                ):
                    if not any(header.startswith(prefix) for header in header_lines):
                        header_lines.append(definition)
                header_lines.append(line)
                return write_vcf(
                    header_lines, map(_annotate_fragment_record, handle), annotated_vcf
                )
            if not line.startswith("##"):
                break
            header_lines.append(line)
    raise ValueError(f"Missing #CHROM header in {input_vcf}")


class BCFtools(External):
    """BCFtools interface for VCF/BCF file processing and filtering.

    Provides low-level callables (path- and string-based) and high-level
    :class:`VcfElement`-returning methods for use in pipeline graphs.

    Low-level / high-level pairing (analogous to :class:`Samtools`):

    +------------------------------+----------------------------------------+
    | Low-level (paths → callable) | High-level (element API)               |
    +==============================+========================================+
    | ``view_vcf``                 | ``view``                               |
    +------------------------------+----------------------------------------+
    | ``filter_vcf``               | ``filter``                             |
    +------------------------------+----------------------------------------+
    | ``index_vcf``                | ``index``                              |
    +------------------------------+----------------------------------------+
    | ``restrict_to_targets_vcf``  | ``restrict_to_targets``                |
    +------------------------------+----------------------------------------+
    | ``hard_filter_vcf``          | ``hard_filter``                        |
    +------------------------------+----------------------------------------+
    | ``count_variants_vcf``       | ``count_variants``                     |
    +------------------------------+----------------------------------------+

    Examples
    --------
    Basic PASS filter and index::

        bt = BCFtools()
        pass_elem = bt.view(mutect_elem, pass_only=True,
                            output_vcf="tumor.pass.vcf.gz")
        idx_elem  = bt.index(pass_elem)

    Full hard-filter pipeline in one call::

        filtered = bt.hard_filter(
            mutect_elem,
            tumor_sample="TUMOR",
            targets_bed=Path("targets.bed"),
            thresholds=HardFilterThresholds(min_dp=10, min_ad=3, min_vaf=0.05),
            output_folder=Path("results/vcf/"),
        )
    """

    def __init__(
        self,
        name: str = "bcftools",
        primary_binary: str = "bcftools",
        version: str | None = None,
        source: str = "https://github.com/samtools/bcftools",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialise BCFtools wrapper.

        Parameters
        ----------
        name : str
            Tool name used for logging and artefact keys, default: "bcftools".
        primary_binary : str
            Executable name on PATH, default: "bcftools".
        version : Optional[str]
            Version string override; detected automatically when None.
        source : str
            Informational URL for the tool.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Set of parameters for invocations. If the tool has subroutines,
            this can be a mapping from subroutine names to parameter sets.
            This will be used for default parameters, validation and
            constructing cli arguments in the ``build_cmd`` function.
        """
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Return the bcftools version string.

        Parameters
        ----------
        fallback : str | None
            Returned when the version cannot be determined.

        Returns
        -------
        str | None
            Version string such as ``"1.19"`` or *fallback*.
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
                check=False,
            )
            output = (cp.stdout or "").strip()
            if output:
                first_line = output.splitlines()[0]
                parts = first_line.split()
                if len(parts) >= 2 and parts[0] == "bcftools":
                    return parts[1]
        except Exception:
            pass

        return fallback

    ###########################################################################
    # Helpers
    ###########################################################################

    ###########################################################################
    # View (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def view(
        self,
        called: VcfElement,
        *,
        targets: Element | Path | str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        pass_only: bool = False,
        biallelic: bool = False,
        samples: str | list[str] | None = None,
        variant_type: str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> VcfElement:
        """Run ``bcftools view`` on an element and return a :class:`VcfElement`.

        Parameters
        ----------
        called : VcfElement
            Input VCF element (uses ``element.artifacts["vcf"]``).
        targets : Element | Path | str | None
            BED file element for target-region restriction.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the input VCF folder,
            ext ``vcf.gz``.
        pass_only : bool
            Retain only PASS-filter records.
        biallelic : bool
            Retain only biallelic variants.
        samples : str | list[str] | None
            Sample(s) to extract.
        variant_type : str | None
            ``"snps"`` or ``"indels"`` etc.
        par : ParType | None
            Additional parameters for vcf view call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        VcfElement
            Element with ``artifacts["vcf"]`` pointing to the output file.
        """
        vcf = called.vcf
        resolved_tag = called.tag.bump(
            method=Method.BCFTOOLS, state=State.FILTER
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=vcf.parent, ext="vcf.gz"), out
        )
        out_vcf = resolved_out.file

        targets_bed = targets.bed if isinstance(targets, Element) else targets
        samples = samples or [called.tag.root]
        if not samples:
            raise ValueError("No samples specified for bcftools view.")

        runner = self.view_vcf(
            input_vcf=vcf,
            output_vcf=out_vcf,
            targets_bed=targets_bed,
            pass_only=pass_only,
            biallelic=biallelic,
            samples=samples,
            variant_type=variant_type,
            par=par,
            cfg=cfg,
        )
        tag_str = None
        if pass_only or biallelic:
            tag_str = "pass" if pass_only else ""
            tag_str += "_biallelic" if biallelic else ""

        suffix = f"{Path(targets_bed).stem}" if targets_bed else None
        key = Element.generate_key(
            resolved_tag, self.version_name, "view", suffix=suffix, params=tag_str
        )
        pres = (called, targets) if isinstance(targets, Element) else (called,)
        determinants = self.signature_determinants(par)
        return VcfElement(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(vcf,),
            artifacts=ArtifactSet(out_vcf, primary_name="vcf"),
            pres=pres,
        )

    @subroutine
    def view_vcf(
        self,
        input_vcf: Path | str,
        output_vcf: Path | str,
        *,
        targets_bed: Path | str | None = None,
        pass_only: bool = False,
        biallelic: bool = False,
        samples: str | list[str] | None = None,
        variant_type: str | list[str] | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run ``bcftools view`` and write compressed VCF output.

        The output is always written with ``-Oz -o <output_vcf>``.

        Parameters
        ----------
        input_vcf : Path | str
            Input VCF/BCF file (may be gzipped).
        output_vcf : Path | str
            Path for the output VCF.gz file.
        pass_only : bool
            When *True* add ``-f PASS`` to retain only PASS-filter records.
        biallelic : bool
            When *True* add ``-m2 -M2`` to retain only biallelic records.
        samples : str | list[str] | None
            Sample(s) to extract with ``-s``; a list is joined with commas.
        variant_type : str | None
            Restrict to variant type via ``-v``; e.g. ``"snps"`` or
            ``"indels"``.
        targets_bed : Path | str | None
            BED file of target regions, passed as ``-R <bed>``.
        par : ParType | None
            Additional parameters for vcf view call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).


        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that executes the command when called.
        """
        input_vcf = Path(input_vcf).absolute()
        output_vcf = Path(output_vcf).absolute()
        ensure(output_vcf.parent)
        arguments: list[str] = ["view"]

        if pass_only:
            arguments += ["-f", "PASS"]
        if biallelic:
            arguments += ["-m2", "-M2"]
        if samples is not None:
            sample_str = samples if isinstance(samples, str) else ",".join(samples)
            arguments += ["-s", sample_str]
        if variant_type is not None:
            arguments += [
                "-v",
                (
                    variant_type
                    if isinstance(variant_type, str)
                    else ",".join(variant_type)
                ),
            ]
        if targets_bed is not None:
            arguments += ["-R", str(Path(targets_bed).absolute())]
        arguments += ["-Oz", "-o", str(output_vcf), str(input_vcf)]
        subcommand = "view"
        return arguments, subcommand, [input_vcf], [output_vcf], None, None, None

    ###########################################################################
    # Normalize (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def normalize(
        self,
        called: VcfElement,
        genome: EnsemblGenome,
        *,
        multiallelic: str = "-any",
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> VcfElement:
        """Normalize a VCF element and preserve pre-split fragment depths.

        Parameters
        ----------
        called : VcfElement
            Input VCF element.
        reference_fasta : Path | str
            Reference FASTA used to left-align and normalize alleles.
        multiallelic : str
            Value passed to bcftools norm -m. Defaults to "-any".
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to a VCF.gz in the
            input VCF folder.
        par : ParType | None
            Additional parameters for the norm call.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        VcfElement
            Element with the normalized VCF as its primary artifact.
        """
        vcf = called.vcf
        reference = Path(genome.fasta).resolve()
        resolved_tag = called.tag.bump(
            method=Method.BCFTOOLS, state=State.NORMALIZED
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=vcf.parent, ext="vcf.gz"), out
        )
        out_vcf = resolved_out.file

        runner = self.normalize_vafs(
            input_vcf=vcf,
            output_vcf=out_vcf,
            reference_fasta=reference,
            multiallelic=multiallelic,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "normalize",
            reference_fasta=reference,
            multiallelic=multiallelic,
        )
        determinants = self.signature_determinants(par)
        return VcfElement(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(vcf, reference),
            artifacts=ArtifactSet(out_vcf, primary_name="vcf"),
            pres=(called,),
        )

    @subroutine
    def normalize_vafs(
        self,
        input_vcf: Path | str,
        output_vcf: Path | str,
        *,
        reference_fasta: Path | str,
        multiallelic: str = "-any",
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Annotate fragment depths and run bcftools norm.

        Each sample's FRAG_VAF uses its own FAD denominator. TOTAL_FAD is
        the sum of complete FAD vectors across samples at the original site.

        Parameters
        ----------
        input_vcf : Path | str
            Input VCF, plain or gzipped.
        output_vcf : Path | str
            Normalized VCF path. A .gz suffix selects compressed output.
        reference_fasta : Path | str
            Reference FASTA passed to -f.
        multiallelic : str
            Value for -m. Defaults to "-any" to split all multiallelic
            records into biallelic rows.
        par : ParType | None
            Additional parameters for the norm call.
        cfg : CfgType | None
            Optional subprocess configuration.

        Returns
        -------
        Callable[[], CompletedProcess]
            Zero-argument callable that annotates and normalizes the VCF.
        """
        input_vcf = Path(input_vcf).absolute()
        output_vcf = Path(output_vcf).absolute()
        reference = Path(reference_fasta).absolute()
        ensure(output_vcf.parent)
        workspace = TemporaryDirectory(prefix="bcftools-norm-")
        annotated_vcf = Path(workspace.name) / "annotated.vcf"

        def annotate() -> None:
            # Keep the temporary directory owned by the runner across retries.
            _annotate_fragment_vafs(
                input_vcf, Path(workspace.name) / annotated_vcf.name
            )

        def cleanup() -> None:
            annotated_vcf.unlink(missing_ok=True)

        arguments = [
            "norm",
            "-f",
            str(reference),
            "-m",
            multiallelic,
            "-Oz" if output_vcf.suffix == ".gz" else "-Ov",
            "-o",
            str(output_vcf),
            str(annotated_vcf),
        ]
        return (
            arguments,
            "norm",
            [input_vcf, reference],
            [output_vcf],
            None,
            annotate,
            cleanup,
        )

    ###########################################################################
    # Filter (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def filter(
        self,
        called: VcfElement,
        *,
        include_expr: str | None = None,
        exclude_expr: str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> VcfElement:
        """
        Run ``bcftools filter`` on an element and return a :class:`VcfElement`.
        Parameters
        ----------
        called : VcfElement | Path | str
            Input VCF element or path.
        include_expr : str | None
            ``-i`` filter expression.
        exclude_expr : str | None
            ``-e`` filter expression.
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the input VCF folder,
            ext ``vcf.gz``.
        par : ParType | None
            Additional parameters for the filter call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        VcfElement
            Element with filtered vcf.
        """
        vcf = called.vcf
        resolved_tag = called.tag.bump(
            method=Method.BCFTOOLS, state=State.FILTER
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=vcf.parent, ext="vcf.gz"), out
        )
        out_vcf = resolved_out.file

        runner = self.filter_vcf(
            input_vcf=vcf,
            output_vcf=out_vcf,
            include_expr=include_expr,
            exclude_expr=exclude_expr,
            par=par,
            cfg=cfg,
        )

        key = Element.generate_key(
            resolved_tag,
            self.version_name,
            "filter",
            include_expr=include_expr,
            exclude_expr=exclude_expr,
        )
        pres = (called,)
        determinants = self.signature_determinants(par)

        return VcfElement(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(vcf,),
            artifacts=ArtifactSet(out_vcf, primary_name="vcf"),
            pres=pres,
        )

    @subroutine
    def filter_vcf(
        self,
        input_vcf: Path | str,
        output_vcf: Path | str,
        *,
        include_expr: str | None = None,
        exclude_expr: str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Run ``bcftools filter`` and write compressed VCF output.

        Exactly one of *include_expr* or *exclude_expr* should be supplied.

        Parameters
        ----------
        input_vcf : Path | str
            Input VCF/BCF file.
        output_vcf : Path | str
            Path for the output VCF.gz file.
        include_expr : str | None
            bcftools ``-i`` (include) filter expression, e.g.
            ``"FMT/DP>=10 && FMT/AD[1]>=3"``.
        exclude_expr : str | None
            bcftools ``-e`` (exclude) filter expression.
        par : ParType | None
            Additional parameters for the filter call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], None]
            Zero-argument callable.
        """
        if include_expr is None and exclude_expr is None:
            raise ValueError("Provide either include_expr or exclude_expr.")

        input_vcf = Path(input_vcf).absolute()
        output_vcf = Path(output_vcf).absolute()
        ensure(output_vcf.parent)

        arguments: list[str] = ["filter"]
        if include_expr is not None:
            arguments += ["-i", include_expr]
        else:
            arguments += ["-e", exclude_expr]
        arguments += ["-Oz", "-o", str(output_vcf), str(input_vcf)]
        subcommand = "filter"
        return arguments, subcommand, [input_vcf], [output_vcf], None, None, None

    ###########################################################################
    # Indexing (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def index(
        self,
        filtered: VcfElement,
        *,
        tag: TagType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Create a tabix index (``.tbi``) for a VCF element.

        Parameters
        ----------
        filtered : VcfElement
            Input VCF element (must point to a ``.vcf.gz`` file).
        tag : TagType | None
            Partial or full Element tag override.
        par : ParType | None
            Additional parameters for the index call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            Element with ``artifacts["tbi"]`` pointing to the ``.tbi`` file.
        """
        vcf = filtered.vcf
        pres = [filtered]
        runner = self.index_vcf(vcf_file=vcf, par=par, cfg=cfg)

        tbi = Path(str(vcf) + ".tbi")
        resolved_tag = filtered.tag.bump(
            stage=Stage.PREP, method=Method.BCFTOOLS, state=State.INDEX
        ).resolve(tag)
        determinants = self.signature_determinants(par)
        key = Element.generate_key(resolved_tag, self.version_name, "index")

        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(vcf,),
            artifacts=ArtifactSet(tbi, primary_name="tbi"),
            pres=tuple(pres),
        )

    @subroutine
    def index_vcf(
        self,
        vcf_file: Path | str,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Create a tabix index (``.tbi``) with ``bcftools index -t``.

        Parameters
        ----------
        vcf_file : Path | str
            Path to the VCF.gz file to index.
        par : ParType | None
            Additional parameters for the index call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], None]
            Zero-argument callable.
        """
        vcf_file = Path(vcf_file).absolute()
        arguments: list[str] = ["index", "-t", str(vcf_file)]

        return arguments, "index", [vcf_file], [], None, None, None

    ###########################################################################
    # Counting (High-level and Low-level wrapper)
    ###########################################################################

    @subroutine
    def count_variants_vcf(
        self,
        input_vcf: Path | str,
        variant_type: str | None = None,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Return a callable that counts variants in a VCF file.

        Runs ``bcftools view -H [-v <type>] <vcf>`` and counts output lines.

        Parameters
        ----------
        input_vcf : Path | str
            VCF file to query.
        variant_type : str | None
            Optional variant type filter passed as ``-v``; e.g. ``"snps"``
            or ``"indels"``.
        par : ParType | None
            Additional parameters for the view call.
        cfg : CfgType | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], int]
            Zero-argument callable that returns the variant count as an int.
        """
        input_vcf = Path(input_vcf).absolute()
        arguments: list[str] = ["view", "-H"]
        if variant_type is not None:
            arguments += ["-v", variant_type]
        arguments.append(str(input_vcf))
        subcommand = "view"
        return arguments, subcommand, [input_vcf], [], None, None, None

    def count_variants_post_filter(
        self,
        vcf_file: Path | str,
        output_json: Path | str,
        sample_name: str | None = None,
        loci_summary_or_mosdepth: Path | str | None = None,
        *,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Callable[[], dict[str, int]]:
        """Return a callable that counts total, SNV, and indel variants.

        The returned callable writes results to *output_json* and returns a
        dict with keys ``"total"``, ``"snps"``, ``"indels"`` (and optionally
        ``"sample"``, ``"mutational_load"``).

        Parameters
        ----------
        vcf_file : Path | str
            VCF file to count variants in.
        output_json : Path | str
            Path for the output JSON file; holds the counts after execution.
        sample_name : str | None
            Optional sample name added to the counts dict.
        loci_summary_or_mosdepth : Path | str | None
            Optional path to a JSON file with callable-bases information
            (used to compute mutational load).
        par : ParType | None
            Additional parameters forwarded to count sub-calls.
        cfg : CfgType | None
            Optional configuration for running the subprocess.

        Returns
        -------
        Callable[[], dict[str, int]]
        """

        def __count(output_json=output_json) -> dict[str, int]:
            "Here we need to combine multiple calls into one callable"
            try:
                ensure(Path(output_json).parent)
                cp_total = self.count_variants_vcf(
                    vcf_file, variant_type=None, par=par, cfg=cfg
                )()
                total = len(cp_total.stdout.splitlines())

                cp_snps = self.count_variants_vcf(
                    vcf_file, variant_type="snps", par=par, cfg=cfg
                )()
                snps = len(cp_snps.stdout.splitlines())
                cp_indels = self.count_variants_vcf(
                    vcf_file, variant_type="indels", par=par, cfg=cfg
                )()
                indels = len(cp_indels.stdout.splitlines())
                counts: dict[str, Any] = {
                    "total": total,
                    "snps": snps,
                    "indels": indels,
                }
                if sample_name:
                    counts["sample"] = sample_name
                if loci_summary_or_mosdepth is not None:
                    callable_mb_data = from_json(Path(loci_summary_or_mosdepth))
                    callable_mb = callable_mb_data.get("callable_mb", None)
                    if not callable_mb:
                        raise ValueError(
                            f"Missing 'callable_mb' in {loci_summary_or_mosdepth}"
                        )
                    callable_mb = float(callable_mb)
                    if callable_mb == 0:
                        raise ValueError(
                            f"Callable megabases is zero in {loci_summary_or_mosdepth}, cannot compute mutational load."  # noqa: E501
                        )
                    mutational_load = total / float(callable_mb)
                    counts["method"] = callable_mb_data.get("method")
                    counts["callable"] = callable_mb_data.get("callable")
                    counts["callable_mb"] = callable_mb
                    counts["mutational_load"] = mutational_load
                else:
                    raise NotImplementedError("we need a mosdepth call for that")

                with open(output_json, "w") as f:
                    json.dump(counts, f, indent=2)
                return counts
            except Exception as e:
                err = PipelineError(
                    e=e,
                    phase="RUN (count_variants_post_filter)",
                )
                raise err

        return __count

    @element
    def count(
        self,
        filtered: VcfElement,
        callable_mb: Element,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Return an element with a JSON file that holds total, SNV, and indel counts.

        Parameters
        ----------
        filtered : VcfElement
            VCF element to count variants in.
        callable_mb : Element
            Element whose ``artifacts["callable"]`` points to a loci-summary
            JSON (used to compute mutational load).
        tag : TagType | None
            Partial or full Element tag override.
        out : OutType | None
            Output specification override. Defaults to the input VCF folder,
            ext ``json``.
        par : ParType | None
            Additional parameters forwarded to count sub-calls.
        cfg : CfgType | None
            Optional configuration for running the subprocess.

        Returns
        -------
        Element
            Element with ``artifacts["json"]`` pointing to the counts file.
        """
        vcf = filtered.vcf
        resolved_tag = filtered.tag.bump(
            stage=Stage.CALL, method=Method.BCFTOOLS, state=State.COUNT
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=vcf.parent, ext="json"), out
        )
        output = resolved_out.file
        key = Element.generate_key(
            resolved_tag, self.version_name, "counts_mutational_load"
        )
        loci_summary_or_mosdepth = callable_mb.artifacts["callable"]

        runner = self.count_variants_post_filter(
            vcf_file=vcf,
            output_json=output,
            sample_name=filtered.name,
            loci_summary_or_mosdepth=loci_summary_or_mosdepth,
            par=par,
            cfg=cfg,
        )
        determinants = self.signature_determinants(par)

        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(vcf,),
            artifacts=ArtifactSet(output, primary_name="json"),
            pres=(filtered, callable_mb),
        )

    ###########################################################################
    # Fill Tags
    ###########################################################################

    @element
    def filltags(
        self,
        vcf_element: Element,
        tags_to_fill: str | None = None,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        index_off=False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """
        filltags fills missing fields in a VCF file using bcftools +fill-tags
        plugin. This is useful for filling in allele frequencies (AF) or other
        annotations that may be missing from the input VCF.

        Parameters
        ----------
        vcf_element : Element
            Element containing the VCF file to be processed. The VCF file
            should be accessible via `vcf_element.artifacts["vcf"]`.
        tags_to_fill : str
            Comma-separated list of tags to fill in the VCF file, e.g. "AF" or
            "AF,AC".
        tag : TagType | None, optional
            Partial or full Element tag override.
        out : OutType | None, optional
            Output specification override. Defaults to the input VCF folder,
            ext ``vcf``.
        index_off : bool, optional
            If True, disables indexing of the output VCF file, by default False
        par : ParType | None, optional
            Additional parameters for the bcftools +fill-tags call, by default
            None
        cfg : CfgType | None, optional
            Optional configuration for running the subprocess, by default None

        Returns
        -------
        Element
            Element containing the VCF file with filled tags.
        """
        par = par or ParType()
        if tags_to_fill is None:
            raise ValueError(
                "tags_to_fill must be provided as a comma-separated string."
            )
        input_vcf = vcf_element.vcf
        resolved_tag = vcf_element.tag.bump(
            stage=Stage.PREP, method=Method.BCFTOOLS, state=State.ANNOTATE
        ).resolve(tag)
        resolved_out = OutSpec.from_tag(resolved_tag).resolve(
            OutSpec(folder=input_vcf.parent, ext="vcf"), out
        )
        output_vcf = resolved_out.file
        artifacts = ArtifactSet(output_vcf, primary_name="vcf")
        if not index_off:
            artifacts = artifacts.with_extra(
                "idx", output_vcf.with_suffix(f"{output_vcf.suffix}.idx")
            )
        runner = self.filltags_vcf(
            input_vcf=input_vcf,
            output_vcf=output_vcf,
            tags_to_fill=tags_to_fill,
            index_off=index_off,
            par=par,
            cfg=cfg,
        )
        key = Element.generate_key(
            resolved_tag, self.version_name, "filltags", fill=tags_to_fill
        )
        determinants = self.signature_determinants(par, "+fill-tags")
        return Element(
            key,
            runner,
            tag=resolved_tag,
            determinants=determinants,
            inputs=(input_vcf,),
            artifacts=artifacts,
            pres=(vcf_element,),
        )

    @subroutine
    def filltags_vcf(
        self,
        input_vcf: Path | str,
        output_vcf: Path | str,
        tags_to_fill: str,
        *,
        index_off=False,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """
        Return a callable that runs bcftools +fill-tags on the input VCF file.

        Parameters
        ----------
        input_vcf : Path | str
            Path to the input VCF file.
        output_vcf : Path | str
            Path to the output VCF file with filled tags.
        tags_to_fill : str
            Comma-separated list of tags to fill in the VCF file, e.g. "AF" or
            "AF,AC".
        index_off : bool, optional
            If True, disables indexing of the output VCF file, by default False
        par : ParType | None, optional
            Additional parameters for the bcftools +fill-tags call, by default
            None
        cfg : CfgType | None, optional
            Optional configuration for running the subprocess, by default None

        Returns
        -------
        Callable[[], None]
            A zero-argument callable that executes the bcftools +fill-tags command when
            called.
        """
        input_vcf = Path(input_vcf).absolute()
        # bcftools +fill-tags .../C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz -- -t AF
        arguments: list[str] = [
            "+fill-tags",
            "-t",
            tags_to_fill,
            "-Oz",
            "-o",
            self.strabs(output_vcf),
            self.strabs(input_vcf),
        ]
        post = GATK().index_feature_file(output_vcf) if not index_off else None
        subcommand = "+fill-tags"
        return arguments, subcommand, [input_vcf], [output_vcf], output_vcf, None, post

    ###########################################################################
    # Convenience methods for common workflows
    ###########################################################################

    def hard_filter(
        self,
        called: VcfElement,
        *,
        targets: Element | Path | str | None = None,
        thresholds: HardFilterThresholds | None = None,
        pass_only: bool = True,
        biallelic: bool = True,
        samples: str | list[str] | None = None,
        variant_type: str | None = "snps,indels",
        tag: TagType | SpecValues[TagType] | None = None,
        out: OutType | SpecValues[OutType] | None = None,
        par: ParType | SpecValues[ParType] | None = None,
        cfg: CfgType | SpecValues[CfgType] | None = None,
    ) -> VcfElement:
        """Run the full hard-filter pipeline on a Mutect2 VCF element.

        Internally chains:

        1. PASS filter (``bcftools view -f PASS``)
        2. Sample extraction + depth/VAF filter (``bcftools filter``)

        Parameters
        ----------
        called : VcfElement
            Input VCF element (Mutect2 filtered output).
        targets : Element | Path | str | None
            Optional BED element or path; restricts the view step to target
            regions.
        thresholds : HardFilterThresholds | None
            Depth / VAF thresholds (defaults to
            ``HardFilterThresholds(min_dp=10, min_ad=3, min_vaf=0.05)``).
        pass_only : bool
            Retain only PASS-filter records in the view step.
        biallelic : bool
            Retain only biallelic variants.
        samples : str | list[str] | None
            Sample(s) to extract; defaults to ``[called.name]``.
        variant_type : str | None
            Variant type filter passed as ``-v``; e.g. ``"snps,indels"``.
        tag : TagType | SpecValues[TagType] | None
            Tag override. For per-step tags, declare
            ``valid_steps=("view", "filter")`` and put values in the
            ``values`` mapping. The filter step defaults to
            ``flag="hard"``.
        out : OutType | SpecValues[OutType] | None
            Output specification override, optionally per step (``view``,
            ``filter``). Defaults to the input VCF folder.
        par : ParType | SpecValues[ParType] | None
            Additional parameters, optionally per step (``view``, ``filter``).
        cfg : CfgType | SpecValues[CfgType] | None
            Run configuration, optionally per step (``view``, ``filter``).

        Returns
        -------
        VcfElement
            Element pointing to the fully filtered VCF, with
            ``artifacts["vcf"]`` set to the final output path.
        """
        thresholds = thresholds or HardFilterThresholds()
        samples = samples or [called.tag.root]
        view = self.view(
            called,
            targets=targets,
            pass_only=pass_only,
            biallelic=biallelic,
            samples=samples,
            variant_type=variant_type,
            tag=pick_spec(tag, "view"),
            out=pick_spec(out, "view"),
            par=pick_spec(par, "view"),
            cfg=pick_spec(cfg, "view"),
        )

        include_expr = (
            "&&".join([thresholds.as_expression(s) for s in samples])
            if len(samples) > 1
            else thresholds.as_expression()
        )
        filtered = self.filter(
            view,
            include_expr=include_expr,
            tag=pick_spec(tag, "filter") or PartialElementTag(flag="hard"),
            out=pick_spec(out, "filter"),
            par=pick_spec(par, "filter"),
            cfg=pick_spec(cfg, "filter"),
        )
        return filtered
