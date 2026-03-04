"""
Module contains a bcftools interface for VCF file processing and filtering.
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any, Callable, Mapping

from mmalignments.models.data import HardFilterThresholds
from mmalignments.models.tasks import Element, VcfElement, element
from mmalignments.services.io import ensure, from_json

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet

logger = logging.getLogger(__name__)


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
        called: VcfElement | Path | str,
        *,
        targets: Element | Path | str | None = None,
        output_vcf: Path | str | None = None,
        pass_only: bool = False,
        biallelic: bool = False,
        samples: str | list[str] | None = None,
        variant_type: str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> VcfElement:
        """Run ``bcftools view`` on an element and return a :class:`VcfElement`.

        Parameters
        ----------
        called : VcfElement | Path | str
            Input VCF element (uses ``element.artifacts["vcf"]``) or a plain
            path to a VCF file.
        targets : Element | Path | str | None
            BED file element for target-region restriction.
        output_vcf : Path | str | None
            Explicit output path; derived from *element* name when *None*.
        pass_only : bool
            Retain only PASS-filter records.
        biallelic : bool
            Retain only biallelic variants.
        samples : str | list[str] | None
            Sample(s) to extract.
        variant_type : str | None
            ``"snps"`` or ``"indels"`` etc.
        params : Params | None
            Additional parameters for vcf view call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        VcfElement
            Element with ``artifacts["vcf"]`` pointing to the output file.
        """
        vcf = called.vcf

        out_vcf = (
            output_vcf or (vcf.parent / f"{vcf.stem}.bcf_prefilter.vcf.gz").absolute()
        )

        targets_bed = targets.bed if targets is not None else None
        samples = samples or [called.name]

        runner = self.view_vcf(
            input_vcf=vcf,
            output_vcf=out_vcf,
            targets_bed=targets_bed,
            pass_only=pass_only,
            biallelic=biallelic,
            samples=samples,
            variant_type=variant_type,
            params=params,
            cfg=cfg,
        )
        tag = "view"
        tag += "_pass" if pass_only else ""
        tag += "_biallelic" if biallelic else ""

        name = f"{called.name}_bcf_{tag}_{self.version_name}"
        key = f"{name}"
        pres = [called]
        if targets is not None:
            key += f"_targets-{Path(targets_bed).stem}"
            pres.append(targets)

        determinants = self.signature_determinants(params)
        return VcfElement(
            name=name,
            key=key,
            run=runner,
            determinants=determinants,
            inputs=[vcf],
            artifacts={"vcf": out_vcf},
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
        include_variants: str | list[str] | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], None]:
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
        include_variants : str | None
            Restrict to variant type via ``-v``; e.g. ``"snps"`` or
            ``"indels"``.
        targets_bed : Path | str | None
            BED file of target regions, passed as ``-R <bed>``.
        params : Params | None
            Additional parameters for vcf view call.
        cfg : ExternalRunConfig | None
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
        if include_variants is not None:
            arguments += [
                "-v",
                (
                    include_variants
                    if isinstance(include_variants, str)
                    else ",".join(include_variants)
                ),
            ]
        if targets_bed is not None:
            arguments += ["-R", str(Path(targets_bed).absolute())]
        arguments += ["-Oz", "-o", str(output_vcf), str(input_vcf)]

        return self.runnable(
            arguments=arguments,
            params=params,
            cfg=cfg,
        )

    ###########################################################################
    # Filter (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def filter(
        self,
        called: VcfElement | Path | str,
        *,
        include_expr: str | None = None,
        exclude_expr: str | None = None,
        output_vcf: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
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
        output_vcf : Path | str | None
            Explicit output path.
        params : Params | None
            Additional parameters for the filter call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        VcfElement
            Element with filtered vcf.
        """
        vcf = called.vcf

        out_vcf = (
            output_vcf or (vcf.parent / f"{vcf.stem}.bcf.filtered.vcf.gz").absolute()
        )

        runner = self.filter_vcf(
            input_vcf=vcf,
            output_vcf=out_vcf,
            include_expr=include_expr,
            exclude_expr=exclude_expr,
            params=params,
            cfg=cfg,
        )
        name = f"{called.name}_filtered"
        key = f"{name}_filtered_{self.version_name}"
        if include_expr is not None:
            key += f"_incl-{include_expr}"
        elif exclude_expr is not None:
            key += f"_excl-{exclude_expr}"
        pres = [called]
        determinants = self.signature_determinants(params)

        return VcfElement(
            name=name,
            key=key,
            run=runner,
            determinants=determinants,
            inputs=[vcf],
            artifacts={"vcf": out_vcf},
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
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], None]:
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
        params : Params | None
            Additional parameters for the filter call.
        cfg : ExternalRunConfig | None
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

        return self.runnable(
            arguments=arguments,
            params=params,
            cfg=cfg,
        )

    ###########################################################################
    # Indexing (High-level and Low-level wrapper)
    ###########################################################################

    @element
    def index(
        self,
        filtered: VcfElement,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Create a tabix index (``.tbi``) for a VCF element.

        Parameters
        ----------
        filtered : VcfElement
            Input VCF element (must point to a ``.vcf.gz`` file).
        params : Params | None
            Additional parameters for the index call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Element
            Element with ``artifacts["tbi"]`` pointing to the ``.tbi`` file.
        """
        vcf = filtered.vcf
        pres = [filtered]
        runner = self.index_vcf(vcf_file=vcf, params=params, cfg=cfg)

        tbi = Path(str(vcf) + ".tbi")
        name = f"{filtered.name}_index"
        determinants = self.signature_determinants(params)

        return Element(
            name=name,
            key=f"{name}_index_{self.version_name}",
            run=runner,
            determinants=determinants,
            inputs=[vcf],
            artifacts={"tbi": tbi},
            pres=pres,
        )

    def index_vcf(
        self,
        vcf_file: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], None]:
        """Create a tabix index (``.tbi``) with ``bcftools index -t``.

        Parameters
        ----------
        vcf_file : Path | str
            Path to the VCF.gz file to index.
        params : Params | None
            Additional parameters for the index call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess (e.g. working
            directory, environment variables).

        Returns
        -------
        Callable[[], None]
            Zero-argument callable.
        """
        vcf_file = Path(vcf_file).absolute()
        arguments: list[str] = ["index", "-t", str(vcf_file)]

        return self.runnable(
            arguments=arguments,
            params=params,
            cfg=cfg,
        )

    ###########################################################################
    # Counting (High-level and Low-level wrapper)
    ###########################################################################

    @subroutine
    def count_variants_vcf(
        self,
        input_vcf: Path | str,
        variant_type: str | None = None,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Callable[[], int]:
        """Return a callable that counts variants in a VCF file.

        Runs ``bcftools view -H [-v <type>] <vcf>`` and counts output lines.

        Parameters
        ----------
        input_vcf : Path | str
            VCF file to query.
        variant_type : str | None
            Optional variant type filter passed as ``-v``; e.g. ``"snps"``
            or ``"indels"``.
        params : Params | None
            Additional parameters for the view call.
        cfg : ExternalRunConfig | None
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

        return self.runnable(
            arguments=arguments,
            params=params,
            cfg=cfg,
        )

    @subroutine
    def count_variants_post_filter(
        self,
        vcf_file: Path | str,
        output_json: Path | str,
        sample_name: str | None = None,
        loci_summary_or_mosdepth: Path | str | None = None,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
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
        params : Params | None
            Additional parameters forwarded to count sub-calls.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess.

        Returns
        -------
        Callable[[], dict[str, int]]
        """

        def __count(output_json=output_json) -> dict[str, int]:
            ensure(Path(output_json).parent)
            total = self.count_variants_vcf(
                vcf_file, variant_type=None, params=params, cfg=cfg
            )()
            snps = self.count_variants_vcf(
                vcf_file, variant_type="snps", params=params, cfg=cfg
            )()
            indels = self.count_variants_vcf(
                vcf_file, variant_type="indels", params=params, cfg=cfg
            )()
            counts: dict[str, Any] = {
                "total": total,
                "snps": snps,
                "indels": indels,
            }
            if sample_name:
                counts["sample"] = sample_name
            if loci_summary_or_mosdepth is not None:
                callable_mb_data = from_json(loci_summary_or_mosdepth)
                callable_mb = callable_mb_data.get("callable_mb")
                mutational_load = total / callable_mb
                counts["method"] = callable_mb_data.get("method")
                counts["callable"] = callable_mb_data.get("callable")
                counts["callable_mb"] = callable_mb
                counts["mutational_load"] = mutational_load

            with open(output_json, "w") as f:
                json.dump(counts, f, indent=2)
            return counts

        return __count

    @element
    def count(
        self,
        filtered: VcfElement,
        callable_mb: Element,
        *,
        output_json: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Return an element with a JSON file that holds total, SNV, and indel counts.

        Parameters
        ----------
        filtered : VcfElement
            VCF element to count variants in.
        callable_mb : Element
            Element whose ``artifacts["callable"]`` points to a loci-summary
            JSON (used to compute mutational load).
        output_json : Path | str | None
            Optional path for the output JSON file; defaults to
            ``<filtered.vcf.parent>/<filtered.vcf.stem>.vcounts.json``.
        params : Params | None
            Additional parameters forwarded to count sub-calls.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess.

        Returns
        -------
        Element
            Element with ``artifacts["json"]`` pointing to the counts file.
        """
        vcf = filtered.vcf
        output = Path(
            output_json or (vcf.parent / f"{vcf.stem}.vcounts.json").absolute()
        )
        name = f"{filtered.name}_mutational_load"
        key = f"{name}_counts_{self.version_name}"
        loci_summary_or_mosdepth = callable_mb.artifacts["callable"]

        runner = self.count_variants_post_filter(
            vcf_file=vcf,
            output_json=output,
            sample_name=filtered.name,
            loci_summary_or_mosdepth=loci_summary_or_mosdepth,
            params=params,
            cfg=cfg,
        )
        determinants = self.signature_determinants(params)

        return Element(
            name=name,
            key=key,
            run=runner,
            determinants=determinants,
            inputs=[vcf],
            artifacts={"json": output},
            pres=[filtered, callable_mb],
        )

    ###########################################################################
    # Convenience methods for common workflows
    ###########################################################################

    def hard_filter(
        self,
        called: VcfElement,
        targets: Element | Path | str | None = None,
        *,
        thresholds: HardFilterThresholds | None = None,
        output_vcf: Path | str | None = None,
        pass_only: bool = True,
        biallelic: bool = True,
        samples: str | list[str] | None = None,
        variant_type: str | None = "snps,indels",
        params_view: Params | None = None,
        params_filter: Params | None = None,
        cfg: ExternalRunConfig | None = None,
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
        output_vcf : Path | str | None
            Explicit output path for the final filtered VCF.
        pass_only : bool
            Retain only PASS-filter records in the view step.
        biallelic : bool
            Retain only biallelic variants.
        samples : str | list[str] | None
            Sample(s) to extract; defaults to ``[called.name]``.
        variant_type : str | None
            Variant type filter passed as ``-v``; e.g. ``"snps,indels"``.
        params_view : Params | None
            Additional parameters for the ``bcftools view`` call.
        params_filter : Params | None
            Additional parameters for the ``bcftools filter`` call.
        cfg : ExternalRunConfig | None
            Optional configuration for running the subprocess.

        Returns
        -------
        VcfElement
            Element pointing to the fully filtered VCF, with
            ``artifacts["vcf"]`` set to the final output path.
        """
        vcf = called.vcf
        thresholds = thresholds or HardFilterThresholds()
        output_vcf = (
            output_vcf or (vcf.parent / f"{vcf.stem}.hard_filtered.vcf.gz").absolute()
        )

        samples = samples or [called.name]
        view = self.view(
            called,
            pass_only=pass_only,
            biallelic=biallelic,
            samples=samples,
            variant_type=variant_type,
            targets=targets,
            params=params_view,
            cfg=cfg,
        )

        filtered = self.filter(
            view,
            include_expr=thresholds.as_expression(),
            output_vcf=output_vcf,
            params=params_filter,
            cfg=cfg,
        )
        return filtered
