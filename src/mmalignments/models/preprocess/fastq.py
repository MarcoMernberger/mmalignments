"""Module contains an interface for fastqgrab (fastqrab), a FASTQ processing tool."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping

import toml

from mmalignments.models.elements import (
    Element,
    NextGenSampleElement,
    element,
    sample_fastqs,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    Omics,
    Stage,
    State,
    PartialElementTag,
    from_prior,
    merge_tag,
)
from mmalignments.services.io import parents

from ..externals import External, ExternalRunConfig, subroutine
from ..parameters import Params, ParamSet

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Barcode type aliases
# ---------------------------------------------------------------------------

# Mapping from barcode_sequence (or "barcode1_barcode2") to sample name
BarcodeMap = Mapping[str, str]


class FastGrab(External):
    """Interface for fastqgrab (fastqrab), a FASTQ processing tool.

    fastqgrab processes FASTQ files with filtering, sampling, slicing,
    demultiplexing, and analysis via a TOML configuration file.

    The typical workflow is:
    1. Create a TOML configuration file (``configure``)
    2. Run fastqgrab on the configuration file (``process``)

    Examples
    --------
    Create a demultiplexing pipeline::

        fg = FastGrab()
        config_element = fg.configure(
            sample=my_sample,
            single_barcodes={"CTGGCA": "sample1", "GGTCGA": "sample2"},
            sample_barcodes={"CTGGCA_CTGGCA": "sample1", "GGTCGA_GGTCGA": "sample2"},
            output_prefix="results/demultiplex/run1",
        )
        process_element = fg.process(config=config_element)
    """

    def __init__(
        self,
        name: str = "fastqgrab",
        primary_binary: str = "fastqrab",
        version: str | None = None,
        folder: Path | None = None,
        source: str | None = None,
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        """Initialise the FastGrab wrapper.

        Parameters
        ----------
        name : str
            Logical tool name, default ``"fastqgrab"``.
        primary_binary : str
            Executable name on ``$PATH``, default ``"fastqrab"``.
        version : str | None
            Optional version string override.
        folder : Path | None
            Optional base output folder.
        source : str | None
            Documentation / source URL.
        parameters : Mapping[str, ParamSet] | ParamSet | None
            Parameter sets for invocations.
        """
        if source is None:
            source = "https://tyberiusprime.github.io/fastqrab/"
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            folder=folder,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        """Return the fastqrab version string.

        Parameters
        ----------
        fallback : str | None
            Value returned when the version cannot be determined.

        Returns
        -------
        str | None
            Version string (e.g. ``"0.3.1"``) or *fallback*.
        """
        import subprocess

        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "version"],
                capture_output=True,
                text=True,
                check=True,
            )
            output = (cp.stdout or cp.stderr).strip()
            if output:
                # last whitespace-separated token is usually the version number
                return output.split()[-1]
        except subprocess.CalledProcessError:
            pass
        return fallback

    # -----------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------

    def default_config_dir(self, sample_name: str) -> Path:
        """Return the default directory for configuration files."""
        return (
            Path("results")
            / "demultiplex"
            / self.version_name
            / sample_name
        )

    def default_output_dir(self, sample_name: str) -> Path:
        """Return the default directory for demultiplexed outputs."""
        return (
            Path("results")
            / "demultiplex"
            / self.version_name
            / sample_name
            / "output"
        )

    # -----------------------------------------------------------------------
    # Configuration (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    @element
    def configure(
        self,
        sample: NextGenSampleElement,
        *,
        sample_barcodes: BarcodeMap,
        barcode_search: list[str] | None = None,
        output_prefix: str | Path | None = None,
        output_dir: Path | None = None,
        config_path: Path | None = None,
        max_mismatches: int = 1,
        max_hamming_distance: int = 1,
        tag: PartialElementTag | ElementTag | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Create an Element that writes a fastqgrab TOML configuration file.

        The generated TOML follows the demultiplex.toml template: it extracts
        barcodes from both reads, applies Hamming correction, concatenates the
        results, demultiplexes by combined barcode, trims at the barcodes and
        stores tags in a TSV table.

        Parameters
        ----------
        sample : NextGenSampleElement
            Input sample with paired FASTQ files.
        sample_barcodes : BarcodeMap
            Mapping of single barcode sequence → sample name
            (e.g. ``{"CTGGCA": "sample1"}``).
        barcode_search : list[str] | None
            Explicit list of barcode sequences to search for.  When *None*,
            the keys of *sample_barcodes* are used.
        output_prefix : str | Path | None
            Prefix for fastqgrab output files (passed to the ``[output]``
            section of the TOML).  When *None*, a default path is derived
            from *output_dir* and the sample name.
        output_dir : Path | None
            Directory in which to write the configuration file and outputs.
        config_path : Path | None
            Explicit path for the written TOML file.  When *None*, a default
            path is derived from *output_dir*.
        max_mismatches : int
            Maximum mismatches allowed for the IUPAC barcode search step.
        max_hamming_distance : int
            Maximum Hamming distance for barcode correction.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        params : Params | None
            Unused; kept for API symmetry.
        cfg : ExternalRunConfig | None
            Unused; kept for API symmetry.

        Returns
        -------
        Element
            Element whose artifact ``"toml"`` is the path to the written TOML
            configuration file.
        """
        fastq_r1, fastq_r2, sample_name, _ = sample_fastqs(sample)

        output_dir = output_dir or self.default_output_dir(sample_name)
        config_dir = output_dir.parent
        config_path = config_path or (config_dir / f"{sample_name}_demultiplex.toml")
        output_prefix = output_prefix or str(output_dir / sample_name)

        # Derive sample_barcodes from single_barcodes if not provided
        if sample_barcodes is None:
            sample_barcodes = {
                f"{bc}_{bc}": sn for bc, sn in single_barcodes.items()
            }

        barcode_search = barcode_search or list(sample_barcodes.keys())

        runner = self.write_config(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            config_path=config_path,
            sample_barcodes=sample_barcodes,
            barcode_search=barcode_search,
            output_prefix=str(output_prefix),
            max_mismatches=max_mismatches,
            max_hamming_distance=max_hamming_distance,
            params=params,
            cfg=cfg,
        )

        default_tag = from_prior(
            sample.tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.DEMULTIPLEX,
            omics=Omics.DNA,
            ext="toml",
        )
        tag = merge_tag(default_tag, tag)

        key = f"{tag.default_name}_config_{self.version_name}"
        artifacts = {"toml": config_path}
        determinants = (
            str(sorted(sample_barcodes.items())),
            str(max_mismatches),
            str(max_hamming_distance),
        )

        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=[fastq_r1] + ([fastq_r2] if fastq_r2 else []),
            pres=[sample],
            name=None,
        )

    @subroutine
    def write_config(
        self,
        fastq_r1: Path | str,
        fastq_r2: Path | str | None,
        config_path: Path | str,
        sample_barcodes: BarcodeMap,
        barcode_search: list[str],
        output_prefix: str,
        *,
        max_mismatches: int = 1,
        max_hamming_distance: int = 1,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple:
        """Write a fastqgrab TOML configuration file and return a subroutine tuple.

        This low-level method constructs the TOML configuration document
        from the supplied parameters and writes it to *config_path*.  Because
        the "command" is purely a Python file-write operation (not a
        subprocess call), the returned tuple uses a Python-callable runner
        rather than a CLI argument list.

        Parameters
        ----------
        fastq_r1 : Path | str
            Absolute path to the R1 FASTQ file.
        fastq_r2 : Path | str | None
            Absolute path to the R2 FASTQ file, or *None* for single-end.
        config_path : Path | str
            Destination path for the written TOML file.
        sample_barcodes : BarcodeMap
            Combined barcode pair → sample mapping for demultiplexing.
        barcode_search : list[str]
            Barcode sequences to search for in each read.
        output_prefix : str
            Prefix for fastqgrab output files.
        max_mismatches : int
            Mismatches allowed in ExtractIUPAC steps.
        max_hamming_distance : int
            Maximum Hamming distance for HammingCorrect steps.
        params : Params | None
            Unused; kept for API symmetry.
        cfg : ExternalRunConfig | None
            Unused; kept for API symmetry.

        Returns
        -------
        tuple
            ``(arguments, paths, output, pre, post)`` consumed by the
            ``@subroutine`` decorator.
        """
        config_path = Path(config_path).absolute()
        fastq_r1 = Path(fastq_r1).absolute()

        toml_doc = _build_demultiplex_toml(
            fastq_r1=fastq_r1,
            fastq_r2=Path(fastq_r2).absolute() if fastq_r2 else None,
            sample_barcodes=dict(sample_barcodes),
            barcode_search=list(barcode_search),
            output_prefix=output_prefix,
            max_mismatches=max_mismatches,
            max_hamming_distance=max_hamming_distance,
        )

        parents(config_path)

        def _write() -> None:
            config_path.write_bytes(toml.dumps(toml_doc).encode())

        # pre = write the file before the (no-op) command
        # We encode the file-write as the pre-hook and pass an empty argument
        # list so the subroutine machinery does not try to call a binary.
        # The actual "command" is represented as a python:// pseudo-URI so
        # that runners can log it meaningfully.
        return (
            ["python://write_toml", str(config_path)],
            [config_path],
            None,
            _write,
            None,
        )

    # -----------------------------------------------------------------------
    # Process (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    @element
    def process(
        self,
        config: Element,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Create an Element that runs ``fastqrab process <config.toml>``.

        Parameters
        ----------
        config : Element
            The configuration Element produced by :meth:`configure`.  Its
            ``"toml"`` artifact is used as the configuration file path.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        params : Params | None
            Additional CLI parameters for the ``process`` sub-command.
        cfg : ExternalRunConfig | None
            Optional subprocess configuration (working directory, threads, …).

        Returns
        -------
        Element
            Element whose artifact ``"config"`` points to the TOML file that
            was processed.  Downstream Elements can declare this as a
            prerequisite.
        """
        toml_path = Path(config.artifacts["toml"]).absolute()
        cfg = cfg or ExternalRunConfig(cwd=toml_path.parent)

        runner = self.run_process(
            config_path=toml_path,
            params=params,
            cfg=cfg,
        )

        tag = from_prior(
            config.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.DEMULTIPLEX,
            ext=None,
        )

        key = f"{tag.default_name}_process_{self.version_name}_{toml_path}"
        determinants = self.signature_determinants(params, subroutine="process")

        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"config": toml_path},
            determinants=determinants,
            inputs=[toml_path],
            pres=[config],
            empty_ok=True,
            name=None,
        )

    @subroutine
    def run_process(
        self,
        config_path: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple:
        """Low-level wrapper for ``fastqrab process <config.toml>``.

        Parameters
        ----------
        config_path : Path | str
            Path to the TOML configuration file.
        params : Params | None
            Extra CLI flags for the ``process`` sub-command.
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        tuple
            ``(arguments, paths, output, pre, post)`` consumed by the
            ``@subroutine`` decorator.
        """
        config_path = Path(config_path).absolute()
        return (
            ["process", str(config_path)],
            [config_path],
            None,
            None,
            None,
        )

    # -----------------------------------------------------------------------
    # Validate (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    @element
    def validate(
        self,
        config: Element,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Create an Element that runs ``fastqrab validate <config.toml>``.

        Validation checks the TOML configuration without producing any
        processed output.  This is useful to verify a configuration before
        submitting expensive processing jobs.

        Parameters
        ----------
        config : Element
            The configuration Element produced by :meth:`configure`.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        params : Params | None
            Additional CLI parameters.
        cfg : ExternalRunConfig | None
            Optional subprocess configuration.

        Returns
        -------
        Element
            Element representing the validation step; its only artifact is
            the path to the validated TOML.
        """
        toml_path = Path(config.artifacts["toml"]).absolute()
        cfg = cfg or ExternalRunConfig(cwd=toml_path.parent)

        runner = self.run_validate(
            config_path=toml_path,
            params=params,
            cfg=cfg,
        )

        default_tag = from_prior(
            config.tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.DEMULTIPLEX,
            ext=None,
        )
        tag = merge_tag(default_tag, tag)

        key = f"{tag.default_name}_validate_{self.version_name}_{toml_path}"
        determinants = self.signature_determinants(params, subroutine="validate")

        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"config": toml_path},
            determinants=determinants,
            inputs=[toml_path],
            pres=[config],
            empty_ok=True,
            name=None,
        )

    @subroutine
    def run_validate(
        self,
        config_path: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple:
        """Low-level wrapper for ``fastqrab validate <config.toml>``.

        Parameters
        ----------
        config_path : Path | str
            Path to the TOML configuration file.
        params : Params | None
            Extra CLI flags.
        cfg : ExternalRunConfig | None
            Subprocess configuration.

        Returns
        -------
        tuple
            ``(arguments, paths, output, pre, post)`` consumed by the
            ``@subroutine`` decorator.
        """
        config_path = Path(config_path).absolute()
        return (
            ["validate", str(config_path)],
            [config_path],
            None,
            None,
            None,
        )

    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------

    def configure_and_process(
        self,
        sample: NextGenSampleElement,
        *,
        sample_barcodes: BarcodeMap,
        barcode_search: list[str] | None = None,
        output_prefix: str | Path | None = None,
        output_dir: Path | None = None,
        max_mismatches: int = 1,
        max_hamming_distance: int = 1,
        params: Mapping[str, Params] | None = None,
        cfgs: Mapping[str, ExternalRunConfig] | None = None,
    ) -> tuple[Element, Element]:
        """Configure and process a demultiplexing run in one step.

        This convenience method calls :meth:`configure` followed by
        :meth:`process` and returns both Elements.

        Parameters
        ----------
        sample : NextGenSampleElement
            Input sample with paired FASTQ files.
        sample_barcodes : BarcodeMap
            Barcode → sample mapping.
        barcode_search : list[str] | None
            Explicit barcode sequences to search for.
        output_prefix : str | Path | None
            Prefix for fastqgrab output files.
        output_dir : Path | None
            Base output directory.
        max_mismatches : int
            Mismatches for IUPAC extraction.
        max_hamming_distance : int
            Hamming distance for correction.
        params : Mapping[str, Params] | None
            Parameter sets keyed by step name (``"configure"``,
            ``"process"``).
        cfgs : Mapping[str, ExternalRunConfig] | None
            Run configurations keyed by step name.

        Returns
        -------
        tuple[Element, Element]
            ``(config_element, process_element)``
        """
        params = params or {}
        cfgs = cfgs or {}

        config_el = self.configure(
            sample=sample,
            sample_barcodes=sample_barcodes,
            barcode_search=barcode_search,
            output_prefix=output_prefix,
            output_dir=output_dir,
            max_mismatches=max_mismatches,
            max_hamming_distance=max_hamming_distance,
            params=params.get("configure"),
            cfg=cfgs.get("configure"),
        )
        process_el = self.process(
            config=config_el,
            params=params.get("process"),
            cfg=cfgs.get("process"),
        )
        return config_el, process_el


# ---------------------------------------------------------------------------
# TOML builder (private helper)
# ---------------------------------------------------------------------------


def _build_demultiplex_toml(
    *,
    fastq_r1: Path,
    fastq_r2: Path | None,
    sample_barcodes: dict[str, str],
    barcode_search: list[str],
    output_prefix: str,
    max_mismatches: int = 1,
    max_hamming_distance: int = 1,
) -> dict:
    """Build a fastqgrab TOML configuration document as a plain dict.

    The structure mirrors the ``demultiplex.toml`` template:

    - ``[input]``  – R1 (and optionally R2) paths
    - ``[barcodes]`` – single and combined barcode → sample maps
    - ``[[step]]`` list – Report, ExtractIUPAC (×2), HammingCorrect (×2),
      StoreTagInComment (×2), ConcatTags, Demultiplex, TrimAtTag (×2),
      StoreTagsInTable
    - ``[output]``

    Parameters
    ----------
    fastq_r1, fastq_r2 : Path | None
        Input FASTQ paths.
    sample_barcodes : dict[str, str]
        Combined barcode → sample name for demultiplexing.
    barcode_search : list[str]
        Barcode sequences to search for in each read.
    output_prefix : str
        Prefix for all fastqgrab output files.
    max_mismatches : int
        Mismatches allowed during barcode extraction.
    max_hamming_distance : int
        Maximum Hamming distance for barcode correction.

    Returns
    -------
    dict
        A nested dict ready to be serialised with ``toml.dumps``.
    """
    doc: dict = {}

    # [input]
    input_section: dict = {"read1": str(fastq_r1)}
    if fastq_r2 is not None:
        input_section["read2"] = str(fastq_r2)
    doc["input"] = input_section

    # [barcodes]
    doc["barcodes"] = {
        "sample_barcodes": sample_barcodes,
    }

    # [[step]] list
    steps: list[dict] = []

    # Report before
    steps.append(
        {
            "action": "Report",
            "name": "report_before",
            "count": True,
            "base_statistics": True,
            "length_distribution": True,
        }
    )

    # ExtractIUPAC for R1 and R2
    for segment, label in [("read1", "barcode_in_R1"), ("read2", "barcode_in_R2")]:
        steps.append(
            {
                "action": "ExtractIUPAC",
                "out_label": label,
                "anchor": "Anywhere",
                "search": barcode_search,
                "segment": segment,
                "max_mismatches": max_mismatches,
                "on_tie": "leftmost",
            }
        )

    # HammingCorrect for R1 and R2
    for in_label, out_label in [
        ("barcode_in_R1", "barcode_in_R1_corrected"),
        ("barcode_in_R2", "barcode_in_R2_corrected"),
    ]:
        steps.append(
            {
                "action": "HammingCorrect",
                "in_label": in_label,
                "out_label": out_label,
                "barcodes": "sample_barcodes",
                "max_hamming_distance": max_hamming_distance,
                "on_no_match": "keep",
                "on_tie": "first",
            }
        )

    # StoreTagInComment for raw barcode tags
    for label in ["barcode_in_R1", "barcode_in_R2"]:
        steps.append({"action": "StoreTagInComment", "in_label": label})

    # ConcatTags
    steps.append(
        {
            "action": "ConcatTags",
            "in_labels": ["barcode_in_R1_corrected", "barcode_in_R2_corrected"],
            "out_label": "combined_barcodes",
            "on_missing": "merge_present",
            "separator": "_",
        }
    )

    # Demultiplex
    steps.append(
        {
            "action": "Demultiplex",
            "in_label": "combined_barcodes",
            "barcodes": "sample_barcodes",
            "output_unmatched": True,
        }
    )

    # TrimAtTag for corrected barcodes in R1 and R2
    for label in ["barcode_in_R1_corrected", "barcode_in_R2_corrected"]:
        steps.append(
            {
                "action": "TrimAtTag",
                "in_label": label,
                "direction": "Start",
                "keep_tag": False,
            }
        )

    # StoreTagsInTable
    steps.append(
        {
            "action": "StoreTagsInTable",
            "infix": "tags",
            "compression": "Raw",
            "region_separator": "_",
            "include_read_name": True,
        }
    )

    doc["step"] = steps

    # [output]
    doc["output"] = {
        "prefix": output_prefix,
        "format": "FASTQ",
        "report_json": False,
        "report_html": True,
    }

    return doc
