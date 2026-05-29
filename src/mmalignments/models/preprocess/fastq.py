"""Module contains an interface for fastqgrab (fastqrab), a FASTQ processing tool."""

from __future__ import annotations

import logging
import re
import gzip
from pathlib import Path
from typing import Mapping

from tomlkit import TOMLDocument, document, dumps, loads, table, nl

from mmalignments.models.elements import (
    Element,
    FileElement,
    NextGenSampleElement,
    element,
    sample_fastqs,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    Omics,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.io import parents
from mmalignments.services.logging import current_call_to_string
from mmalignments.services.toml import indent_regular_tables
from mmalignments.services.dependencies import depends
from ..externals import External, ExternalRunConfig, Runnable, subroutine, SubroutineIn
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
                return output.split()[0]
        except subprocess.CalledProcessError:
            pass
        return fallback

    # -----------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------

    def default_config_dir(self, sample_name: str) -> Path:
        """Return the default directory for configuration files."""
        return (
            Path("cache") / "demultiplex" / self.version_name / "config" / sample_name
        )

    def default_output_dir(self, sample_name: str) -> Path:
        """Return the default directory for demultiplexed outputs."""
        return (
            Path("cache") / "demultiplex" / self.version_name / sample_name
        )

    # -----------------------------------------------------------------------
    # Configuration (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    @element
    def configure(
        self,
        sample: NextGenSampleElement,
        template: FileElement,
        sample_barcodes: Mapping[str, str],
        *,
        primers: tuple[str, str] | None = None,
        demultiplex_dir: Path | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
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
        params = params or Params()
        fastq_r1, fastq_r2, _, _ = sample_fastqs(sample)
        tag = from_prior(
            sample.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.CONFIG,
            omics=Omics.DNA,
            ext="toml",
        )

        output_dir = Path(outdir or self.default_config_dir(sample.root))
        config_file = filename or tag.default_output
        config_path = output_dir / config_file
        demultiplex_destination = demultiplex_dir or self.default_output_dir(
            sample.root
        )
        demultiplex_destination = demultiplex_destination / sample.root
        forward_primer, reverse_primer = primers or (None, None)
        runner = self.write_config(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            config_path=config_path,
            sample_barcodes=sample_barcodes,
            template_path=template.toml,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            barcode_hamming_distance=params.get("barcode_hamming_distance", 1),
            primer_hamming_distance=params.get("primer_hamming_distance", 3),
            barcode_anchor_distance=params.get("barcode_anchor_distance", 12),
            output_prefix=str(demultiplex_destination),
        )
        key, name = self.build_element_name(tag, "configure")
        artifacts = {"toml": config_path, "prefix": str(demultiplex_destination)}
        determinants = tuple(
            [str(sorted(sample_barcodes.items())),
            str(params.get("barcode_hamming_distance", 1)),
            str(params.get("primer_hamming_distance", 3)),
            str(params.get("barcode_anchor_distance", 1)),
            forward_primer,
            reverse_primer,
            ] + list(sample_barcodes.keys())
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,),
            pres=(sample, template),
            name=name,
        )

    def write_config(
        self,
        fastq_r1: Path | str,
        fastq_r2: Path | str | None,
        config_path: Path | str,
        sample_barcodes: BarcodeMap,
        template_path: Path | str,
        *,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        barcode_hamming_distance: int = 1,
        primer_hamming_distance: int = 3,
        barcode_anchor_distance: int = 12,
        output_prefix: str = "output",
    ) -> Runnable:
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
            Mapping of sample/replicate name → barcode sequence.
        template_path : Path | str
            Path to the TOML template file with ``${variable}`` placeholders.
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
            template_path=Path(template_path).absolute(),
            fastq_r1=fastq_r1,
            fastq_r2=Path(fastq_r2).absolute() if fastq_r2 else None,
            sample_barcodes=dict(sample_barcodes),
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            barcode_hamming_distance=barcode_hamming_distance,
            primer_hamming_distance=primer_hamming_distance,
            barcode_anchor_distance=barcode_anchor_distance,
            output_prefix=output_prefix,
        )

        parents(config_path)

        @depends(_build_demultiplex_toml)
        def _write():
            text = indent_regular_tables(dumps(toml_doc))
            config_path.write_bytes(text.encode())

        return Runnable(
            _write,
            cmd=[current_call_to_string()],
            display="write_config",
        )

        # pre = write the file before the (no-op) command
        # We encode the file-write as the pre-hook and pass an empty argument
        # list so the subroutine machinery does not try to call a binary.
        # The actual "command" is represented as a python:// pseudo-URI so
        # that runners can log it meaningfully.
        # return (
        #     ["python://write_toml", str(config_path)],
        #     [config_path],
        #     None,
        #     _write,
        #     None,
        # )

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

        Process will run whatever is defined in the configuration. As main
        artifact the path to the report html and json files are assumed.
        Although, with demultiplexing, this will also create a bunch of FASTQ 
        files as intermediate outputs, which are not declared as artifacts here 
        since they are not fixed outputs of the process 
        (e.g. the number and names depend on the barcodes in the config).
        
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
        runner = self.process_config(
            config_path=config.toml
        )

        tag = from_prior(
            config.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.REPORT,
            ext=None,
        )
        key, name = self.build_element_name(tag, "process")

        determinants = self.signature_determinants(params, subroutine="process")

        return Element(
            key,
            runner,
            tag=tag,
            artifacts={"html": Path(config.prefix + "_report.html"), "json": Path(config.prefix + "_report.json")},
            determinants=determinants,
            inputs=(config.toml,),
            pres=(config,),
            empty_ok=True,
            name=name,
        )

    @subroutine
    def process_config(
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
        pre = self.run_validate(
            config_path=config_path
        )
        return (
            ["process", str(config_path), "--allow-overwrite"],
            "process",
            [config_path],
            [],
            None,
            pre,
            None
        )

    # -----------------------------------------------------------------------
    # Validate (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    # @element
    # def validate(
    #     self,
    #     config: Element,
    #     *,
    #     tag: PartialElementTag | ElementTag | None = None,
    # ) -> Element:
    #     """Create an Element that runs ``fastqrab validate <config.toml>``.

    #     Validation checks the TOML configuration without producing any
    #     processed output.  This is useful to verify a configuration before
    #     submitting expensive processing jobs.

    #     Parameters
    #     ----------
    #     config : Element
    #         The configuration Element produced by :meth:`configure`.
    #     tag : PartialElementTag | ElementTag | None
    #         Optional tag override.
    #     params : Params | None
    #         Additional CLI parameters.
    #     cfg : ExternalRunConfig | None
    #         Optional subprocess configuration.

    #     Returns
    #     -------
    #     Element
    #         Element representing the validation step; its only artifact is
    #         the path to the validated TOML.
    #     """
    #     toml_path = Path(config.toml).absolute()

    #     runner = self.run_validate(
    #         config_path=toml_path,
    #     )

    #     tag = from_prior(
    #         config.tag,
    #         tag,
    #         stage=Stage.PREP,
    #         method=Method.FASTQGRAB,
    #         state=State.DEMULTIPLEX,
    #         ext=None,
    #     )

    #     key, name = self.build_element_name(tag, "validate")

    #     return Element(
    #         key,
    #         runner,
    #         tag=tag,
    #         artifacts={"config": toml_path},
    #         inputs=[toml_path],
    #         pres=[config],
    #         empty_ok=True,
    #         name=name,
    #     )

    @subroutine
    def run_validate(
        self,
        config_path: Path | str,
        *,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for ``fastqrab validate <config.toml>``.

        Parameters
        ----------
        config_path : Path | str
            Path to the TOML configuration file.

        Returns
        -------
        tuple
            ``(arguments, subcommand, inpaths, outpaths, pipeoutput, pre, post)``
            consumed by the ``@subroutine`` decorator.
        """
        config_path = Path(config_path).absolute()
        return (
            ["validate", str(config_path)],
            "validate",
            [config_path],
            [],
            None,
            None,
            None
        )

    # -----------------------------------------------------------------------
    # Consolidate (high-level @element + low-level subroutine)
    # -----------------------------------------------------------------------

    @element
    def consolidate(
        self,
        process: Element,
        config: Element,
        sample_barcodes: Mapping[str, str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
    ) -> Element:
        """Create an Element that consolidates demultiplexed FASTQ files.
        This is done after demultiplexing to combine undetermined fastqs into 
        a single file and rename files in a consistent way. 

        1. Merges all ``barcode_unambiguous=false`` and ``nobarcode`` FASTQs
           into ``undetermined_R1.fq`` and ``undetermined_R2.fq``, adding
           reason tags to read names.
        2. Renames ``barcode_unambiguous=true`` FASTQs from
           ``{prefix}_barcode_unambiguous=true_{sample}_read{N}.fq``
           to ``{prefix}_{sample}_R{N}.fq``.
        3. Removes the original intermediate files.

        Parameters
        ----------
        process : Element
            The process Element that produced the demultiplexed FASTQs.
        sample_barcodes : Mapping[str, str]
            Barcode → sample mapping (same as used in configure).
        tag : PartialElementTag | ElementTag | None
            Optional tag override.

        Returns
        -------
        Element
            Element with artifacts for each consolidated sample FASTQ and
            the undetermined FASTQs.
        """
        prefix = Path(config.prefix)
        sample_names = list(sample_barcodes.keys())


        tag = from_prior(
            process.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.REPORT,
            ext="fq",
        )

        key, name = self.build_element_name(tag, "consolidate")
        # Build artifacts: one per sample (R1 + R2) plus undetermined
        artifacts = {}
        for sample_name in sample_names:
            artifacts[f"{sample_name}_R1"] = prefix.parent / f"{prefix.name}_{sample_name}_R1.fq"
            artifacts[f"{sample_name}_R2"] = prefix.parent / f"{prefix.name}_{sample_name}_R2.fq"
        artifacts["undetermined_R1"] = prefix.parent / f"{prefix.name}_undetermined_R1.fq"
        artifacts["undetermined_R2"] = prefix.parent / f"{prefix.name}_undetermined_R2.fq"

        # Input files (all the process-generated FASTQs)
        input_pattern = str(prefix / f"{prefix.name}_barcode_unambiguous=*")
        
        runner = self.run_consolidate(
            prefix=prefix,
            sample_names=sample_names,
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=tuple(sorted(sample_names)),
            inputs=(process.html, process.json),  # symbolic dependency
            pres=(process,),
            name=name,
        )

    def run_consolidate(
        self,
        prefix: Path,
        sample_names: list[str],
    ) -> Runnable:
        """Consolidate demultiplexed FASTQ files.

        Parameters
        ----------
        prefix : Path
            Output prefix used by fastqgrab (directory + base name).
        sample_names : list[str]
            List of sample names from the barcode mapping.

        Returns
        -------
        Runnable
            Runnable that performs the file consolidation.
        """

        prefix = Path(prefix).absolute()
        prefix_dir = prefix.parent
        prefix_name = prefix.name

        def _consolidate():
            """Consolidate FASTQ files after demultiplexing."""
            # Collect files
            all_fqs = list(prefix_dir.glob(f"{prefix_name}_barcode_unambiguous=*"))
            
            undetermined_r1 = prefix_dir / f"{prefix_name}_undetermined_R1.fq"
            undetermined_r2 = prefix_dir / f"{prefix_name}_undetermined_R2.fq"
            consolidated_log = prefix_dir / f"{prefix_name}_consolidation.log"
            with open(consolidated_log, 'w') as log_fh:
                log_fh.write(f"Consolidating FASTQ files for prefix {prefix_name}\n")
                log_fh.write(f"Found {len(all_fqs)} FASTQ files:\n")
                for fq in all_fqs:
                    log_fh.write(f"  {fq.name}\n")
                # Process undetermined (false + nobarcode)
                with open(undetermined_r1, 'w') as out_r1, open(undetermined_r2, 'w') as out_r2:
                    for fq_path in all_fqs:
                        fq_name = fq_path.name
                        
                        # Skip if unambiguous=true (these get renamed, not merged)
                        if "barcode_unambiguous=true" in fq_name:
                            continue
                        
                        # Determine reason tag
                        if "nobarcode" in fq_name:
                            reason = "reason:nobarcode"
                        elif "barcode_unambiguous=false" in fq_name:
                            reason = "reason:ambiguous_barcode"
                        else:
                            continue
                        
                        # Determine read direction
                        if "_read1.fq" in fq_name:
                            log_fh.write(f"Adding {fq_name} to undetermined_R1.fq with tag {reason}\n")
                            out_fh = out_r1
                        elif "_read2.fq" in fq_name:
                            log_fh.write(f"Adding {fq_name} to undetermined_R2.fq with tag {reason}\n")
                            out_fh = out_r2
                        else:
                            continue
                        
                        # Copy reads with modified names
                        with open(fq_path, 'r') as in_fh:
                            lines = in_fh.readlines()
                            for i in range(0, len(lines), 4):
                                if i + 3 < len(lines):
                                    # Modify read name (first line)
                                    header = lines[i].strip()
                                    if header.startswith('@'):
                                        header = f"{header} {reason}\n"
                                    else:
                                        header = lines[i]
                                    out_fh.write(header)
                                    out_fh.write(lines[i + 1])  # sequence
                                    out_fh.write(lines[i + 2])  # plus
                                    out_fh.write(lines[i + 3])  # quality

                # Rename unambiguous=true files
                for sample_name in sample_names:
                    for read_num, read_suffix in [(1, "read1"), (2, "read2")]:
                        old_name = prefix_dir / f"{prefix_name}_barcode_unambiguous=true_{sample_name}_{read_suffix}.fq"
                        new_name = prefix_dir / f"{prefix_name}_{sample_name}_R{read_num}.fq"
                        
                        log_fh.write(f"Renaming {old_name} to {new_name}\n")
                        if old_name.exists():
                            old_name.rename(new_name)

                # Clean up intermediate files
                for fq_path in all_fqs:
                    if fq_path.exists():
                        log_fh.write(f"Removing intermediate file {fq_path}\n")
                        fq_path.unlink()
                
        return Runnable(
            _consolidate,
            cmd=[current_call_to_string()],
            display="consolidate",
        )

    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------

    def grabfast(
        self,
        sample: NextGenSampleElement,
        template: FileElement,
        sample_barcodes: Mapping[str, str],
        *,
        primers: tuple[str, str] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        configdir: Path | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[Element, Element, Element]:
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
        tuple[Element, Element, Element]
            ``(consolidate_element, process_element, config_element)``
        """
        config_el = self.configure(
            sample=sample,
            template=template,
            sample_barcodes=sample_barcodes,
            primers=primers,
            demultiplex_dir=outdir,
            tag=tag,
            outdir=configdir,
            filename=filename,
            params=params,
            cfg=cfg,
        )
        # validate_el = self.validate(
        #     config=config_el,
        # )  # this is now a pre call before process, not a separate element
        process_el = self.process(
            config=config_el,
        )
        consolidate_el = self.consolidate(
            process=process_el,
            config=config_el,
            sample_barcodes=sample_barcodes,
        )
        return consolidate_el, process_el, config_el


    def demultiplex(
        self,
        sample: NextGenSampleElement,
        template: FileElement,
        sample_barcodes: Mapping[str, str],
        *,
        primers: tuple[str, str] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        configdir: Path | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> dict[str, NextGenSampleElement]:
        """Configure and process a demultiplexing run in one step.

        In addition, a mapping of sample name → demultiplexed 
        NextGenSampleElement is returned, containing the paths to the 
        consolidated FASTQ files for each demultiplexed sample.

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
            Parameter for the ``"configure"`` step.
        cfgs : Mapping[str, ExternalRunConfig] | None
            Run configuration for the configure step.

        Returns
        -------
        dict[str, NextGenSampleElement]
            Mapping of sample name to demultiplexed NextGenSampleElement.
        """
        consolidate_el, _, _ = self.grabfast(
            sample=sample,
            template=template,
            sample_barcodes=sample_barcodes,
            primers=primers,
            tag=tag,
            outdir=outdir,
            configdir=configdir,
            filename=filename,
            params=params,
            cfg=cfg,
        )
        return self.samples(consolidate_el)

    def samples(self, consolidate_el: Element) -> tuple[dict[str, NextGenSampleElement], Element]:
        """Create NextGenSampleElements for each demultiplexed sample.

        This method takes the consolidation Element (which has artifacts for
        each consolidated FASTQ) and creates a NextGenSampleElement for each
        sample, using the consolidated FASTQs as input.

        Parameters
        ----------
        consolidate_el : Element
            The consolidation Element produced by :meth:`consolidate`.

        Returns
        -------
        dict[str, NextGenSampleElement]
            Mapping of sample name to NextGenSampleElement with FASTQ paths.
        """
        samples = {}
        for artifact_name, fq_path in consolidate_el.artifacts.items():
            if artifact_name.startswith("undetermined"):
                continue  # skip undetermined files
            if fq_path.suffix not in {".fq", ".fastq", ".fq.gz", ".fastq.gz"}:
                continue  # skip non-FASTQ artifacts
            
            sample_name = f"{consolidate_el.root}_{artifact_name[:-3]}"
            if sample_name not in samples:
                samples[sample_name] = {}
            if artifact_name.endswith("_R1"):
                samples[sample_name]["fastq_r1"] = fq_path

            elif artifact_name.endswith("_R2"):
                samples[sample_name]["fastq_r2"] = fq_path
            else:
                continue  # skip artifacts that don't match expected pattern
            # Extract sample name and read number from artifact name
        
        # Create NextGenSampleElement for each sample
        sample_elements = {}
        for sample_name, paths in samples.items():
            if paths.get("fastq_r1") is None:
                logger.warning(f"Sample {sample_name} is missing R1 FASTQ; skipping")
                continue
            tag =from_prior(
                    consolidate_el.tag,
                    root=sample_name,
                    method=Method.FASTQGRAB,
                    state=State.DEMULTIPLEX,
                    omics=Omics.DNA,
                    ext="fq",
                )
            sample_elements[sample_name] = NextGenSampleElement(
                path=paths,
                tag=tag,
                pres=(consolidate_el,),
            )
        
        return sample_elements, consolidate_el
    
    def undetermined(self, consolidate_el: Element) -> NextGenSampleElement:
        """Create a NextGenSampleElement for the undetermined reads."""
        undetermined_r1 = consolidate_el.artifacts.get("undetermined_R1")
        undetermined_r2 = consolidate_el.artifacts.get("undetermined_R2")
        if not undetermined_r1 or not undetermined_r2:
            logger.warning("Undetermined FASTQ files are missing; returning None")
            return None
        root=f"{consolidate_el.root}_Undetermined"
        
        return NextGenSampleElement(
            path={"fastq_r1": undetermined_r1, "fastq_r2": undetermined_r2},
            root=root,
            tag=from_prior(
                consolidate_el.tag,
                root=root,
                method=Method.FASTQGRAB,
                state=State.DEMULTIPLEX,
                omics=Omics.DNA,
                ext="fq",
            ),
            pres=(consolidate_el,),
        )
    

# ---------------------------------------------------------------------------
# TOML builder (private helper)
# ---------------------------------------------------------------------------

def _build_demultiplex_toml(
    *,
    template_path: Path,
    fastq_r1: Path,
    fastq_r2: Path | None,
    sample_barcodes: dict[str, str],
    forward_primer: str | None = None,
    reverse_primer: str | None = None,
    barcode_hamming_distance: int = 1,
    primer_hamming_distance: int = 3,
    barcode_anchor_distance: int = 12,
    output_prefix: str = "output",
) -> TOMLDocument:
    """Build a fastqgrab TOML configuration document.

    The ``[input]``, ``[barcodes]``, and ``[output]`` sections are constructed
    programmatically.  The ``[[step]]`` blocks are read from *template_path*
    after substituting all ``${variable}`` placeholders with the supplied
    values.

    Parameters
    ----------
    template_path : Path
        TOML template file whose ``[[step]]`` blocks define the processing
        pipeline.  Variables of the form ``${name}`` inside those blocks are
        replaced before parsing.
    fastq_r1 : Path
        Absolute path to the R1 FASTQ file.
    fastq_r2 : Path | None
        Absolute path to the R2 FASTQ file, or *None* for single-end.
    sample_barcodes : dict[str, str]
        Mapping of sample/replicate name → barcode sequence.
    forward_primer : str | None
        Forward primer sequence (IUPAC), or *None* if no primer-based
        orientation detection is needed.
    reverse_primer : str | None
        Reverse primer sequence (IUPAC).
    barcode_hamming_distance : int
        Mismatches allowed for barcode extraction and correction.
    primer_hamming_distance : int
        Mismatches allowed when searching for primers.
    barcode_anchor_distance : int
        Maximum anchor distance for barcode search.
    output_prefix : str
        Prefix for fastqgrab output files.

    Returns
    -------
    TOMLDocument
        Complete TOML document ready to be written to disk.
    """
    # substitute variables inside the step text
    text = template_path.read_text()

    variables: dict[str, object] = {
        "fwd_primer": forward_primer or "",
        "rev_primer": reverse_primer or "",
        "primer_hamming_distance": primer_hamming_distance,
        "barcode_hamming_distance": barcode_hamming_distance,
        "barcode_anchor_distance": barcode_anchor_distance,
        "sample_barcodes": list(sample_barcodes.values()),
    }
    text = _substitute_template_vars(text, variables)
    # extract steps
    steps = loads(text).get("step", [])
    
    # build the document
    doc = document()

    # [input]
    inputs = table()
    inputs.add("read1", str(fastq_r1))
    if fastq_r2 is not None:
        inputs.add("read2", str(fastq_r2))
    doc.add("input", inputs)
    doc.add(nl())
    doc.add(nl())

    # [barcodes]
    sample_barcodes_table = table()
    for sample_name, barcode in sample_barcodes.items():
        sample_barcodes_table.add(barcode, sample_name)
    barcodes_section = table()
    barcodes_section.add("sample_barcodes", sample_barcodes_table)
    doc.add("barcodes", barcodes_section)
    doc.add(nl())
    doc.add(nl())

    # [[step]] taken from template
    doc.add("step", steps)

    # [output]
    output_section = table()
    output_section.add("prefix", output_prefix)
    output_section.add("format", "FASTQ")
    output_section.add("output", ["read1", "read2"])
    output_section.add("report_json", True)
    output_section.add("report_html", True)
    doc.add("output", output_section)
    return doc


_TEMPLATE_VAR_RE = re.compile(r"""(['"])\$\{(\w+):(int|str|list)\}(['"])""")


def _substitute_template_vars(text: str, variables: dict) -> str:
    """Replace ``${varname:type}`` placeholders in *text*.

    - ``int``  – TOML-Quotes werden entfernt, Wert wird ein bare integer.
    - ``str``  – TOML-Quotes bleiben, nur der Platzhalter-Inhalt wird ersetzt.
    - ``list`` – TOML-Quotes werden entfernt, Wert wird ein Inline-Array.
    """

    def _replace(m: re.Match) -> str:
        open_q, varname, typ, close_q = m.groups()
        if varname not in variables:
            return m.group(0)
        value = variables[varname]
        if typ == "int":
            return str(int(value))  # type: ignore[arg-type]
        if typ == "str":
            return f"{open_q}{value}{close_q}"
        if typ == "list":
            return "[" + ", ".join(f'"{v}"' for v in value) + "]"  # type: ignore[union-attr]
        return m.group(0)

    return _TEMPLATE_VAR_RE.sub(_replace, text)
