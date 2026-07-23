"""Module contains an interface for fastqgrab (fastqrab), a FASTQ processing tool."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterator
from io import TextIOWrapper
from pathlib import Path
from typing import Callable, Literal, Mapping, Sequence

import pandas as pd
from tomlkit import (  # type: ignore[import]
    TOMLDocument,
    document,
    dumps,
    inline_table,
    loads,
    nl,
    table,
)

from mmalignments.models.artifacts import ArtifactSet, FastqArtifact, OutputSpec
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FastqSource,
    FileSource,
    NextGenSample,
    TableElement,
    element,
    # sample_fastqs,
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
from mmalignments.services.dependencies import depends
from mmalignments.services.io import concat_fastq, parents, write_fasta
from mmalignments.services.logging import current_call_to_string
from mmalignments.services.toml import indent_regular_tables

from ..externals import (
    External,
    ExternalRunConfig,
    Runnable,
    SubroutineIn,
    subroutine,
)
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
        return Path("results") / "demultiplex" / self.version_name / sample_name

    def default_output(self, sample_name: str) -> OutputSpec:
        """Return the default output specification for demultiplexed outputs."""
        return OutputSpec(
            sample_name,
            self.default_output_dir(sample_name),
        )

    # -----------------------------------------------------------------------
    # Configuration (high-level @element + low-level @subroutine)
    # -----------------------------------------------------------------------

    @element
    def configure(
        self,
        sample: NextGenSample,
        template: FileSource | Sequence[FileSource],
        barcodes: Callable | Element | Mapping[str, Element],
        *,
        demultiplex_dir: Path | None = None,
        config_dir: Path | None = None,
        compression: Literal["Raw", "Gzip"] = "Gzip",
        tag: PartialElementTag | ElementTag | None = None,
        outspec: OutputSpec | None = None,
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
        sample : NextGenSample
            Input sample with paired FASTQ files.
        template : FileSource | Sequence[FileSource]
            TOML template(s) with steps to be resolved by fastgrab.
        barcodes : Callable | Element | Mapping[str, Element]
            Mapping of barcode set name → Element with FASTA artifact or Callable
            returning such a mapping. If it's just a single demultipolexing step,
            this can be a single Element with FASTA artifact.
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
        compression: Literal["Raw", "Gzip"] = "Gzip"
            Compression format for the output files.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outspec : OutputSpec | None
            Optional output specification for the configuration file.
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
        fastq_r1, fastq_r2 = sample.r1, sample.r2
        tag = from_prior(
            sample.tag,
            tag,
            root=sample.root,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.CONFIG,
            omics=Omics.DNA,
        )
        spec = OutputSpec(
            tag.default_name,
            config_dir or self.default_config_dir(sample.root),
            ext="toml",
        ).merge(outspec)
        # output_dir = Path(outdir or self.default_config_dir(sample.root))
        # config_file = filename or tag.default_output
        # config_path = output_dir / config_file
        config_path = spec.path()
        demultiplex_destination = demultiplex_dir or self.default_output_dir(
            sample.root
        )
        # demultiplex_destination = Path(demultiplex_destination) / sample.root
        # forward_primer, reverse_primer = primers or (None, None)
        # Allow `template` to be a single FileElement or a sequence thereof
        pres = sample.pres
        template_sources = [template] if isinstance(template, FileSource) else template
        template_paths = []
        for temp in template_sources:
            if not isinstance(temp, FileSource):
                raise ValueError(
                    f"template must be a FileSource or a sequence of FileSource objects, was {type(template)}."  # noqa: E501
                )
            template_paths.append(temp.artifacts.primary.resolve())
            pres += temp.pres
        determinants = None
        barcodes_map: Mapping[str, Path] | Mapping[str, dict[str, str]] = {}
        if callable(barcodes):
            barcodes_map = barcodes()
            determinants = tuple(str(sorted(barcodes_map.items())))
        elif isinstance(barcodes, Element):
            # we assume the element has a Callable Artifact that returns the barcode map
            barcodes_map[barcodes.root] = barcodes.fasta
            pres += (barcodes,)
        else:
            for name, elem in barcodes.items():
                barcodes_map[name] = elem.primary.resolve()
            pres += tuple(barcodes.values())
        runner = self.write_config(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            config_path=config_path,
            barcodes_map=barcodes_map,
            template_path=template_paths,
            variables=params.to_dict(),
            compression=compression,
            output_prefix=str(demultiplex_destination / sample.root),
        )
        key, name = self.build_element_name(tag, "configure")
        artifacts = {
            "toml": config_path,
            "prefix": demultiplex_destination,
            "library_name": sample.root,
            "suffix": "fq" if compression == "Raw" else "fq.gz",
        }
        inputs = (
            tuple(template_paths) + (fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,)
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            pres=pres,
            inputs=inputs,
            name=name,
        )

    def write_config(
        self,
        fastq_r1: Path | str,
        fastq_r2: Path | str | None,
        config_path: Path | str,
        barcodes_map: Mapping[str, Path] | Mapping[str, dict[str, str]],
        template_path: Path | str | Sequence[Path | str],
        *,
        variables: dict[str, object] | None = None,
        compression: Literal["Raw", "Gzip"] = "Gzip",
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
        barcodes_map : Mapping[str, Path] |  Mapping[str, dict[str, str]]
            Mapping of sample/replicate name → path to the barcode FASTA file. Or it can
            be a mapping of sample/replicate name → dict of barcode sequence → sample
            name.
        template_path : Path | str | Sequence[Path | str]
            Path to the TOML template file with ``${variable}`` placeholders, or a
            sequence of such paths.
        variables : dict[str, object] | None
            Optional mapping of variable names → values to substitute into the
            template(s).
        compression : Literal["Raw", "Gzip"]
            Compression format for the output files.
        output_prefix : str
            Prefix for fastqgrab output files.

        Returns
        -------
        tuple
            ``(arguments, paths, output, pre, post)`` consumed by the
            ``@subroutine`` decorator.
        """
        config_path = Path(config_path).resolve()
        fastq_r1 = Path(fastq_r1).resolve()

        # Normalize template_path to a list of absolute Paths
        if isinstance(template_path, (str, Path)):
            template_paths = [Path(template_path).resolve()]
        else:
            template_paths = [Path(p).resolve() for p in template_path]

        toml_doc = _build_demultiplex_toml(
            template_paths=template_paths,
            fastq_r1=fastq_r1,
            fastq_r2=Path(fastq_r2).resolve() if fastq_r2 else None,
            barcodes_map=barcodes_map,
            variables=variables or {},
            compression=compression,
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
        runner = self.process_config(config_path=config.toml)

        tag = from_prior(
            config.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.DEMULTIPLEX,
            ext="html",
        )
        spec = OutputSpec(
            config.library_name, config.prefix, ext="json", exts=("html",)
        )
        key, name = self.build_element_name(tag, "process")
        determinants = self.signature_determinants(params, subroutine="process")

        artifacts = ArtifactSet(
            spec.files["json"], primary_name="json", html=spec.files["html"]
        )
        # this does not include the fastqs themselves, because its hidden in the toml
        return Element(
            key,
            runner,
            tag=tag,
            determinants=determinants,
            artifacts=artifacts,
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
        pre = self.run_validate(config_path=config_path)
        return (
            ["process", str(config_path), "--allow-overwrite"],
            "process",
            [config_path],
            [],
            None,
            pre,
            None,
        )

    @subroutine
    def run_validate(
        self,
        config_path: Path | str,
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
            None,
        )

    # -----------------------------------------------------------------------
    # Consolidate (high-level @element + low-level subroutine)
    # -----------------------------------------------------------------------

    @element
    def consolidate(
        self,
        process: Element,
        config: Element,
        *,
        rename: Callable[[str], str] | None = None,
        names: Iterator[str] | list[str] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
    ) -> Element:
        """Create an Element that consolidates demultiplexed FASTQ files.
        This is done after demultiplexing to combine undetermined fastqs into
        a single file and rename files in a consistent way.

        1. Merges all ``barcode_unambiguous=false`` and ``nobarcode`` FASTQs
           into ``undetermined_R1.fq`` and ``undetermined_R2.fq``, adding
           reason tags to read names.
        2. Renames ``barcode_unambiguousundetermined_r1true`` FASTQs from
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
        tag = from_prior(
            process.tag,
            tag,
            stage=Stage.PREP,
            method=Method.FASTQGRAB,
            state=State.DEMULTIPLEX,
        )

        key, name = self.build_element_name(tag, "consolidate")
        # Build artifacts: one per sample (R1 + R2) plus undetermined
        log_file = config.prefix / f"{config.prefix.name}_consolidation.log"
        artifacts = ArtifactSet(log_file, primary_name="log")
        fastqs = {}

        def check_output(sample_name: str, paired: bool = True) -> FastqArtifact:
            r1 = config.prefix / f"{sample_name}_R1.{config.suffix}"
            r2 = config.prefix / f"{sample_name}_R2.{config.suffix}" if paired else None
            return FastqArtifact(r1, r2)

        for sample_name in names:  # type: ignore
            full_name = f"{config.library_name}_{sample_name}"
            fastqs[full_name] = check_output(full_name)
        undetermined_name = f"{config.library_name}_undetermined"
        fastqs["Undetermined"] = check_output(undetermined_name)
        artifacts = artifacts.with_extras(fastqs)
        # Input files (all the process-generated FASTQs)
        # input_pattern = str(prefix / f"{prefix.name}_barcode_unambiguous=*")
        runner = self.run_consolidate(
            paths=tuple(artifacts.iterfiles()),
            prefix=config.prefix,
            rename=rename,
            ext=config.suffix,
        )
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            inputs=(),
            pres=(process,),
            name=name,
            empty_ok=True,
        )

    def run_consolidate(
        self,
        paths: tuple[Path, ...],
        prefix: Path,
        rename: Callable[[str], str] | None = None,
        ext: str = "fq",
    ) -> Runnable:
        """Consolidate demultiplexed FASTQ files.

        Parameters
        ----------
        prefix : Path
            Output prefix used by fastqgrab (directory + base name).
        ext : str, default "fq"
            File extension for the FASTQ files.

        Returns
        -------
        Runnable
            Runnable that performs the file consolidation.
        """

        def default_rename(filename: str) -> str:
            return (
                filename.replace("barcode_unambiguous=true_", "")
                .replace("_read2", "_R2")
                .replace("_read1", "_R1")
            )

        rename = rename or default_rename
        return Runnable(
            self.__consolidate(paths, rename, prefix, prefix.name, ext),
            cmd=[current_call_to_string()],
            display="consolidate",
        )

    def __consolidate(
        self,
        paths: tuple[Path, ...],
        rename: Callable[[str], str],
        prefix_dir: Path,
        prefix_name: str,
        ext: str,
    ) -> Callable[[], None]:
        """Consolidate FASTQ files after demultiplexing."""

        def __call():
            all_fqs = list(prefix_dir.resolve().glob(f"{prefix_name}*.{ext}"))
            is_paired = any(
                (re.search(r"_R2\.", fq.name) or re.search(r"_read2\.", fq.name))
                for fq in all_fqs
            )
            undetermined_r1 = prefix_dir / f"{prefix_name}_undetermined_R1.{ext}"
            undetermined_r2 = (
                prefix_dir / f"{prefix_name}_undetermined_R2.{ext}"
                if is_paired
                else None
            )
            remaining = [x for x in all_fqs if x not in paths]

            log_file = prefix_dir / f"{prefix_name}_consolidation.log"
            # raise ValueError()
            if remaining:
                with open(log_file, "w") as log_fh:
                    log_fh.write(
                        f"Consolidating FASTQ files for prefix {prefix_name}\n"
                    )
                    log_fh.write(f"Found {len(all_fqs)} FASTQ files:\n")

                    for fq in all_fqs:
                        log_fh.write(f"  {fq.name}\n")

                    successful, temporary = _merge_undetermined(
                        remaining,
                        undetermined_r1,
                        undetermined_r2,
                        log_fh,
                    )
                    # collect the old files
                    ren = _rename_successful(
                        successful,
                        rename=rename,
                        prefix_dir=prefix_dir,
                        log_fh=log_fh,
                    )
                    temporary.extend(ren)

                    _cleanup(
                        temporary,
                        log_fh,
                    )

        return __call

    # -----------------------------------------------------------------------
    # create barcode fasta
    # -----------------------------------------------------------------------

    @element
    def barcodefasta(
        self,
        barcodes: TableElement | FileSource,
        *,
        name_column: str = "sample",
        barcode_column: str = "barcode",
        outspec: OutputSpec | None = None,
        tag: PartialElementTag | ElementTag | None = None,
    ) -> Element:

        tag = from_prior(
            barcodes.tag,
            state=State.PREPROCESS,
            method=Method.FASTQGRAB,
            ext="fasta",
            param=barcode_column,
        )
        outspec = OutputSpec(
            stem=tag.default_name,
            outdir=self.default_config_dir(barcodes.root),
            ext="fasta",
        ).merge(outspec)
        barcodefile = barcodes.artifacts.primary.resolve()
        outfile = outspec.path(tag.default_name)
        key = f"{tag.default_name}_flanks_fasta"
        runner = self.write_barcode_fasta(
            sensorfile=barcodefile,
            outfile=outfile,
            name_column=name_column,
            barcode_column=barcode_column,
        )
        artifacts = {
            "fasta": outfile,
        }
        return Element(
            key,
            runner,
            tag=tag,
            inputs=(barcodefile,),
            artifacts=artifacts,
            pres=(barcodes,) if isinstance(barcodes, Element) else (),
            name=key,
        )

    def write_barcode_fasta(
        self, sensorfile: Path, outfile: Path, name_column: str, barcode_column: str
    ) -> Runnable:

        def __run():
            df_constructs = pd.read_csv(sensorfile, sep="\t")
            sample_barcodes = dict(
                zip(df_constructs[name_column], df_constructs[barcode_column])
            )
            write_fasta(outfile, sample_barcodes)

        callspec = CallSpec(
            path=("write_barcode_fasta",),
            kwargs={
                "sensorfile": sensorfile,
                "outfile": outfile,
                "name_column": name_column,
                "barcode_column": barcode_column,
            },
        ).render()
        return Runnable(
            __run,
            display=callspec,
        )

    # we also need a "trim and mark" element here to extract the flanking regions and
    # mark the reads accordingly, also to trim off the scaffold parts

    # -----------------------------------------------------------------------
    # Convenience
    # -----------------------------------------------------------------------

    def grabfast(
        self,
        sample: NextGenSample,
        template: FileSource | Sequence[FileSource],
        barcodes: Callable | Element | Mapping[str, Element],
        *,
        rename: Callable[[str], str] | None = None,
        names: Iterator[str] | list[str] | None = None,
        compression: Literal["Raw", "Gzip"] = "Gzip",
        tag: PartialElementTag | ElementTag | None = None,
        # outspec: OutputSpec | None = None,
        demultiplex_dir: Path | None = None,
        config_dir: Path | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[Element, Element, Element]:
        """Configure and process a demultiplexing run in one step.

        This convenience method calls :meth:`configure` followed by
        :meth:`process` and returns both Elements.

        Parameters
        ----------
        sample : NextGenSample
            Input sample with paired FASTQ files.
        template : FileSource | Sequence[FileSource]
            TOML template(s) with steps to be resolved by fastgrab.
        barcodes : Callable | Element | Mapping[str, Element]
            Mapping of barcode set name → Element with FASTA artifact or Callable
            returning such a mapping. If it's just a single demultiplexing step,
            this can be a single Element with FASTA artifact.
        rename : Callable[[str], str] | None
            Optional function to rename consolidated FASTQ files. If *None*, a default
            renaming scheme is applied that removes the "barcode_unambiguous=true_"
            prefix and replaces "_read1" and "_read2" with "_R1" and "_R2".
        names : Iterator[str] | list[str] | None
            Optional list of sample names to expect in the consolidated output. If
            *None*, all samples found in the process output are included.
        tag : PartialElementTag | ElementTag | None
            Optional tag override for all created Elements.
        compression : Literal["Raw", "Gzip"]
            Compression method for the consolidated FASTQ files.  If *Gzip*, the files
            are compressed with gzip; if *Raw*, they are left uncompressed.
        demultiplex_dir : Path | None
            Base output directory for the demultiplexed FASTQs.  If *None*, a default
            directory is derived from the sample name.
        config_dir : Path | None
            Directory in which to write the configuration file.  If *None*, a default
            directory is derived from the sample name.
        filename : Path | str | None
        outspec : OutputSpec | None
            Output specification for the demultiplexed FASTQs.  If *None*, a default
            specification is derived from the sample name.
        outdir : Path | str | None
            Base output directory for the demultiplexed FASTQs.  If *None*, a default
            directory is derived from the sample name.
        configdir : Path | None
            Directory in which to write the configuration file.  If *None*, a default
            directory is derived from the sample name.
        filename : Path | str | None
            Filename for the configuration file.  If *None*, a default name is derived
            from the sample name.
        params : Params | None
            Parameter sets keyed by step name (``"configure"``,
            ``"process"``).
        cfg : ExternalRunConfig | None
            Run configurations keyed by step name.

        Returns
        -------
        tuple[Element, Element, Element]
            ``(consolidate_element, process_element, config_element)``
        """
        config_el = self.configure(
            sample=sample,
            template=template,
            barcodes=barcodes,
            demultiplex_dir=demultiplex_dir,
            config_dir=config_dir,
            compression=compression,
            tag=tag,
            params=params,
            cfg=cfg,
        )
        process_el = self.process(config=config_el, tag=tag)
        consolidate_el = self.consolidate(
            process=process_el, config=config_el, rename=rename, names=names, tag=tag
        )
        return consolidate_el, process_el, config_el

    def demultiplex(
        self,
        sample: NextGenSample,
        template: FileSource | Sequence[FileSource],
        barcodes: Callable | Element | Mapping[str, Element],
        *,
        rename: Callable[[str], str] | None = None,
        names: Iterator[str] | list[str] | None = None,
        compression: Literal["Raw", "Gzip"] = "Gzip",
        tag: PartialElementTag | ElementTag | None = None,
        demultiplex_dir: Path | None = None,
        config_dir: Path | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> tuple[dict[str, NextGenSample], NextGenSample | None]:
        """Configure and process a demultiplexing run in one step.

        In addition, a mapping of sample name → demultiplexed
        NextGenSampleElement is returned, containing the paths to the
        consolidated FASTQ files for each demultiplexed sample.

        Parameters
        ----------
        sample : NextGenSample
            Input sample with paired FASTQ files.
        template : FileSource | Sequence[FileSource]
            TOML template(s) with steps to be resolved by fastgrab.
        barcodes : Callable | Element | Mapping[str, Element]
            Mapping of barcode set name → Element with FASTA artifact or Callable
            returning such a mapping. If it's just a single demultipolexing step,
            this can be a single Element with FASTA artifact.
        rename : Callable[[str], str] | None
            Optional function to rename consolidated FASTQ files. If *None*, a default
            renaming scheme is applied that removes the "barcode_unambiguous=true_"
            prefix and replaces "_read1" and "_read2" with "_R1" and "_R2".
        names : Iterator[str] | list[str] | None
            Optional list of sample names to expect in the consolidated output. If
            *None*, all samples found in the process output are included.
        compression : Literal["Raw", "Gzip"]
            Compression method for the consolidated FASTQ files.  If *Gzip*, the files
            are compressed with gzip; if *Raw*, they are left uncompressed.
        tag : PartialElementTag | ElementTag | None
            Optional tag override for all created Elements.
        demultiplex_dir : Path | None
            Base output directory for the demultiplexed FASTQs.  If *None*, a default
            directory is derived from the sample name.
        config_dir : Path | None
            Directory in which to write the configuration file.  If *None*, a default
            directory is derived from the sample name.
        filename : Path | str | None
            Filename for the configuration file.  If *None*, a default name is derived
            from the sample name.
        params : Mapping[str, Params] | None
            Parameter for the ``"configure"`` step.
        cfg : ExternalRunConfig | None
            Run configuration for the configure step.

        Returns
        -------
        dict[str, NextGenSample]
            Mapping of sample name to demultiplexed NextGenSample.
        """
        consolidate_el, _, _ = self.grabfast(
            sample=sample,
            template=template,
            barcodes=barcodes,
            rename=rename,
            names=names,
            compression=compression,
            tag=tag,
            demultiplex_dir=demultiplex_dir,
            config_dir=config_dir,
            filename=filename,
            params=params,
            cfg=cfg,
        )
        samples, undetermined = self.samples(consolidate_el)
        return samples, undetermined

    def samples(
        self, consolidate_el: Element
    ) -> tuple[dict[str, NextGenSample], NextGenSample | None]:
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
        dict[str, NextGenSample]
            Mapping of sample name to NextGenSample with FASTQ paths.
        """
        samples = {}
        sample_elements: dict[str, NextGenSample] = {}
        for artifact_name, artifact in consolidate_el.artifacts.items():
            if artifact_name.startswith("Undetermined"):
                continue  # skip undetermined files
            if not isinstance(artifact, FastqArtifact):
                continue  # skip non-FASTQ artifacts
            samples[artifact_name] = artifact
            if artifact.paired:
                logger.warning(f"Sample {artifact_name} is paired")
            else:
                logger.warning(f"Sample {artifact_name} is single-end")

            tag = from_prior(
                consolidate_el.tag,
                root=artifact_name,
                method=Method.FASTQGRAB,
                state=State.DEMULTIPLEX,
                omics=Omics.DNA,
            )
            source = FastqSource(artifact, consolidate_el)
            sample_elements[artifact_name] = NextGenSample(
                artifact_name,
                source,
                tag=tag,
            )
        undetermined = self.undetermined(consolidate_el)
        return sample_elements, undetermined

    def undetermined(self, consolidate_el: Element) -> NextGenSample | None:
        """Create a NextGenSampleElement for the undetermined reads."""
        undetermined = consolidate_el.artifacts.get("Undetermined")
        if not undetermined:
            logger.warning("Undetermined FASTQ files are missing; returning None")
            return None
        root = f"{consolidate_el.root}_Undetermined"
        source = FastqSource(
            FastqArtifact(
                r1=undetermined.r1,
                r2=undetermined.r2,
            ),
            consolidate_el,
        )

        return NextGenSample(
            "Undetermined",
            source,
            tag=from_prior(
                consolidate_el.tag,
                root=root,
                method=Method.FASTQGRAB,
                state=State.DEMULTIPLEX,
                omics=Omics.DNA,
                ext="fq",
            ),
        )


# ---------------------------------------------------------------------------
# TOML builder (private helper)
# ---------------------------------------------------------------------------


def _build_demultiplex_toml(
    *,
    template_paths: Sequence[Path],
    fastq_r1: Path,
    fastq_r2: Path | None,
    barcodes_map: Mapping[str, Path] | Mapping[str, dict[str, str]],
    variables: dict[str, object] | None = None,
    compression: Literal["Raw", "Gzip"] = "Gzip",
    output_prefix: str = "output",
) -> TOMLDocument:
    """Build a fastqgrab TOML configuration document.

    The ``[input]``, ``[barcodes]``, and ``[output]`` sections are constructed
    after substituting all ``${variable}`` placeholders with the supplied
    values.  The templates are applied in the order they appear in the
    *template_paths* sequence.
    values.

    Parameters
    ----------
    template_paths : Sequence[Path]
        TOML template files whose ``[[step]]`` blocks define the processing
        pipeline.  Variables of the form ``${name}`` inside those blocks are
        replaced before parsing.
    fastq_r1 : Path
        Absolute path to the R1 FASTQ file.
    fastq_r2 : Path | None
        Absolute path to the R2 FASTQ file, or *None* for single-end.
    barcodes_map : Mapping[str, Path] | Mapping[str, dict[str, str]]
        Mapping of sample/replicate name → barcode fasta. Or it can be a mapping of
        sample/replicate name → dict of barcode sequence → sample name.
    variables : dict[str, object] | None
        Optional mapping of variable names → values to substitute into the
        template(s).
    compression : Literal["Raw", "Gzip"]
        Compression method for the output FASTQ files.
    output_prefix : str
        Prefix for fastqgrab output files.

    Returns
    -------
    TOMLDocument
        Complete TOML document ready to be written to disk.
    """

    steps: list = []
    for tp in template_paths:
        text = Path(tp).read_text()
        if variables:
            text = _substitute_template_vars(text, variables)
        steps_part = loads(text).get("step", [])
        if steps_part:
            # extend in-order so templates append sequentially
            steps.extend(steps_part)

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
    barcodes_section = table()
    for map_name, sample_barcodes in barcodes_map.items():
        subtable = table()
        if isinstance(sample_barcodes, Path):
            # if it is a file pointer
            from_file = inline_table()
            from_file["filename"] = str(sample_barcodes)
            subtable["from_file"] = from_file
            barcodes_section[map_name] = subtable
        else:
            # if it is a direct mapping
            sample_barcodes_table = table()
            for sample_name, barcode in sample_barcodes.items():
                sample_barcodes_table.add(barcode, sample_name)
                barcodes_section = table()
                barcodes_section.add(map_name, sample_barcodes_table)

    doc.add("barcodes", barcodes_section)
    doc.add(nl())
    doc.add(nl())
    # [[step]] taken from template
    doc.add("step", steps)

    # [output]
    output_section = table()
    output_section.add("prefix", output_prefix)
    output_section.add("format", "FASTQ")
    # output_section.add("output", ["read1", "read2"])
    suffix = "fq" if compression == "Raw" else "fq.gz"
    output_section.add("suffix", suffix)
    output_section.add("compression", compression)
    output_section.add("report_json", True)
    output_section.add("report_html", True)
    doc.add("output", output_section)
    return doc


_TEMPLATE_VAR_RE = re.compile(r"""(['"])\$\{(\w+):(int|str|list)\}(['"])""")


def _substitute_template_vars(text: str, variables: dict[str, object]) -> str:
    """Replace ``${varname:type}`` placeholders in *text*.

    - ``int``  – remove TOML-Quotes, value will be a bare integer.
    - ``str``  – TOML-Quotes remain, only the placeholder content is replaced.
    - ``list`` – remove TOML-Quotes, value will be an inline array.
    """

    def _replace(m: re.Match) -> str:
        open_q, varname, typ, close_q = m.groups()
        if varname not in variables:
            raise ValueError(
                "Substitution pattern found for undefined variable: " + m.group(0)
            )
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


def generate_sample_iterator(
    library_name: str, replicates: list[str] = ["Rep1", "Rep2", "Rep3"]
) -> Iterator[str]:

    def sample_names():
        for rep in replicates:
            yield f"{library_name}-{rep}"

    return sample_names()


########################################################################################
# consolidate helper
########################################################################################


def _reason_from_name(name: str) -> str:
    if "nobarcode" in name:
        return "reason:nobarcode"
    if "barcode_unambiguous=false" in name:
        return "reason:ambiguous_barcode"
    raise ValueError(f"Unexpected undetermined FASTQ file name: {name}")


def _get_file_handle(
    name: str,
    undetermined_r1: Path,
    undetermined_r2: Path | None,
) -> Path:
    if undetermined_r2 is None:
        return undetermined_r1
    if "_read1." in name or "_R1" in name:
        return undetermined_r1
    if "_read2." in name:
        return undetermined_r2
    raise ValueError(
        f"Unexpected FASTQ file name (cannot determine read direction): {name}"
    )


def _merge_undetermined(
    files: list[Path],
    undetermined_r1: Path,
    undetermined_r2: Path | None,
    log_fh: TextIOWrapper,
) -> tuple[list[Path], list[Path]]:
    def is_undetermined(name: str) -> bool:
        return "barcode_unambiguous=false" in name or "nobarcode" in name

    successful = []
    temporary = []

    for fq_path in files:
        name = fq_path.name

        if not is_undetermined(name):
            successful.append(fq_path)
            continue

        reason = _reason_from_name(name)
        out_fh = _get_file_handle(name, undetermined_r1, undetermined_r2)

        log_fh.write(f"Adding {name} to undetermined: {reason}\n")
        concat_fastq(inputs=(fq_path,), output=out_fh)
        temporary.append(fq_path)

    return successful, temporary


def _rename_successful(
    files: list[Path],
    rename: Callable[[str], str],
    prefix_dir: Path,
    log_fh: TextIOWrapper,
) -> list[Path]:

    renamed = []
    for fq_path in files:
        new_name = rename(fq_path.name)
        # if "barcode_unambiguous=true" not in fq_path.name:
        #     continue

        new_path = prefix_dir / new_name

        log_fh.write(f"Renaming {fq_path} to {new_name}\n")

        if fq_path.exists():
            fq_path.rename(new_path)
        renamed.append(fq_path)

    return renamed


def _cleanup(
    files: list[Path],
    log_fh: TextIOWrapper,
):
    for fq_path in files:
        if fq_path.exists():
            log_fh.write(f"Removing intermediate file {fq_path}\n")
            fq_path.unlink()
