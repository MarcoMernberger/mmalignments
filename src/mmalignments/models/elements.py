from __future__ import annotations

import subprocess
import traceback
from functools import cached_property, wraps
from pathlib import Path
from subprocess import CompletedProcess
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    Mapping,
    ParamSpec,
    TypeAlias,
    TypeVar,
    cast,
    overload,
    runtime_checkable,
)

from mmalignments.models.data import Pairing
from mmalignments.services.dependencies import (
    Signifiable,
    combined_signature,
    file_signature,
    stable_hash,
)
from mmalignments.services.io import parents, paths_exists

from .artifacts import (
    Artifact,
    ArtifactSet,
    FastqArtifact,
    TableArtifact,
)
from .registry import current_element_registry
from .specs import CallSpec, Runnable, ValidationPolicy  # type: ignore[import]
from .tags import (
    Method,
    Omics,
    PartialElementTag,
    Stage,
    State,
)
from mmalignments.models.overlay import (
    ElementTag,
)


def generate_element_key_name(
    tag: ElementTag | None,
    tool_name: str,
    tool_version: str | None = None,
    subcommand: str | None = None,
    suffix: str | None = None,
    **params_str: Any,
) -> tuple[str, str]:
    """Generate deterministic element key/name fragments for caching and provenance."""
    tag_name = getattr(tag, "default_name", "default") if tag is not None else "default"
    parts = [tag_name, tool_name]
    if subcommand:
        parts.append(subcommand)
    short = "_".join(parts)
    if tool_version:
        parts.append(tool_version)
    if suffix:
        parts.append(suffix)
    for k, v in params_str.items():
        if v is not None:
            parts.append(f"{k}-{v}")
    return "_".join(parts), short


################################################################################
# Aliases
################################################################################

RunType: TypeAlias = (
    Callable[[], CompletedProcess | None | bool | Any] | Runnable
)  # noqa: E501
# FilesSourceType: TypeAlias = Path | str | Mapping[str, Path | str] | RunType


################################################################################
# Sources
################################################################################


@runtime_checkable
class Prerequisite(Signifiable):
    """Source for input files"""

    @cached_property
    def key(self) -> str: ...

    """a unique key for this source."""

    @cached_property
    def provenance(self) -> str: ...

    """the provenance string for this source: which steps came before?"""


@runtime_checkable
class Source(Prerequisite):
    """Source for input files, not necessarily an Element."""

    @property
    def tag(self) -> ElementTag: ...

    """A description tag used by downstream Elements."""

    @property
    def artifacts(self) -> ArtifactSet: ...

    """The set of artifacts (files) associated with this source."""

    @property
    def primary(self) -> Artifact: ...

    """The primary artifact (e.g. FastqArtifact) associated with this source."""


class FileSource(Source):  # the new FileElement for existing files
    """
    A Source that represents a set of pre-existing files, typically used as
    input for an Element.
    """

    def __init__(
        self,
        path: Path | str | TableArtifact | Mapping[str, Path | TableArtifact],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        is_prefix: bool = False,
    ):
        """
        Initialize a FileSource instance.

        Parameters
        ----------
        path : Path | str | TableArtifact | Mapping[str, Path  |  TableArtifact]
            The path to the file(s) or a mapping of artifact names to paths.
        tag : PartialElementTag | ElementTag | None, optional
            The tag associated with this FileSource, defaults to a standard tag
            if not provided.
        is_prefix : bool, optional
            Whether the provided path is a prefix for multiple files, by default
            False
        """
        self._artifacts = self.normalize(path, is_prefix)
        self._tag = ElementTag(
            root=self._artifacts.primary.stem,
            level=0,
            omics=Omics.DNA,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
        ).patch(tag)
        self._key = f"FileSource:{self._tag}" + "::".join(
            sorted(self._artifacts.keys())
        )

    ############################################################################
    # Properties
    ############################################################################

    @property
    def root(self) -> str:
        """The root name of this FileSource, derived from its tag."""
        return self.tag.root

    # Prerequisite protocol –––––––––––––––––––––––––––––––––––––––––––––––––––

    @cached_property
    def signature(self) -> str:
        """The signature of the source used for invalidation."""
        return self.artifacts.primary.signature()

    @cached_property
    def key(self) -> str:
        """The unique key for this Source, used for registry identification."""
        return self._key

    @cached_property
    def provenance(self) -> str:
        """The provenance string for this Source, describing its lineage."""
        return self.tag.default_name

    # Source protocol ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

    @property
    def tag(self) -> ElementTag:
        """A description tag used by downstream Elements."""
        return self._tag

    @property
    def primary(self) -> Artifact:
        """The primary artifact associated with this Source."""
        return self.artifacts.primary  # (e.g. FastqArtifact)

    @property
    def artifacts(self) -> ArtifactSet:
        """The set of artifacts (files) associated with this Source."""
        return self._artifacts

    ############################################################################
    # Initalization
    ############################################################################

    @classmethod
    def artifacts_from_prefix(
        cls, prefix: str | Path, primary_key: str | None = None
    ) -> ArtifactSet:
        """
        Generate artifacts from a directory prefix.

        This method scans the directory containing the given prefix and collects
        all files whose names start with the prefix. It returns an ArtifactSet
        containing these files.

        Parameters
        ----------
        prefix : str | Path
            The directory prefix to scan for files.
        primary_key : str | None, optional
            The key of the primary artifact, by default None.

        Returns
        -------
        ArtifactSet
            An ArtifactSet containing the files matching the prefix.

        Raises
        ------
        ValueError
            If no files are found with the given prefix.
        """
        artifacts = {}
        p = Path(prefix)
        file_dir = p.parent
        file_prefix = p.stem
        for p in file_dir.iterdir():
            if p.stem.startswith(file_prefix):
                artifacts[p.stem] = p.resolve()
        if len(artifacts) == 0:
            raise ValueError(f"No files found with prefix {prefix}")
        primary_key = primary_key or next(iter(artifacts.keys()))
        return ArtifactSet.from_any(artifacts, primary_key)

    @classmethod
    def normalize(
        cls,
        path: Path | str | TableArtifact | Mapping[str, Path | TableArtifact],
        is_prefix: bool,
        primary_key: str | None = None,
    ) -> ArtifactSet:
        """
        Normalize a path into an ArtifactSet.

        This method handles different types of input paths, including directory
        prefixes, individual files, and mappings of artifact names to paths. It
        returns an ArtifactSet containing the resolved artifacts.

        Parameters
        ----------
        path : Path | str | TableArtifact | Mapping[str, Path  |  TableArtifact]
            The path to normalize. Can be a directory prefix, individual file,
            or a mapping of artifact names to paths.
        is_prefix : bool
            Whether the given path is a directory prefix.
        primary_key : str | None, optional
            The key of the primary artifact, by default None.

        Returns
        -------
        ArtifactSet
            An ArtifactSet containing the resolved artifacts.
        """
        if is_prefix and isinstance(path, (str, Path)):
            artifacts = cls.artifacts_from_prefix(path, primary_key)
        else:
            if isinstance(path, Mapping):
                primary_key = primary_key or next(iter(path.keys()))
                artifacts = ArtifactSet.from_any(
                    path,
                    primary_key,
                )
            elif isinstance(path, Artifact):
                artifacts = ArtifactSet(
                    path,
                    primary_name=primary_key
                    or path.path.suffix.lstrip("."),  # noqa: E501
                )
            else:
                path = Path(path).resolve()
                primary_key = primary_key or path.suffix.lstrip(".")
                artifacts = ArtifactSet.from_any(
                    {primary_key: path}, primary_key
                )  # noqa: E501

        return artifacts

    def __getattr__(self, name: str) -> Any:
        """
        Allow convenient access to artifacts by name.

        Example:
            source.toml
            source.parquet

        resolves to:
            source.artifacts["toml"]
            source.artifacts["parquet"]
        """
        artifacts = object.__getattribute__(self, "_artifacts")

        if name in artifacts:
            return artifacts[name]

        raise AttributeError(
            f"{type(self).__name__!s} has no attribute {name!r} "
            f"and no artifact with that name exists"
        )


class FastqSource(Source):
    """Source class for handling FASTQ files."""

    def __init__(
        self,
        fastqs: FastqArtifact | Path | str | Mapping[str, Path],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        is_prefix: bool = False,
    ):
        """
        Initialize a FastqSource instance.

        This constructor sets up the FastqSource with the given FASTQ files,
        tag, and prefix information. It normalizes the input FASTQ files into
        an ArtifactSet and constructs the element tag.

        Parameters
        ----------
        fastqs : FastqArtifact | Path | str | Mapping[str, Path]
            The FASTQ files to normalize. Can be a directory prefix, individual
            file, or a mapping of artifact names to paths.
        tag : PartialElementTag | ElementTag | None, optional
            The element tag for the FastqSource, by default None.
        is_prefix : bool, optional
            Whether the given fastqs path is a directory prefix, by default
            False.
        """
        self._artifacts = self.normalize(fastqs, is_prefix)
        tag = ElementTag(
            root=self._artifacts.primary.stem,
            level=0,
            omics=Omics.DNA,
            stage=Stage.INPUT,
            method=Method.CHECK,
            state=State.RAW,
        ).merge(tag)
        self._tag = tag
        self._key = f"FastqSource:{self._tag}" + "::".join(
            ff for ff in self.artifacts.primary
        )

    ############################################################################
    # Properties
    ############################################################################

    @property
    def root(self) -> str:
        """The root of the element tag. A short handle for the source."""
        return self.tag.root

    # Prerequisite protocol ––––––––––––––––––––––––––––––––––––––––––––––––––––

    @cached_property
    def key(self) -> str:
        """The unique key for this Source, used for registry identification."""
        return self._key

    @cached_property
    def provenance(self) -> str:
        """The provenance string for this Source, describing its lineage."""
        return self.key

    @cached_property
    def signature(self) -> str:
        """The signature of the source used for invalidation."""
        return self.artifacts.primary.signature()

    # Source protocol ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

    @property
    def tag(self) -> ElementTag:
        """The description tag used by downstream Elements."""
        return self._tag

    @property
    def primary(self) -> FastqArtifact:
        """The primary artifact associated with this Source."""
        return self._artifacts.primary  # (e.g. FastqArtifact)

    @property
    def artifacts(self) -> ArtifactSet:
        """The set of artifacts (files) associated with this Source."""
        return self._artifacts

    ############################################################################
    # Initalization
    ############################################################################

    def normalize(
        self,
        fastqs: FastqArtifact | Path | str | Mapping[str, Path],
        is_prefix: bool,
    ) -> ArtifactSet:
        """
        Normalize FASTQ files into an ArtifactSet.

        This method takes various forms of FASTQ input (single files, directory
        prefixes, or mappings) and normalizes them into a consistent ArtifactSet
        structure.

        Parameters
        ----------
        fastqs : FastqArtifact | Path | str | Mapping[str, Path]
            The FASTQ input to normalize.
        is_prefix : bool
            Whether the input is a directory prefix.

        Returns
        -------
        ArtifactSet
            The normalized ArtifactSet containing the FASTQ files.

        Raises
        ------
        ValueError
            If no R1 files are found with the given prefix.
        ValueError
            If multiple FASTQ files per lane are found.
        NotImplementedError
            If the input type is not supported.
        TypeError
            If the input type is incorrect.
        """

        if isinstance(fastqs, FastqArtifact):
            return ArtifactSet(
                fastqs,
                primary_name="fastq",
            )

        if is_prefix and isinstance(fastqs, (str, Path)):
            r1 = []
            r2 = []
            p = Path(fastqs)
            file_prefix = p.stem
            file_dir = p.parent
            file_prefix = p.stem
            for p in file_dir.iterdir():
                if p.stem.startswith(file_prefix):
                    if "_R1_" in p.stem or "_1" in p.stem:
                        r1.append(p.resolve())
                    elif "_R2_" in p.stem or "_2" in p.stem:
                        r2.append(p.resolve())
                    else:
                        raise ValueError(
                            f"File {p} does not match expected R1/R2 naming convention"  # noqa: E501
                        )
            if len(r1) == 0:
                raise ValueError(f"No files found with prefix {file_prefix}")
            elif len(r1) > 1:
                print("Found: ", r1)
                raise NotImplementedError(
                    "Multiple fastq files per lane are discouraged. Please use FastqNormalize first."  # noqa: E501
                )
            else:
                fastq = FastqArtifact(r1[0], r2[0] if r2 else None)
                primary_name = "fastq"
            return ArtifactSet(fastq, primary_name=primary_name)

        if isinstance(fastqs, (str, Path)):
            return ArtifactSet(
                FastqArtifact(Path(fastqs)),
                primary_name="fastq",
            )

        raise TypeError(
            f"Unsupported type for fastqs: {type(fastqs)}. Must be FastqArtifact, Path, str, or Mapping[str, Path]."  # noqa: E501
        )


class NextGenSample(Source):
    """A Source that represents a next-generation sequencing sample."""

    def __init__(
        self,
        root: str,  # short sample name
        source: Source,  # An Element producing a FastqArtifact or a FastqSource
        *,
        tag: PartialElementTag | ElementTag | None = None,
        read_group: str | None = None,
        reverse_reads: bool = False,
        cache_dir: Path | None = None,
        result_dir: Path | None = None,
    ):
        """
        Initialize a NextGenSample instance.

        This class represents a next-generation sequencing sample, encapsulating
        the source of the sample, associated artifacts, and metadata such as
        read group and pairing information.

        Parameters
        ----------
        root : str
            The short sample name, used as a handle for the source.
        tag : PartialElementTag | ElementTag | None, optional
            The tag associated with the sample, by default None
        read_group : str | None, optional
            The read group identifier for the sample, by default None
        reverse_reads : bool, optional
            Whether the reads are reversed, by default False
        cache_dir : Path | None, optional
            The directory to use for caching, by default None
        result_dir : Path | None, optional
            The directory to use for storing results, by default None
        """
        self.source = source
        self._artifacts = source.artifacts
        self._key = f"{self._tag.default_name}" + "::".join(
            [str(filepath) for filepath in self._artifacts.primary]
        )
        self._tag = self.source.tag.merge(tag).merge(
            PartialElementTag(root=root)
        )  # noqa: E501
        self.cache_dir = (
            cache_dir or Path(f"cache/samples/{self._tag.root}").resolve()
        )  # noqa: E501
        self.result_dir = (
            result_dir or Path(f"results/samples/{self._tag.root}").resolve()
        )
        self.read_group = read_group or self._tag.root
        self.reverse_reads = reverse_reads

    ############################################################################
    # Properties
    ############################################################################

    @property
    def pres(self) -> tuple[Prerequisite, ...]:
        """
        Return the prerequisites for this NextGenSample if there is one. If
        the source is an Element, it is included as a prerequisite; otherwise,
        an empty tuple is returned.
        """
        if isinstance(self.source, Element):
            return (self.source,)
        return ()

    # Prerequisite protocol –––––––––––––––––––––––––––––––––––––––––––––––––––

    @cached_property
    def signature(self) -> str:
        """The signature of the source used for invalidation."""
        return self.source.signature

    @cached_property
    def key(self) -> str:
        """The unique key for this Source, used for registry identification."""
        return self._key

    @cached_property
    def provenance(self) -> str:
        """
        A human-readable provenance string that describes the lineage of this
        Element.
        """
        return self.source.provenance + "->" + self.root

    # Source protocol –––––––––––––––––––––––––––––––––––––––––––––––––––––––––

    @property
    def tag(self) -> ElementTag:
        """A description tag used by downstream Elements."""
        return self._tag

    @property
    def primary(self) -> FastqArtifact:
        """The primary artifact associated with this Source."""
        return self.source.artifacts.primary

    @property
    def artifacts(self) -> ArtifactSet:
        """The set of artifacts (files) associated with this Source."""
        return self._artifacts

    @cached_property
    def root(self) -> str:
        """The root name of this NextGenSample, derived from its tag."""
        return self._tag.root

    @property
    def pairing(self) -> Pairing:
        """The pairing type of the primary artifact."""
        return Pairing.PAIRED if self.primary.paired else Pairing.SINGLE

    @property
    def r1(self) -> Path:
        """The path to the R1 FASTQ file."""
        return self.primary.r1

    @property
    def r2(self) -> Path | None:
        """The path to the R2 FASTQ file, if it exists."""
        return self.primary.r2

    def fastqs(self) -> Iterator[Path]:
        """An iterator over the FASTQ files associated with this Source."""
        return iter(self.primary)

    @cached_property
    def files(self) -> tuple[Path, ...]:
        """A tuple of all output files for this NextGenSample."""
        return self.primary.files


###############################################################################
# Elements
###############################################################################


class Element(Prerequisite):
    """
    An Element represents a unit of work in a computational pipeline,
    encapsulating its key, run, tag, artifacts, and other relevant metadata.
    This is the key class for defining and managing computational steps in a
    reproducible manner.
    """

    REQUIRED = ["key", "run", "tag"]

    def __init__(
        self,
        key: str,
        run: RunType,
        tag: ElementTag,
        artifacts: ArtifactSet,
        *,
        determinants: tuple[str, ...] | None = None,
        inputs: tuple[Path, ...] | None = None,
        validator: Callable[[], tuple[bool, str]] | None = None,
        pres: tuple[Prerequisite, ...] | None = None,
        empty_ok: bool = False,
    ):
        """
        Initialize an Element.

        This method sets up the Element with its key, run, tag, artifacts, and
        other optional parameters.

        Parameters
        ----------
        key : str
            A unique key to identify this Element in the registry. Identical
            keys will raise Exceptions.
        run : RunType
            A parameterless callable function that executes the Element's work
            and creates new artifacts, then returns them as a result.
        tag : ElementTag
            The description tag associated with this Element.
        artifacts : ArtifactSet
            The set of artifacts (files) associated with this Element.
        determinants : tuple[str, ...] | None, optional
            The determinants (key factors) for this Element, by default None.
            This should be a list of parameters that may affect the output of
            the Element. These determinants are used to compute the signature of
            the Element and determine if it needs to be re-run.
        inputs : tuple[Path, ...] | None, optional
            The input files required by this Element, by default None. These
            also contribute to the signature of the Element and are used to
            determine if the Element needs to be re-run.
        validator : Callable[[], tuple[bool, str]] | None, optional
            A callable that validates the Element, returning a tuple of a
            boolean indicating success and a string message, by default None.
        pres : tuple[Element, ...] | None, optional
            The prerequisite Elements that must be completed before this
            Element, by default None. Only the immediate predecessors are
            required, not the full transitive closure.
        empty_ok : bool, optional
            Whether it is acceptable for this Element to have no output files,
            by default False.
        """
        self._key = key
        self.run = run
        self.tag = tag
        self.name = tag.default_name
        self.artifacts = ArtifactSet.from_any(artifacts)
        self.pres = tuple(pres or ())
        self.validator = validator
        self.determinants = tuple(str(det) for det in (determinants or ()))
        if self.determinants:
            self._key = f"{self._key}:" + self.determinants_as_str()
        self.inputs = tuple(sorted(inputs or [], key=lambda p: str(p)))
        self._signature: str | None = None  # cache
        self._validation_policy = ValidationPolicy.CHECK
        self.validate_fields()
        self.build_provenance()
        self.empty_ok = empty_ok
        # self.output_files  # trigger validation of output file paths
        # for error stack traces, we keep the creation trace
        self.creation_trace = "".join(traceback.format_stack()[:-1])
        # self.threads = getattr(
        #     run, "threads", None
        # )  # optional attribute for external runners
        # self.name = (
        #     name or self.default_name
        # )  # supplied name overrided tag-based default name

    def validate_fields(self) -> None:
        """
        Validate that required fields are present and that output files exist.

        Raises
        ------
        ValueError
            If any required field is missing.
        AssertionError
            If any output file does not exist or is not a valid path.
        AssertionError
            If any input file does not exist or is not a valid path.
        """
        required_fields = ["key", "run", "tag", "artifacts"]
        for rfield in required_fields:
            if getattr(self, rfield) is None:
                raise ValueError(
                    f"Element '{self.key}' is missing required field: {rfield}."
                )
        if self.files is not None:
            for path in self.files:
                try:
                    assert isinstance(path, Path) or isinstance(path, str)
                except AssertionError:
                    raise AssertionError(
                        f"Output file '{path}' for {self.name} does not exist or is not a valid path."  # noqa: E501
                    )
        if self.inputs is not None:
            for path in self.inputs:
                try:
                    assert isinstance(path, Path) or (isinstance(path, str))
                except AssertionError:
                    raise AssertionError(
                        f"Input file '{path}' does not exist or is not a valid path."  # noqa: E501
                    )

    ############################################################################
    # Properties
    ############################################################################

    @cached_property
    def key(self) -> str:
        """
        The unique Key to identify this Element in the registry. Identical keys
        will raise Exceptions.
        """
        return self._key

    @cached_property
    def provenance(self) -> str:
        """
        A human-readable provenance string that describes the lineage of this
        Element.
        """
        return self.build_provenance()

    @cached_property
    def signature(self) -> str:
        "The signature used to invalidate this Element and triggerr a re-run."
        return combined_signature(*self.signature_data.values())

    @cached_property
    def signature_data(self) -> Mapping[str, Any]:
        "The data used to compute the signature of this Element."
        return self.collect_sig_data()

    @property
    def validation_policy(self) -> ValidationPolicy:
        return self._validation_policy

    @validation_policy.setter
    def validation_policy(self, value: ValidationPolicy) -> None:
        self._validation_policy = value

    ############################################################################
    # Convenience
    ############################################################################
    @property
    def root(self) -> str:
        """
        The root name of this Element, derived from its tag. A short sample id,
        for example.
        """
        return self.tag.root

    @cached_property
    def short_name(self) -> str:
        """
        Generate a short name from the Element's key by taking the tag and
        tool_name from the key.

        Returns
        -------
        str
            The generated short name.
        """
        return "::".join(self.key.split("::")[:1])

    @cached_property
    def files(self) -> tuple[Path, ...]:
        """
        Return a tuple of all output files for this Element, or None if there
        are no output files. This is a convenience method that collects all
        files from the artifacts of this Element.

        Returns
        -------
        tuple[Path, ...] | None
            A tuple of all output files for this Element, or None if there are
            no output files.
        """
        files: tuple[Path, ...] = ()
        for _, v in sorted(self.artifacts.items()):
            if hasattr(v, "files"):  # covers TableArtifact, FileArtifact
                files += v.files
            elif isinstance(v, Path):  # covers Path
                files += (v.resolve(),)
        return files if files else ()

    @cached_property
    def iterfiles(self) -> Iterable[Path]:
        """
        Return an iterator over all output files for this Element.

        This is a convenience method that iterates over the files collected
        from the artifacts of this Element.

        Returns
        -------
        Iterable[Path]
            An iterable of all output files for this Element.

        Yields
        ------
        Iterator[Path]
            An iterator over all output files for this Element.
        """
        for filepath in self.files or ():
            yield filepath

    @property
    def file(self) -> Path | None:
        "Short for the primary file path, assuming there is such a thing."
        return self.artifacts["primary"].resolve()

    @property
    def primary(self) -> Any:
        "Short for the primary artifact, assuming there is such a thing."
        return self.artifacts.primary

    def __getattr__(self, name: str) -> Any:
        """
        Convenient attribute lookup that allows access to artifacts as if they
        were attributes of the Element.

        Example: element.bam will return the artifact named "bam" if it exists.

        Parameters
        ----------
        name : str
            Name of the attribute to retrieve.

        Returns
        -------
        Any
            The artifact corresponding to the given name.

        Raises
        ------
        AttributeError
            Raised if the attribute is not found among the artifacts.
        """
        artifacts = self.__dict__.get("artifacts", {})
        if name in artifacts:
            return artifacts[name]
        raise AttributeError(
            f"{self.__class__.__name__} has no attribute '{name}'"
        )  # noqa: E501

    ############################################################################
    # builder - run once
    ############################################################################

    def build_provenance(self) -> str:
        """
        Build a human-readable provenance string for this Element.

        Returns
        -------
        str
            A string that describes the lineage of this Element, including its
            name and the provenance of its prerequisites.
        """
        if not self.pres:
            prefix = ""
        elif len(self.pres) == 1:
            prefix = self.pres[0].provenance + "->"
        else:
            prefix = "(" + ",".join(p.provenance for p in self.pres) + ")->"
        return prefix + f"{self.name}"

    def collect_sig_data(self) -> Mapping[str, Any]:
        """
        Collect the data used to compute the signature of this Element. This
        includes the signature of the artifacts, the determinants, the inputs,
        the key, and the signatures of the prerequisites.

        As dicts preserve the insertion order, the returned dict is sorted by
        key to ensure a consistent order.

        Returns
        -------
        Mapping[str, Any]
            A dictionary containing the data used to compute the signature of
            this Element.
        """
        sig_data = {
            "key": self.key,
            "artifacts": self.artifacts.signature,
            "determinants": self.determinants_as_str(),
            "inputs": [file_signature(p) for p in self.inputs],
            "pre_sigs": [
                pre.signature for pre in sorted(self.pres, key=lambda p: p.key)
            ],
        }
        if isinstance(self.run, Runnable):
            sig_data["run_hash"] = self.run.signature
        return dict(sorted(sig_data.items()))

    @classmethod
    def generate_key(
        cls,
        tag: ElementTag,
        tool_name: str,
        subcommand: str | None = None,
        tool_version: str | None = None,
        suffix: str | None = None,
        **params_str: Any,
    ) -> str:
        """
        generate_element_key generates a unique key for an Element based on its tag
        and other parameters. We want the keys to be somewhat informative for
        debugging and provenance, but also unique and stable for caching.

        Parameters
        ----------
        tag : ElementTag
            The tag associated with the element.
        tool_name : str
            The name of the tool.
        tool_version : str | None, optional
            The version of the tool, by default None
        subcommand : str | None, optional
            The subcommand used, by default None
        suffix : str | None, optional
            The suffix to append to the key, by default None
        params : str | None, optional
            Additional parameters to include in the key, by default None

        Returns
        -------
        str
            The generated element key.
        """
        parts = [tag.default_name, tool_name]
        if subcommand:
            parts.append(subcommand)
        if tool_version:
            parts.append(tool_version)
        if suffix:
            parts.append(suffix)
        for k, v in params_str.items():
            if v is not None:
                parts.append(f"{k}-{v}")
        key = "::".join(parts)
        return key

    ############################################################################
    # Pipeline integration
    ############################################################################

    def __call__(self) -> Any:
        "Run the Element's run function and return the result."
        return self.run()

    def force(self) -> Element:
        "Forces a re-run of this Element, regardless of signature."
        self.validation_policy = ValidationPolicy.FORCE_RUN
        return self

    # add skip

    def is_done(
        self,
        cached_signature: str | None = None,
        cached_sig_data: dict[str, Any] | None = None,
    ) -> tuple[bool, str]:
        """
        Determine if the Element is invalidated and triggers a (re-)run.
        This is the main invalidation check that is used by the pipeline to
        decide whether to skip or run.

        This method checks the validation policy, cached signature, and output
        files to determine if the Element can be skipped.

        Parameters
        ----------
        cached_signature : str | None, optional
            Cached signature from the previous run, by default None
        cached_sig_data : dict[str, Any] | None, optional
            Cached signature data from the previous run, by default None

        Returns
        -------
        tuple[bool, str]
            Tuple where the first element is a boolean indicating if the Element
            is up-to-date, and the second element is a string explaining the
            reason for invalidating.
        """
        # TODO turn this into Validation Policy
        if self.validation_policy == ValidationPolicy.FORCE_RUN:
            return False, "Validation policy forces run"

        elif self.validation_policy == ValidationPolicy.FORCE_SKIP:
            return True, "Validation policy forces skip"

        # Has the Element been run before?
        if cached_signature is None:
            return False, "First run"

        # Check if the signature has changed
        if cached_signature != self.signature:
            return False, self._explain_signature_diff(cached_sig_data)

        # Are the output files present and non-empty?
        check_output, reason = self.outputs_ok()
        if not check_output:
            return check_output, reason

        # Is there a custom validator?
        if self.validator:
            return self.validator()

        # if no-one complains, accept the result as valid
        return True, "No relevant changes detected"

    def _explain_signature_diff(self, cached: Mapping[str, Any] | None) -> str:
        """Compare current sig_data against a previously stored one and return
        a human-readable diff.  When *cached_sig_data* is ``None`` or the two
        dicts share no keys the method falls back to a simple mismatch message.
        """
        current = self.signature_data

        if not cached:
            return "Cached signature does not match (no cached sig_data available)"  # noqa: E501
        cached_keys = list(cached.keys())
        current_keys = list(current.keys())
        missing_common = [k for k in current_keys if k not in cached_keys]
        missing_current = [k for k in cached_keys if k not in current_keys]
        if missing_common:
            return f"Cached signature does not match (missing keys in cached sig_data: {', '.join(sorted(missing_common))})"  # noqa: E501
        if missing_current:
            return f"Cached signature does not match (missing keys in current sig_data: {', '.join(sorted(missing_current))})"  # noqa: E501

        lines: list[str] = []
        for k in current_keys:
            old_v = cached.get(k, "<missing>")
            new_v = current.get(k, "<missing>")
            if stable_hash(old_v) != stable_hash(new_v):
                lines.append(f"  {k}: {old_v!r} → {new_v!r}")
        if lines:
            return "Signature mismatch in:\n" + "\n".join(lines)
        # hashes of individual keys matched but hash differs (ordering/encoding)
        return "Cached signature does not match (unknown diff)"

    def outputs_ok(self) -> tuple[bool, str]:
        """
        Check if the output files have been generated correctly. If not, this
        Element will be invalidated. This method checks for the existence and
        non-emptiness of the output files.

        Returns
        -------
        tuple[bool, str]
            Tuple containing a boolean indicating if the outputs are OK and a
            message.
        """
        missing = []
        empty = []
        check = True
        for p in self.files or ():
            if not p.exists():
                missing.append(str(p))
                check = False
        if not check:
            return check, "Missing output files: " + ", ".join(missing)
        elif not self.empty_ok:
            for p in self.files:
                if p.stat().st_size == 0 and not self.empty_ok:
                    empty.append(str(p))
                    check = False
            if not check:
                return check, "Empty output files: " + ", ".join(empty)
        return True, "Output files are OK"

    ############################################################################
    # Print
    ############################################################################

    def print_sig_data(self) -> None:
        "Print the signature data for this Element."
        print(self.signature_data)

    def determinants_as_str(self) -> str:
        """
        Return a comma-separated string of determinants. This is useful for
        creating a unique key for the Element based on its determinants and also
        calculate signatures.

        Returns
        -------
        str
            Comma-separated string of determinants.
        """
        return ",".join(str(d) for d in self.determinants)

    def describe(self) -> str:
        """
        Return a string description of the Element.

        The description includes key, threads, determinants, inputs,
        artifacts, and pres.

        Returns
        -------
        str
            A string description of the Element.
        """
        artifactlist = [f"{k}: {v}" for k, v in self.artifacts.items()]
        return (
            f"{self.__class__.__name__}:\n"
            f"  key:          {self.key}\n"
            f"  threads:      {self.threads}\n"
            f"  determinants: {', '.join(self.determinants)}\n"
            f"  inputs:       {[str(p) for p in self.inputs]}\n"
            f"  artifacts:    {{ {', '.join(artifactlist)}}}\n"
            f"  pres:         {[p.key for p in self.pres]}\n"
        )

    def __repr__(self) -> str:
        "Return a string representation of the Element."
        return (
            f"{self.__class__.__name__}("
            f"key='{self.key}', "
            f"provenience={self.provenience}, "
            f"inputs={len(self.inputs)}, "
            f"artifacts={list(self.artifacts.keys())}"
            f")"
        )

    def __hash__(self) -> int:
        "Make it hashable."
        return hash(self.key)

    def __eq__(self, other: object) -> bool:
        "Check equality based on the key."
        return isinstance(other, Element) and self.key == other.key

    # def _init_store_attributes(self):
    #     reserved = {
    #         "name",
    #         "key",
    #         "run",
    #         "validator",
    #         "determinants",
    #         "inputs",
    #         "artifacts",
    #         "pres",
    #     }  # do not overwrite by accident
    #     for attr, value in self._store_attributes.items():
    #         if attr in reserved:
    #             raise ValueError(
    #                 f"store_attributes key '{attr}' collides with reserved attribute" # noqa: E501
    #             )
    #         setattr(self, attr, value)

    # def _artifact_sig(self) -> dict[str, Any]:
    #     sig: dict[str, Any] = {}
    #     for k, v in sorted(self.artifacts.items()):
    #         if isinstance(v, TransientArtifact):
    #             continue
    #         if isinstance(v, Path):
    #             sig[k] = str(v)
    #         elif isinstance(v, (str, int, float, bool)) or v is None:
    #             sig[k] = v
    #         elif isinstance(v, Artifact):
    #             sig[k] = v.signature()
    #         else:
    #             raise TypeError(
    #                 f"artifact '{k}' has unsupported type for signature: {type(v)}" # noqa: E501
    #             )
    #     return sig


class TableElement(Element):
    """Compatibility alias for table-producing Elements used by downstream modules."""


class MappedElement(Element):

    @property
    def bam(self) -> Path:
        v = self.artifacts["bam"]
        if not isinstance(v, Path):
            raise TypeError("artifact 'bam' must be a Path")
        return v

    # @property
    # def bai(self) -> Path:
    #     v = self.artifacts["bai"]
    #     if not isinstance(v, Path):
    #         raise TypeError("artifact 'bai' must be a Path")
    #     return v


class VcfElement(Element):
    """An Element that wraps a VCF/BCF file."""

    @property
    def vcf(self) -> Path:
        v = self.artifacts["vcf"]
        if not isinstance(v, Path):
            raise TypeError("artifact 'vcf' must be a Path")
        return v

    # This does not necessarily exist
    # @property
    # def tbi(self) -> Optional[Path]:
    #     """Path to the tabix index (.tbi) if present in artifacts."""
    #     v = self.artifacts.get("tbi")
    #     return Path(v) if v is not None else None
    #     return Path(v) if v is not None else None


################################################################################
# Fastq
################################################################################


class FastqSelector:
    """
    A class to select FASTQ files from a given folder based on sample names and
    read identifiers.
    """

    def __init__(self, folder: Path | str):
        """
        Initialize the FastqSelector with the given folder.

        Parameters
        ----------
        folder : Path | str
            The folder containing FASTQ files.
        """
        self.folder = Path(folder) if isinstance(folder, str) else folder

    @classmethod
    def remove_suffixes(cls, path: Path) -> tuple[str, str]:
        """Remove all suffixes from a file path to get the base name."""
        if path.name.endswith(".gz"):
            suffixes = "".join(path.suffixes[-2:])
        else:
            suffixes = "".join(path.suffixes)
        stem = path.name.removesuffix(suffixes)
        return stem, suffixes

    def is_sample(self, suffix: str) -> bool:
        """
        Check if the given suffix corresponds to a FASTQ file.

        Parameters
        ----------
        suffix : str
            The suffix of the file to check.

        Returns
        -------
        bool
            True if the suffix corresponds to a FASTQ file, False otherwise.
        """
        return suffix in {".fastq", ".fastq.gz", ".fq.gz", ".fq"}

    def select(self, sample: str) -> Mapping[str, list[Path]]:
        """
        Select FASTQ files for the given sample.

        This method searches the folder for FASTQ files that match the sample
        name and read identifiers (R1/R2).

        Parameters
        ----------
        sample : str
            The name of the sample to select FASTQ files for.

        Returns
        -------
        Mapping[str, list[Path]]
            A mapping with keys "R1" and "R2" containing lists of Paths to the
            corresponding FASTQ files.

        Raises
        ------
        FileNotFoundError
            If no FASTQ files are found for the given sample in the folder.
        """
        files = {"R1": [], "R2": []}
        for filepath in self.folder.iterdir():
            stem, suffixes = FastqSelector.remove_suffixes(filepath)
            if filepath.is_file() and self.is_sample(suffixes):
                read, check = self.match_read(stem, sample)
                if check:
                    files[read].append(filepath.resolve())
        if len(files["R1"]) == 0 and len(files["R2"]) == 0:
            raise FileNotFoundError(
                f"No {self.__class__.__name__.capitalize()} files found for sample '{sample}' in folder '{self.folder}'"  # noqa: E501
            )
        return files

    def match_read(self, stem: str, sample: str) -> tuple[str, bool]:
        """
        Match the read identifier (R1/R2) for the given sample.

        This method checks if the stem of a FASTQ file corresponds to the given
        sample and determines whether it is an R1 or R2 read.

        Parameters
        ----------
        stem : str
            The stem of the FASTQ file (filename without suffixes).
        sample : str
            The name of the sample to match against the stem.

        Returns
        -------
        tuple[str, bool]
            A tuple containing the read identifier ("R1" or "R2") and a boolean
            indicating whether the stem matches the sample. If no match is found,
            returns ("", False).
        """
        if stem.startswith(sample):
            if "R1" in stem:
                return "R1", True
            elif "R2" in stem:
                return "R2", True
        return "", False


class NovogeneSelector(FastqSelector):
    """
    A FastqSelector for Novogene-style FASTQ files, which typically have
    filenames ending with '_1' for R1 and '_2' for R2 reads.
    """

    def match_read(self, stem: str, sample: str) -> tuple[str, bool]:
        """
        Match the read identifier (R1/R2) for Novogene-style FASTQ files.

        Parameters
        ----------
        stem : str
            The stem of the FASTQ file (filename without suffixes).
        sample : str
            The name of the sample to match against the stem.

        Returns
        -------
        tuple[str, bool]
            A tuple containing the read identifier ("R1" or "R2") and a boolean
            indicating whether the stem matches the sample. If no match is
            found, returns ("", False).
        """
        if stem.startswith(sample):
            if stem.endswith("_1"):
                return "R1", True
            elif stem.endswith("_2"):
                return "R2", True
        return "", False


class IlluminaSelector(FastqSelector):
    """
    An IlluminaSelector for Illumina-style FASTQ files, which typically have
    filenames containing '_R1_' for R1 and '_R2_' for R2 reads.
    """

    def match_read(self, stem: str, sample: str) -> tuple[str, bool]:
        """
        Match the read identifier (R1/R2) for Illumina-style FASTQ files.

        Parameters
        ----------
        stem : str
            The stem of the FASTQ file (filename without suffixes).
        sample : str
            The name of the sample to match against the stem.

        Returns
        -------
        tuple[str, bool]
            A tuple containing the read identifier ("R1" or "R2") and a boolean
            indicating whether the stem matches the sample. If no match is
            found, returns ("", False).
        """
        if stem.startswith(sample):
            if "_R1_" in stem:
                return "R1", True
            elif "_R2_" in stem:
                return "R2", True
        return "", False


class UndeterminedSelector(FastqSelector):
    """
    A FastqSelector for undetermined reads, typically labeled as "Undetermined"
    in the sample name.
    """

    def match_read(self, stem: str, sample: str = "Undetermined") -> tuple[str, bool]:
        """
        Match the read identifier (R1/R2) for undetermined FASTQ files.

        Parameters
        ----------
        stem : str
            The stem of the FASTQ file (filename without suffixes).
        sample : str, optional
            The name of the sample to match against the stem, by default
            "Undetermined".

        Returns
        -------
        tuple[str, bool]
            A tuple containing the read identifier ("R1" or "R2") and a boolean
            indicating whether the stem matches the sample. If no match is
            found, returns ("", False).
        """
        if sample in stem:
            if "_R1_" in stem:
                return "R1", True
            elif "_R2_" in stem:
                return "R2", True
        return "", False


class FastqConcat(Element):
    """
    An Element that concatenates multiple FASTQ files for a given sample into
    single R1 and R2 files. This is useful for samples that have been sequenced
    across multiple lanes or runs, where the resulting FASTQ files need to be
    combined for downstream analysis.
    """

    def __init__(
        self,
        name: str,
        folder: Path | str,
        output_folder: Path = Path("cache/fastq"),
        *,
        selector: type[FastqSelector] = NovogeneSelector,
        tag: ElementTag | PartialElementTag | None = None,
    ):
        """
        Initialize a FastqConcat element.

        This constructor sets up the FastqConcat element, which concatenates
        multiple FASTQ files for a given sample into single R1 and R2 files.

        Parameters
        ----------
        name : str
            The name of the FastqConcat element.
        folder : Path | str
            The folder containing the input FASTQ files.
        output_folder : Path, optional
            The folder to store the concatenated FASTQ files, by default
            Path("cache/fastq").
        selector : type[FastqSelector], optional
            The selector class to use for identifying FASTQ files, by default
            NovogeneSelector
        tag : ElementTag | PartialElementTag | None, optional
            An optional tag for the element, by default None
        """
        # Creates FastqArtifacts by concatenation
        self.path = Path(folder).resolve()
        self.selector = selector(self.path)
        self.output_folder = output_folder
        tag = ElementTag(
            root=name,
            level=0,
            omics=Omics.DNA,
            stage=Stage.INPUT,
            method=Method.CUSTOM,
            state=State.RAW,
        ).merge(tag)
        key = Element.generate_key(tag, "FastqConcat")
        normalized, files_to_merge = self.setup_normalization()
        artifacts = ArtifactSet(FastqArtifact(normalized["R1"], normalized.get("R2")))
        run = self.normalize(files_to_merge)

        super().__init__(
            key,
            run,
            tag,
            artifacts,
            pres=(),
        )  # Element

    def resolve_output_filename(self, read: str, suffix: str) -> str:
        if self.type == "illumina":
            return f"{self.name}_{read}_001.{suffix}"
        elif self.type == "novogene":
            return f"{self.name}_{read}.{suffix}"
        else:
            raise ValueError(f"Unsupported type: {self.type}")

    def r1(self) -> Path:
        return self.artifacts["r1"]

    def r2(self) -> Path | None:
        return self.artifacts.get("r2", None)

    def normalize(self, files_to_merge: Mapping[Path, list[Path]]) -> Runnable:
        def __run():
            for output_path, input_files in files_to_merge.items():
                FastqConcat.concat(
                    input_files=input_files,
                    output=output_path,
                )
            return True

        display = CallSpec(
            path=("FastqSource", "normalize", "__run"),
            kwargs={"name": self.name, "folder": self.path, "type": self.type},
        ).render()

        return Runnable(__run, display=display)

    def setup_normalization(self):
        files_to_concat = self.selector.select(self.name)
        if not files_to_concat:
            raise FileNotFoundError(f"No files found for sample '{self.name}'")
        normalized: dict[str, Path] = {}
        files_to_merge: dict[Path, list[Path]] = {}
        for key, tuple_of_files in files_to_concat.items():
            if not tuple_of_files:
                raise FileNotFoundError(
                    f"No files found for key '{key}' in sample '{self.name}'"
                )
            if len(tuple_of_files) == 1:
                normalized[key] = tuple_of_files[0]
            else:
                output_filename = self.resolve_output_filename(
                    read=key, suffix=tuple_of_files[0].suffix.lstrip(".")
                )
                output_path = self.output_folder / output_filename
                files_to_merge[output_path] = list(tuple_of_files)
                normalized[key] = output_path
        return normalized, files_to_merge

    @classmethod
    def concat(cls, input_files: Iterable[Path], output: Path):

        command = ["cat", *map(str, input_files)]
        parents(output)
        with output.open("wb") as out:
            subprocess.run(
                command,
                stdout=out,
                check=True,
            )


def __check_pre_sigs(current: dict[str, Any], stored: dict[str, Any]) -> str:
    cur_pre = sorted(current.get("pre_sigs", []))
    sto_pre = sorted(stored.get("pre_sigs", []))
    if cur_pre != sto_pre:
        only_current = set(cur_pre) - set(sto_pre)
        only_stored = set(sto_pre) - set(cur_pre)
        detail_lines = ["Predecessor signatures differ:"]
        if only_current:
            detail_lines.append(f"    only in current: {sorted(only_current)}")
        if only_stored:
            detail_lines.append(f"    only in stored:  {sorted(only_stored)}")
        return "\n".join(detail_lines)
    return ""


def __check_key(current: dict[str, Any], stored: dict[str, Any]) -> str:
    if current.get("key") != stored.get("key"):
        return (
            f"Keys differ:\n"
            f"    current: {current.get('key')!r}\n"
            f"    stored:  {stored.get('key')!r}"
        )
    return ""


def __check_determinants(current: dict[str, Any], stored: dict[str, Any]) -> str:
    cur_det = list(current.get("determinants", []))
    sto_det = list(stored.get("determinants", []))
    if cur_det != sto_det:
        return (
            f"Determinants differ:\n"
            f"    current: {cur_det}\n"
            f"    stored:  {sto_det}"
        )
    return ""


def __check_inputs(current: dict[str, Any], stored: dict[str, Any]) -> str:
    cur_inputs = current.get("inputs", [])
    sto_inputs = stored.get("inputs", [])
    if cur_inputs != sto_inputs:
        detail_lines = [
            f"Inputs differ: ({len(cur_inputs)} current vs {len(sto_inputs)} stored)"
        ]  # noqa: E501
        all_paths = {inp.get("path") for inp in cur_inputs + sto_inputs}
        for path in sorted(all_paths):
            cur_entry = next((i for i in cur_inputs if i.get("path") == path), None)
            sto_entry = next((i for i in sto_inputs if i.get("path") == path), None)
            if cur_entry != sto_entry:
                detail_lines.append(f"    path: {path!r}")
                detail_lines.append(f"      current: {cur_entry}")
                detail_lines.append(f"      stored:  {sto_entry}")
        return "\n".join(detail_lines)
    return ""


def __check_artifacts(current: dict[str, Any], stored: dict[str, Any]) -> str:
    cur_art = current.get("artifacts", {})
    sto_art = stored.get("artifacts", {})
    if cur_art != sto_art:
        all_artifact_keys = set(cur_art) | set(sto_art)
        detail_lines = ["Artifacts differ:"]
        for k in sorted(all_artifact_keys):
            if cur_art.get(k) != sto_art.get(k):
                detail_lines.append(f"    [{k!r}]")
                detail_lines.append(f"      current: {cur_art.get(k)!r}")
                detail_lines.append(f"      stored:  {sto_art.get(k)!r}")
        return "\n".join(detail_lines)
    return ""


def explain_signature_diff(
    element,
    store: dict[str, Any],
) -> tuple[bool, str]:
    """
    Compares the sig_data of an Element with the stored entry.

    Returns True if identical, False if different.
    Provides a detailed explanation in case of differences.
    """
    key = element.key
    current = element.sig_data()

    if key not in store:
        raise ValueError(f"Key not found in store: {key!r}")

    stored = store[key]

    diffs: list[str] = []

    diff = __check_key(current, stored)
    if diff:
        diffs.append(diff)

    # --- determinants ---
    diff = __check_determinants(current, stored)
    if diff:
        diffs.append(diff)

    # --- inputs ---
    diff = __check_inputs(current, stored)
    if diff:
        diffs.append(diff)

    # --- artifacts ---
    diff = __check_artifacts(current, stored)
    if diff:
        diffs.append(diff)

    # --- pre_sigs ---

    diff = __check_pre_sigs(current, stored)
    if diff:
        diffs.append(diff)

    # --- Result ---
    if not diffs:
        result = True
        msg = f"[✓] No differences in sig_data for {key!r}"
    else:
        result = False
        msg = f"[✗] Sig_data differs for {key!r}:"
        for diff in diffs:
            msg += f"\n{diff}"

    return result, msg


def exists(paths: tuple[Path, ...]) -> Runnable:
    def __run():
        return paths_exists(*paths)

    spec = CallSpec(path=("io", "paths_exists"), args=paths).render()
    return Runnable(__run, display=spec)


P = ParamSpec("P")
R = TypeVar("R", bound=Element)


def _as_path(x: Any) -> Path | None:
    if isinstance(x, Path):
        return x
    if isinstance(x, str) and x.strip():
        return Path(x)
    return None


def _looks_like_filepath(p: Path) -> bool:
    s = str(p)
    # relativ oder absolut mit "/" oder Windows "\" -> likely a path
    if ("/" in s) or ("\\" in s):
        return True
    # oder hat eine Endung -> likely a file
    if p.suffix:
        return True
    return False


def get_candidates(
    arts: Mapping[str, Any],
    outputs: str | Iterable[str] | None = None,
    output_files: Iterable[Path] | None = None,
    auto_outputs: bool = True,
) -> list[Any]:
    if outputs is not None:
        keys = [outputs] if isinstance(outputs, str) else list(outputs)
        candidates = [arts.get(k) for k in keys]
    elif output_files and auto_outputs:
        candidates = list(
            output_files
        )  # your Element.output_files uses artifacts Paths
    else:
        candidates = []
    return candidates


@overload
def element(fn: Callable[P, R]) -> Callable[P, R]: ...


@overload
def element(
    *,
    outputs: str | Iterable[str] | None = None,
    auto_outputs: bool = True,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def element(
    fn: Callable[P, R] | None = None,
    *,
    outputs: str | Iterable[str] | None = None,
    auto_outputs: bool = True,
) -> Callable[P, R] | Callable[[Callable[P, R]], Callable[P, R]]:
    """
        Usable as:
          @element
          def builder(...): ...
    Callable[[], CompletedProcess | None | bool | Any] | Runnable
        or:
          @element(outputs="bam")
          def builder(...): ...
    """

    def deco(func: Callable[P, R]) -> Callable[P, R]:
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            e = func(*args, **kwargs)

            # intern (optional)
            reg = current_element_registry()
            if reg is not None:
                if isinstance(e, Element):
                    e = cast(R, reg.intern(e))
                elif isinstance(e, Source) and e.producer is not None:
                    e = cast(R, reg.intern(e.producer))
            # mkdir parents for outputs (independent of registry)
            arts = e.artifacts
            output_files = e.files

            candidates = get_candidates(arts, outputs, output_files, auto_outputs)

            out_paths: list[Path] = []
            for c in candidates:
                p = _as_path(c)
                if p is None:
                    continue
                if _looks_like_filepath(p):
                    out_paths.append(p)

            if out_paths:
                parents(*out_paths)

            return e

        return wrapper

    # called as @element
    if fn is not None:
        return deco(fn)

    # called as @element(...)
    return deco


@element
def register(element: Element):
    return element
