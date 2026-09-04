from typing import Protocol, Mapping
from pathlib import Path
from mmalignments.models.elements import (
    Element,
)
from mmalignments.models.overlay import (
    ElementTag,
    PartialElementTag,
)
from mmalignments.models.artifacts import (
    ArtifactSet,
    FastqArtifact,
)
from mmalignments.models.tags import (
    Omics,
    Stage,
    Method,
    State,
)
from mmalignments.models.specs import CallSpec, Runnable
from typing import Iterable
import subprocess


def remove_suffixes(path: Path) -> tuple[str, str]:
    """Remove all suffixes from a file path to get the base name.

    Shared by all selectors so each one doesn't need its own copy.
    """
    suffixes = "".join(path.suffixes)
    stem = path.name.removesuffix(suffixes)
    return stem, suffixes


class FastqSelector(Protocol):
    def __call__(self, folder: Path) -> Mapping[str, list[Path]]: ...

    def is_sample(self, suffix: str) -> bool:
        return suffix in {".fastq", ".fastq.gz", ".fq.gz", ".fq"}

    def output_filename(self, name: str, read: str, suffix: str) -> str:
        """Default naming convention for merged output files.

        Selectors that need a different convention (e.g. Illumina's
        `_001` suffix) override this.
        """
        return f"{name}_{read}.{suffix}"


class NovogeneSelector(FastqSelector):
    def __call__(self, folder: Path) -> Mapping[str, list[Path]]:
        files: dict[str, list[Path]] = {"R1": [], "R2": []}
        if not folder.exists():
            raise FileNotFoundError(f"Sample path {folder} does not exist.")
        # NOTE: iterdir() does not recurse. If novogene deliveries are
        # ever nested in subfolders, switch this to folder.rglob("*").
        for filepath in folder.iterdir():
            if filepath.is_file():
                stem, suffixes = remove_suffixes(filepath)
                if self.is_sample(suffixes):
                    if stem.endswith("_1"):
                        files["R1"].append(filepath.resolve())
                    elif stem.endswith("_2"):
                        files["R2"].append(filepath.resolve())
        return files


class IlluminaSelector(FastqSelector):
    def __call__(self, folder: Path) -> Mapping[str, list[Path]]:
        files: dict[str, list[Path]] = {"R1": [], "R2": []}
        for filepath in folder.iterdir():
            stem, suffixes = remove_suffixes(filepath)
            if filepath.is_file() and self.is_sample(suffixes):
                if "_R1_" in stem:
                    files["R1"].append(filepath.resolve())
                elif "_R2_" in stem:
                    files["R2"].append(filepath.resolve())
        return files

    def output_filename(self, name: str, read: str, suffix: str) -> str:
        return f"{name}_{read}_001.{suffix}"


class UndeterminedSelector(FastqSelector):
    def __call__(self, folder: Path) -> Mapping[str, list[Path]]:
        files: dict[str, list[Path]] = {"R1": [], "R2": []}
        for filepath in folder.iterdir():
            stem, suffixes = remove_suffixes(filepath)
            if filepath.is_file() and self.is_sample(suffixes):
                if "_R1_" in stem and "Undetermined" in stem:
                    files["R1"].append(filepath.resolve())
                elif "_R2_" in stem and "Undetermined" in stem:
                    files["R2"].append(filepath.resolve())
        return files

    def output_filename(self, name: str, read: str, suffix: str) -> str:
        return f"{name}_{read}_001.{suffix}"


class FastqConcat(Element):

    def __init__(
        self,
        root: str,
        folder: Path | str,
        output_folder: Path = Path("cache/fastq"),
        *,
        selector: FastqSelector = NovogeneSelector(),
        tag: ElementTag | PartialElementTag | None = None,
    ):
        self.path = Path(folder).resolve()
        self.selector = selector
        self.output_folder = output_folder
        fulltag = ElementTag(
            root=root,
            level=0,
            omics=Omics.DNA,
            stage=Stage.INPUT,
            method=Method.CUSTOM,
            state=State.RAW,
        ).patch(tag)
        key = Element.generate_key(fulltag, "FastqSource", subcommand="concat")
        normalized, files_to_merge = self.setup_normalization()
        artifacts = ArtifactSet(FastqArtifact(normalized["R1"], normalized.get("R2")))
        self.normalized = normalized
        self.files_to_merge = files_to_merge

        super().__init__(
            key,
            self.concat_runnable(),
            fulltag,
            artifacts=artifacts,
            pres=(),
        )

    def r1(self) -> Path:
        return self.artifacts["r1"]

    def r2(self) -> Path | None:
        return self.artifacts.get("r2", None)

    def setup_normalization(self):
        files_to_merge: dict[Path, list[Path]] = {}
        normalized: dict[str, Path] = {}

        input_files_dict = self.selector(self.path)
        if not input_files_dict:
            raise FileNotFoundError(
                f"No files found for sample '{self.key}' in folder '{self.path}' "
                f"using {type(self.selector).__name__}"
            )
        for key, tuple_of_files in input_files_dict.items():
            if not tuple_of_files:
                raise FileNotFoundError(
                    f"No files found for key '{key}' in sample '{self.key}'"
                )
            if len(tuple_of_files) == 1:
                normalized[key] = tuple_of_files[0]
            else:
                # NOTE: Path.suffix only returns the LAST suffix, so for
                # "sample.fastq.gz" this yields "gz", not "fastq.gz".
                # Pre-existing behavior, carried over as-is; flagging in
                # case downstream naming depends on the full suffix.
                output_filename = self.selector.output_filename(
                    self.name, key, tuple_of_files[0].suffix.lstrip(".")
                )
                output_path = self.output_folder / output_filename
                files_to_merge[output_path] = list(tuple_of_files)
                normalized[key] = output_path
        return normalized, files_to_merge

    def concat_runnable(self) -> Runnable:
        def __run():
            for output_path, input_files in self.files_to_merge.items():
                self._concat(input_files, output_path)
            return True

        display = CallSpec(
            path=("FastqSource", "concat", "__run"),
            kwargs={
                "root": self.root,
                "folder": self.path,
                "selector": type(self.selector).__name__,
            },
        ).render()

        return Runnable(__run, display=display)

    @staticmethod
    def _concat(input_files: Iterable[Path], output: Path) -> None:
        command = ["cat", *map(str, input_files)]
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("wb") as out:
            subprocess.run(command, stdout=out, check=True)
