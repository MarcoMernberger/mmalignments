"""A module to take care of I/O related services in mmalignments."""

import gzip
import json
from pathlib import Path
from typing import IO, Any, Callable

import pandas as pd  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]


def ensure(*files: Path | str) -> bool:
    """
    ensure_dir ensures that the directory for the given path exists.

    Parameters
    ----------
    file : Path | str
        The folder path to ensure it exists.

    Returns
    -------
    bool
        True if the directory was successfully ensured, False otherwise.
    """
    ret = True
    for path in files:
        ret &= ensure_path(path)
    return ret


def parents(*files: Path | str) -> bool:
    """
    ensure_dir ensures that the directory for the given path exists.

    Parameters
    ----------
    file : Path | str
        The folder path to ensure it exists.

    Returns
    -------
    bool
        True if the directory was successfully ensured, False otherwise.
    """
    return ensure(*(Path(f).parent for f in files))


def ensure_path(path: Path | str) -> bool:
    """
    ensure_dir ensures that the directory for the given path exists.

    Parameters
    ----------
    path : Path | str
        The folder path to ensure it exists.

    Returns
    -------
    bool
        True if the directory was successfully ensured, False otherwise.
    """
    path = Path(path)
    if not (isinstance(path, Path) or isinstance(path, str)):
        raise ValueError(f"Expected Path or str, got {type(path)}")
    try:
        path.mkdir(parents=True, exist_ok=True)
        return True
    except Exception:
        return False


def open_fastq(path):
    """Open FASTQ file (gzipped or plain)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def write_json(path: Path, results: dict[str, str]) -> None:
    with open(path, "w") as f:
        json.dump(results, f, indent=2)


def write_fastq_check_results(
    path: Path, sample_name: str, pairing: str, results: dict
) -> None:
    with open(path, "w") as f:
        f.write(f"FASTQ Input Check for {sample_name}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Status: {results['status']}\n")
        f.write(f"Pairing: {pairing}\n\n")

        f.write("Checks:\n")
        for key, value in results["checks"].items():
            f.write(f"  {key}: {value}\n")

        if results["errors"]:
            f.write("\nErrors:\n")
            for err in results["errors"]:
                f.write(f"  - {err}\n")

        if results["warnings"]:
            f.write("\nWarnings:\n")
            for warn in results["warnings"]:
                f.write(f"  - {warn}\n")


def from_json(infile: Path, encoding="utf-8") -> dict[str, str | int | float]:
    if not infile.exists():
        raise FileNotFoundError(f"JSON file not found: {infile}")

    with open(infile, "r", encoding=encoding) as f:
        data = json.load(f)
        return data


def load_param_json(path: Path) -> dict[str, Any]:
    try:
        obj = from_json(path)
        if not isinstance(obj, dict):
            raise ValueError(
                "Top-level JSON must be an object: {subroutine: {param: spec}}"
            )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Param spec JSON not found: {path}") from e
    return obj


def open_target(
    target: Path | None | IO = None, *, append: bool, encoding: str = "utf-8"
):
    if target is None:
        return None
    if isinstance(target, Path):
        parents(target)
        mode = "a" if append else "w"
        return open(target, mode, encoding=encoding)

    return target


def absolutize(*paths: Path | str) -> tuple[Path, ...]:
    """Convert one or more paths to absolute Path objects.

    Parameters
    ----------
    *paths : Path | str
        One or more paths to convert.

    Returns
    -------
    tuple[Path, ...]
        List of absolute Path objects corresponding to the input paths.
    """
    return tuple(Path(p).absolute() for p in paths)


def paths_exists(*paths: Path | str) -> Callable[[], bool]:
    """Check if all given paths exist."""

    def check():
        return all(Path(p).exists() for p in paths)

    return check


def paths_exists_raise(*paths: Path | str) -> Callable[[], None]:
    """Check if all given paths exist."""

    def check():
        for p in paths:
            if not Path(p).exists():
                raise FileNotFoundError(f"Required path does not exist: {p}")

    return check


def exists(path: Path | str) -> Callable[[], bool]:
    """Check if a file or directory exists at the given path."""

    def check():
        return Path(path).exists()

    return check


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    parents(path)
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


def write_frame(df: DataFrame, path: Path, **kwargs) -> None:
    ext = path.suffix.lower()
    parents(path)
    params = {"sep": "\t", "index": False}
    params.update(kwargs)
    if ext in (".tsv", ".csv", ".txt"):
        df.to_csv(path, **params)
    elif ext == ".parquet":
        params.pop("sep", None)
        df.to_parquet(path, **params)
    elif ext in (".xlsx", ".xls"):
        params.pop("sep", None)
        df.to_excel(path, **params)
    else:
        raise ValueError(f"Unsupported file format: {ext}")


def read_frame(path: Path, **kwargs) -> DataFrame:

    if path.suffix == ".tsv":
        return pd.read_csv(path, sep="\t", **kwargs)
    elif path.suffix == ".parquet":
        return pd.read_parquet(path, **kwargs)
    elif path.suffix in (".csv", ".txt"):
        return pd.read_csv(path, **kwargs)
    elif path.suffix in (".xlsx", ".xls"):
        return pd.read_excel(path, **kwargs)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
