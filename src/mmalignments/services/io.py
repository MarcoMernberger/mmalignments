"""A module to take care of I/O related services in mmalignments."""

import gzip
import json
import logging
from datetime import datetime
from logging import Logger
from pathlib import Path
from typing import Any, Callable


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


def open_target(target, *, append: bool):
    if target is None:
        return None
    if isinstance(target, Path):
        parents(target)
        mode = "a" if append else "w"
        return open(target, mode, encoding="utf-8")
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


def exists(path: Path | str) -> Callable[[], bool]:
    """Check if a file or directory exists at the given path."""

    def check():
        return Path(path).exists()

    return check
