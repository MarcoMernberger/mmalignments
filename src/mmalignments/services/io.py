"""A module to take care of I/O related services in mmalignments."""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Callable, Iterable, Mapping

import pandas as pd  # type: ignore[import]
import pyarrow.parquet as pq  # type: ignore[import]
from matplotlib.figure import Figure  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

if TYPE_CHECKING:
    from artifacts import TableArtifact  # type: ignore[import]


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
    """
    Load a JSON file and return its contents as a dictionary.

    Parameters
    ----------
    infile : Path
        The path to the JSON file to load.
    encoding : str, optional
        The encoding to use when reading the file, by default "utf-8"

    Returns
    -------
    dict[str, str | int | float]
        The contents of the JSON file as a dictionary.

    Raises
    ------
    FileNotFoundError
        If the JSON file does not exist.
    """
    if not infile.exists():
        raise FileNotFoundError(f"JSON file not found: {infile}")

    with open(infile, "r", encoding=encoding) as f:
        data = json.load(f)
        return data


def load_param_json(path: Path) -> dict[str, Any]:
    """
    Load a parameter JSON file and return its contents as a dictionary.

    Parameters
    ----------
    path : Path
        The path to the parameter JSON file.

    Returns
    -------
    dict[str, Any]
        The contents of the parameter JSON file as a dictionary.

    Raises
    ------
    ValueError
        If the top-level JSON is not an object.
    FileNotFoundError
        If the parameter JSON file does not exist.
    """
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
    """
    Open a target file or return an existing file-like object.

    Parameters
    ----------
    append : bool
        Whether to append to the file if it exists.
    target : Path | None | IO, optional
        The target file path or file-like object, by default None
    encoding : str, optional
        The encoding to use when opening the file, by default "utf-8"
    """
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
    return tuple(Path(p).resolve() for p in paths)


def paths_exists(*paths: Path | TableArtifact) -> Callable[[], bool]:
    """Check if all given paths exist."""

    def check():
        return all(p.resolve().exists() for p in paths)

    return check


def paths_exists_raise(*paths: Path | TableArtifact) -> Callable[[], None]:
    """Check if all given paths exist."""

    def check():
        for p in paths:
            if not p.resolve().exists():
                raise FileNotFoundError(f"Required path does not exist: {p}")

    return check


def exists(path: Path | str) -> Callable[[], bool]:
    """Check if a file or directory exists at the given path."""

    def check():
        return Path(path).resolve().exists()

    return check


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    """
    Write sequences to a FASTA file.

    Parameters
    ----------
    path : Path
        The path to the output FASTA file.
    sequences : dict[str, str]
        A dictionary mapping sequence names to sequences.
    """
    parents(path)
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


def write_fasta_from_list(path: Path, names: list[str], sequences: list[str]) -> None:
    """
    Write sequences to a FASTA file.

    Parameters
    ----------
    path : Path
        The path to the output FASTA file.
    names : list[str]
        A list of sequence names.
    sequences : list[str]
        A list of sequences corresponding to the names.
    """
    parents(path)
    with open(path, "w") as f:
        for name, seq in zip(names, sequences):
            f.write(f">{name}\n{seq}\n")


def write_frames(
    df: DataFrame,
    paths: Iterable[Path],
    **kwargs,
) -> None:
    """
    Write a DataFrame to a file.

    Parameters
    ----------
    df : DataFrame
        The DataFrame to write.
    path : Path
        The path to the file to write.
    mode : Literal["tsv", "parquet", "both"], optional
        The mode to write the file, by default "both".

    Raises
    ------
    ValueError
        If the file format is unsupported.
    """
    for path in paths:
        write_frame(df, path, **kwargs)


def write_frame(df: DataFrame, path: Path, **kwargs) -> None:
    """
    Write a DataFrame to a file.

    Parameters
    ----------
    df : DataFrame
        The DataFrame to write.
    path : Path
        The path to the file to write.

    Raises
    ------
    ValueError
        If the file format is unsupported.
    """
    ext = path.suffix.lower()
    parents(path)
    params = {"sep": "\t", "index": True}
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


def read_frame(path: Path, drop_unnamed_columns: bool = True, **kwargs) -> DataFrame:
    """
    Read a DataFrame from a file.

    Parameters
    ----------
    path : Path
        The path to the file to read.

    Returns
    -------
    DataFrame
        The DataFrame read from the file.

    Raises
    ------
    ValueError
        If the file format is unsupported.
    """
    if drop_unnamed_columns:
        kwargs.setdefault("usecols", lambda c: not c.startswith("Unnamed"))
    if path.suffix in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t", **kwargs)
    elif path.suffix == ".parquet":
        return pd.read_parquet(path, **kwargs)
    elif path.suffix in (".csv",):
        return pd.read_csv(path, **kwargs)
    elif path.suffix in (".xlsx", ".xls"):
        return pd.read_excel(path, **kwargs)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


def read_schema(path: Path, **kwargs) -> list[str]:
    """
    Read a DataFrame schema from a file.

    Parameters
    ----------
    path : Path
        The path to the file to read.

    Returns
    -------
    list[str]
        The list of column names from the DataFrame.

    Raises
    ------
    ValueError
        If the file format is unsupported.
    """
    if path.suffix in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t", nrows=0, **kwargs).columns.tolist()
    elif path.suffix == ".parquet":
        return pq.read_schema(path).names
    elif path.suffix in (".csv",):
        return pd.read_csv(path, nrows=0, **kwargs).columns.tolist()
    elif path.suffix in (".xlsx", ".xls"):
        return pd.read_excel(path, nrows=0, **kwargs).columns.tolist()
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


def concat_files(output_file: Path, *input_files: Path) -> None:
    """
    Concatenate multiple input files into a single output file.

    Parameters
    ----------
    output_file : Path
        The path to the output file.
    input_files : Path
        One or more paths to the input files to concatenate.
    """
    with open(output_file, "w") as out_f:
        for input_file in input_files:
            with open(input_file, "r") as in_f:
                for line in in_f:
                    out_f.write(line)


def get_paths_from_prefix_path(path: str | Path) -> Mapping[str, Path]:
    """
    Get all paths in the same directory as the given path that start with the
    same prefix.

    Parameters
    ----------
    path : str | Path
        The path to use as the prefix.

    Returns
    -------
    Mapping[str, Path]
        A mapping of file stem to absolute Path for all matching files.
    """
    paths = {}
    p = Path(path)
    file_dir = p.parent
    prefix = p.stem
    for p in file_dir.iterdir():
        if p.stem.startswith(prefix):
            paths[p.stem] = p.absolute()
    return paths


def concat_fastq(inputs: tuple[Path, ...], output: Path) -> None:
    """
    Concatenate multiple FASTQ files into a single output FASTQ file.

    Parameters
    ----------
    inputs : tuple[Path, ...]
        The input FASTQ files to concatenate.
    output : Path
        The output FASTQ file.

    """
    output.parent.mkdir(parents=True, exist_ok=True)

    files = [str(p) for p in inputs]

    with output.open("wb") as out:
        subprocess.run(
            ["cat", *files],
            stdout=out,
            check=True,
        )


def save_figure(fig: Figure, outfile: Path) -> None:
    """
    Save a matplotlib figure to a file.

    Parameters
    ----------
    fig : Figure
        The matplotlib figure to save.
    outfile : Path
        The path to the output file.
    """
    parents(outfile)

    fig.savefig(
        str(outfile),
        format=outfile.suffix[1:],
        bbox_inches="tight",
    )
