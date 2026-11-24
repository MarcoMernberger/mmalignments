"""
MM Alignments – demultiplexing and read processing tools.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mmalignments")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.0"
