"""
Job wrappers for running MM Alignments inside pipelines
(e.g. pypipegraph, anysnake, Snakemake wrappers, etc.).
"""

from .aligner_jobs import write_lengths, write_path, create_file

__all__ = [
    "write_lengths",
    "write_path",
    "create_file",
]
