"""
MM Alignments – demultiplexing and read processing tools.

This package provides:
- core
- models
- services
- jobs
- api
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mmalignments")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.0"

from .core import ensemblmap  # noqa: F401
from .jobs.aligner_jobs import write_lengths  # noqa: F401
from .jobs.aligner_jobs import index_genome, jobify  # noqa: F401
from .models.aligners import (BCFtools, Bedtools, BWAMem2,  # noqa: F401
                              Samtools)
from .models.callers.gatk import GATK  # noqa: F401
from .models.data import Genome, Sample  # noqa: F401
from .models.executor import Executor  # noqa: F401
from .models.qc import (FastQC, MultiQC, post_mapping_qc_with_multiqc,
                        pre_alignment_qc)
from .models.reports.report import MutationalLoadReport  # noqa: F401
from .models.tasks import Element, FileElement, ValidationPolicy  # noqa: F401
from .services import ensure, initlog  # noqa: F401

__all__ = [
    "__version__",
    "Bedtools",
    "BCFtools",
    "BWAMem2",
    "Element",
    "ensure",
    "ensure_dir",
    "ensure_dirs",
    "Executor",
    "FastQC",
    "FileElement",
    "ensemblmap",
    "GATK",
    "Genome",
    "index_genome",
    "initlog",
    "jobify",
    "MultiQC",
    "MutationalLoadReport",
    "pre_alignment_qc",
    "post_mapping_qc_with_multiqc",
    "Sample",
    "Samtools",
    "ValidationPolicy",
    "write_lengths",
    "Samtools",
    "ValidationPolicy",
    "write_lengths",
    "ValidationPolicy",
    "write_lengths",
]
    "ValidationPolicy",
    "write_lengths",
]
