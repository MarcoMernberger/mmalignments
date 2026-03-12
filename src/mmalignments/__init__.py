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

from mmalignments.jobs.aligner_jobs import (  # noqa: F401
    index_genome,
    jobify,
    write_lengths,  # noqa: F401
)
from mmalignments.models.aligners import (
    BCFtools,
    Bedtools,  # noqa: F401
    BWAMem2,
    Samtools,
)
from mmalignments.models.callers.gatk import GATK  # noqa: F401
from mmalignments.models.data import Genome, Sample, Pairing  # noqa: F401
from mmalignments.models.elements import (
    Element,
    FileElement,  # noqa: F401
    FilesElement,
    NextGenSampleElement,  # noqa: F401
    ValidationPolicy,
)
from mmalignments.models.executor import Executor  # noqa: F401
from mmalignments.models.qc import (
    FastQC,
    MultiQC,
    post_mapping_qc_with_multiqc,
    pre_alignment_qc,
)
from mmalignments.models.reports.report import MutationalLoadReport  # noqa: F401
from mmalignments.services import ensure, initlog  # noqa: F401

from .core import ensemblmap  # noqa: F401
from mmalignments.models.parameters import Params  # noqa: F401
from mmalignments.models.externals import ExternalRunConfig  # noqa: F401
from mmalignments.models.tags import ElementTag  # noqa: F401

__all__ = [
    "Bedtools",
    "BCFtools",
    "BWAMem2",
    "Element",
    "ElementTag",
    "ExternalRunConfig",
    "FilesElement",
    "FileElement",
    "ensure",
    "Executor",
    "FastQC",
    "FilesElement",
    "ensemblmap",
    "GATK",
    "Genome",
    "index_genome",
    "initlog",
    "jobify",
    "MultiQC",
    "MutationalLoadReport",
    "NextGenSampleElement",
    "Params",
    "pre_alignment_qc",
    "post_mapping_qc_with_multiqc",
    "Sample",
    "Samtools",
    "ValidationPolicy",
    "__version__",
    "write_lengths",
]
