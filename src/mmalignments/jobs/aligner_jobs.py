"""A module to serve as an interface to pypipegraph for alignment jobs."""

import logging
import inspect
from pathlib import Path
from typing import Mapping, Any

from pypipegraph2 import (  # type: ignore
    FileGeneratingJob,
    Job,
    MultiFileGeneratingJob,
    ParameterInvariant,
    FileInvariant,
)

from mmalignments.models.elements import Element, FilesElement
from mmalignments.services.io import parents

logger = logging.getLogger(__name__)


def jobify(
    elements: Element | tuple[Element, ...],
    job_cache: Mapping[str, Job] | None = None,
    dependencies: list[Job] | Job | None = None,
) -> list[Job]:
    """
    Create (or reuse) a pypipegraph Job for an Element. Dependencies are
    derived from runner.pres.
    extra_dependencies are optional additional job deps.

    Parameters
    ----------
    elements : Element | tuple[Element, ...]
        An Element object or a tuple f Element objects that run the task.
    job_cache : Optional[Dict[str, Job]], optional
        A dictionary to cache created jobs, by default None
    extra_dependencies : Optional[Union[List[Job]], optional
        A list of further job dependencies, by default None

    Returns
    -------
    Job
        A pypipegraph Job representing the task.

    Raises
    ------
    ValueError
        If the Element has no output files, a ValueError is raised.
    """
    elements_for_job = elements if isinstance(elements, tuple) else (elements,)

    return jobify_multi(
        elements_for_job, job_cache=job_cache, dependencies=dependencies
    )


def jobify_multi(
    elements: tuple[Element, ...],
    job_cache: Mapping[str, Job] | None = None,
    dependencies: list[Job] | Job | None = None,
) -> list[Job]:
    """
    Jobify multiple Elements and return a list of Jobs.

    Parameters
    ----------
    elements : list[Element]
        A list of Element objects to be jobified.
    job_cache : Mapping[str, Job] | None, optional
        A dictionary to cache created jobs, by default None
    dependencies : list[Job] | Job | None, optional
        A list of further job dependencies, by default None

    Returns
    -------
    list[Job]
        A list of pypipegraph Jobs representing the tasks.
    """
    jobs = []
    for element in elements:
        jobs.extend(
            jobify_single(element, job_cache=job_cache, dependencies=dependencies)
        )
    return jobs


def jobify_single(
    runner: Element,
    job_cache: Mapping[str, Job] | None = None,
    dependencies: list[Job] | Job | None = None,
) -> list[Job]:
    """
    Create (or reuse) a pypipegraph Job for an Element. Dependencies are
    derived from runner.pres.
    dependencies are optional additional job deps.

    Parameters
    ----------
    runner : Element
        An Element object that a runner that runs the task.
    job_cache : Mapping[str, Job] | None, optional
        A dictionary to cache created jobs, by default None
    dependencies : list[Job] | Job | None, optional
        A list of further job dependencies, by default None

    Returns
    -------
    Job
        A pypipegraph Job representing the task.

    Raises
    ------
    ValueError
        If the Element has no output files, a ValueError is raised.
    """
    logger.info("Jobyfing %s\n ...", runner)
    skip, reason = runner.skip()
    if skip:
        logger.info(f"Skipping %s because {reason}", runner)
        return []
    job_cache = job_cache or {}

    # reuse existing job
    if runner.key in job_cache:
        job = job_cache[runner.key]
        if dependencies is not None:
            job.depends_on(check_dependencies(dependencies))
        return job

    # build prerequisite jobs first
    pre_jobs = []
    for pre in runner.pres:
        skip, _ = pre.skip()
        if not skip:
            pre_job = jobify(pre, job_cache=job_cache)
            pre_jobs.extend(pre_job)
    jobclass = _select_appropriate_job_class(runner)
    args = _build_class_arguments(runner, jobclass)

    try:
        print(args)
        job = jobclass(*args)
        print(runner.name, jobclass.__name__, runner.skip(), runner.artifacts)
        print("")
    except TypeError as e:
        raise TypeError(
            f"Error creating job for {runner.name} with output files {runner.output_files}: {e}"
        ) from e
    if not (jobclass is FileInvariant):
        # signature invariant (skip logic in ppg form)
        signature_dep = ParameterInvariant(runner.key, runner.signature)
        job.depends_on(signature_dep)

    # depends on prerequisites
    if pre_jobs:
        job.depends_on(pre_jobs)

    # additional deps
    if dependencies is not None:
        job.depends_on(check_dependencies(dependencies))
    # cores
    _set_cores(job, runner)
    # register
    job_cache[runner.key] = job
    return job


def _select_appropriate_job_class(runner: Element) -> type:
    """Select the appropriate Job class based on the number of output files."""
    if isinstance(runner, FilesElement):
        return FileInvariant
    elif runner.output_files:
        if len(runner.output_files) == 0:
            raise NotImplementedError(
                f"Element {runner.key} has no output_files; cannot create FileGeneratingJob"
            )
        elif len(runner.output_files) == 1:
            return FileGeneratingJob
        else:
            return MultiFileGeneratingJob
    else:
        raise NotImplementedError(
            f"Element {runner.key} has no output_files; cannot determine job class"
        )


def _set_cores(job: Job, runner: Element) -> None:
    """Set cores_needed for a job based on multi_core flag."""
    job.cores_needed = -1 if getattr(runner, "multi_core", True) else 1


def check_dependencies(
    dependencies: list[Job] | Job | None = None,
) -> list[Job]:
    """
    check_dependencies ensures the dependencies list is valid.

    Parameters
    ----------
    dependencies : Optional[Union[List[Job], Job]], optional
        A list of pypipegraph Jobs that this job depends on, by default None

    Returns
    -------
    List[Job]
        A valid list of dependencies.
    """
    if dependencies is None:
        dependencies = []
    elif isinstance(dependencies, Job):
        dependencies = [dependencies]
    else:
        dependencies = list(dependencies)
    for dep in dependencies:
        if not isinstance(dep, Job):
            raise TypeError(f"Dependency {dep} is not a pypipegraph Job.")
    return dependencies


def index_genome(genome, aligner) -> Job:
    """
    index_genome creates a pypipegraph Job for genome indexing.

    Parameters
    ----------
    genome : Genome
        A Genome object representing the reference genome to be indexed.
    index_runner : Callable
        A zero-argument callable that performs the indexing when invoked.

    Returns
    -------
    Job
        A pypipegraph Job representing the genome indexing task.
    """

    def __build_index(outfile, genome=genome, aligner=aligner):
        fasta_file = genome.find_file("genome.fasta")
        out_prefix = (
            genome.prebuild_manager.prebuilt_path
            / genome.prebuild_manager.hostname
            / genome.prebuild_prefix
            / aligner.name
            / aligner.version_name
            / genome.name
            / "index"
            / "bwa_mem2_index"
        )
        aligner.index(
            fasta_file=fasta_file,
            output_prefix=out_prefix,
            cwd=None,
        )()

    return FileGeneratingJob(f"{genome.name}_index_done.txt", __build_index).depends_on(
        genome.download_genome()
    )


def create_file():
    def __write(outfile):
        with open(outfile, "w") as f:
            f.write("This is a test file.\n")

    return FileGeneratingJob("test.txt", __write)


def write_path(genome, aligner):
    def __dump(outfile, genome=genome, aligner=aligner):
        with open(outfile, "w") as f:
            f.write(str(genome.find_file("genome.fasta")) + "\n")
            f.write(str(genome.find_file("genes.gtf")) + "\n")
            f.write(
                str(
                    genome.prebuild_manager.prebuilt_path
                    / genome.prebuild_manager.hostname
                    / genome.prebuild_prefix
                    / aligner.name
                    / "genome"
                )
                + "\n"
            )
            fasta_file = genome.find_file("genome.fasta")
            out_prefix = (
                genome.prebuild_manager.prebuilt_path
                / genome.prebuild_manager.hostname
                / genome.prebuild_prefix
                / aligner.name
                / aligner.version_name
                / genome.name
                / "index"
                / "bwa_mem2_index"
            )
            f.write(str(fasta_file) + "\n")
            f.write(str(out_prefix) + "\n")

    return FileGeneratingJob("genome_path.txt", __dump).depends_on(
        genome.download_genome()
    )


def write_lengths(
    mbfgenome,
    cache_chromosome_lengths_file="cache/genome.chrom.sizes",
    canocical_chromosomes_only=True,
):
    def _write_lengths(
        cache_chromosome_lengths_file,
        mbfgenome=mbfgenome,
        canocical_chromosomes_only=canocical_chromosomes_only,
    ):
        with open(cache_chromosome_lengths_file, "w") as f:
            canocical_chromosomes = mbfgenome.get_true_chromosomes()
            for chrom, length in mbfgenome.get_chromosome_lengths().items():
                if canocical_chromosomes_only and chrom not in canocical_chromosomes:
                    continue
                f.write(f"{chrom}\t{length}\n")

    return FileGeneratingJob(cache_chromosome_lengths_file, _write_lengths).depends_on(
        mbfgenome.download_genome()
    )


def get_job_func(runner: Element) -> Any:
    """Get the function to run for a given Element."""
    if hasattr(runner, "run") and callable(runner.run):

        def func(*args, runner=runner, **kwargs):
            return runner.run()

        return func
    else:
        raise NotImplementedError(f"Element {runner.key} has no callable 'run' method.")


def _build_class_arguments(runner: Element, jobclass: type) -> Any:
    """Convert runner.output_files to relative paths."""
    if jobclass is FileGeneratingJob:
        if runner.output_files is None:
            raise ValueError(
                f"Element {runner.key} has no output_files; cannot create FileGeneratingJob"
            )
        args = _relative_path(runner.output_files[0]), get_job_func(runner)

    elif jobclass is MultiFileGeneratingJob:
        if runner.output_files is None:
            raise ValueError(
                f"Element {runner.key} has no output_files; cannot create MultiFileGeneratingJob"
            )
        args = [_relative_path(f) for f in runner.output_files], get_job_func(runner)
    elif jobclass is FileInvariant:
        return runner.inputs[0]
    else:
        raise NotImplementedError(f"Unsupported job class: {jobclass}")
    return args


def _relative_path(filepath: Path) -> Path:
    filepath = filepath.resolve()
    try:
        res = filepath.relative_to(Path.cwd())
    except ValueError:
        cwd = Path.cwd()
        res = Path(cwd.name) / ".." / filepath.relative_to(cwd.anchor)
    return res
