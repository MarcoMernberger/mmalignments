"""Test samtools wrapper functionality."""

import subprocess
from pathlib import Path
from mmalignments.models.aligners.samtools import Samtools


def test_samtools_sort_callable_signature(monkeypatch):
    """Verify that sort() returns a zero-argument callable."""
    st = Samtools()

    # Mock External.run to return a fake callable
    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(Samtools, "run", fake_run)

    sort_callable = st.sort(input_bam="in.bam", output_bam="out.bam")
    assert callable(sort_callable)

    # Invoke should not raise
    sort_callable()


def test_samtools_index_callable_signature(monkeypatch):
    """Verify that index() returns a zero-argument callable."""
    st = Samtools()

    # Mock External.run to return a fake callable
    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(Samtools, "run", fake_run)

    index_callable = st.index(bam_file="sorted.bam")
    assert callable(index_callable)

    # Invoke should not raise
    index_callable()


def test_samtools_get_version():
    """Check that get_version returns a string or fallback."""
    st = Samtools()
    version = st.get_version(fallback="unknown")
    # If samtools is installed, version should be a string; otherwise fallback
    assert isinstance(version, (str, type(None)))


def test_chained_post_with_bwamem2_and_samtools(monkeypatch):
    """Verify that sort and index callables work in a post list."""
    from mmalignments.models.aligners.bwamem2 import BWAMem2

    calls = []

    def fake_bwa_run(*args, **kwargs):
        # BWAMem2.align wraps post list into a single post(cp) callable
        post_wrapper = kwargs.get("post")
        if post_wrapper:
            post_wrapper(subprocess.CompletedProcess(args=["cmd"], returncode=0))
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    def fake_samtools_run(*args, **kwargs):
        # Samtools.sort and .index wrap External.run into a zero-arg callable
        def runner():
            calls.append(kwargs.get("extra_arguments", ["unknown"])[0])
            return subprocess.CompletedProcess(args=["cmd"], returncode=0)

        return runner

    monkeypatch.setattr(BWAMem2, "run", fake_bwa_run)
    monkeypatch.setattr(Samtools, "run", fake_samtools_run)

    aligner = BWAMem2()
    st = Samtools()

    sort_callable = st.sort(input_bam="aligned.bam", output_bam="sorted.bam")
    index_callable = st.index(bam_file="sorted.bam")

    # Pass both callables to align as post list
    runner = aligner.align(
        index_prefix="/tmp/idx",
        fastq_r1="r1.fq",
        output_bam="aligned.bam",
        post=[sort_callable, index_callable],
    )

    # Verify that both samtools commands were invoked in order
    assert calls == ["sort", "index"]
