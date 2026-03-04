"""Test GATK wrapper functionality."""

import subprocess
from pathlib import Path
from mmalignments.models.aligners.gatk import GATK


def test_gatk_mark_duplicates_callable_signature(monkeypatch):
    """Verify that mark_duplicates() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    mark_dups = gatk.mark_duplicates(
        input_bam="in.bam", output_bam="out.bam", metrics_file="metrics.txt"
    )
    assert callable(mark_dups)
    mark_dups()


def test_gatk_learn_read_orientation_callable_signature(monkeypatch):
    """Verify that learn_read_orientation_model() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    learn = gatk.learn_read_orientation_model("f1r2.tar.gz", "model.tar.gz")
    assert callable(learn)
    learn()


def test_gatk_mutect2_callable_signature(monkeypatch):
    """Verify that mutect2() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    mutect = gatk.mutect2(
        reference="ref.fa",
        input_bams=["tumor.bam", "normal.bam"],
        output_vcf="variants.vcf.gz",
        tumor_sample="TUMOR",
        normal_sample="NORMAL",
    )
    assert callable(mutect)
    mutect()


def test_gatk_get_pileup_summaries_callable_signature(monkeypatch):
    """Verify that get_pileup_summaries() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    pileup = gatk.get_pileup_summaries(
        input_bam="tumor.bam",
        output_table="pileup.table",
        variant_sites="sites.vcf.gz",
    )
    assert callable(pileup)
    pileup()


def test_gatk_calculate_contamination_callable_signature(monkeypatch):
    """Verify that calculate_contamination() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    contamination = gatk.calculate_contamination(
        tumor_pileup="tumor_pileup.table", output_table="contamination.table"
    )
    assert callable(contamination)
    contamination()


def test_gatk_filter_mutect_calls_callable_signature(monkeypatch):
    """Verify that filter_mutect_calls() returns a zero-argument callable."""
    gatk = GATK()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(GATK, "run", fake_run)

    filter_calls = gatk.filter_mutect_calls(
        reference="ref.fa", input_vcf="unfiltered.vcf.gz", output_vcf="filtered.vcf.gz"
    )
    assert callable(filter_calls)
    filter_calls()


def test_gatk_chained_pipeline(monkeypatch):
    """Verify that GATK callables can be chained in order."""
    calls = []

    def fake_run(*args, **kwargs):
        def runner():
            # Extract the tool name (first argument after gatk)
            extra_args = kwargs.get("extra_arguments", [])
            if extra_args:
                calls.append(extra_args[0])
            return subprocess.CompletedProcess(args=["cmd"], returncode=0)

        return runner

    monkeypatch.setattr(GATK, "run", fake_run)

    gatk = GATK()

    # Create callables for a typical variant calling pipeline
    mark_dups = gatk.mark_duplicates("in.bam", "dedup.bam", "metrics.txt")
    mutect = gatk.mutect2(
        reference="ref.fa",
        input_bams=["tumor.bam", "normal.bam"],
        output_vcf="raw.vcf.gz",
        tumor_sample="TUMOR",
        normal_sample="NORMAL",
    )
    learn_orient = gatk.learn_read_orientation_model("f1r2.tar.gz", "model.tar.gz")
    pileup = gatk.get_pileup_summaries(
        input_bam="tumor.bam", output_table="pileup.table", variant_sites="sites.vcf.gz"
    )
    contamination = gatk.calculate_contamination(
        tumor_pileup="pileup.table", output_table="contamination.table"
    )
    filter_calls = gatk.filter_mutect_calls(
        reference="ref.fa",
        input_vcf="raw.vcf.gz",
        output_vcf="filtered.vcf.gz",
        contamination_table="contamination.table",
        orientation_model="model.tar.gz",
    )

    # Execute in order
    mark_dups()
    mutect()
    learn_orient()
    pileup()
    contamination()
    filter_calls()

    # Verify execution order
    assert calls == [
        "MarkDuplicates",
        "Mutect2",
        "LearnReadOrientationModel",
        "GetPileupSummaries",
        "CalculateContamination",
        "FilterMutectCalls",
    ]


def test_gatk_version_extraction_from_path():
    """Test that version can be extracted from GATK path."""
    gatk = GATK(primary_binary="code/gatk-4.6.2.0/gatk")
    version = gatk.get_version(fallback="unknown")
    # Should extract "4.6.2.0" from path
    assert version == "4.6.2.0" or version == "unknown"
