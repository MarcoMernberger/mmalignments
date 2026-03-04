"""Test FastQC and MultiQC wrapper functionality."""

import subprocess
from pathlib import Path

from mmalignments.models.qc.fastqc import FastQC
from mmalignments.models.qc.multiqc import MultiQC


def test_fastqc_run_qc_single_file(monkeypatch):
    """Verify that run_qc() returns a zero-argument callable for single file."""
    fastqc = FastQC()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(FastQC, "run", fake_run)

    qc = fastqc.run_qc(input_file="sample.bam", output_dir="qc_results")
    assert callable(qc)
    qc()


def test_fastqc_run_qc_multiple_files(monkeypatch):
    """Verify that run_qc() handles multiple files."""
    fastqc = FastQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["extra_args"] = kwargs.get("extra_arguments", [])
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(FastQC, "run", fake_run)

    qc = fastqc.run_qc(
        input_file=["sample1.bam", "sample2.bam"], output_dir="qc_results"
    )
    assert callable(qc)
    qc()

    # Check that both files were added to command
    extra_args = captured.get("extra_args", [])
    assert "sample1.bam" in str(extra_args)
    assert "sample2.bam" in str(extra_args)


def test_fastqc_run_qc_with_extract(monkeypatch):
    """Verify that extract flag is added when requested."""
    fastqc = FastQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["extra_args"] = kwargs.get("extra_arguments", [])
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(FastQC, "run", fake_run)

    qc = fastqc.run_qc(input_file="sample.bam", output_dir="qc", extract=True)
    qc()

    # Check that --extract flag was added
    extra_args = captured.get("extra_args", [])
    assert "--extract" in extra_args


def test_fastqc_run_qc_with_threads(monkeypatch):
    """Verify that threads parameter is passed correctly."""
    fastqc = FastQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["params"] = kwargs
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(FastQC, "run", fake_run)

    qc = fastqc.run_qc(input_file="sample.bam", output_dir="qc", threads=8)
    qc()

    # Check that --threads parameter was set
    params = captured.get("params", {})
    assert params.get("--threads") == 8


def test_multiqc_run_single_dir(monkeypatch):
    """Verify that run() returns a zero-argument callable for single directory."""
    multiqc = MultiQC()

    def fake_run(*args, **kwargs):
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(MultiQC, "run", fake_run)

    report = multiqc.run(analysis_dir="qc_results", output_dir="reports")
    assert callable(report)
    report()


def test_multiqc_run_multiple_dirs(monkeypatch):
    """Verify that run() handles multiple analysis directories."""
    multiqc = MultiQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["extra_args"] = kwargs.get("extra_arguments", [])
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(MultiQC, "run", fake_run)

    report = multiqc.run(analysis_dir=["qc_results", "other_qc"], output_dir="reports")
    assert callable(report)
    report()

    # Check that both directories were added
    extra_args = captured.get("extra_args", [])
    assert "qc_results" in str(extra_args)
    assert "other_qc" in str(extra_args)


def test_multiqc_run_with_report_name(monkeypatch):
    """Verify that report name is passed correctly."""
    multiqc = MultiQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["extra_args"] = kwargs.get("extra_arguments", [])
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(MultiQC, "run", fake_run)

    report = multiqc.run(
        analysis_dir="qc_results", output_dir="reports", report_name="my_report"
    )
    report()

    # Check that -n and report name were added
    extra_args = captured.get("extra_args", [])
    assert "-n" in extra_args
    assert "my_report" in extra_args


def test_multiqc_run_with_force(monkeypatch):
    """Verify that force flag is added when requested."""
    multiqc = MultiQC()

    captured = {}

    def fake_run(*args, **kwargs):
        captured["extra_args"] = kwargs.get("extra_arguments", [])
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(MultiQC, "run", fake_run)

    report = multiqc.run(analysis_dir="qc_results", output_dir="reports", force=True)
    report()

    # Check that --force flag was added
    extra_args = captured.get("extra_args", [])
    assert "--force" in extra_args


def test_qc_pipeline_integration(monkeypatch):
    """Verify that FastQC and MultiQC can be chained."""
    calls = []

    def fake_fastqc_run(*args, **kwargs):
        def runner():
            calls.append("fastqc")
            return subprocess.CompletedProcess(args=["cmd"], returncode=0)

        return runner

    def fake_multiqc_run(*args, **kwargs):
        def runner():
            calls.append("multiqc")
            return subprocess.CompletedProcess(args=["cmd"], returncode=0)

        return runner

    monkeypatch.setattr(FastQC, "run", fake_fastqc_run)
    monkeypatch.setattr(MultiQC, "run", fake_multiqc_run)

    fastqc = FastQC()
    multiqc = MultiQC()

    # Create callables for QC pipeline
    qc1 = fastqc.run_qc("sample1.bam", "qc_results")
    qc2 = fastqc.run_qc("sample2.bam", "qc_results")
    report = multiqc.run("qc_results", "reports")

    # Execute in order
    qc1()
    qc2()
    report()

    # Verify execution order
    assert calls == ["fastqc", "fastqc", "multiqc"]


def _write_fastq(file_path: Path, records: list[tuple[str, str, str]]) -> None:
    """Helper to create a simple FASTQ file with given ID, seq and qual."""
    with open(file_path, "w") as fh:
        for read_id, seq, qual in records:
            fh.write(f"@{read_id}\n")
            fh.write(f"{seq}\n")
            fh.write("+\n")
            fh.write(f"{qual}\n")


def test_build_input_check_single(tmp_path):
    """Input check should process single-end FASTQ correctly."""
    from mmalignments.models.data import Sample
    from mmalignments.models.qc.qc import build_input_check

    r1 = tmp_path / "r1.fastq"
    _write_fastq(r1, [("read1", "ACTG", "!!!!"), ("read2", "GGTA", "####")])

    sample = Sample(name="s1", pairing="single", fastq_r1_path=r1)
    elem = build_input_check(sample, out_dir=tmp_path / "qc")
    elem.run()

    json_file = tmp_path / "qc" / "s1_input_check.json"
    txt_file = tmp_path / "qc" / "s1_input_check.txt"
    assert json_file.exists()
    assert txt_file.exists()

    data = json.loads(json_file.read_text())
    assert data["status"] == "PASS"
    assert data["pairing"] == "single"
    assert data["checks"]["r1_records_checked"] == 2
    assert data["checks"]["r1_format_ok"]


def test_build_input_check_paired_matching(tmp_path):
    """Paired-end check should validate matching ID pairs."""
    from mmalignments.models.data import Sample
    from mmalignments.models.qc.qc import build_input_check

    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    records = [("A/1", "AAA", "!!!"), ("B/1", "CCC", "@@@")]
    _write_fastq(r1, records)
    _write_fastq(r2, records)

    sample = Sample(name="s2", pairing="paired", fastq_r1_path=r1, fastq_r2_path=r2)
    elem = build_input_check(sample, out_dir=tmp_path / "qc")
    elem.run()

    data = json.loads((tmp_path / "qc" / "s2_input_check.json").read_text())
    assert data["status"] == "PASS"
    assert data["checks"]["ids_match"]
    assert data["checks"]["id_mismatches"] == 0
    assert not data["warnings"]


def test_build_input_check_paired_mismatch(tmp_path):
    """Check should report warnings when read IDs differ between R1 and R2."""
    from mmalignments.models.data import Sample
    from mmalignments.models.qc.qc import build_input_check

    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, [("X", "AAA", "!!!")])
    _write_fastq(r2, [("Y", "AAA", "!!!")])

    sample = Sample(name="s3", pairing="paired", fastq_r1_path=r1, fastq_r2_path=r2)
    elem = build_input_check(sample, out_dir=tmp_path / "qc")
    elem.run()

    data = json.loads((tmp_path / "qc" / "s3_input_check.json").read_text())
    assert data["status"] == "PASS"  # still succeeds, but warns
    assert not data["checks"]["ids_match"]
    assert data["checks"]["id_mismatches"] == 1
    assert any("ID mismatch" in w for w in data["warnings"])
    assert any("ID mismatch" in w for w in data["warnings"])
