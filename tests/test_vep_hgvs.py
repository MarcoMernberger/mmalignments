"""Offline tests for the VEP annotation element."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tarfile
from io import BytesIO
from pathlib import Path
from typing import Callable

import pandas as pd
import pytest

from mmalignments.models.externals import External
from mmalignments.models.genome import EnsemblGenome
from mmalignments.models.mutations import vep_hgvs
from mmalignments.models.mutations.vep_hgvs import VEPAnnotationError, VEPAnnotator
from mmalignments.models.overlay import Par, PartialOutSpec


@pytest.fixture
def annotate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Callable[..., tuple[pd.DataFrame, list[list[str]], dict]]:
    """Run the annotation DAG with mocked VEP output."""
    genome = EnsemblGenome("Homo_sapiens", 116, base_dir=tmp_path)
    fasta = genome.fasta
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text(">1\nCAAAA\nAGTCA\nTATAT\nGCC\n")
    genome.fai.write_text("1\t18\t3\t5\t6\n")

    def run(
        frame: pd.DataFrame,
        output_rows: str = "",
        policy: str = "all",
        header: str = "Feature\tFeature_type\tConsequence\tExtra",
        unsupported_spdi: bool = False,
        skipped_message: str | None = None,
        download_cache: bool = False,
        refseq: bool = False,
        one_letter: bool = True,
        extension: str | None = None,
    ) -> tuple[pd.DataFrame, list[list[str]], dict]:
        vcf_records: list[list[str]] = []
        species = "homo_sapiens_refseq" if refseq else "homo_sapiens"
        cache = tmp_path / species / "116_GRCh38"
        cache.mkdir(parents=True, exist_ok=True)
        (cache / "info.txt").write_text("assembly\tGRCh38\n")

        def probe(
            command: list[str], **kwargs: object
        ) -> subprocess.CompletedProcess[str]:
            assert "--help" in command
            if unsupported_spdi and "--spdi" in command:
                return subprocess.CompletedProcess(
                    command, 1, "", "Unknown option: spdi"
                )
            return subprocess.CompletedProcess(command, 0, "ensembl-vep : 117.0\n", "")

        def runnable(
            self: External, *, arguments: list[str], params: Par, **kwargs: object
        ) -> Callable[[], None]:
            def annotate_vcf() -> None:
                assert {
                    "--offline",
                    "--cache",
                    "--check_ref",
                    "--transcript_version",
                    "--hgvsg_use_accession",
                    "--mane",
                } <= set(arguments)
                assert ("--refseq" in self.build_cmd(arguments, params)) == refseq
                assert ("--spdi" in arguments) != unsupported_spdi
                assert arguments[arguments.index("--shift_hgvs") + 1] == "1"
                vcf = Path(arguments[arguments.index("--input_file") + 1])
                lines = vcf.read_text().splitlines()
                assert lines[:3] == [
                    "##fileformat=VCFv4.2",
                    "##contig=<ID=1,length=18>",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                ]
                vcf_records.extend(line.split("\t") for line in lines[3:])
                output = Path(arguments[arguments.index("--output_file") + 1])
                output.write_text(
                    f"## VEP metadata\n#Uploaded_variation\t{header}\n{output_rows}"
                )
                if skipped_message:
                    skipped = Path(
                        arguments[arguments.index("--skipped_variants_file") + 1]
                    )
                    skipped.write_text(skipped_message)

            return annotate_vcf

        monkeypatch.setattr(subprocess, "run", probe)
        monkeypatch.setattr(External, "runnable", runnable)
        source = tmp_path / "input.parquet"
        frame.to_parquet(source)
        tool = VEPAnnotator(
            genome=genome,
            cache_dir=tmp_path,
            executable=sys.executable,
        )
        assert isinstance(tool, External)
        element, process, prepare = tool.annotate(
            source,
            transcript_policy=policy,
            out=PartialOutSpec(folder=tmp_path / "output", ext=extension),
        )
        if not one_letter:
            element = tool.report(
                prepare.pres[0],
                prepare,
                process,
                element.file,
                transcript_policy=policy,
                one_letter=False,
            )
        assert prepare.pres[0].primary.path == source
        prepared_step, cache_step = process.pres
        assert prepared_step is prepare
        info = cache / "info.txt"
        assert cache_step.file == info
        assert tool.is_cache_downloaded(refseq=refseq)
        assert not tool.is_cache_downloaded(refseq=not refseq)
        if download_cache:
            payload = info.read_bytes()
            info.unlink()
            cache.rmdir()
            compressed = BytesIO()
            with tarfile.open(fileobj=compressed, mode="w:gz") as archive:
                member = tarfile.TarInfo(f"{species}/116_GRCh38/info.txt")
                member.size = len(payload)
                archive.addfile(member, BytesIO(payload))

            def fetch(url: str, *, timeout: int) -> BytesIO:
                assert url.endswith(
                    f"/release-116/variation/indexed_vep_cache/{species}_vep_116_GRCh38.tar.gz"
                )
                return BytesIO(compressed.getvalue())

            monkeypatch.setattr(vep_hgvs, "urlopen", fetch)
            assert not tool.is_cache_downloaded(refseq=refseq)
        cache_step()
        assert tool.is_cache_downloaded(refseq=refseq)
        assert cache_step.outputs_ok()[0]
        assert element.pres == (prepare.pres[0], prepare, process)
        for step in (prepare, process, element):
            assert not step.file.exists()
            step()
            assert step.outputs_ok()[0]
        metadata = json.loads(element.artifacts["metadata"].read_text())
        assert (
            metadata["reference_sha256"]
            == hashlib.sha256(fasta.read_bytes()).hexdigest()
        )
        if extension is None:
            result = pd.read_parquet(element.file)
        else:
            suffix = "." + extension.removeprefix(".")
            assert element.file.suffix == suffix
            assert not element.file.name.endswith("." + suffix)
            result = pd.read_csv(
                element.file, sep="\t" if suffix == ".tsv" else ",", index_col=0
            )
        return result, vcf_records, metadata

    return run


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"pos": 0, "ref": "C", "alt": "A"}, "Invalid chromosome"),
        ({"pos": 1, "ref": "G", "alt": "A"}, "Reference mismatch"),
        ({"pos": 1, "ref": "C", "alt": "A,<DEL>"}, "Invalid chromosome"),
        ({"pos": 1, "ref": "C", "alt": "C"}, "Invalid chromosome"),
        (
            {"pos": 1, "ref": "C", "alt": "T", "transcript_id": "NP_1"},
            "Invalid transcript",
        ),
        ({"pos": 19, "ref": "C", "alt": "T"}, "outside the FASTA"),
        ({"pos": 1, "ref": "C"}, "Input requires"),
    ],
)
def test_invalid_input(annotate: Callable, values: dict, message: str) -> None:
    """Reject malformed variants and reference mismatches."""
    with pytest.raises(VEPAnnotationError, match=message):
        annotate(pd.DataFrame([{"chrom": "1", **values}]))


@pytest.mark.parametrize("one_letter", [True, False])
@pytest.mark.parametrize("extension", [None, ".tsv", "csv"])
def test_biallelic_ids_join_duplicates_and_missing(
    annotate: Callable, one_letter: bool, extension: str | None
) -> None:
    """Join each consequence to the original row and ALT."""
    frame = pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "position": [1, 1],
            "ref": ["C", "C"],
            "alt": ["T,G", "T,G"],
            "label": ["same", "same"],
        },
        index=pd.Index(["duplicate", "duplicate"], name="input_id"),
    )
    tab = (
        "row1_alt0\tENST2.5\tTranscript\tmissense_variant\t"
        "HGVSc=ENST2.5:c.1C>T;HGVSp=ENSP1.1:p.(Arg1Cys);HGVSg=NC_000001.11:g.1C>T;SPDI=NC_000001.11:0:C:T\n"
        "row0_alt1\tENST1.3\tTranscript\tmissense_variant\tHGVSc=ENST1.3:c.1C>G;HGVSp=-\n"
        "row0_alt1\tENST2.5\tTranscript\tsynonymous_variant\tHGVSc=ENST2.5:c.1C>G\n"
    )
    result, vcf, metadata = annotate(
        frame, tab, one_letter=one_letter, extension=extension
    )
    assert [line[:5] for line in vcf] == [
        ["1", "1", "row0_alt0", "C", "T"],
        ["1", "1", "row0_alt1", "C", "G"],
        ["1", "1", "row1_alt0", "C", "T"],
        ["1", "1", "row1_alt1", "C", "G"],
    ]
    assert result["vep_variant_id"].tolist() == [
        "row0_alt0",
        "row0_alt1",
        "row0_alt1",
        "row1_alt0",
        "row1_alt1",
    ]
    pd.testing.assert_frame_equal(result[frame.columns], frame.iloc[[0, 0, 0, 1, 1]])
    assert result["vep_row"].tolist() == [0, 0, 0, 1, 1]
    values = result.astype(object).where(result.notna(), None)
    assert values["vep_transcript_id"].tolist() == [
        None,
        "ENST1.3",
        "ENST2.5",
        "ENST2.5",
        None,
    ]
    assert values["hgvs_c"].tolist() == [
        None,
        "ENST1.3:c.1C>G",
        "ENST2.5:c.1C>G",
        "ENST2.5:c.1C>T",
        None,
    ]
    protein = "ENSP1.1:p.(R1C)" if one_letter else "ENSP1.1:p.(Arg1Cys)"
    assert values["hgvs_p"].tolist() == [None, None, None, protein, None]
    assert result["hgvs_g"].iloc[3] == "NC_000001.11:g.1C>T"
    assert result["hgvs_g"].iloc[0] == result["hgvs_g"].iloc[3]
    assert result["spdi"].iloc[3] == "NC_000001.11:0:C:T"
    assert metadata["ensembl_cache_release"] == 116
    assert "independently unverified" in metadata["normalization"]["hgvs_c"]


@pytest.mark.parametrize(
    ("policy", "requested", "expected"),
    [
        ("all", None, ["ENST1.3", "ENST2.5", "ENST3.2"]),
        ("canonical", None, ["ENST1.3", "ENST3.2"]),
        ("mane", None, ["ENST2.5"]),
        ("ENST2", None, ["ENST2.5"]),
        ("ENST2.5", None, ["ENST2.5"]),
        ("ENST2.4", None, []),
        ("ENST2.5", "ENST1.3", ["ENST2.5"]),
        ("mane", "ENST1.3", ["ENST1.3"]),
        ("canonical", "ENST4", [None]),
    ],
)
def test_transcript_selection(
    annotate: Callable, policy: str, requested: str | None, expected: list[str | None]
) -> None:
    """Apply the requested transcript policy and version matching."""
    frame = pd.DataFrame(
        {
            "chrom": ["1"],
            "pos": [1],
            "ref": ["C"],
            "alt": ["T"],
            "transcript_id": [requested],
        }
    )
    header = "Feature\tFeature_type\tConsequence\tCANONICAL\tMANE_SELECT"
    tab = (
        "row0_alt0\tENST1.3\tTranscript\tmissense_variant\tYES\t-\n"
        "row0_alt0\tENST2.5\tTranscript\tmissense_variant\t-\tNM_001.2\n"
        "row0_alt0\tENST3.2\tTranscript\tintron_variant\tYES\t-\n"
    )
    result, _, _ = annotate(frame, tab, policy, header)
    assert result["vep_transcript_id"].tolist() == expected


def test_genomic_three_prime_placement_and_vcf_left_alignment(
    annotate: Callable,
) -> None:
    """Check alignment in both directions for equivalent repeat alleles."""
    frame = pd.DataFrame(
        [
            ("1", 1, "CA", "C"),
            ("1", 5, "AA", "A"),  # delete an A
            ("1", 1, "C", "CA"),
            ("1", 5, "A", "AA"),  # duplicate an A
            ("1", 9, "CAT", "C"),
            ("1", 12, "ATA", "A"),  # delete AT/TA in a repeat
            ("1", 9, "C", "CAT"),
            ("1", 10, "A", "ATA"),  # insert rotated AT/TA
            ("1", 7, "GT", "TC"),  # preserve the MNV
        ],
        columns=["chrom", "pos", "ref", "alt"],
    )
    result, vcf, metadata = annotate(frame, unsupported_spdi=True)
    by_id = {line[2]: (line[1], *line[3:5]) for line in vcf}
    assert by_id["row0_alt0"] == by_id["row1_alt0"] == ("1", "CA", "C")
    assert by_id["row2_alt0"] == by_id["row3_alt0"] == ("1", "C", "CA")
    assert by_id["row4_alt0"] == by_id["row5_alt0"] == ("9", "CAT", "C")
    assert by_id["row6_alt0"] == by_id["row7_alt0"] == ("9", "C", "CAT")
    assert by_id["row8_alt0"] == ("7", "GT", "TC")
    assert result["hgvs_g"].tolist() == [
        "1:g.6del",
        "1:g.6del",
        "1:g.6dup",
        "1:g.6dup",
        "1:g.14_15del",
        "1:g.14_15del",
        "1:g.14_15dup",
        "1:g.14_15dup",
        "1:g.7_8delinsTC",
    ]
    assert result["spdi"].isna().all()
    assert "--spdi" not in metadata["vep_flags"]


def test_malformed_vep_id(annotate: Callable) -> None:
    """Reject a VEP record with an unknown input ID."""
    frame = pd.DataFrame({"chrom": ["1"], "pos": [1], "ref": ["C"], "alt": ["T"]})
    with pytest.raises(VEPAnnotationError, match="unknown variant ID"):
        annotate(frame, "other\tENST1.3\tTranscript\tmissense_variant\t-\n")


def test_vep_reference_mismatch(annotate: Callable) -> None:
    """Surface reference mismatches reported in the VEP skip log."""
    frame = pd.DataFrame({"chrom": ["1"], "pos": [1], "ref": ["C"], "alt": ["T"]})
    with pytest.raises(VEPAnnotationError, match="VEP reference mismatch"):
        annotate(
            frame,
            skipped_message="Specified reference allele C does not match Ensembl reference allele",
        )


def test_mane_falls_back_to_available_transcripts(annotate: Callable) -> None:
    """Keep available transcripts when MANE Select is absent."""
    frame = pd.DataFrame({"chrom": ["1"], "pos": [1], "ref": ["C"], "alt": ["T"]})
    tab = (
        "row0_alt0\tENST1.3\tTranscript\tmissense_variant\tMANE_SELECT=-\n"
        "row0_alt0\tENST2.5\tTranscript\tmissense_variant\tMANE_PLUS_CLINICAL=NM_001.2\n"
    )
    result, _, _ = annotate(frame, tab, "mane")
    assert result["vep_transcript_id"].tolist() == ["ENST1.3", "ENST2.5"]


def test_empty_table_completes_all_steps(annotate: Callable) -> None:
    """Publish a valid empty report without invoking VEP."""
    frame = pd.DataFrame(columns=["chrom", "pos", "ref", "alt"])
    result, vcf, metadata = annotate(frame)
    pd.testing.assert_frame_equal(result[frame.columns], frame)
    assert result.empty
    assert {"hgvs_g", "hgvs_c", "hgvs_p", "spdi"} <= set(result.columns)
    assert not vcf
    assert metadata["assembly"] == "GRCh38"


@pytest.mark.parametrize("refseq", [False, True])
def test_download_cache_for_genome_revision(annotate: Callable, refseq: bool) -> None:
    """Install the reference release and use it for offline annotation."""
    frame = pd.DataFrame({"chrom": ["1"], "pos": [1], "ref": ["C"], "alt": ["A"]})
    policy = "NM_000546.6" if refseq else "ENST00000269305.9"
    rows = (
        f"row0_alt0\t{policy}\tTranscript\tmissense_variant\t-\n"
        f"row0_alt0\t{policy.partition('.')[0]}.1\tTranscript\tmissense_variant\t-\n"
    )
    result, _, metadata = annotate(
        frame, rows, policy=policy, download_cache=True, refseq=refseq
    )
    assert result["vep_transcript_id"].tolist() == [policy]
    assert ("--refseq" in metadata["vep_flags"]) == refseq
    assert metadata["ensembl_cache_release"] == 116
