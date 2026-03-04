import sys
import pathlib
import pytest
import io

# from pypipegraph.testing.fixtures import (  # noqa:F401
#    new_pipegraph,
#    pytest_runtest_makereport,
# )
# from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401
# from pypipegraph.testing import force_load
import pandas as pd
from pathlib import Path

from sklearn.metrics import r2_score
from mbf.align import Sample

# from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401
from pypipegraph.testing.fixtures import new_pipegraph  # noqa:F401

root = pathlib.Path(__file__).parent.parent
sys.path.append(str(root / "src"))
from mmdemultiplex.util import Read, Fragment
from mmdemultiplex.samples import DemultiplexInputSample

R1 = """@A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT
NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFF:FFFFFF
@A01284:56:HNNKWDRXY:1:2101:2248:1000 1:N:0:TAGCTT
NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFF,FFFFFFFFFF:FFFF
"""

R2 = """@A01284:56:HNNKWDRXY:1:2101:1524:1000 2:N:0:TAGCTT
NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC
+
#FF,FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A01284:56:HNNKWDRXY:1:2101:2248:1000 2:N:0:TAGCTT
NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC
+
#FF:FFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF
"""

"""
r2 GCCACC
r1 TCGACC > GGTGGC
r2 AAGTGC
r1 CACAGT > GCACTT
Exaamples below:
#TCA_CTGGCA_TGCCCAGGGTCCGGAGGC_TTTCCC_ATCGATCG_GGGCCC_GGGTGGTTGTCAGTGGCCCTCC_CTTTTG_CAGCTA  # CTGGCA:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer R1
#AGT_GGTCGA_TGCCCAGGGTCCGGAGGC_TTTCCC_ATCGATCG_TTTGG_GGGTGGTTGTCAGTGGCCCTCC_GGTGGC_ATGTAC  # fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer R1
#CTAT_ACTGTG_TGCCCAGGGTCCGGAGGC_TTTCCC_ATCGATCG_CCCT_GGGTGGTTGTCAGTGGCCCTCC_GCACTT_GCTCT  # fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer R1
#TAGCTG_CAAAAG_GGAGGGCCACTGACAACCACCC_GGGTTT_CGATCGAT_GGGCCC_GCCTCCGGACCCTGGGCA_TGCCAG_TGA  # rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_rev-comp-primer R2
#GTACAT_GCCACC_GGAGGGCCACTGACAACCACCC_CCAAA_CGATCGAT_GGGAAA_GCCTCCGGACCCTGGGCA_TCGACC_ACT  # rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_rev-comp-primer R2
#AGAGC_AAGTGC_GGAGGGCCACTGACAACCACCC_AGGG_CGATCGAT_GGGAAA_GCCTCCGGACCCTGGGCA_CACAGT_ATAG  # rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_rev-comp-primer R2
"""
test_reads_pe = {
    "Sample1": (
        (
            "@Sample1:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:001",
            "TCACTGGCATGCCCAGGGTCCGGAGGCTTTCCCATCGATCGGGGCCCGGGTGGTTGTCAGTGGCCCTCCCTTTTGCAGCTA",
        ),
        (
            "@Sample1:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:001",
            "TAGCTGCAAAAGGGAGGGCCACTGACAACCACCCGGGCCCCGATCGATGGGAAAGCCTCCGGACCCTGGGCATGCCAGTGA",
        ),
    ),
    "Sample2": (
        (
            "@Sample2:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:002",
            "AGTGGTCGATGCCCAGGGTCCGGAGGCTTTCCCATCGATCGTTTGGGGGTGGTTGTCAGTGGCCCTCCGGTGGCATGTAC",
        ),
        (
            "@Sample2:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:002",
            "GTACATGCCACCGGAGGGCCACTGACAACCACCCCCAAACGATCGATGGGAAAGCCTCCGGACCCTGGGCATCGACCACT",
        ),
    ),
    "Sample3": (
        (
            "@Sample3:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:003",
            "CTATACTGTGTGCCCAGGGTCCGGAGGCTTTCCCATCGATCGCCCTGGGTGGTTGTCAGTGGCCCTCCGCACTTGCTCT",
        ),
        (
            "@Sample3:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:003",
            "AGAGCAAGTGCGGAGGGCCACTGACAACCACCCAGGGCGATCGATGGGAAAGCCTCCGGACCCTGGGCACACAGTATAG",
        ),
    ),
    "Sample-discard": (
        (
            "@Sample-discard:r1:004",
            "CTATACTATCTAGGATATTCTTCTTATATTAAACTCTTCTTATATATTCGGGCGATCGATCGATGCTAGCTAGCTACGA",
        ),
        (
            "@Sample-discard:r2:004",
            "CTATACTATCTAGGATATTCTTCTTATATTAAACTCTTCTTATATATTCGGGCGATCGATCGATGCTAGCTAGCTACGA",
        ),
    ),
    "Sample2-empty": (
        (
            "@Sample2-empty:fw-overhang_fw-barcode_fw-constant:r1:005",
            "AGTGGTCGATGCCCAGGGTCCGGAGGC",
        ),
        (
            "@Sample2_empty:rev-overhang_rev-barcode_rev-constant:r2:005",
            "GTACATGCCACCGGAGGGCCACTGACAACCACCC",
        ),
    ),
    "Sample3-empty": (
        (
            "@Sample3-empty:fw-overhang_fw-barcode_fw-constant:r1:006",
            "CTATACTGTGTGCCCAGGGTCCGGAGGC",
        ),
        (
            "@Sample3-empty:rev-overhang_rev-barcode_rev-constant_:r2:006",
            "AGAGCAAGTGCGGAGGGCCACTGACAACCACCC",
        ),
    ),
    "Sample2-1mismatch": (
        (
            "@Sample2-1mismatch:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:007",
            "AGTGGTCCATGCCCAGGGTCCGGAGGCTTTCCCATCGATCGTTTGGGGGTGGTTGTCAGTGGCCCTCCCGTGGCATGTAC",
        ),
        (
            "@Sample2-1mismatch:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:007",
            "GTACATGCCTCCGGAGGGCCACTGACAACCACCCCCAAACGATCGATGGGAAAGCCTCCGGACCCTGGGCATCGACCACT",
        ),
    ),
    "Sample2-empty-truncated": (
        (
            "@Sample2-empty-truncated:fw-overhang_fw-barcode_fw-constant:r1:008",
            "AGTGGTCGATGCCCAGGGTCC",
        ),
        (
            "@Sample2-empty-truncated:rev-overhang_rev-barcode_rev-constant_:r2:008",
            "GTACATGCCACCGGAGGGCCACTGACAA",
        ),
    ),
    "Sample2-reversed": (
        (
            "@Sample2-reversed:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:009",
            "GTACATGCCACCGGAGGGCCACTGACAACCACCCCCAAACGATCGATGGGAAAGCCTCCGGACCCTGGGCATCGACCACT",
        ),
        (
            "@Sample2-reversed:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:009",
            "AGTGGTCGATGCCCAGGGTCCGGAGGCTTTCCCATCGATCGTTTGGGGGTGGTTGTCAGTGGCCCTCCGGTGGCATGTAC",
        ),
    ),
    "Sample2-no-adapter": (
        (
            "@Sample2-no-adapter:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:010",
            "AGTCCCCGATGCCCAGGGTCCGGAGGCTTTCCCATCGATCGTTTGGGGGTGGTTGTCAGTGGCCCTCCGGTGGCATGTAC",
        ),
        (
            "@Sample2-no-adapter:rev-overhang_rev-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:010",
            "GTACATCCCCGAGGAGGGCCACTGACAACCACCCCCAAACCTTCGATGGGAAAGCCTCCGGACCCTGGGCATCGACCACT",
        ),
    ),
    "Sample4-fw-rev-same": (
        (
            "@Sample4:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_fw-barcode:r1:011",
            "CTAT_ACTGTG_TGCCCAGGGTCCGGAGGC_TTTCCC_ATCGATCG_CCCT_GGGTGGTTGTCAGTGGCCCTCC_CACAGT_GCTCT",
        ),
        (
            "@Sample4:rev-overhang_fw-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:011",
            "AGAGC_ACTGTG_GGAGGGCCACTGACAACCACCC_AGGG_CGATCGAT_GGGAAA_GCCTCCGGACCCTGGGCA_CACAGT_ATAG",
        ),
    ),
    "Sample4-fw-rev-mismatch": (
        (
            "@Sample4:fw-overhang_fw-barcode_fw-constant_fw-adapter_amplicon_space_fw-barcode:r1:012",
            "CTAT_ACTTTG_TGCCCAGGGTCCGGAGGC_TTTCCC_ATCGATCG_CCCT_GGGTGGTTGTCAGTGGCCCTCC_CAGAGT_GCTCT",
        ),
        (
            "@Sample4:rev-overhang_fw-barcode_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-primer:r2:012",
            "AGAGC_ACTGCG_GGAGGGCCACTGACAACCACCC_AGGG_CGATCGAT_GGGAAA_GCCTCCGGACCCTGGGCA_CACAGT_ATAG",
        ),
    ),
    "Sample3-partial-and-error": (
        (
            "@Sample3:fw-barcode_trunc_fw-constant_fw-adapter_amplicon_space_rev-comp-primer:r1:013",
            "CTGTGTGCCCAGGGTCCGGAGGCTTTCCCATCGATCGCCCTGGGTGGTTGTCAGTGGCCCTCCGCACTTGCTCT",
        ),
        (
            "@Sample3:rev-barcode-trunc_rev-constant_rev-adapter_reverse-comp-amplicon_space_fw-comp-barcode-trunc-error:r2:013",
            "AGTGCGGAGGGCCACTGACAACCACCCAGGGCGATCGATGGGAAAGCCTCCGGACCCTGGGCACACAC",
        ),
    ),
}


class MockBlockedFileAdapter:
    def __init__(self, filename):
        self.filename = filename
        if "_R2_" in self.filename:
            self.file_content = io.BytesIO(R2.encode())
        else:
            self.file_content = io.BytesIO(R1.encode())
        self.file_iterator = None

    def iterator(self):
        for line in self.file_content.readlines():
            if line == "+\n":
                continue
            yield line

    def readline(self):
        if self.file_iterator is None:
            self.file_iterator = self.iterator()
        return next(self.file_iterator)


class DummySample:
    def __init__(self, name):
        self.name = name

    def get_aligner_input_filenames(self):
        return self.filenames


class DummySamplePE(DummySample):
    def __init__(self, name):
        super().__init__(name)
        self.filenames = (f"{self.name}_R1_.fastq", f"{self.name}_R2_.fastq")
        self.is_paired = True


class DummySampleSE(DummySample):
    def __init__(self, name):
        super().__init__(name)
        self.filenames = (f"{self.name}_R1_.fastq",)
        self.is_paired = False


class DummyDemultiplexInputSample:
    def __init__(self, name):
        self.name = name
        self.filenames = [(f"{self.name}_R1_.fastq", f"{self.name}_R2_.fastq")]
        self.is_paired = True

    def get_aligner_input_filenames(self):
        return self.filenames


@pytest.fixture
def se_sample():
    return DummySampleSE("SampleSE")


@pytest.fixture
def pe_sample():
    return DummySamplePE("SamplePE")


@pytest.fixture
def pe_sample_demultiplex():
    return DummyDemultiplexInputSample("SamplePE_DemultiplexInputSample")


@pytest.fixture
def paired_fragments():
    fragments = {}
    for name in test_reads_pe:
        r1_tup, r2_tup = test_reads_pe[name]
        r1 = Read(r1_tup[0], r1_tup[1], "F" * len(r1_tup[1]))
        r2 = Read(r2_tup[0], r2_tup[1], "F" * len(r2_tup[1]))
        fragment = Fragment(r1, r2)
        fragments[name] = fragment
    return fragments


def barcode_df_full_callback():
    df = pd.DataFrame(
        {
            "key": ["Sample1", "Sample2", "Sample3", "Sample4"],
            "start_barcode": ["CTGGCA", "GGTCGA", "ACTGTG", "ACTGTG"],
            "end_barcode": ["CAAAAG", "GCCACC", "AAGTGC", "ACTGTG"],
            "trim_after_start": [24, 24, 24, 27],
            "trim_before_end": [28, 27, 26, 29],
            "maximal_errors_start": [0, 1, 1, 0],
            "maximal_errors_end": [0, 1, 0, 0],
            "minimal_overlap_start": [6, 6, 5, 6],
            "minimal_overlap_end": [6, 6, 5, 6],
        }
    )
    df = df.set_index("key")
    return df


def barcode_df_callback():
    df = pd.DataFrame(
        {
            "key": ["Sample1", "Sample2", "Sample3"],
            "start_barcode": ["CTGGCA", "GGTCGA", "ACTGTG"],
            "end_barcode": ["CAAAAG", "GCCACC", "AAGTGC"],
        }
    )
    df = df.set_index("key")
    return df
