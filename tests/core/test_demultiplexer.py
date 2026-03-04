# 7  # -*- coding: utf-8 -*-

# import pytest
# import pandas as pd
# import pypipegraph as ppg
# import mbf.align
# from mmdemultiplex import (
#     Demultiplexer,
#     DemultiplexStrategy,
#     PE_Decide_On_Start_Trim_Start_End,
#     Fragment,
# )
# from mmdemultiplex.util import get_fastq_iterator
# from pathlib import Path
# from unittest.mock import patch
# from pypipegraph import Job
# from conftest import (
#     MockBlockedFileAdapter,
#     DummyDemultiplexInputSample,
# )

# __author__ = "MarcoMernberger"
# __copyright__ = "MarcoMernberger"
# __license__ = "mit"


# def barcode_df_callback():
#     df = pd.DataFrame(
#         {"key": ["Sample1"], "start_barcode": ["ATCG"], "end_barcode": ["TTCG"]}
#     )
#     df = df.set_index("key")
#     return df


# class MockDecisionCallback:
#     def __init__(self):
#         self.name = "MockStrategy"

#     def match_and_trim(self, fragment):
#         if fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT":
#             return fragment
#         else:
#             return False


# def test_init(tmp_path, pe_sample, se_sample):
#     barcode_df = barcode_df_callback()
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path, prefix="DM_"
#     )
#     assert demultiplexer.name == f"DM_{pe_sample.name}"
#     assert hasattr(demultiplexer, "barcode_df")
#     assert isinstance(demultiplexer.barcode_df, pd.DataFrame)
#     assert demultiplexer.library_name == pe_sample.name
#     assert demultiplexer.output_folder == tmp_path / demultiplexer.name
#     assert demultiplexer.output_folder.exists()
#     assert demultiplexer.input_sample == pe_sample
#     assert isinstance(demultiplexer.input_files, list)
#     assert isinstance(demultiplexer.input_files[0], tuple)
#     assert demultiplexer.input_files == [pe_sample.get_aligner_input_filenames()]
#     assert issubclass(demultiplexer.strategy, DemultiplexStrategy)
#     assert demultiplexer.is_paired
#     assert hasattr(demultiplexer, "decision_callbacks")
#     assert len(demultiplexer.decision_callbacks) == barcode_df.shape[0]
#     assert isinstance(
#         demultiplexer.decision_callbacks["Sample1"], PE_Decide_On_Start_Trim_Start_End
#     )
#     # with SE
#     demultiplexer = Demultiplexer(
#         se_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     assert isinstance(demultiplexer.input_files, list)
#     assert isinstance(demultiplexer.input_files[0], tuple)
#     assert demultiplexer.name == f"{se_sample.name}"
#     assert not demultiplexer.is_paired
#     # with default output folder
#     with patch("pathlib.Path.mkdir", return_value=None):
#         demultiplexer = Demultiplexer(pe_sample, barcode_df_callback)
#         assert demultiplexer.output_folder == (Path("cache") / demultiplexer.name)


# def test_with_DemultiplexInputSample(pe_sample_demultiplex, tmp_path):
#     demultiplexer = Demultiplexer(
#         pe_sample_demultiplex, barcode_df_callback, output_folder=tmp_path
#     )
#     assert demultiplexer.name == f"{pe_sample_demultiplex.name}"
#     assert demultiplexer.is_paired
#     assert isinstance(demultiplexer.input_files, list)
#     assert isinstance(demultiplexer.input_files[0], tuple)
#     sample_with_lists = DummyDemultiplexInputSample("SamplePE_DemultiplexInputSample")
#     sample_with_lists.filenames = ["R1", "R2"]
#     demultiplexer = Demultiplexer(
#         sample_with_lists, barcode_df_callback, output_folder=tmp_path
#     )
#     print(demultiplexer.input_files)
#     assert isinstance(demultiplexer.input_files, list)
#     assert isinstance(demultiplexer.input_files[0], tuple)


# @pytest.mark.usefixtures("new_pipegraph")
# def test_get_dependencies(tmp_path, pe_sample):
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     for job in demultiplexer.get_dependencies():
#         assert isinstance(job, Job) or isinstance(job, str)
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     pe_sample.prepare_input = lambda: "prepare_input"
#     pe_sample.dependencies = ["own_deps"]
#     deps = demultiplexer.get_dependencies()
#     assert "own_deps" in deps
#     assert "prepare_input" in deps


# def test_parameters(tmp_path, pe_sample):
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     parameters = demultiplexer.parameters
#     assert parameters[0] == demultiplexer.name
#     assert parameters[1] == pe_sample.name
#     assert parameters[2] == "PE_Decide_On_Start_Trim_Start_End"
#     for item in parameters:
#         assert hasattr(item, "__hash__")


# def test_parameter_string(tmp_path, pe_sample):
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     assert isinstance(demultiplexer.parameter_string(), str)


# def test_fastq_iterator(tmp_path, pe_sample, se_sample):
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer = Demultiplexer(
#             pe_sample, barcode_df_callback, output_folder=tmp_path
#         )
#         iterator = demultiplexer.get_fastq_iterator()
#         assert callable(iterator)
#         for files_tuple in demultiplexer.input_files:
#             for fragment in iterator(files_tuple):
#                 assert isinstance(fragment, Fragment)
#                 assert hasattr(fragment, "Read1")
#                 assert hasattr(fragment, "Read2")
#                 assert hasattr(fragment, "reads")
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer = Demultiplexer(
#             se_sample, barcode_df_callback, output_folder=tmp_path
#         )
#         iterator = demultiplexer.get_fastq_iterator()
#         assert callable(iterator)
#         for files_tuple in demultiplexer.input_files:
#             for fragment in iterator(files_tuple):
#                 assert isinstance(fragment, Fragment)
#                 assert hasattr(fragment, "Read1")
#                 assert not hasattr(fragment, "Read2")
#                 assert hasattr(fragment, "reads")


# def test_demultiplexer_prefix(tmp_path, pe_sample):
#     prefix = "DM_"
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path, prefix=prefix
#     )
#     assert demultiplexer.name == f"{prefix}{pe_sample.name}"
#     assert demultiplexer.output_folder == tmp_path / f"{prefix}{pe_sample.name}"


# @pytest.mark.usefixtures("new_pipegraph")
# def test_do_demultiplex_pe(tmp_path, pe_sample):
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     first_read_sample_name = f"{pe_sample.name}_first_read"
#     discarded_sample_name = f"{pe_sample.name}_discarded"
#     files_created = {
#         "first_read": (
#             demultiplexer.output_folder
#             / first_read_sample_name
#             / f"{first_read_sample_name}_R1_.fastq",
#             demultiplexer.output_folder
#             / first_read_sample_name
#             / f"{first_read_sample_name}_R2_.fastq",
#         ),
#         "discarded": (
#             demultiplexer.output_folder
#             / discarded_sample_name
#             / f"{discarded_sample_name}_R1_.fastq",
#             demultiplexer.output_folder
#             / discarded_sample_name
#             / f"{discarded_sample_name}_R2_.fastq",
#         ),
#     }
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         job = demultiplexer.do_demultiplex()
#         ppg.run_pipegraph()
#         sentinel = demultiplexer.output_folder / "done.txt"
#         assert sentinel.exists()
#         for filepath in job.filenames:
#             assert Path(filepath).exists()

#     fastq_iterator = get_fastq_iterator(pe_sample.is_paired)
#     for fragment in fastq_iterator(files_created["first_read"]):
#         assert (
#             fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
#         )
#         assert (
#             fragment.Read2.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 2:N:0:TAGCTT"
#         )
#         assert (
#             fragment.Read1.Sequence
#             == "NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
#         )
#         assert (
#             fragment.Read2.Sequence
#             == "NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC"
#         )
#         assert (
#             fragment.Read1.Quality
#             == "#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFF:FFFFFF"
#         )
#         assert (
#             fragment.Read2.Quality
#             == "#FF,FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
#         )
#         break
#     for fragment in fastq_iterator(files_created["discarded"]):
#         assert (
#             fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:2248:1000 1:N:0:TAGCTT"
#         )
#         assert (
#             fragment.Read2.Name == "A01284:56:HNNKWDRXY:1:2101:2248:1000 2:N:0:TAGCTT"
#         )
#         assert (
#             fragment.Read1.Sequence
#             == "NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
#         )
#         assert (
#             fragment.Read2.Sequence
#             == "NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC"
#         )
#         assert (
#             fragment.Read1.Quality
#             == "#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFF,FFFFFFFFFF:FFFF"
#         )
#         assert (
#             fragment.Read2.Quality
#             == "#FF:FFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF"
#         )
#         break


# @pytest.mark.usefixtures("new_pipegraph")
# def test_do_demultiplex_se(tmp_path, se_sample):
#     first_read_sample_name = f"{se_sample.name}_first_read"
#     discarded_sample_name = f"{se_sample.name}_discarded"
#     files_created = {
#         "first_read": (
#             tmp_path
#             / se_sample.name
#             / first_read_sample_name
#             / f"{first_read_sample_name}_R1_.fastq",
#         ),
#         "discarded": (
#             tmp_path
#             / se_sample.name
#             / discarded_sample_name
#             / f"{discarded_sample_name}_R1_.fastq",
#         ),
#     }
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer = Demultiplexer(
#             se_sample, barcode_df_callback, output_folder=tmp_path
#         )
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         job = demultiplexer.do_demultiplex()
#         ppg.run_pipegraph()
#         for filepath in job.filenames:
#             assert Path(filepath).exists()
#         sentinel = demultiplexer.output_folder / "done.txt"
#         assert sentinel.exists()
#     fastq_iterator = get_fastq_iterator(se_sample.is_paired)
#     for fragment in fastq_iterator(files_created["first_read"]):
#         assert (
#             fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
#         )
#         assert not hasattr(fragment, "Read2")
#         break
#     for fragment in fastq_iterator(files_created["discarded"]):
#         assert (
#             fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:2248:1000 1:N:0:TAGCTT"
#         )
#         break


# def test_decide_on_barcode(tmp_path, pe_sample):
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         fragments = []
#         demultiplexer = Demultiplexer(
#             pe_sample, barcode_df_callback, output_folder=tmp_path
#         )
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         iterator = demultiplexer.get_fastq_iterator()
#         for files_tuple in demultiplexer.input_files:
#             for fragment in iterator(files_tuple):
#                 fragments.append(fragment)
#         first_result = demultiplexer._decide_on_barcode(fragments[0])
#         assert first_result[0] == "first_read"
#         assert first_result[1] == fragments[0]
#         second_result = demultiplexer._decide_on_barcode(fragments[1])
#         assert second_result[0] == "discarded"
#         assert second_result[1] == fragments[1]


# @pytest.mark.usefixtures("new_pipegraph")
# def test__make_samples(pe_sample, tmp_path):
#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         job = demultiplexer.do_demultiplex()
#         demultiplexed_samples = demultiplexer._make_samples()
#         ppg.run_pipegraph()
#         assert isinstance(demultiplexed_samples, dict)
#         a_sample = list(demultiplexed_samples.values())[0]
#         assert isinstance(a_sample, mbf.align.raw.Sample)


# @pytest.mark.usefixtures("new_pipegraph")
# def test_get_samples(pe_sample, tmp_path):
#     def mock_make():

#         raise ValueError("this should not be called")

#     demultiplexer = Demultiplexer(
#         pe_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         demultiplexer.do_demultiplex()
#         demultiplexed_samples = demultiplexer.get_samples()
#         assert isinstance(demultiplexed_samples, dict)
#         a_sample = list(demultiplexed_samples.values())[0]
#         assert a_sample.pairing == "paired"
#         assert not a_sample.reverse_reads
#         assert isinstance(a_sample.fastq_processor, mbf.align.fastq2.Straight)
#         assert isinstance(a_sample, mbf.align.raw.Sample)

#         with patch("mmdemultiplex.Demultiplexer._make_samples", mock_make):
#             demultiplexed_samples2 = demultiplexer.get_samples()
#         assert demultiplexed_samples2 == demultiplexed_samples
#         ppg.run_pipegraph()


# @pytest.mark.usefixtures("new_pipegraph")
# def test_get_samples_se(se_sample, tmp_path):
#     demultiplexer = Demultiplexer(
#         se_sample, barcode_df_callback, output_folder=tmp_path
#     )
#     with patch("mbf.align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
#         demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
#         demultiplexer.do_demultiplex()
#         demultiplexed_samples = demultiplexer.get_samples()
#         a_sample = list(demultiplexed_samples.values())[0]
#         assert a_sample.pairing == "single"
#         ppg.run_pipegraph()


# def get_samples(self):
#     if not hasattr(self, "raw_samples"):
#         self.raw_samples = self._make_samples()
#     return self.raw_samples


# @pytest.mark.usefixtures("new_pipegraph")
# def test_divide_reads(tmp_path, pe_sample):
#     demultiplexer = Demultiplexer(pe_sample, barcode_df_callback)
#     demultiplexer.divide_reads(
#         input_files=[
#             (
#                 Path.cwd().parent.parent / "data" / "one_R1_.fastq",
#                 Path.cwd().parent.parent / "data" / "one_R2_.fastq",
#             ),
#         ],
#         new_output_folder=tmp_path / "test",
#     )

#     ppg.run_pipegraph()
#     outfiles = (
#         tmp_path / "test" / "one_R1_.fastq",
#         tmp_path / "test" / "one_R2_.fastq",
#     )
#     for fragment in demultiplexer.get_fastq_iterator()(outfiles):
#         assert len(fragment.Read1.Sequence) >= len(fragment.Read2.Sequence)
