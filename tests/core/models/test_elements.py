"""Test class for elements.py"""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import MappingProxyType

import pytest

from mmalignments.models.elements import (
    Element,
    FileElement,
    FilesElement,
    MappedElement,
    NextGenSampleElement,
    Sample,
    ValidationPolicy,
    VcfElement,
    _as_path,
    _looks_like_filepath,
    element,
    file_sig,
    generate_element_key_name,
    get_candidates,
    sample_fastqs,
    short_hash,
    stable_hash,
)
from mmalignments.models.registry import ElementRegistry, element_build_context
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
)

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def make_tag(root: str = "sample", level: int = 1) -> ElementTag:
    return ElementTag(
        root=root,
        level=level,
        omics=None,
        stage=Stage.ALIGN,
        method=Method.BWAMEM2,
        state=State.RAW,
        ext="bam",
    )


def make_element(
    key: str = "mykey",
    *,
    run=None,
    tag: ElementTag | None = None,
    artifacts: dict | None = None,
    inputs: tuple | None = None,
    determinants: tuple | None = None,
    pres: tuple | None = None,
    validator=None,
    empty_ok: bool = False,
    name: str | None = None,
) -> Element:
    if tag is None:
        tag = make_tag()
    if run is None:
        run = lambda: None  # noqa: E731
    return Element(
        key=key,
        run=run,
        tag=tag,
        artifacts=artifacts,
        inputs=inputs,
        determinants=determinants,
        pres=pres,
        validator=validator,
        empty_ok=empty_ok,
        name=name,
    )


# ---------------------------------------------------------------------------
# file_sig
# ---------------------------------------------------------------------------


def test_file_sig_missing(tmp_path):
    p = tmp_path / "no_such_file.txt"
    result = file_sig(p)
    assert result == {"path": str(p), "missing": True}


def test_file_sig_existing(tmp_path):
    p = tmp_path / "data.txt"
    p.write_bytes(b"hello world")
    result = file_sig(p)
    assert result["path"] == str(p)
    assert result["size"] == len(b"hello world")
    assert "head_sha256" in result
    # verify the hash value
    expected = hashlib.sha256(b"hello world").hexdigest()
    assert result["head_sha256"] == expected


# ---------------------------------------------------------------------------
# stable_hash / short_hash
# ---------------------------------------------------------------------------


def test_stable_hash_deterministic():
    obj = {"a": 1, "b": [2, 3]}
    assert stable_hash(obj) == stable_hash(obj)


def test_stable_hash_different_for_different_objects():
    assert stable_hash({"a": 1}) != stable_hash({"a": 2})


def test_stable_hash_sort_keys():
    assert stable_hash({"a": 1, "b": 2}) == stable_hash({"b": 2, "a": 1})


def test_short_hash():
    h = "abcdef1234567890"
    assert short_hash(h) == "abcdef12"
    assert short_hash(h, n=4) == "abcd"


# ---------------------------------------------------------------------------
# ValidationPolicy
# ---------------------------------------------------------------------------


def test_validation_policy_values():
    assert ValidationPolicy.CHECK.value == "check"
    assert ValidationPolicy.FORCE_RUN.value == "force_run"
    assert ValidationPolicy.FORCE_SKIP.value == "force_skip"


# ---------------------------------------------------------------------------
# Element – basic construction
# ---------------------------------------------------------------------------


def test_element_init_minimal():
    tag = make_tag()
    e = make_element(tag=tag)
    assert e.key == "mykey"
    assert e.tag is tag
    assert e.artifacts == MappingProxyType({})
    assert e.pres == ()
    assert e.inputs == ()
    assert e.determinants == ()
    assert e.validation_policy == ValidationPolicy.CHECK


def test_element_init_with_name():
    e = make_element(name="custom_name")
    assert e.name == "custom_name"


def test_element_default_name_from_tag():
    tag = make_tag(root="sample", level=1)
    e = make_element(tag=tag)
    assert e.name == tag.default_name


def test_element_tag_none_raises():
    with pytest.raises(ValueError, match="tag must be provided"):
        Element(key="k", run=lambda: None, tag=None)


def test_element_init_with_determinants():
    e = make_element(determinants=("v1.0", "param=foo"))
    assert "v1.0" in e.key
    assert "param=foo" in e.key
    assert e.determinants == ("v1.0", "param=foo")


def test_element_determinants_as_str():
    e = make_element(determinants=("a", "b", "c"))
    assert e.determinants_as_str() == "a,b,c"


def test_element_init_with_inputs(tmp_path):
    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    p1.touch()
    p2.touch()
    e = make_element(inputs=(p1, p2))
    assert set(e.inputs) == {p1, p2}


def test_element_inputs_sorted(tmp_path):
    p_z = tmp_path / "z.txt"
    p_a = tmp_path / "a.txt"
    p_z.touch()
    p_a.touch()
    e = make_element(inputs=(p_z, p_a))
    assert e.inputs[0] == p_a
    assert e.inputs[1] == p_z


def test_element_with_path_artifacts(tmp_path):
    bam = tmp_path / "out.bam"
    e = make_element(artifacts={"bam": bam})
    assert e.artifacts["bam"] == bam


def test_element_output_files_from_artifacts(tmp_path):
    bam = tmp_path / "out.bam"
    e = make_element(artifacts={"bam": bam, "label": "hello"})
    assert bam in e.output_files


def test_element_output_files_none_when_no_paths():
    e = make_element(artifacts={"label": "hello"})
    assert e.output_files is None


def test_element_output_files_none_when_no_artifacts():
    e = make_element()
    assert e.output_files is None


# ---------------------------------------------------------------------------
# Element – validate_fields
# ---------------------------------------------------------------------------


def test_validate_fields_invalid_output_path():
    tag = make_tag()
    with pytest.raises(AssertionError):
        Element(
            key="k",
            run=lambda: None,
            tag=tag,
            artifacts={"out": 12345},  # integer, not a Path or str
        )


def test_validate_fields_invalid_input_path():
    tag = make_tag()
    with pytest.raises(AssertionError):
        Element(
            key="k",
            run=lambda: None,
            tag=tag,
            inputs=(42,),  # not a Path
        )


# ---------------------------------------------------------------------------
# Element – validation policy & force
# ---------------------------------------------------------------------------


def test_element_force_sets_policy():
    e = make_element()
    returned = e.force()
    assert e.validation_policy == ValidationPolicy.FORCE_RUN
    assert returned is e


def test_element_validation_policy_setter():
    e = make_element()
    e.validation_policy = ValidationPolicy.FORCE_SKIP
    assert e.validation_policy == ValidationPolicy.FORCE_SKIP


# ---------------------------------------------------------------------------
# Element – signature
# ---------------------------------------------------------------------------


def test_element_signature_is_string():
    e = make_element()
    sig = e.signature
    assert isinstance(sig, str)
    assert len(sig) == 64  # sha256 hex


def test_element_signature_stable():
    e = make_element(key="stable_key")
    assert e.signature == e.signature


def test_element_signature_changes_with_key():
    e1 = make_element(key="key_a")
    e2 = make_element(key="key_b")
    assert e1.signature != e2.signature


# ---------------------------------------------------------------------------
# Element – sig_data / print_sig_data
# ---------------------------------------------------------------------------


def test_element_sig_data_structure():
    e = make_element(key="k", determinants=("d1",))
    data = e.sig_data()
    assert "key" in data
    assert "determinants" in data
    assert "inputs" in data
    assert "artifacts" in data
    assert "pre_sigs" in data


def test_element_print_sig_data(capsys):
    e = make_element()
    e.print_sig_data()
    captured = capsys.readouterr()
    assert "key" in captured.out


# ---------------------------------------------------------------------------
# Element – artifact signature types
# ---------------------------------------------------------------------------


def test_artifact_sig_with_path(tmp_path):
    p = tmp_path / "out.bam"
    e = make_element(artifacts={"bam": p})
    sig = e._artifact_sig()
    assert sig["bam"] == str(p)


def test_artifact_sig_with_primitives():
    e = make_element(
        artifacts={
            "count": 5,
            "flag": True,
            "label": "x",
            "ratio": 1.5,
            "none_val": None,
        }
    )
    sig = e._artifact_sig()
    assert sig["count"] == 5
    assert sig["flag"] is True
    assert sig["label"] == "x"
    assert sig["ratio"] == 1.5
    assert sig["none_val"] is None


def test_artifact_sig_unsupported_type():
    e = make_element(artifacts={"bad": [1, 2, 3]})
    with pytest.raises(TypeError, match="unsupported type for signature"):
        e._artifact_sig()


# ---------------------------------------------------------------------------
# Element – outputs_ok
# ---------------------------------------------------------------------------


def test_outputs_ok_no_output_files():
    e = make_element()
    ok, reason = e.outputs_ok()
    assert ok is True


def test_outputs_ok_file_exists(tmp_path):
    p = tmp_path / "out.bam"
    p.write_bytes(b"data")
    e = make_element(artifacts={"bam": p})
    ok, reason = e.outputs_ok()
    assert ok is True


def test_outputs_ok_file_missing(tmp_path):
    p = tmp_path / "missing.bam"
    e = make_element(artifacts={"bam": p})
    ok, reason = e.outputs_ok()
    assert ok is False
    assert "Missing" in reason


def test_outputs_ok_empty_file(tmp_path):
    p = tmp_path / "empty.bam"
    p.touch()  # 0 bytes
    e = make_element(artifacts={"bam": p})
    ok, reason = e.outputs_ok()
    assert ok is False
    assert "Empty" in reason


def test_outputs_ok_empty_file_allowed(tmp_path):
    p = tmp_path / "empty.bam"
    p.touch()
    e = make_element(artifacts={"bam": p}, empty_ok=True)
    ok, reason = e.outputs_ok()
    assert ok is True


# ---------------------------------------------------------------------------
# Element – skip logic
# ---------------------------------------------------------------------------


def test_skip_force_run():
    e = make_element()
    e.validation_policy = ValidationPolicy.FORCE_RUN
    do_skip, reason = e.skip()
    assert do_skip is False
    assert "forces run" in reason


def test_skip_force_skip():
    e = make_element()
    e.validation_policy = ValidationPolicy.FORCE_SKIP
    do_skip, reason = e.skip()
    assert do_skip is True
    assert "forces skip" in reason


def test_skip_no_cached_signature():
    e = make_element()
    do_skip, reason = e.skip(cached_signature=None)
    assert do_skip is False
    assert "First run" in reason


def test_skip_mismatched_signature():
    e = make_element()
    do_skip, reason = e.skip(cached_signature="wrong_sig")
    assert do_skip is False
    assert "Cached signature" in reason


def test_skip_matching_signature_outputs_ok(tmp_path):
    p = tmp_path / "out.bam"
    p.write_bytes(b"content")
    e = make_element(artifacts={"bam": p})
    sig = e.signature
    do_skip, reason = e.skip(cached_signature=sig)
    assert do_skip is True


def test_skip_matching_signature_outputs_missing(tmp_path):
    p = tmp_path / "missing.bam"
    e = make_element(artifacts={"bam": p})
    sig = e.signature
    do_skip, reason = e.skip(cached_signature=sig)
    assert do_skip is False


def test_skip_with_custom_validator_returning_false(tmp_path):
    p = tmp_path / "out.bam"
    p.write_bytes(b"content")
    e = make_element(
        artifacts={"bam": p}, validator=lambda: (False, "validator says no")
    )
    sig = e.signature
    do_skip, reason = e.skip(cached_signature=sig)
    assert do_skip is False
    assert "validator says no" in reason


def test_skip_with_custom_validator_returning_true(tmp_path):
    p = tmp_path / "out.bam"
    p.write_bytes(b"content")
    e = make_element(
        artifacts={"bam": p}, validator=lambda: (True, "validator says yes")
    )
    sig = e.signature
    do_skip, reason = e.skip(cached_signature=sig)
    assert do_skip is True


# ---------------------------------------------------------------------------
# Element – __call__
# ---------------------------------------------------------------------------


def test_element_call_invokes_run():
    called = []

    def my_run():
        called.append(1)
        return True

    e = make_element(run=my_run)
    result = e()
    assert called == [1]
    assert result is True


# ---------------------------------------------------------------------------
# Element – provenance
# ---------------------------------------------------------------------------


def test_provenance_no_pres():
    e = make_element()
    assert "->" not in e.provenance or e.provenance.endswith(e.name)


def test_provenance_single_pre():
    pre = make_element(key="pre_key")
    e = make_element(key="child_key", pres=(pre,))
    assert pre.name in e.provenance
    assert "->" in e.provenance


def test_provenance_multiple_pres():
    pre1 = make_element(key="pre1_key")
    pre2 = make_element(key="pre2_key")
    e = make_element(key="child_key", pres=(pre1, pre2))
    assert "(" in e.provenance
    assert "," in e.provenance


# ---------------------------------------------------------------------------
# Element – __repr__ and describe
# ---------------------------------------------------------------------------


def test_element_repr():
    e = make_element()
    r = repr(e)
    assert "Element" in r
    assert "key=" in r


def test_element_describe():
    e = make_element()
    d = e.describe()
    assert "key" in d
    assert "signature" in d


# ---------------------------------------------------------------------------
# Element – __getattr__ delegation to artifacts
# ---------------------------------------------------------------------------


def test_getattr_from_artifacts(tmp_path):
    p = tmp_path / "out.bam"
    e = make_element(artifacts={"bam": p})
    assert e.bam == p


def test_getattr_missing_raises():
    e = make_element()
    with pytest.raises(AttributeError):
        _ = e.nonexistent_attribute


# ---------------------------------------------------------------------------
# Element – __hash__ and __eq__
# ---------------------------------------------------------------------------


def test_element_hash():
    e = make_element(key="hashkey")
    assert hash(e) == hash("hashkey")


def test_element_eq_same_key():
    e1 = make_element(key="same")
    e2 = make_element(key="same")
    assert e1 == e2


def test_element_eq_different_key():
    e1 = make_element(key="a")
    e2 = make_element(key="b")
    assert e1 != e2


def test_element_eq_non_element():
    e = make_element(key="k")
    assert e != "not_an_element"


# ---------------------------------------------------------------------------
# Element – root attribute
# ---------------------------------------------------------------------------


def test_element_root():
    tag = make_tag(root="kidney_1")
    e = make_element(tag=tag)
    assert e.root == "kidney_1"


# ---------------------------------------------------------------------------
# MappedElement
# ---------------------------------------------------------------------------


def test_mapped_element_bam(tmp_path):
    bam = tmp_path / "out.bam"
    tag = make_tag()
    e = MappedElement(key="mapped", run=lambda: None, tag=tag, artifacts={"bam": bam})
    assert e.bam == bam


def test_mapped_element_bam_wrong_type():
    tag = make_tag()
    e = MappedElement(
        key="mapped", run=lambda: None, tag=tag, artifacts={"bam": "not_a_path"}
    )
    with pytest.raises(TypeError, match="must be a Path"):
        _ = e.bam


# ---------------------------------------------------------------------------
# VcfElement
# ---------------------------------------------------------------------------


def test_vcf_element_vcf(tmp_path):
    vcf = tmp_path / "out.vcf"
    tag = make_tag()
    e = VcfElement(key="vcf_el", run=lambda: None, tag=tag, artifacts={"vcf": vcf})
    assert e.vcf == vcf


def test_vcf_element_vcf_wrong_type():
    tag = make_tag()
    e = VcfElement(
        key="vcf_el", run=lambda: None, tag=tag, artifacts={"vcf": "not_a_path"}
    )
    with pytest.raises(TypeError, match="must be a Path"):
        _ = e.vcf


# ---------------------------------------------------------------------------
# _as_path
# ---------------------------------------------------------------------------


def test_as_path_from_path(tmp_path):
    p = tmp_path / "f.txt"
    assert _as_path(p) == p


def test_as_path_from_str():
    result = _as_path("/some/path.bam")
    assert isinstance(result, Path)


def test_as_path_from_empty_str():
    assert _as_path("  ") is None


def test_as_path_from_none():
    assert _as_path(None) is None


def test_as_path_from_int():
    assert _as_path(42) is None


# ---------------------------------------------------------------------------
# _looks_like_filepath
# ---------------------------------------------------------------------------


def test_looks_like_filepath_with_slash():
    assert _looks_like_filepath(Path("/some/path")) is True


def test_looks_like_filepath_with_extension():
    assert _looks_like_filepath(Path("file.bam")) is True


def test_looks_like_filepath_bare_name():
    assert _looks_like_filepath(Path("justname")) is False


# ---------------------------------------------------------------------------
# get_candidates
# ---------------------------------------------------------------------------


def test_get_candidates_with_outputs_key(tmp_path):
    p = tmp_path / "out.bam"
    arts = {"bam": p, "label": "x"}
    candidates = get_candidates(arts, outputs="bam")
    assert p in candidates


def test_get_candidates_with_outputs_list(tmp_path):
    p1 = tmp_path / "a.bam"
    p2 = tmp_path / "b.bam"
    arts = {"a": p1, "b": p2}
    candidates = get_candidates(arts, outputs=["a", "b"])
    assert p1 in candidates
    assert p2 in candidates


def test_get_candidates_with_output_files(tmp_path):
    p = tmp_path / "out.bam"
    candidates = get_candidates({}, output_files=[p], auto_outputs=True)
    assert p in candidates


def test_get_candidates_auto_outputs_false(tmp_path):
    p = tmp_path / "out.bam"
    candidates = get_candidates({}, output_files=[p], auto_outputs=False)
    assert candidates == []


def test_get_candidates_no_outputs():
    candidates = get_candidates({})
    assert candidates == []


# ---------------------------------------------------------------------------
# @element decorator
# ---------------------------------------------------------------------------


def test_element_decorator_basic():
    tag = make_tag()

    @element
    def build_element():
        return Element(key="deco_key", run=lambda: None, tag=tag)

    e = build_element()
    assert isinstance(e, Element)
    assert e.key == "deco_key"


def test_element_decorator_with_outputs(tmp_path):
    tag = make_tag()
    bam = tmp_path / "subdir" / "out.bam"

    @element(outputs="bam")
    def build_element():
        return Element(
            key="deco_out", run=lambda: None, tag=tag, artifacts={"bam": bam}
        )

    e = build_element()
    # parent directory should have been created
    assert bam.parent.exists()


def test_element_decorator_with_registry():
    tag = make_tag()
    registry = ElementRegistry()

    @element
    def build_element():
        return Element(key="reg_key", run=lambda: None, tag=tag)

    with element_build_context(registry):
        e = build_element()

    assert registry.get("reg_key") is e


def test_element_decorator_auto_outputs_false(tmp_path):
    tag = make_tag()
    bam = tmp_path / "nodir" / "out.bam"

    @element(auto_outputs=False)
    def build_element():
        return Element(key="no_out", run=lambda: None, tag=tag, artifacts={"bam": bam})

    build_element()
    # directory should NOT be created because auto_outputs=False
    assert not bam.parent.exists()


# ---------------------------------------------------------------------------
# FilesElement
# ---------------------------------------------------------------------------


def test_files_element_from_single_path(tmp_path):
    p = tmp_path / "data.bam"
    p.write_bytes(b"data")
    fe = FilesElement(p)
    assert p.absolute() in fe.files


def test_files_element_from_mapping(tmp_path):
    p1 = tmp_path / "r1.fastq"
    p2 = tmp_path / "r2.fastq"
    p1.write_bytes(b"data")
    p2.write_bytes(b"data")
    fe = FilesElement({"r1": p1, "r2": p2})
    assert len(fe.files) == 2


def test_files_element_output_files_is_none(tmp_path):
    p = tmp_path / "data.bam"
    p.write_bytes(b"data")
    fe = FilesElement(p)
    assert fe.output_files is None


def test_files_element_validate(tmp_path):
    p = tmp_path / "data.bam"
    p.write_bytes(b"data")
    fe = FilesElement(p)
    ok, reason = fe.validate()
    assert ok is True


def test_files_element_md5sum(tmp_path):
    p = tmp_path / "data.bam"
    p.write_bytes(b"hello")
    fe = FilesElement(p)
    md5 = fe.md5sum
    assert md5 is not None
    assert len(md5) > 0


def test_files_element_calc_md5sum_missing(tmp_path):
    p = tmp_path / "missing.bam"
    fe = FilesElement(str(p))
    # file does not exist → empty string
    assert fe.calc_md5sum() == ""


def test_files_element_with_tag(tmp_path):
    p = tmp_path / "data.vcf"
    p.write_bytes(b"data")
    tag = PartialElementTag(stage=Stage.CALL, method=Method.GATK, state=State.RAW)
    fe = FilesElement(p, tag=tag)
    assert fe.ext == "vcf"


def test_files_element_with_root(tmp_path):
    p = tmp_path / "data.bam"
    p.write_bytes(b"data")
    fe = FilesElement(p, root="my_root")
    assert "my_root" in fe.key


# ---------------------------------------------------------------------------
# FileElement
# ---------------------------------------------------------------------------


def test_file_element(tmp_path):
    p = tmp_path / "sample.bam"
    p.write_bytes(b"content")
    fe = FileElement(p)
    assert fe.ext == "bam"


def test_file_element_with_root(tmp_path):
    p = tmp_path / "sample.bam"
    p.write_bytes(b"content")
    fe = FileElement(p, root="custom_root")
    assert "custom_root" in fe.key


# ---------------------------------------------------------------------------
# Sample (FilesElement subclass)
# ---------------------------------------------------------------------------


def test_sample_element(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"@read\nACGT\n+\nIIII\n")
    s = Sample(p)
    assert p.absolute() in s.files


def test_sample_element_with_name(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"data")
    s = Sample(p, name="my_sample")
    # name parameter is accepted (passed to super via kwargs)
    assert s is not None


# ---------------------------------------------------------------------------
# NextGenSampleElement
# ---------------------------------------------------------------------------


def test_next_gen_sample_single(tmp_path):
    p = tmp_path / "sample_R1.fastq"
    p.write_bytes(b"data")
    ns = NextGenSampleElement(p, root="sample1")
    assert ns.pairing.value == "single"
    assert len(ns.input_files) == 1


def test_next_gen_sample_paired(tmp_path):
    p1 = tmp_path / "r1.fastq"
    p2 = tmp_path / "r2.fastq"
    p1.write_bytes(b"data")
    p2.write_bytes(b"data")
    ns = NextGenSampleElement({"r1": p1, "r2": p2}, root="sample1")
    assert ns.pairing.value == "paired"
    assert len(ns.input_files) == 2


def test_next_gen_sample_default_dirs(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"data")
    ns = NextGenSampleElement(p, root="s1")
    assert "s1" in str(ns.cache_dir) or "cache" in str(ns.cache_dir)


def test_next_gen_sample_custom_dirs(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"data")
    cache = tmp_path / "my_cache"
    result = tmp_path / "my_result"
    ns = NextGenSampleElement(p, root="s1", cache_dir=cache, result_dir=result)
    assert ns.cache_dir == cache
    assert ns.result_dir == result


def test_next_gen_sample_read_group(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"data")
    ns = NextGenSampleElement(p, root="s1", read_group="@RG\tID:1")
    assert ns.read_group == "@RG\tID:1"


def test_next_gen_sample_reverse_reads(tmp_path):
    p = tmp_path / "sample.fastq"
    p.write_bytes(b"data")
    ns = NextGenSampleElement(p, root="s1", reverse_reads=True)
    assert ns.reverse_reads is True


# ---------------------------------------------------------------------------
# sample_fastqs
# ---------------------------------------------------------------------------


def test_sample_fastqs_from_next_gen_single(tmp_path):
    p = tmp_path / "r1.fastq"
    p.write_bytes(b"data")
    ns = NextGenSampleElement(p, root="s1")
    r1, r2, name, rg = sample_fastqs(ns)
    assert r1 == p.absolute()
    assert r2 is None
    assert name == ns.name


def test_sample_fastqs_from_next_gen_paired(tmp_path):
    p1 = tmp_path / "r1.fastq"
    p2 = tmp_path / "r2.fastq"
    p1.write_bytes(b"data")
    p2.write_bytes(b"data")
    ns = NextGenSampleElement({"r1": p1, "r2": p2}, root="s1")
    r1, r2, name, rg = sample_fastqs(ns)
    assert r1 is not None
    assert r2 is not None


def test_sample_fastqs_from_plain_element(tmp_path):
    p1 = tmp_path / "r1.fastq"
    p2 = tmp_path / "r2.fastq"
    p1.write_bytes(b"data")
    p2.write_bytes(b"data")
    tag = make_tag()
    e = Element(
        key="plain_el",
        run=lambda: None,
        tag=tag,
        artifacts={
            "fastq_r1": p1,
            "fastq_r2": p2,
            "sample_name": "my_sample",
            "read_group": "rg1",
        },
    )
    r1, r2, name, rg = sample_fastqs(e)
    assert r1 == p1.absolute()
    assert r2 == p2.absolute()
    assert name == "my_sample"
    assert rg == "rg1"


def test_sample_fastqs_from_plain_element_no_r2(tmp_path):
    p1 = tmp_path / "r1.fastq"
    p1.write_bytes(b"data")
    tag = make_tag()
    e = Element(
        key="plain_no_r2",
        run=lambda: None,
        tag=tag,
        artifacts={"fastq_r1": p1},
    )
    r1, r2, name, rg = sample_fastqs(e)
    assert r1 == p1.absolute()
    assert r2 is None


# ---------------------------------------------------------------------------
# generate_element_key_name
# ---------------------------------------------------------------------------


def test_generate_element_key_name_minimal():
    tag = make_tag(root="sample", level=1)
    key, short = generate_element_key_name(tag, "bwa", "2.2.1")
    assert "bwa" in key
    assert "2.2.1" in key
    assert "bwa" in short


def test_generate_element_key_name_with_subcommand():
    tag = make_tag()
    key, short = generate_element_key_name(
        tag, "gatk", "4.6", subcommand="HaplotypeCaller"
    )
    assert "HaplotypeCaller" in key
    assert "HaplotypeCaller" in short


def test_generate_element_key_name_with_suffix():
    tag = make_tag()
    key, short = generate_element_key_name(tag, "tool", "1.0", suffix="extra")
    assert "extra" in key


def test_generate_element_key_name_with_params():
    tag = make_tag()
    key, short = generate_element_key_name(tag, "tool", "1.0", min_mapq="30")
    assert "min_mapq-30" in key


def test_generate_element_key_name_params_none_excluded():
    tag = make_tag()
    key, short = generate_element_key_name(tag, "tool", "1.0", min_mapq=None)
    assert "min_mapq" not in key


# ---------------------------------------------------------------------------
# Element with pres – signature depends on pres
# ---------------------------------------------------------------------------


def test_element_sig_includes_pre_sigs():
    pre = make_element(key="pre")
    e = make_element(key="child", pres=(pre,))
    data = e.sig_data()
    assert pre.signature in data["pre_sigs"]


# ---------------------------------------------------------------------------
# element_build_context with registry intern (strict key uniqueness)
# ---------------------------------------------------------------------------


def test_registry_intern_same_key_different_element_raises():
    tag = make_tag()
    registry = ElementRegistry()

    @element
    def build():
        return Element(key="shared_key", run=lambda: None, tag=tag)

    with pytest.raises(ValueError, match="Duplicate key collision"):
        with element_build_context(registry):
            build()
            build()


def test_registry_intern_key_collision_different_sig():
    tag1 = make_tag(root="r1")
    tag2 = make_tag(root="r2")
    registry = ElementRegistry()

    with pytest.raises(ValueError, match="Duplicate key collision"):
        with element_build_context(registry):
            e1 = Element(key="collision", run=lambda: None, tag=tag1)
            registry.intern(e1)
            e2 = Element(key="collision", run=lambda: None, tag=tag2)
            registry.intern(e2)
