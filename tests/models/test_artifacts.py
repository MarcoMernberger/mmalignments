from __future__ import annotations

import importlib
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
MMALIGNMENTS_ROOT = SRC_ROOT / "mmalignments"


def _ensure_namespace_package(name: str, path: Path) -> None:
    """Create a lightweight namespace package in sys.modules if missing."""
    if name in sys.modules:
        return

    module = types.ModuleType(name)
    module.__path__ = [str(path)]
    sys.modules[name] = module


_ensure_namespace_package("mmalignments", MMALIGNMENTS_ROOT)
_ensure_namespace_package("mmalignments.models", MMALIGNMENTS_ROOT / "models")
_ensure_namespace_package("mmalignments.services", MMALIGNMENTS_ROOT / "services")

artifacts_module = importlib.import_module("mmalignments.models.artifacts")
dependencies_module = importlib.import_module("mmalignments.services.dependencies")

Artifact = artifacts_module.Artifact
ArtifactSet = artifacts_module.ArtifactSet
FastqArtifact = artifacts_module.FastqArtifact
FileArtifact = artifacts_module.FileArtifact
TableArtifact = artifacts_module.TableArtifact

combined_signature = dependencies_module.combined_signature
file_signature = dependencies_module.file_signature
stable_hash = dependencies_module.stable_hash


class DummyArtifact(Artifact):
    """Concrete Artifact subclass used to exercise the abstract base type."""


class BareFileArtifact(FileArtifact):
    """Deliberately incomplete file artifact to test base-class recursion behavior."""


class ResolvablePath:
    """Simple helper object exposing resolve() so _infer_primary_name can inspect it."""

    def __init__(self, path: Path):
        self.path = path

    def resolve(self) -> Path:
        return self.path


class ResolvableNonPath:
    """Helper object whose resolve() result is not a Path."""

    def resolve(self) -> str:
        return "not-a-path"


class HasFiles:
    """Container with a files attribute to mimic file-like artifacts in ArtifactSet."""

    def __init__(self, *paths: Path):
        self.files = tuple(path.resolve() for path in paths)


def test_artifact_base_class_can_be_subclassed():
    """Validate that the Artifact base class is usable via normal subclassing."""
    # This test intentionally keeps scope narrow: it checks that the abstract base
    # class is structurally sound and can be instantiated through a concrete child.
    artifact = DummyArtifact()

    # We assert both the child and base relationships so coverage includes both the
    # concrete helper and the base class behavior in one easy-to-read assertion set.
    assert isinstance(artifact, DummyArtifact)
    assert isinstance(artifact, Artifact)


def test_fileartifact_default_protocol_is_recursively_defined():
    """Exercise the default FileArtifact implementation and document recursion failure."""
    # The default FileArtifact.files implementation iterates over self, while
    # __iter__ iterates over self.files; this causes recursion unless subclasses
    # override at least one of the two. We validate this behavior explicitly.
    artifact = BareFileArtifact()

    # Accessing files on the incomplete subclass should recurse until Python raises.
    with pytest.raises(RecursionError):
        _ = artifact.files

    # Signature depends on files, so it should fail with the same recursion issue.
    with pytest.raises(RecursionError):
        _ = artifact.signature

    # Explicit iteration also follows the same recursion path and should therefore
    # produce the same error mode.
    with pytest.raises(RecursionError):
        next(iter(artifact))


@pytest.mark.parametrize(
    ("filename", "expected_stem"),
    [
        ("sample.fq", "sample"),
        ("sample.fastq", "sample"),
        ("sample.fq.gz", "sample"),
        ("sample.fastq.gz", "sample"),
        ("sample.custom", "sample.custom"),
    ],
)
def test_fastqartifact_stem_and_pairing_logic(
    tmp_path: Path, filename: str, expected_stem: str
):
    """Verify FASTQ stem stripping across all configured suffixes and non-matching names."""
    # We write a tiny dummy file so that all path operations run on real filesystem
    # entries, which mirrors actual production usage.
    r1 = tmp_path / filename
    r1.write_text("@read\nACGT\n+\n!!!!\n", encoding="utf-8")

    # Single-end artifacts should report unpaired and expose exactly one file.
    artifact = FastqArtifact(r1=r1)
    assert artifact.paired is False
    assert artifact.files == (r1.resolve(),)
    assert artifact.stem == expected_stem


def test_fastqartifact_paired_files_iter_signature_and_str(tmp_path: Path):
    """Cover paired-end file resolution, iteration, signature hashing, and string rendering."""
    # Two FASTQ files are created so that all paired-end code paths can execute with
    # concrete inputs and deterministic signatures.
    r1 = tmp_path / "pair_R1.fastq"
    r2 = tmp_path / "pair_R2.fastq.gz"
    r1.write_text("@r1\nAAAA\n+\n####\n", encoding="utf-8")
    r2.write_text("@r2\nCCCC\n+\n####\n", encoding="utf-8")

    artifact = FastqArtifact(r1=r1, r2=r2)

    # Paired status and resolved files are central invariants for this type.
    assert artifact.paired is True
    assert artifact.files == (r1.resolve(), r2.resolve())

    # Iteration is inherited from FileArtifact and must yield the same resolved files.
    assert tuple(iter(artifact)) == artifact.files

    # Signature should match the exact hashing contract from FileArtifact.signature.
    expected_signature = stable_hash(
        tuple(file_signature(path) for path in artifact.files)
    )
    assert artifact.signature == expected_signature

    # String representation includes both file references and follows the documented format.
    assert str(artifact) == f"FastqArtifact({r1}\n{r2}\n)"


def test_tableartifact_signature_files_stem_and_resolve(tmp_path: Path):
    """Verify the artifact protocol methods for tabular file wrappers."""
    # This test focuses on deterministic protocol-level properties that downstream
    # components use to manage dependencies and output identities.
    csv_path = tmp_path / "counts.tsv"
    csv_path.write_text("sample\tvalue\nA\t1\n", encoding="utf-8")

    artifact = TableArtifact(path=csv_path)

    # All protocol values should map directly to the wrapped path as implemented.
    assert artifact.signature == file_signature(csv_path)
    assert artifact.files == (csv_path.resolve(),)
    assert artifact.stem == "counts"
    assert artifact.resolve() == csv_path.resolve()


def test_tableartifact_frame_property_reads_table_and_drops_unnamed_columns(
    tmp_path: Path,
):
    """Ensure frame uses read_frame defaults and therefore removes unnamed index columns."""
    # Writing with index=True intentionally creates an unnamed first column so that
    # read_frame's default usecols filter has a visible effect.
    csv_path = tmp_path / "table.csv"
    pd.DataFrame({"sample": ["s1", "s2"], "value": [10, 20]}).to_csv(
        csv_path,
        index=True,
    )

    artifact = TableArtifact(path=csv_path)
    frame = artifact.frame

    # The unnamed column should be removed while real data columns remain intact.
    assert "Unnamed: 0" not in frame.columns
    assert list(frame.columns) == ["sample", "value"]


def test_tableartifact_read_passes_index_col_for_csv_and_tsv(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Confirm read forwards index_col for CSV/TSV inputs when an index column is requested."""
    # We monkeypatch read_frame to capture call kwargs, which makes the behavior
    # under test explicit without relying on pandas internals.
    captured: dict[str, object] = {}

    def fake_read_frame(path: Path, drop_unnamed_columns: bool = True, **kwargs):
        captured["path"] = path
        captured["drop_unnamed_columns"] = drop_unnamed_columns
        captured["kwargs"] = kwargs
        return pd.DataFrame(
            {"value": [1, 2]}, index=pd.Index(["a", "b"], name="sample")
        )

    monkeypatch.setattr(artifacts_module, "read_frame", fake_read_frame)

    csv_path = tmp_path / "table.csv"
    csv_path.write_text("sample,value\na,1\nb,2\n", encoding="utf-8")
    artifact = TableArtifact(path=csv_path)

    frame = artifact.read(drop_unnamed_columns=False, index_column="sample")

    # For CSV/TSV, index_col must be forwarded so pandas can build the index directly.
    assert captured["path"] == csv_path.resolve()
    assert captured["drop_unnamed_columns"] is False
    assert captured["kwargs"] == {"index_col": "sample"}

    # Because sample is already the index in returned frame, no second set_index call occurs.
    assert frame.index.name == "sample"
    assert "sample" not in frame.columns


def test_tableartifact_read_sets_index_post_read_for_non_csv_suffix(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Confirm read applies set_index after load for non-CSV/TSV formats when needed."""
    # This branch exists because index_col is only passed for CSV/TSV, so for other
    # suffixes the method must set the index manually if the column is still present.
    captured: dict[str, object] = {}

    def fake_read_frame(path: Path, drop_unnamed_columns: bool = True, **kwargs):
        captured["path"] = path
        captured["drop_unnamed_columns"] = drop_unnamed_columns
        captured["kwargs"] = kwargs
        return pd.DataFrame({"sample": ["x", "y"], "value": [5, 6]})

    monkeypatch.setattr(artifacts_module, "read_frame", fake_read_frame)

    parquet_path = tmp_path / "table.parquet"
    parquet_path.write_text("not-used", encoding="utf-8")
    artifact = TableArtifact(path=parquet_path)

    frame = artifact.read(index_column="sample")

    # Non-CSV/TSV means no index_col kwarg should be forwarded during read_frame call.
    assert captured["path"] == parquet_path.resolve()
    assert captured["drop_unnamed_columns"] is True
    assert captured["kwargs"] == {}

    # Index should be moved from a normal column to the DataFrame index by TableArtifact.read.
    assert frame.index.name == "sample"
    assert "sample" not in frame.columns


def test_tableartifact_frame_cached_property_reads_only_once(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Demonstrate cached_property semantics for frame to avoid duplicate I/O work."""
    # Caching is a key performance property, so we count invocations directly.
    calls = {"count": 0}

    def fake_read_frame(path: Path, drop_unnamed_columns: bool = True, **kwargs):
        calls["count"] += 1
        return pd.DataFrame({"value": [42]})

    monkeypatch.setattr(artifacts_module, "read_frame", fake_read_frame)

    path = tmp_path / "cached.csv"
    path.write_text("value\n42\n", encoding="utf-8")
    artifact = TableArtifact(path=path)

    first = artifact.frame
    second = artifact.frame

    # Both accesses should return equivalent content with exactly one backend read.
    assert first.equals(second)
    assert calls["count"] == 1


def test_artifactset_init_rejects_reserved_primary_key():
    """Ensure constructor forbids using reserved name primary inside extras."""
    # The reserved-key contract prevents ambiguous lookups and must fail early.
    with pytest.raises(KeyError, match="reserved"):
        ArtifactSet(Path("a.bam"), primary=Path("x"))


def test_artifactset_init_rejects_primary_name_collision_with_extras():
    """Ensure inferred or provided primary_name cannot collide with extras."""
    # Here we force a collision to verify the protective constructor check.
    with pytest.raises(KeyError, match="collides"):
        ArtifactSet(Path("a.bam"), primary_name="bam", bam=Path("other.bam"))


def test_artifactset_basic_mapping_behavior_properties_and_repr(tmp_path: Path):
    """Cover mapping protocol methods, aliases, properties, and both repr branches."""
    # We prepare small but realistic artifacts so files/aliases/contains can be checked
    # in the same assertion block with readable intent.
    primary = tmp_path / "sample.bam"
    secondary = tmp_path / "sample.tsv"
    primary.write_text("bam", encoding="utf-8")
    secondary.write_text("id\tvalue\n1\t2\n", encoding="utf-8")

    aset = ArtifactSet(primary, tsv=secondary)

    # Primary name should be inferred from path suffix and work with both aliases.
    assert aset.primary == primary
    assert aset.primary_name == "bam"
    assert aset["primary"] == primary
    assert aset["bam"] == primary
    assert aset["tsv"] == secondary

    # Mapping behavior: iteration order (primary then extras), length, and contains logic.
    assert list(aset) == ["bam", "tsv"]
    assert len(aset) == 2
    assert "primary" in aset
    assert "bam" in aset
    assert "tsv" in aset
    assert "missing" not in aset

    # extras property should expose the immutable mapping and include only extras.
    assert dict(aset.extras) == {"tsv": secondary}

    # Representation should include extras when present and omit them otherwise.
    assert repr(aset) == f"ArtifactSet(bam={primary!r}, tsv={secondary!r})"
    assert repr(ArtifactSet(primary)) == f"ArtifactSet(bam={primary!r})"


def test_artifactset_signature_data_signature_and_calculation(tmp_path: Path):
    """Validate signature_data sorting plus both signature computation entry points."""
    # The order of signature_data matters because deterministic hashing depends on
    # stable key sorting across runs.
    primary = tmp_path / "x.bam"
    extra_a = tmp_path / "a.txt"
    extra_z = tmp_path / "z.txt"
    aset = ArtifactSet(primary, zeta=extra_z, alpha=extra_a)

    signature_data = aset.signature_data

    # Keys should be sorted lexicographically by dict(sorted(self.items())).
    assert list(signature_data.keys()) == ["alpha", "bam", "zeta"]

    expected = combined_signature(*signature_data.values())

    # calculate_signature and cached signature property must agree bit-for-bit.
    assert aset.calculate_signature() == expected
    assert aset.signature == expected


def test_artifactset_files_and_iterfiles_collect_paths_from_supported_types(
    tmp_path: Path,
):
    """Check that files aggregates Path values and objects exposing a files tuple."""
    # This test covers both collection branches and verifies unsupported values are ignored.
    p1 = tmp_path / "a.fastq"
    p2 = tmp_path / "b.fastq"
    p3 = tmp_path / "c.tsv"
    for path in (p1, p2, p3):
        path.write_text("x", encoding="utf-8")

    file_like = HasFiles(p1, p2)
    aset = ArtifactSet(p3, file_like=file_like, ignored=123)

    expected_files = (p3.resolve(), p1.resolve(), p2.resolve())

    # files should include resolved paths from primary Path and from the files-enabled helper.
    assert aset.files == expected_files

    # iterfiles is a straightforward generator wrapper and should yield in the same order.
    assert tuple(aset.iterfiles()) == expected_files


def test_infer_primary_name_covers_path_and_resolve_variants(tmp_path: Path):
    """Exercise all _infer_primary_name branches including fallback behavior."""
    # Direct Path values with suffix should yield the suffix minus the leading dot.
    assert ArtifactSet._infer_primary_name(tmp_path / "sample.parquet") == "parquet"

    # Objects with resolve() returning Path should be treated exactly like Paths.
    assert (
        ArtifactSet._infer_primary_name(ResolvablePath(tmp_path / "reads.fastq"))
        == "fastq"
    )

    # Paths with no suffix and non-Path resolve() values should both fall back to primary.
    assert ArtifactSet._infer_primary_name(tmp_path / "README") == "primary"
    assert ArtifactSet._infer_primary_name(ResolvableNonPath()) == "primary"

    # Plain objects without resolve() should also follow the fallback branch.
    assert ArtifactSet._infer_primary_name(123) == "primary"


def test_from_any_identity_explicit_primary_candidate_and_fallback_selection(
    tmp_path: Path,
):
    """Cover all from_any selection strategies and identity behavior."""
    # Identity path: passing an ArtifactSet should return the same object unchanged.
    original = ArtifactSet(tmp_path / "sample.bam", txt=tmp_path / "meta.txt")
    assert ArtifactSet.from_any(original) is original

    # Explicit primary key path: custom key is used as both primary and primary_name.
    explicit = ArtifactSet.from_any({"main": 1, "aux": 2}, primary_key="main")
    assert explicit.primary_name == "main"
    assert explicit.primary == 1
    assert explicit["aux"] == 2

    # Candidate detection path: bam should be preferred when no explicit primary key exists.
    candidate = ArtifactSet.from_any({"foo": 1, "bam": 2, "bar": 3})
    assert candidate.primary_name == "bam"
    assert candidate.primary == 2

    # Later-candidate detection: ensure the scan continues beyond the first choices.
    late_candidate = ArtifactSet.from_any({"foo": 1, "fastq": 7, "bar": 3})
    assert late_candidate.primary_name == "fastq"
    assert late_candidate.primary == 7

    # Fallback path: if no candidate names exist, the first mapping item becomes primary.
    fallback = ArtifactSet.from_any({"first": "A", "second": "B"})
    assert fallback.primary_name == "first"
    assert fallback.primary == "A"


def test_from_any_rejects_non_mapping_inputs():
    """Ensure from_any raises a clear TypeError for unsupported input types."""
    # The method contract accepts only ArtifactSet or Mapping, so a plain integer
    # should trigger a deterministic type error message.
    with pytest.raises(TypeError, match="Cannot create ArtifactSet"):
        ArtifactSet.from_any(5)  # type: ignore[arg-type]


def test_with_extra_success_and_conflict_detection(tmp_path: Path):
    """Verify with_extra can append new keys and rejects existing names."""
    # Start from a minimal set and add one valid extra.
    primary = tmp_path / "sample.bam"
    aset = ArtifactSet(primary)
    extended = aset.with_extra("stats", tmp_path / "stats.tsv")
    assert "stats" in extended
    assert extended["stats"] == tmp_path / "stats.tsv"

    # Attempting to add an existing key should fail to prevent accidental overwrite.
    with pytest.raises(KeyError, match="already exists"):
        _ = extended.with_extra("stats", "dup")


def test_with_extras_success_and_conflict_detection(tmp_path: Path):
    """Verify with_extras handles batch additions and conflict reporting."""
    # Successful branch with two new keys in one call.
    primary = tmp_path / "sample.bam"
    aset = ArtifactSet(primary)
    enriched = aset.with_extras({"a": 1, "b": 2})
    assert enriched["a"] == 1
    assert enriched["b"] == 2

    # Conflict branch should list clashing keys in the raised KeyError.
    with pytest.raises(KeyError, match="already exist"):
        _ = enriched.with_extras({"b": 3, "c": 4})


def test_merge_covers_mapping_conversion_primary_preservation_and_equal_primary_name(
    tmp_path: Path,
):
    """Cover merge behavior for mapping conversion and primary-name equality branch."""
    # Base set keeps its primary by contract regardless of merged content.
    base = ArtifactSet(tmp_path / "base.bam", left=tmp_path / "left.txt")

    # Merge from a plain mapping to exercise ArtifactSet.from_any conversion inside merge.
    merged_from_mapping = base.merge({"parquet": tmp_path / "out.parquet", "qc": "ok"})
    assert merged_from_mapping.primary == base.primary
    assert merged_from_mapping.primary_name == base.primary_name
    assert merged_from_mapping["parquet"] == tmp_path / "out.parquet"
    assert merged_from_mapping["qc"] == "ok"

    # Equal primary names should avoid re-adding other.primary as an extra artifact.
    other_same_primary_name = ArtifactSet(tmp_path / "second.bam", extra=99)
    merged_same_name = base.merge(other_same_primary_name)
    assert "bam" not in dict(merged_same_name.extras)
    assert merged_same_name["extra"] == 99


def test_merge_detects_primary_name_and_extra_conflicts(tmp_path: Path):
    """Exercise both merge conflict branches for robust collision handling."""
    # First conflict branch: other.primary_name collides with an existing extra key.
    base = ArtifactSet(tmp_path / "base.bam", parquet=tmp_path / "already.parquet")
    other_primary_collision = ArtifactSet(
        tmp_path / "new.parquet",
        primary_name="parquet",
        fresh=1,
    )
    with pytest.raises(KeyError, match="already exists"):
        _ = base.merge(other_primary_collision)

    # Second conflict branch: extras overlap after primary-handling logic.
    base2 = ArtifactSet(tmp_path / "base2.bam", shared="left")
    other_extra_collision = ArtifactSet(
        tmp_path / "other.fastq",
        primary_name="fastq",
        shared="right",
    )
    with pytest.raises(KeyError, match="already exist"):
        _ = base2.merge(other_extra_collision)
