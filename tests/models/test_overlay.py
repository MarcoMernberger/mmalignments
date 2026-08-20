from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import pytest

from mmalignments.models.overlay import (
    E,
    O,
    P,
    T,
    ElementTag,
    ExternalRunConfig,
    FileSpec,
    OutputSpec,
    Overlay,
    Params,
    PartialElementTag,
    PartialExternalRunConfig,
    PartialOutputSpec,
)
from mmalignments.models.tags import Method, Omics, Stage, State


@dataclass(frozen=True)
class DemoOverlay(Overlay):
    a: int | None = None
    b: str | None = None

    def resolve(self, *patches: Overlay | None) -> DemoOverlay:
        result: DemoOverlay = self
        for patch in patches:
            result = result.patch(patch)
        return result


@dataclass(frozen=True)
class NoResolveOverlay(Overlay):
    value: int = 1


def _tag(
    *,
    root: str = "sample",
    level: int = 1,
    stage: Stage = Stage.PREP,
    method: Method = Method.CUSTOM,
    state: State = State.RAW,
    omics: Omics | None = Omics.DNA,
    param: str | None = "flag",
) -> ElementTag:
    return ElementTag(
        root=root,
        level=level,
        stage=stage,
        method=method,
        state=state,
        omics=omics,
        param=param,
    )


def test_overlay_dict_and_mapping_helpers() -> None:
    overlay = DemoOverlay(a=1, b=None)

    assert overlay.to_dict() == {"a": 1}
    assert overlay.to_dict(drop_none=False) == {"a": 1, "b": None}
    assert list(overlay.keys()) == ["a", "b"]
    assert list(overlay.values()) == [1, None]
    assert list(overlay.items()) == [("a", 1), ("b", None)]
    assert overlay["a"] == 1
    assert list(iter(overlay)) == ["a", "b"]
    assert "a" in overlay
    assert "missing" not in overlay


def test_overlay_update_patch_and_resolve_order() -> None:
    base = DemoOverlay(a=1, b="x")

    assert base.update(a=2) == DemoOverlay(a=2, b="x")
    assert base.patch(None) is base

    # None values from patch are dropped by to_dict(drop_none=True)
    patched = base.patch(DemoOverlay(a=None, b="y"))
    assert patched == DemoOverlay(a=1, b="y")

    resolved = base.resolve(DemoOverlay(a=5), None, DemoOverlay(b="z"))
    assert resolved == DemoOverlay(a=5, b="z")


def test_overlay_resolve_base_raises_not_implemented() -> None:
    with pytest.raises(NotImplementedError):
        NoResolveOverlay().resolve()


def test_element_tag_validation_and_default_name_format() -> None:
    with pytest.raises(ValueError, match="root"):
        ElementTag(  # type: ignore[arg-type]
            root=None,
            level=1,
            stage=Stage.PREP,
            method=Method.CUSTOM,
            state=State.RAW,
        )

    tag = _tag()
    expected = ".".join(
        [
            tag.root,
            tag.format_level(tag.level),
            tag.omics,
            tag.stage,
            tag.method,
            tag.state,
            tag.param,
        ]
    )
    assert tag.format_level(3) == "S03"
    assert tag.default_name == expected

    tag_without_optional = _tag(omics=None, param=None)
    expected_no_optional = ".".join(
        [
            tag_without_optional.root,
            tag_without_optional.format_level(tag_without_optional.level),
            tag_without_optional.stage,
            tag_without_optional.method,
            tag_without_optional.state,
        ]
    )
    assert tag_without_optional.default_name == expected_no_optional


def test_element_tag_resolve_and_from_prior_variants() -> None:
    prior = _tag(level=4, stage=Stage.ALIGN, method=Method.BWAMEM2, state=State.MAP)

    resolved = prior.resolve(None, PartialElementTag(state=State.SORT))
    assert resolved.state == State.SORT
    assert resolved.level == prior.level

    derived = ElementTag.from_prior(
        prior,
        tag=PartialElementTag(state=State.FILTER),
    )
    assert derived.root == prior.root
    assert derived.level == prior.level + 1
    assert derived.stage == prior.stage
    assert derived.method == prior.method
    assert derived.omics == prior.omics
    assert derived.state == State.FILTER
    assert derived.param == "flag"  # no override by None!

    derived_overridden = ElementTag.from_prior(
        prior,
        tag=None,
        state=State.COUNT,
        level=9,
        method=Method.CUSTOM,
    )
    assert derived_overridden.state == State.COUNT
    assert derived_overridden.level == 9
    assert derived_overridden.method == Method.CUSTOM


def test_partial_element_tag_inherits_overlay_patch_semantics() -> None:
    base = PartialElementTag(root="r", level=1)
    patch = PartialElementTag(root=None, stage=Stage.PREP)

    patched = base.patch(patch)
    assert patched.root == "r"
    assert patched.level == 1
    assert patched.stage == Stage.PREP


def test_params_mapping_and_patching() -> None:
    params = Params(a=1, b=None)

    assert params["a"] == 1
    assert list(iter(params)) == ["a", "b"]
    assert len(params) == 2
    assert params.to_dict() == {"a": 1}
    assert params.to_dict(drop_none=False) == {"a": 1, "b": None}
    assert params.get("a") == 1
    assert params.get("missing", "fallback") == "fallback"
    assert params.patch(None) is params

    patched = params.patch(Params(a=3, c=4))
    assert patched.to_dict(drop_none=False) == {"a": 3, "b": None, "c": 4}

    resolved = params.resolve(None, Params(a=5), Params(c=None))
    assert resolved.to_dict(drop_none=False) == {"a": 5, "b": None, "c": None}
    assert resolved.to_dict() == {"a": 5}
    assert isinstance(params, Mapping)


def test_output_specs_and_filespec_resolution() -> None:
    file_spec = FileSpec(name="results", ext="tsv")
    assert file_spec.name == "results"
    assert file_spec.ext == "tsv"

    base = OutputSpec(
        stem="base",
        outdir=Path("/tmp"),
        prefix="pre_",
        suffix="_suf",
        ext="parquet",
        exts=("parquet", "tsv"),
        additional_output={"tsv": file_spec},
    )
    resolved = base.resolve(None, PartialOutputSpec(prefix="new_", ext="csv"))

    assert resolved.stem == "base"
    assert resolved.prefix == "new_"
    assert resolved.ext == "csv"
    assert resolved.additional_output["tsv"] == file_spec


def test_external_run_config_validation_and_resolve() -> None:
    with pytest.raises(ValueError, match="stdout and stderr"):
        ExternalRunConfig(stdout=Path("same.log"), stderr=Path("same.log"))

    with pytest.raises(ValueError, match="log_dir must be a directory"):
        ExternalRunConfig(log_dir="not-a-path")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="threads must be a positive integer"):
        ExternalRunConfig(threads=0)

    config = ExternalRunConfig()
    assert config.threads >= 1

    base = ExternalRunConfig(cwd=Path("/tmp"), threads=4, multi=True)
    resolved = base.resolve(None, PartialExternalRunConfig(threads=2, append=True))

    assert resolved.cwd == Path("/tmp")
    assert resolved.threads == 2
    assert resolved.multi is True
    assert resolved.append is True


def test_shortcut_aliases_point_to_expected_classes() -> None:
    assert T is PartialElementTag
    assert O is PartialOutputSpec
    assert P is Params
    assert E is PartialExternalRunConfig
