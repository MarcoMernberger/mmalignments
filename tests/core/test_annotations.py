# tests/test_annotations.py

import pytest

from mmalignments.core.annotations import (
    ColumnTag,
    View,
    _legacy_to_tag,
    describe_metric,
)


def test_describe_metric_known():
    assert describe_metric("gu") == "Gene unstranded count"


def test_describe_metric_unknown():
    assert describe_metric("foo") == "(unknown metric segment: 'foo')"


@pytest.mark.parametrize(
    "colname, expected",
    [
        (
            "gene_stable_id",
            ColumnTag(
                metric=("gene_stable_id",),
                sample_id=None,
                pipeline="META",
            ),
        ),
        (
            "chr",
            ColumnTag(
                metric=("chr",),
                sample_id=None,
                pipeline="META",
            ),
        ),
    ],
)
def test_legacy_meta_columns(colname, expected):
    assert _legacy_to_tag(colname, []) == expected


def test_legacy_count_gene_unstranded():
    tag = _legacy_to_tag(
        "gene unstranded tag count sample1",
        ["sample1"],
    )

    assert tag == ColumnTag(
        metric=("Count", "GU"),
        sample_id="sample1",
        pipeline=None,
    )


def test_legacy_count_gene_unstranded_dedup():
    tag = _legacy_to_tag(
        "gene unstranded tag count sample1_STAR_DEDUP",
        ["sample1"],
    )

    assert tag.metric == ("Count", "GU", "Dedup")
    assert tag.sample_id == "sample1"


def test_legacy_count_exon_stranded():
    tag = _legacy_to_tag(
        "exon, protein coding, stranded tag count S1",
        ["S1"],
    )

    assert tag.metric == ("Count", "ES")


def test_legacy_count_exon_stranded_dedup():
    tag = _legacy_to_tag(
        "exon, protein coding, stranded tag count S1_dedup",
        ["S1"],
    )

    assert tag.metric == ("Count", "ES", "Dedup")


def test_legacy_unknown_column_returns_none():
    assert (
        _legacy_to_tag(
            "something completely different",
            ["S1"],
        )
        is None
    )


def test_legacy_unknown_sample_raises():
    with pytest.raises(ValueError):
        _legacy_to_tag(
            "gene unstranded tag count unknown",
            ["S1"],
        )


def test_legacy_unknown_category_raises():
    with pytest.raises(ValueError):
        _legacy_to_tag(
            "unknown category tag count S1",
            ["S1"],
        )


@pytest.mark.parametrize(
    "column, expected",
    [
        (
            "Count.GU (S1)",
            ColumnTag(
                metric=("Count", "GU"),
                sample_id="S1",
                pipeline=None,
            ),
        ),
        (
            "Count.GU (S1) [STAR]",
            ColumnTag(
                metric=("Count", "GU"),
                sample_id="S1",
                pipeline="STAR",
            ),
        ),
        (
            "VST",
            ColumnTag(
                metric=("VST",),
                sample_id=None,
                pipeline=None,
            ),
        ),
    ],
)
def test_decode(column, expected):
    assert ColumnTag.decode(column) == expected


def test_decode_strips_whitespace():
    tag = ColumnTag.decode(" Count.GU ( S1 ) [ STAR ] ")

    assert tag.metric == ("Count", "GU")
    assert tag.sample_id == "S1"
    assert tag.pipeline == "STAR"


def test_decode_fallback():
    tag = ColumnTag.decode("not canonical")

    assert tag.metric == ("not canonical",)
    assert tag.sample_id is None


def test_encode():
    tag = ColumnTag(
        metric=("Count", "GU"),
        sample_id="S1",
        pipeline="STAR",
    )

    assert tag.encode() == "Count.GU (S1) [STAR]"


def test_encode_without_optional_fields():
    tag = ColumnTag(metric=("VST",))

    assert tag.encode() == "VST"


def test_generate_new_tag_append():
    old = ColumnTag(
        metric=("Count", "GU"),
        sample_id="S1",
        pipeline="STAR",
    )

    view = View(metric=("Normalized",))

    result = ColumnTag.generate_new_tag(
        old,
        view,
    )

    assert result.metric == (
        "Count",
        "GU",
        "Normalized",
    )


def test_generate_new_tag_replace():
    old = ColumnTag(
        metric=("Count",),
        sample_id="S1",
    )

    result = ColumnTag.generate_new_tag(
        old,
        View(metric=("VST",)),
        append=False,
    )

    assert result.metric == ("VST",)


def test_generate_new_tag_metric_none():
    old = ColumnTag(metric=("Count",))

    result = ColumnTag.generate_new_tag(
        old,
        View(metric=None),
    )

    assert result.metric == ("Count",)


def test_generate_new_tag_sample_override():
    old = ColumnTag(
        metric=("Count",),
        sample_id="A",
    )

    result = ColumnTag.generate_new_tag(
        old,
        View(sample_id="B"),
    )

    assert result.sample_id == "B"


def test_generate_new_tag_pipeline_override():
    old = ColumnTag(
        metric=("Count",),
        pipeline="OLD",
    )

    result = ColumnTag.generate_new_tag(
        old,
        View(pipeline="NEW"),
    )

    assert result.pipeline == "NEW"


def test_rename_columns():
    result = ColumnTag.rename_columns(
        [
            "Count.GU (S1)",
            "VST",
        ],
        View(metric=("X",)),
    )

    assert result == {
        "Count.GU (S1)": "Count.GU.X (S1)",
        "VST": "VST.X",
    }


def test_matches_view_empty_view():
    tag = ColumnTag(
        metric=("Count",),
        sample_id="S1",
        pipeline="STAR",
    )

    assert tag.matches_view(View())


def test_matches_view_metric():
    tag = ColumnTag(metric=("Count",))

    assert tag.matches_view(View(metric=("Count",)))

    assert not tag.matches_view(View(metric=("VST",)))


def test_matches_view_pipeline():
    tag = ColumnTag(
        metric=("Count",),
        pipeline="STAR",
    )

    assert tag.matches_view(View(pipeline="STAR"))

    assert not tag.matches_view(View(pipeline="OTHER"))


def test_matches_view_sample_string():
    tag = ColumnTag(
        metric=("Count",),
        sample_id="S1",
    )

    assert tag.matches_view(View(sample_id="S1"))

    assert not tag.matches_view(View(sample_id="S2"))


def test_matches_view_sample_tuple():
    tag = ColumnTag(
        metric=("Count",),
        sample_id="S1",
    )

    assert tag.matches_view(View(sample_id=("S1", "S2")))


def test_select_from_view_none():
    cols = [
        "Count.GU (S1)",
        "VST",
    ]

    assert (
        ColumnTag.select_from_view(
            cols,
            None,
        )
        == cols
    )


def test_select_from_view():
    cols = [
        "Count.GU (S1)",
        "Count.GU (S2)",
        "VST",
    ]

    result = ColumnTag.select_from_view(
        cols,
        View(
            metric=("Count", "GU"),
            sample_id="S1",
        ),
    )

    assert result == ["Count.GU (S1)"]
