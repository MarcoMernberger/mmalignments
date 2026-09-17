import pandas as pd
import pytest

from mmalignments.models.vep.vep import annotate_samples


def test_annotate_samples_keeps_wt_total_from_raw_data_even_when_below_min_count():
    dominant_variants = pd.DataFrame(
        [
            {
                "Sample": "S1",
                "Gene": "G1",
                "variant": "A",
                "count": 10,
                "frequency": 0.6667,
                "vep_consequence": "missense",
            }
        ]
    )

    dominant_haplotypes = pd.DataFrame(
        [
            {
                "Sample": "S1",
                "Gene": "G1",
                "haplotype": "A",
                "count": 10,
                "frequency": 0.6667,
                "functional_class": "missense",
                "haplotype_category": "MAJOR",
            }
        ]
    )

    raw_df = pd.DataFrame(
        [
            {"Sample": "S1", "Gene": "G1", "Count": 5, "variants": ""},
            {"Sample": "S1", "Gene": "G1", "Count": 10, "variants": "A"},
        ]
    )

    result = annotate_samples(
        dominant_variants,
        dominant_haplotypes,
        raw_df=raw_df,
        min_frequency=0.01,
        major_haplo_threshold=0.05,
    )

    assert result.loc[0, "wt_count"] == 5
    assert result.loc[0, "wt_frequency"] == pytest.approx(5 / 15)
    assert result.loc[0, "wt_interpretation"] == "33.3%"


def test_annotate_samples_counts_major_and_minor_haplotypes_case_insensitive():
    dominant_variants = pd.DataFrame(
        [
            {
                "Sample": "S1",
                "Gene": "G1",
                "variant": "A",
                "count": 12,
                "frequency": 0.60,
                "vep_consequence": "missense",
            },
            {
                "Sample": "S1",
                "Gene": "G1",
                "variant": "B",
                "count": 6,
                "frequency": 0.30,
                "vep_consequence": "synonymous",
            },
        ]
    )

    dominant_haplotypes = pd.DataFrame(
        [
            {
                "Sample": "S1",
                "Gene": "G1",
                "haplotype": "A",
                "count": 12,
                "frequency": 0.60,
                "functional_class": "missense",
                "haplotype_category": "major",
            },
            {
                "Sample": "S1",
                "Gene": "G1",
                "haplotype": "B",
                "count": 6,
                "frequency": 0.30,
                "functional_class": "synonymous",
                "haplotype_category": "minor",
            },
            {
                "Sample": "S1",
                "Gene": "G1",
                "haplotype": "WT",
                "count": 2,
                "frequency": 0.10,
                "functional_class": "WT",
                "haplotype_category": "WT",
            },
        ]
    )

    result = annotate_samples(
        dominant_variants,
        dominant_haplotypes,
        min_frequency=0.01,
        major_haplo_threshold=0.05,
    )

    assert result.loc[0, "n_major_haplotypes"] == 1
    assert result.loc[0, "n_minor_haplotypes"] == 1
    assert result.loc[0, "major_mutant_count"] == 12
    assert result.loc[0, "major_mutant_frequency"] == pytest.approx(0.60)
