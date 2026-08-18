from __future__ import annotations

import json
import math
import re

# ============================================================================
# ALIGNER
# ============================================================================
# def create_aligner(
#     match_score: float = 2.0,
#     mismatch_score: float = -2.0,
#     gap_open: float = -5.0,
#     gap_extend: float = -1.0,
# ) -> PairwiseAligner:
#     """Create a local sequence aligner."""
#     aligner = PairwiseAligner()
#     aligner.mode = "local"
#     aligner.match_score = match_score
#     aligner.mismatch_score = mismatch_score
#     aligner.open_gap_score = gap_open
#     aligner.extend_gap_score = gap_extend
#     return aligner
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import pysam
from Bio.Align import PairwiseAligner  # type: ignore[import]


@dataclass
class Mutation:
    type: Literal["SNV", "DEL", "INS"]
    position: int
    end: int
    ref: str
    alt: str
    description: str


@dataclass
class MutationProfile:
    mutations: tuple[Mutation, ...]  # sorted!

    def __str__(self):
        return ",".join([mut.description for mut in self.mutations])

    def profile_to_string(self) -> str:
        return json.dumps(asdict(self))

    @classmethod
    def profile_from_string(cls, s: str) -> MutationProfile:
        data = json.loads(s)

        return MutationProfile(
            mutations=tuple(Mutation(**mutation) for mutation in data["mutations"])
        )

    def __contains__(self, mutation: Mutation) -> bool:
        return mutation in self.mutations

    def __iter__(self):
        return iter(self.mutations)

    def __len__(self) -> int:
        return len(self.mutations)


def extract_mutations(
    alignment,
    wt: str,
    read: str,
) -> MutationProfile:
    """
    Extract SNVs and true internal indels from a local alignment.

    Regions not included in the local alignment are ignored.

    Therefore a terminal truncation or a region skipped by local alignment
    is NOT classified as an indel.

    Return sorted tuple (by position).
    """

    mutations = []

    wt_blocks, read_blocks = alignment.aligned

    for i, (
        (wt_start, wt_end),
        (read_start, read_end),
    ) in enumerate(zip(wt_blocks, read_blocks)):

        # ---------------------------------------------------------
        # SNVs within aligned block
        # ---------------------------------------------------------

        wt_block = wt[wt_start:wt_end]
        read_block = read[read_start:read_end]

        if len(wt_block) != len(read_block):
            raise ValueError(
                "Aligned block lengths differ: " f"{len(wt_block)} vs {len(read_block)}"
            )

        for offset, (wt_base, read_base) in enumerate(zip(wt_block, read_block)):

            if wt_base != read_base:

                position = int(wt_start + offset + 1)

                mutations.append(
                    Mutation(
                        type="SNV",
                        position=position,
                        end=position,
                        ref=wt_base,
                        alt=read_base,
                        description=(f"{position}{wt_base}>{read_base}"),
                    )
                )

        # ---------------------------------------------------------
        # Gap to next alignment block
        # ---------------------------------------------------------

        if i + 1 >= len(wt_blocks):
            continue

        next_wt_start = wt_blocks[i + 1][0]
        next_read_start = read_blocks[i + 1][0]

        wt_gap = next_wt_start - wt_end
        read_gap = next_read_start - read_end

        # ---------------------------------------------------------
        # True deletion
        # ---------------------------------------------------------

        if wt_gap > 0 and read_gap == 0:

            deleted = wt[wt_end:next_wt_start]

            start = wt_end + 1
            end = next_wt_start

            mutations.append(
                Mutation(
                    type="DEL",
                    position=start,
                    end=end,
                    ref=deleted,
                    alt="-",
                    description=f"del_{start}-{end}",
                )
            )

        # ---------------------------------------------------------
        # True insertion
        # ---------------------------------------------------------

        elif wt_gap == 0 and read_gap > 0:

            inserted = read[read_end:next_read_start]

            # Position convention:
            #
            # insertion after WT position `wt_end`
            #
            position = int(wt_end)

            mutations.append(
                Mutation(
                    type="INS",
                    position=position,
                    end=position,
                    ref="-",
                    alt=inserted,
                    description=(f"ins_{position}_{inserted}"),
                )
            )

        # ---------------------------------------------------------
        # Both sequences advance
        # ---------------------------------------------------------
        #
        # This is a skipped region in a local alignment.
        # It is NOT an indel.
        #
        elif wt_gap > 0 and read_gap > 0:
            pass

        elif wt_gap == 0 and read_gap == 0:
            pass

    mutations = tuple(sorted(mutations, key=lambda x: x.position))
    return MutationProfile(mutations)


########################################################################################
# Create Aligner
########################################################################################


def create_aligner(
    match_score: float = 2,
    mismatch_score: float = -2,
    gap_open: float = -5,
    gap_extend: float = -1,
    mode: str = "local",
) -> PairwiseAligner:
    """Create a pairwise sequence aligner."""

    aligner = PairwiseAligner()

    aligner.mode = mode

    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend

    return aligner


########################################################################################
# ALign and extract a Mutation Profile
########################################################################################


def get_mutation_profile(
    sequence: str,
    wt: str,
    aligner: PairwiseAligner,
) -> MutationProfile:
    """
    Determine the complete mutation profile of `sequence` relative to WT.

    The profile contains all:
        - SNVs
        - deletions
        - insertions

    This is the same mutation extraction logic used for reads.

    Parameters
    ----------
    sequence:
        Sequence to compare against WT.

    wt:
        Full-length WT sequence.

    aligner:
        PairwiseAligner used for the comparison.

        For predefined full-length sequences this should normally be
        a global aligner.

        For reads with possible truncations this can be a local aligner.

    Returns
    -------
    tuple[str, ...]
        Sorted mutation descriptions.

    Examples
    --------
    ("42C>T",)

    ("42C>T", "87C>G")

    ("54G>T", "55A>G", "del_98-100")

    ("ins_123_A", "157A>G")
    """

    alignments = aligner.align(wt, sequence)

    if len(alignments) == 0:
        raise ValueError("Could not align sequence against WT.")

    alignment = alignments[0]

    profile = extract_mutations(
        alignment=alignment,
        wt=wt,
        read=sequence,
    )
    return profile


########################################################################################
# Analyze read with predefined sequences (main API)
########################################################################################


def analyze_read_with_predefined(
    read: str,
    wt: str,
    predefined_sequences: dict[str, str],
    aligner: PairwiseAligner | None = None,
    min_coverage: float = 0.80,
    predefined_profile_strings: dict[str, str] | None = None,
) -> dict[str, Any]:
    """
    Convenience wrapper.

    Predefined sequences are expected as:

        {
            "Mut_001": "ACTG...",
            "Mut_002": "ACTA...",
        }

    Exact matches are identified first.
    All non-exact reads are aligned against WT and subsequently
    compared against the predefined mutation profiles.
    """

    prealigner = aligner or create_aligner(mode="global")
    if predefined_profile_strings:
        predefined_mutation_profiles = {}
        for key, str_rep in predefined_profile_strings:
            predefined_mutation_profiles[key] = MutationProfile.profile_from_string(
                str_rep
            )
    else:
        # build new
        predefined_mutation_profiles = build_mutation_profiles(
            predefined_sequences=predefined_sequences,
            wt=wt,
        )

    return analyze_read(
        read=read,
        wt=wt,
        aligner=prealigner,
        min_coverage=min_coverage,
        predefined_sequences=predefined_sequences,
        predefined_mutation_profiles=predefined_mutation_profiles,
    )


def analyze_read(
    read: str,
    wt: str,
    aligner: PairwiseAligner | None = None,
    min_coverage: float = 0.80,
    predefined_sequences: dict[str, str] | None = None,
    predefined_mutation_profiles: dict[str, MutationProfile] | None = None,
) -> dict[str, Any]:
    """
    Analyze one read against WT.

    Workflow
    --------
    1. Exact predefined match
    2. Local WT alignment
    3. Coverage
    4. SNVs / indels
    5. Truncation
    6. Best predefined mutation profile

    Important
    ---------
    The best predefined sequence is searched for regardless of
    coverage/truncation.

    A truncated read can therefore still receive a best-predefined
    assignment if the observed mutations provide sufficient evidence.

    Mutation profiles contain all mutation types:

        - SNV
        - DEL
        - INS

    The same mutation representation is used for reads and predefined
    sequences.
    """

    readaligner = aligner or create_aligner(mode="local")

    # ====================================================================
    # EXACT PREDEFINED MATCH
    # ====================================================================

    exact_match = None

    if predefined_sequences is not None:
        exact_match = find_exact_predefined(
            read=read,
            predefined_sequences=predefined_sequences,
        )

    # ====================================================================
    # ALIGNMENT
    # ====================================================================

    alignments = readaligner.align(wt, read)

    if len(alignments) == 0:

        return {
            "read": read,
            "wt": wt,
            "coverage": 0.0,
            "covered_bases": 0,
            "covered_intervals": [],
            "n_snv": 0,
            "n_indel": 0,
            "mutations": [],
            "mutation_profile": (),
            "classification": (
                "exact_predefined" if exact_match is not None else "ambiguous"
            ),
            "wt_aligned": "",
            "read_aligned": "",
            "alignment_score": None,
            "wt_start": None,
            "wt_end": None,
            "read_start": None,
            "read_end": None,
            "truncated_left": False,
            "truncated_right": False,
            "matched_sequence": exact_match,
            "best_predefined": exact_match,
            "best_predefined_score": None,
            "predefined_ambiguous": False,
            "predefined_matches": [],
            "predefined_missing": [],
            "predefined_mismatches": [],
            "predefined_extra": [],
            "predefined_uncovered_expected": [],
            "predefined_candidates": [],
        }

    alignment = alignments[0]

    # ====================================================================
    # COORDINATES / COVERAGE
    # ====================================================================

    coordinates = get_alignment_coordinates(alignment)

    wt_start = coordinates["wt_start"]
    wt_end = coordinates["wt_end"]

    read_start = coordinates["read_start"]
    read_end = coordinates["read_end"]

    covered_intervals = coordinates["covered_intervals"]
    covered_bases = coordinates["covered_bases"]

    coverage = covered_bases / len(wt) if len(wt) > 0 else 0.0

    # --------------------------------------------------------------------
    # Truncation
    # --------------------------------------------------------------------

    truncated_left = wt_start > 0
    truncated_right = wt_end < len(wt)

    # ====================================================================
    # Extract Mutations for this read
    # ====================================================================

    # IMPORTANT:
    # This is deliberately generated from the exact same mutation
    # dictionaries used for predefined sequences.
    read_profile: MutationProfile = extract_mutations(
        alignment=alignment,
        wt=wt,
        read=read,
    )

    # mutation_profile = tuple(sorted(mutation["description"] for mutation in mutations))

    # ====================================================================
    # Create gapped alignment strings for visualization
    # ====================================================================

    wt_aligned, read_aligned = alignment_to_strings(
        alignment=alignment,
        wt=wt,
        read=read,
    )

    # ====================================================================
    # Find the best match among the predefined sequnces
    # ====================================================================
    #
    # IMPORTANT:
    #
    # This is deliberately performed BEFORE classification.
    #
    # Therefore a truncated read can still receive a best-predefined
    # assignment.
    #
    # Example:
    #
    # predefined:
    #
    #     Mut_A = 54G>T
    #     Mut_B = 87C>G
    #
    # read:
    #
    #     54G>T
    #     right side truncated
    #
    # => Mut_A can still be the best predefined candidate.
    #
    # find_best_predefined() must distinguish:
    #
    #     missing
    #         Expected mutation is inside covered sequence but was not
    #         observed.
    #
    #     uncovered_expected
    #         Expected mutation lies outside the covered region.
    #
    # The latter is NOT evidence against the predefined sequence.
    # ====================================================================

    predefined_result = None

    best_predefined = None
    best_predefined_score = None
    predefined_ambiguous = False

    if predefined_mutation_profiles:
        predefined_result = find_best_predefined(
            read_mutation_profile=read_profile,
            predefined_mutation_profiles=predefined_mutation_profiles,
            covered_intervals=covered_intervals,
        )
        best_predefined = predefined_result["best_match"]
        best_predefined_score = predefined_result["best_score"]
        best_predefined_profile = predefined_result["best_profile"]
        predefined_ambiguous = predefined_result["ambiguous"]

    # ====================================================================
    # Why is this possibly the predefined source?
    # ====================================================================

    predefined_matches: list[str] = []
    predefined_missing: list[str] = []
    predefined_mismatches: list[dict[str, Any]] = []
    predefined_extra: list[str] = []
    predefined_uncovered_expected: list[str] = []

    if predefined_result is not None and predefined_result["candidates"]:
        best_candidate = predefined_result["candidates"][0]
        predefined_matches = best_candidate["matches"]
        predefined_missing = best_candidate["missing"]
        predefined_mismatches = best_candidate["conflicts"]
        predefined_extra = best_candidate["extra"]
        predefined_uncovered_expected = best_candidate["uncovered_expected"]

    # ====================================================================
    # Classify the read
    # ====================================================================
    snvs = [x for x in read_profile.mutations if x.type == "SNV"]
    indels = [x for x in read_profile.mutations if x.type in ("DEL", "INS")]
    classification = __classify(
        exact_match=exact_match,
        coverage=coverage,
        min_coverage=min_coverage,
        snvs=snvs,
        indels=indels,
        truncated_left=truncated_left,
        truncated_right=truncated_right,
        best_predefined=best_predefined,
        predefined_result=predefined_result,
        predefined_ambiguous=predefined_ambiguous,
    )

    # ====================================================================
    # RESULT
    # ====================================================================
    candidates = (
        [
            (x["name"], x["profile"], x["score"])
            for x in predefined_result["candidates"][:3]
        ]
        if predefined_result is not None
        else []
    )

    return {
        "read": read,
        "wt": wt,
        "coverage": coverage,
        "covered_bases": covered_bases,
        "covered_intervals": covered_intervals,
        "wt_start": (wt_start + 1 if wt_start is not None else None),
        "wt_end": wt_end,
        "read_start": (read_start + 1 if read_start is not None else None),
        "read_end": read_end,
        "truncated_left": truncated_left,
        "truncated_right": truncated_right,
        "n_snv": len(snvs),
        "n_indel": len(indels),
        # "mutations": str(read_profile),
        "mutation_profile": read_profile,
        # ---------------------------------------------------------------
        # Alignment
        # ---------------------------------------------------------------
        "wt_aligned": wt_aligned,
        "read_aligned": read_aligned,
        "alignment_score": alignment.score,
        # ---------------------------------------------------------------
        # Classification
        # ---------------------------------------------------------------
        "classification": classification,
        # ---------------------------------------------------------------
        # Exact predefined
        # ---------------------------------------------------------------
        "matched_sequence": exact_match,
        # ---------------------------------------------------------------
        # Best predefined
        # ---------------------------------------------------------------
        "best_predefined": best_predefined,
        "best_predefined_score": best_predefined_score,
        "best_predefined_profile": best_predefined_profile,
        "predefined_ambiguous": predefined_ambiguous,
        # ---------------------------------------------------------------
        # Evidence
        # ---------------------------------------------------------------
        "predefined_matches": predefined_matches,
        "predefined_missing": predefined_missing,
        "predefined_mismatches": predefined_mismatches,
        "predefined_extra": predefined_extra,
        "predefined_uncovered_expected": (predefined_uncovered_expected),
        # ---------------------------------------------------------------
        # All candidates
        # ---------------------------------------------------------------
        "predefined_candidates": (candidates),
    }


# # ============================================================================
# # PREDEFINED SEQUENCES / MUTATION PROFILES
# # ============================================================================


# def build_mutation_lookup(
#     predefined_sequences: dict[str, str],
#     wt: str,
#     aligner: PairwiseAligner | None = None,
# ) -> dict[MutationProfile, list[str]]:
#     """
#     Build a reverse lookup:

#         mutation profile -> predefined sequence name

#     The mutation profile contains all mutations relative to WT,
#     including:

#         - SNVs
#         - deletions
#         - insertions

#     Example
#     -------
#     {
#         ("42C>T",): "Mut_001",
#         ("87C>G",): "Mut_002",
#         ("54G>T", "55A>G", "56G>A"): "Mut_003",
#         ("42C>T", "del_100-105"): "Mut_004",
#     }

#     Parameters
#     ----------
#     predefined_sequences:
#         Mapping from predefined sequence name to full predefined
#         sequence.

#     wt:
#         Full-length WT sequence.

#     aligner:
#         Aligner used to determine the mutation profile.

#         For predefined sequences this should normally be a global
#         aligner. If omitted, a global aligner is created.

#     Returns
#     -------
#     dict[tuple[str, ...], str]
#         Mapping from complete mutation profile to predefined name.

#     Raises
#     ------
#     ValueError
#         If two predefined sequences have the same mutation profile.
#     """

#     if aligner is None:
#         aligner = create_aligner(mode="global")

#     lookup: dict[tuple[str, ...], str] = {}

#     for name, sequence in predefined_sequences.items():

#         profile = get_mutation_profile(
#             sequence=sequence,
#             wt=wt,
#             aligner=aligner,
#         )

#         key = tuple(sorted(profile))

#         if key in lookup:
#             raise ValueError(
#                 "Duplicate mutation profile:\n"
#                 f"  profile: {key}\n"
#                 f"  sequence 1: {lookup[key]!r}\n"
#                 f"  sequence 2: {name!r}"
#             )

#         lookup[key] = name

#     return lookup


def build_mutation_profiles(
    predefined_sequences: dict[str, str],
    wt: str,
    aligner: PairwiseAligner | None = None,
) -> dict[str, MutationProfile]:
    """
    Build:

        predefined sequence name -> complete mutation profile

    The mutation profile contains all differences relative to WT,
    including:

        - SNVs
        - deletions
        - insertions

    Example
    -------
    {
        "Mut_001": ("42C>T",),
        "Mut_002": ("87C>G",),
        "Mut_003": (
            "54G>T",
            "55A>G",
            "56G>A",
        ),
        "Mut_004": (
            "42C>T",
            "del_100-105",
        ),
    }

    Parameters
    ----------
    predefined_sequences:
        Mapping from predefined sequence name to full sequence.

    wt:
        Full-length WT sequence.

    aligner:
        PairwiseAligner used for the comparison.

        For predefined sequences this should normally be a global
        aligner. If omitted, a global aligner is created.

    Returns
    -------
    dict[str, tuple[str, ...]]
        Mapping from predefined sequence name to its complete
        mutation profile.
    """

    if aligner is None:
        aligner = create_aligner(mode="global")

    profiles: dict[str, MutationProfile] = {}

    for name, sequence in predefined_sequences.items():
        profile = get_mutation_profile(
            sequence=sequence,
            wt=wt,
            aligner=aligner,
        )

        profiles[name] = profile

    return profiles


# ============================================================================
# ALIGNMENT COORDINATES
# ============================================================================


def get_alignment_coordinates(
    alignment,
) -> dict[str, Any]:
    """
    Extract coordinates and true WT coverage from a local alignment.

    Coordinates returned here are 0-based internally.

    covered_intervals contains half-open intervals:

        (start, end)

    corresponding to WT positions:

        start ... end-1
    """

    wt_blocks, read_blocks = alignment.aligned

    if len(wt_blocks) == 0:
        return {
            "wt_start": None,
            "wt_end": None,
            "read_start": None,
            "read_end": None,
            "covered_intervals": [],
            "covered_bases": 0,
        }

    covered_intervals = [(int(start), int(end)) for start, end in wt_blocks]

    covered_bases = sum(end - start for start, end in covered_intervals)

    return {
        "wt_start": int(wt_blocks[0][0]),
        "wt_end": int(wt_blocks[-1][1]),
        "read_start": int(read_blocks[0][0]),
        "read_end": int(read_blocks[-1][1]),
        "covered_intervals": covered_intervals,
        "covered_bases": covered_bases,
    }


def position_is_covered(
    position: int,
    covered_intervals: list[tuple[int, int]],
) -> bool:
    """
    Check whether a 1-based WT position is covered by the alignment.
    """

    position0 = position - 1

    return any(start <= position0 < end for start, end in covered_intervals)


def calculate_coverage(
    alignment,
    wt_length: int,
) -> tuple[int, float]:
    """
    Calculate true WT coverage from alignment blocks.

    This is intentionally based on alignment.aligned and NOT on
    wt_aligned/read_aligned strings.
    """

    wt_blocks, _ = alignment.aligned

    covered_bases = sum(int(end) - int(start) for start, end in wt_blocks)

    coverage = covered_bases / wt_length if wt_length > 0 else 0.0

    return covered_bases, coverage


# ============================================================================
# ALIGNMENT -> GAPPED STRINGS
# ============================================================================


def alignment_to_strings(
    alignment,
    wt: str,
    read: str,
) -> tuple[str, str]:
    """
    Convert a local Bio.Align alignment into gapped WT/read strings.

    Only gaps that are part of the actual alignment are represented.

    Regions skipped by the local alignment are NOT interpreted as
    insertions/deletions.
    """

    wt_blocks, read_blocks = alignment.aligned

    wt_aligned: list[str] = []
    read_aligned: list[str] = []

    for i, (
        (wt_start, wt_end),
        (read_start, read_end),
    ) in enumerate(zip(wt_blocks, read_blocks)):

        # ---------------------------------------------------------
        # Aligned block
        # ---------------------------------------------------------

        wt_block = wt[wt_start:wt_end]
        read_block = read[read_start:read_end]

        if len(wt_block) != len(read_block):
            raise ValueError(
                "Aligned block lengths differ: " f"{len(wt_block)} vs {len(read_block)}"
            )

        wt_aligned.append(wt_block)
        read_aligned.append(read_block)

        # ---------------------------------------------------------
        # Gap to next alignment block
        # ---------------------------------------------------------

        if i + 1 >= len(wt_blocks):
            continue

        next_wt_start = wt_blocks[i + 1][0]
        next_read_start = read_blocks[i + 1][0]

        wt_gap = next_wt_start - wt_end
        read_gap = next_read_start - read_end

        # True deletion:
        #
        # WT advances, read does not.
        #
        # WT:   ABCDEFG
        # READ: ABC---G
        #
        if wt_gap > 0 and read_gap == 0:

            wt_aligned.append(wt[wt_end:next_wt_start])

            read_aligned.append("-" * wt_gap)

        # True insertion:
        #
        # WT:   ABC---DE
        # READ: ABCXYZDE
        #
        elif wt_gap == 0 and read_gap > 0:

            wt_aligned.append("-" * read_gap)

            read_aligned.append(read[read_end:next_read_start])

        # Both sequences advance:
        #
        # This is a skipped region caused by the local alignment.
        #
        # Do NOT interpret it as an indel.
        #
        elif wt_gap > 0 and read_gap > 0:
            pass

        elif wt_gap == 0 and read_gap == 0:
            pass

        else:
            raise ValueError(
                f"Invalid alignment gap: " f"WT gap={wt_gap}, READ gap={read_gap}"
            )

    return (
        "".join(wt_aligned),
        "".join(read_aligned),
    )


# ============================================================================
# MUTATION EXTRACTION
# ============================================================================


# ============================================================================
# READ SNV PROFILE
# ============================================================================


def get_read_profile(
    mutations: list[dict[str, Any]],
) -> tuple[str, ...]:
    """
    Return the SNV profile observed in a read.
    """

    read_profile = tuple(sorted(mutation["description"] for mutation in mutations))
    return read_profile


########################################################################################
# find the most similar predefined sequence to this read
########################################################################################


def find_best_predefined(
    read_mutation_profile: MutationProfile,
    predefined_mutation_profiles: dict[str, MutationProfile],
    covered_intervals: list[tuple[int, int]],
    min_score: float = 0.0,
    min_expected_match_fraction: float = 0.5,
    ambiguity_delta: float = 0.0,
) -> dict[str, Any]:
    """
    Find the predefined sequence that best explains a read.

    `ambiguity_delta=0` means that only equal-scoring candidates are
    considered ambiguous.
    """
    candidates: list[dict[str, Any]] = []

    for name, profile in predefined_mutation_profiles.items():
        score_result = score_against_predefined(
            read_mutation_profile=read_mutation_profile,
            predefined_mutation_profile=profile,
            covered_intervals=covered_intervals,
        )
        candidates.append(
            {
                "name": name,
                "profile": profile,
                **score_result,
            }
        )

    if not candidates:

        return {
            "best_match": None,
            "best_score": None,
            "best_profile": None,
            "ambiguous": False,
            "candidates": [],
        }

    # ---------------------------------------------------------
    # Sort by explanatory power
    # ---------------------------------------------------------

    candidates.sort(
        key=lambda x: (
            x["score"],
            x["n_matches"],
            -x["n_conflicts"],
            -x["n_extra"],
        ),
        reverse=True,
    )

    # ---------------------------------------------------------
    # Filter implausible candidates
    # ---------------------------------------------------------
    plausible = []

    for candidate in candidates:

        if candidate["score"] < min_score:
            continue

        fraction = candidate["expected_match_fraction"]

        if fraction is not None and fraction < min_expected_match_fraction:
            continue

        plausible.append(candidate)
    if not plausible:

        return {
            "best_match": None,
            "best_score": None,
            "best_profile": None,
            "ambiguous": False,
            "candidates": candidates,
        }

    best = plausible[0]

    ambiguous = False

    if len(plausible) > 1:

        second = plausible[1]

        if best["score"] - second["score"] <= ambiguity_delta:
            ambiguous = True

    return {
        "best_match": best["name"],
        "best_score": best["score"],
        "best_profile": best["profile"],
        "ambiguous": ambiguous,
        "candidates": plausible,
    }


########################################################################################
# Comparison score
########################################################################################


def score_against_predefined(
    read_mutation_profile: MutationProfile,
    predefined_mutation_profile: MutationProfile,
    covered_intervals: list[tuple[int, int]],
) -> dict[str, Any]:
    """
    Compare observed read mutations against one predefined
    mutation profile.

    The comparison supports:

        - SNVs
        - deletions
        - insertions

    Important
    ---------
    A predefined mutation outside the read coverage is NOT
    considered missing.

    It is instead reported as:

        uncovered_expected

    This allows truncated reads to still receive a plausible
    predefined assignment.

    Example
    -------
    Predefined:

        54G>T
        87C>G

    Read:

        54G>T

    with coverage ending before position 87.

    Result:

        matches = ["54G>T"]
        uncovered_expected = ["87C>G"]

    and the second mutation receives no penalty.
    """
    matches: list[str] = []
    missing: list[str] = []
    mismatches: list[dict[str, Any]] = []
    uncovered_expected: list[str] = []
    conflicts: list[Mutation] = []
    extra: list[str] = []

    for expected_mutation in predefined_mutation_profile:
        # is it covered?
        if not mutation_is_covered(
            expected_mutation,
            covered_intervals,
        ):
            uncovered_expected.append(expected_mutation.description)
            continue
        # is the mutation also in the read's profile?
        if expected_mutation in read_mutation_profile:
            matches.append(expected_mutation.description)
        else:
            missing.append(expected_mutation.description)

    # remaining mutations in read
    remaining = [
        mutation
        for mutation in read_mutation_profile
        if mutation not in predefined_mutation_profile
    ]

    for observed_mutation in remaining:
        if mutation_conflicts(observed_mutation, predefined_mutation_profile):
            conflicts.append(observed_mutation)
            mismatches.append(
                {
                    "position": expected_mutation.position,
                    "expected": expected_mutation.description,
                    "observed": observed_mutation.description,
                    "expected_mutation": expected_mutation,
                    "observed_mutation": observed_mutation,
                }
            )

        elif not mutation_is_covered(
            observed_mutation,
            covered_intervals,
        ):
            continue
        else:
            extra.append(observed_mutation.description)

    score = similarity_score(
        matches=len(matches),
        missing=len(missing),
        extra=len(extra),
        # uncovered_expected=len(uncovered_expected),
        conflicts=len(mismatches),
    )
    n_observable_expected = len(matches) + len(missing) + len(mismatches)

    if n_observable_expected > 0:

        expected_match_fraction = len(matches) / n_observable_expected

    else:

        expected_match_fraction = None

    # ---------------------------------------------------------
    # Evidence strength
    # ---------------------------------------------------------

    n_expected = len(predefined_mutation_profile)

    n_uncovered = len(uncovered_expected)

    if n_expected > 0:

        observable_fraction = n_observable_expected / n_expected

    else:

        observable_fraction = None

    # ---------------------------------------------------------
    # Return
    # ---------------------------------------------------------

    return {
        "score": score,
        "matches": matches,
        "missing": missing,
        "conflicts": conflicts,
        "mismatches": mismatches,
        "uncovered_expected": uncovered_expected,
        "extra": extra,
        "n_matches": len(matches),
        "n_missing": len(missing),
        "n_conflicts": len(conflicts),
        "n_uncovered_expected": n_uncovered,
        "n_extra": len(extra),
        "n_expected": n_expected,
        "n_observable_expected": n_observable_expected,
        "expected_match_fraction": expected_match_fraction,
        "observable_fraction": observable_fraction,
        "n_read_mutations": len(read_mutation_profile),
        "n_predefined_mutations": len(predefined_mutation_profile),
    }


def similarity_score(
    matches: int,
    missing: int,
    extra: int,
    # uncovered_expected: int,
    conflicts: int,
    conflict_weight: float = 2.0,
) -> float:

    denominator = matches + missing + extra

    if denominator == 0:
        return 1.0 if conflicts == 0 else 0.0

    jaccard = matches / denominator

    conflict_penalty = math.exp(-conflict_weight * conflicts)

    return jaccard * conflict_penalty


def mutation_conflicts(
    observed_mutation: Mutation, predefined_mutation_profile: MutationProfile
) -> bool:
    for expected_mutation in predefined_mutation_profile:
        if conflicts(observed_mutation, expected_mutation):
            return True
    return False


def conflicts(a: Mutation, b: Mutation) -> bool:
    # Same mutation is obviously not a conflict
    if a == b:
        return False

    a_is_ins = a.type == "INS"
    b_is_ins = b.type == "INS"

    # Two insertions
    if a_is_ins and b_is_ins:
        return a.position == b.position

    # Insertion vs SNV/DEL
    if a_is_ins:
        return insertion_conflicts(a, b)

    if b_is_ins:
        return insertion_conflicts(b, a)

    # SNV vs SNV, SNV vs DEL, DEL vs DEL
    return intervals_overlap(a, b)


def intervals_overlap(a: Mutation, b: Mutation) -> bool:
    return a.position <= b.end and b.position <= a.end


def insertion_conflicts(ins: Mutation, other: Mutation) -> bool:
    # An insertion at a position inside the reference interval
    # conflicts with a mutation affecting that reference base.
    return other.position <= ins.position <= other.end


def mutation_overlaps(
    mutation_a: dict[str, Any],
    mutation_b: dict[str, Any],
) -> bool:
    """
    Determine whether two mutation events affect the same WT region.

    This is used to distinguish a true mismatch from an extra mutation.

    Examples
    --------
    42C>T vs 42C>G
        -> True

    del_98-103 vs del_98-103
        -> True

    del_98-103 vs del_100-105
        -> True

    42C>T vs 87C>G
        -> False

    42C>T vs del_42-50
        -> True
    """

    a_start = mutation_a["position"]
    a_end = mutation_a["end"]

    b_start = mutation_b["position"]
    b_end = mutation_b["end"]

    return not (a_end < b_start or b_end < a_start)


def mutation_is_covered(
    mutation: Mutation,
    covered_intervals: list[tuple[int, int]],
) -> bool:
    """
    Return True if the complete WT-coordinate range affected by
    the mutation is covered by the read.

    Parameters
    ----------
    mutation
        Mutation dictionary produced by extract_mutations().

    covered_intervals
        WT intervals covered by the read.

        IMPORTANT:
        These intervals are expected to use 1-based inclusive
        coordinates.

    Examples
    --------
    SNV at 42:

        covered interval (1, 100)
        -> True

    deletion 40-50:

        covered interval (1, 45)
        -> False

    deletion 40-50:

        covered interval (1, 100)
        -> True
    """

    for interval_start, interval_end in covered_intervals:

        if interval_start <= mutation.position and mutation.end <= interval_end:
            return True

    return False


# def mutation_position_range(
#     mutation: dict[str, Any],
# ) -> tuple[int, int]:
#     """
#     Return the WT-coordinate range affected by a mutation.

#     Coordinates are 1-based.

#     SNV:
#         42C>T -> (42, 42)

#     DEL:
#         del_100-105 -> (100, 105)

#     INS:
#         ins_100_ACG -> (100, 100)

#     Insertions are associated with the position immediately
#     preceding the insertion.
#     """

#     mutation_type = mutation["type"]

#     start = int(mutation["position"])

#     if mutation_type == "SNV":
#         return start, start

#     if mutation_type == "DEL":
#         end = int(mutation["end"])
#         return start, end

#     if mutation_type == "INS":
#         return start, start

#     raise ValueError(f"Unknown mutation type: {mutation_type!r}")


########################################################################################
# Helper
########################################################################################


def mutations_to_profile(
    mutations: list[dict[str, Any]],
) -> tuple[tuple, ...]:
    """
    Convert a list of mutation dictionaries into a normalized,
    hashable mutation profile.

    This is useful for comparing predefined sequences and reads
    using exactly the same representation.

    The profile is sorted to make its order independent.

    Example
    -------
    [
        {"type": "SNV", "position": 42, ...},
        {"type": "DEL", "position": 100, ...},
    ]

    becomes something like:

    (
        ("DEL", 100, 105, "ACGT"),
        ("SNV", 42, "C", "T"),
    )
    """

    return tuple(sorted(mutation_to_key(mutation) for mutation in mutations))


def mutation_to_key(mutation: dict[str, Any]) -> tuple:
    """
    Convert a mutation dictionary into a hashable comparison key.

    The key contains all information relevant for comparing mutations.

    Examples
    --------
    SNV:
        {
            "type": "SNV",
            "position": 42,
            "wt": "C",
            "read": "T",
        }

        -> ("SNV", 42, "C", "T")

    DEL:
        {
            "type": "DEL",
            "position": 100,
            "end": 105,
            "wt": "ACGTAC",
            "read": "-",
        }

        -> ("DEL", 100, 105, "ACGTAC")

    INS:
        {
            "type": "INS",
            "position": 100,
            "wt": "-",
            "read": "ACG",
        }

        -> ("INS", 100, "ACG")
    """

    mutation_type = mutation["type"]

    if mutation_type == "SNV":

        return (
            "SNV",
            int(mutation["position"]),
            mutation["wt"],
            mutation["read"],
        )

    elif mutation_type == "DEL":

        return (
            "DEL",
            int(mutation["position"]),
            int(mutation["end"]),
            mutation["wt"],
        )

    elif mutation_type == "INS":

        return (
            "INS",
            int(mutation["position"]),
            mutation["read"],
        )

    else:
        raise ValueError(f"Unknown mutation type: {mutation_type!r}")


########################################################################################
# not done yet
########################################################################################


def parse_mutation(description: str) -> dict[str, Any]:
    """
    Parse a mutation description.

    Supported formats
    ------------------
    SNV:
        42C>T

    DEL:
        del_98-103

    INS:
        ins_98_A
        ins_98_TGG

    Returns
    -------
    dict

    Examples
    --------
    "42C>T"

        {
            "type": "SNV",
            "position": 42,
            "end": 42,
            "description": "42C>T",
        }

    "del_98-103"

        {
            "type": "DEL",
            "position": 98,
            "end": 103,
            "description": "del_98-103",
        }

    "ins_98_A"

        {
            "type": "INS",
            "position": 98,
            "end": 98,
            "description": "ins_98_A",
        }
    """

    # ================================================================
    # SNV
    # ================================================================

    match = re.fullmatch(
        r"(\d+)([ACGTN])>([ACGTN])",
        description,
    )

    if match:

        position = int(match.group(1))

        return {
            "type": "SNV",
            "position": position,
            "end": position,
            "description": description,
        }

    # ================================================================
    # DEL
    # ================================================================

    match = re.fullmatch(
        r"del_(\d+)-(\d+)",
        description,
    )

    if match:

        start = int(match.group(1))
        end = int(match.group(2))

        if end < start:
            raise ValueError(f"Invalid deletion coordinates: {description!r}")

        return {
            "type": "DEL",
            "position": start,
            "end": end,
            "description": description,
        }

    # ================================================================
    # INS
    # ================================================================

    match = re.fullmatch(
        r"ins_(\d+)_(.+)",
        description,
    )

    if match:

        position = int(match.group(1))
        inserted = match.group(2)

        return {
            "type": "INS",
            "position": position,
            "end": position,
            "sequence": inserted,
            "description": description,
        }

    raise ValueError(f"Could not parse mutation description: {description!r}")


# ============================================================================
# EXACT PREDEFINED MATCH
# ============================================================================


def find_exact_predefined(
    read: str,
    predefined_sequences: dict[str, str],
) -> str | None:
    """
    Return the predefined sequence NAME if the read is an exact
    sequence match.

    Returns None otherwise.
    """

    for name, sequence in predefined_sequences.items():

        if read == sequence:
            return name

    return None


# # ============================================================================
# # BASIC CLASSIFICATION
# # ============================================================================


# def classify_alignment(
#     alignment,
#     wt: str,
#     read: str,
#     min_coverage: float = 0.80,
# ) -> dict[str, Any]:
#     """
#     Basic classification of a read against WT.

#     This function is kept as a separate helper for compatibility,
#     while `analyze_read()` provides the complete analysis.
#     """

#     coordinates = get_alignment_coordinates(alignment)

#     covered_bases = coordinates["covered_bases"]

#     coverage = covered_bases / len(wt) if len(wt) > 0 else 0.0

#     mutations = extract_mutations(
#         alignment=alignment,
#         wt=wt,
#         read=read,
#     )

#     snvs = [m for m in mutations if m["type"] == "SNV"]

#     indels = [m for m in mutations if m["type"] in {"DEL", "INS"}]

#     if coverage < 0.5:

#         classification = "low_coverage"

#     elif coverage < min_coverage:

#         classification = "truncated"

#     elif len(indels) > 0:

#         classification = "indel"

#     elif len(snvs) == 0:

#         classification = "WT"

#     elif len(snvs) == 1:

#         classification = "single_snv_novel"

#     else:

#         classification = "multi_snv"

#     return {
#         "coverage": coverage,
#         "covered_bases": covered_bases,
#         "n_snv": len(snvs),
#         "n_indel": len(indels),
#         "mutations": mutations,
#         "classification": classification,
#     }


# ============================================================================
# COMPLETE READ ANALYSIS
# ============================================================================


def __classify(
    exact_match: str | None,
    coverage: float,
    min_coverage: float,
    snvs: list[Mutation],
    indels: list[Mutation],
    truncated_left: bool,
    truncated_right: bool,
    best_predefined: str | None,
    predefined_result: dict[str, Any] | None,
    predefined_ambiguous: bool,
) -> str:

    # perfect match
    if exact_match is not None:
        return "exact_predefined"

    # If several predefined sequences are equally plausible, don't
    # arbitrarily assign one.
    #
    # We still distinguish this from a completely novel sequence.

    if predefined_ambiguous:
        return __resolve_ambiguous(coverage=coverage, min_coverage=min_coverage)

    # ====================================================================
    # Is there a closest predefined sequence?
    # ====================================================================
    #
    # IMPORTANT:
    #
    # This comes BEFORE the coverage checks.
    #
    # Therefore a truncated read can still be classified as belonging
    # to a predefined variant.
    # ====================================================================

    if best_predefined is not None and predefined_result is not None:

        return __resolve_best_predefined(
            predefined_result=predefined_result,
            coverage=coverage,
            min_coverage=min_coverage,
            truncated_left=truncated_left,
            truncated_right=truncated_right,
            n_indel=len(indels),
        )

    # ====================================================================
    # 4. If we have no best predefined match, that is some other artifact
    # ====================================================================

    if coverage < 0.50:
        return "low_coverage, no match"

    if coverage < min_coverage:
        return "truncated, no match"

    if len(indels) > 0:
        return "indel"

    if len(snvs) == 0:
        return "WT"

    if len(snvs) == 1:
        return "single_snv_novel"

    return "multi_snv"


def __resolve_ambiguous(coverage: float, min_coverage: float) -> str:
    """
    Resolve the classification of a read for which several predefined
    sequences are equally plausible.

    This is deliberately kept separate from the best_predefined
    resolution, because it is a different situation.
    """
    if coverage < 0.50:
        return "ambiguous_predefined_low_coverage"

    if coverage < min_coverage:
        return "ambiguous_predefined_truncated"

    return "ambiguous_predefined"


def __resolve_best_predefined(
    predefined_result: dict[str, Any],
    coverage: float,
    min_coverage: float,
    truncated_left: bool,
    truncated_right: bool,
    n_indel: int,
) -> str:
    """
    Resolve the classification of a read for which a best predefined
    sequence has been identified.

    Important distinction:

        predefined_truncated
            The observed mutations are compatible with the predefined
            sequence, but the read does not cover the complete WT.

        predefined_with_indel
            The read contains an indel, but its SNV profile still points
            to a predefined sequence.

        predefined_with_extra_snv
            The read contains the expected predefined mutations plus
            additional SNVs.

        predefined_conflict
            The observed mutations conflict with the predefined profile.
    """

    candidates = predefined_result.get("candidates", [])

    if not candidates:
        return "predefined_conflict, no candidates"

    candidate = candidates[0]

    n_conflicts = candidate.get("n_conflicts", 0)
    n_missing = candidate.get("n_missing", 0)
    n_extra = candidate.get("n_extra", 0)

    n_uncovered_expected = len(candidate.get("uncovered_expected", []))

    # ====================================================================
    # HARD CONFLICT
    # ====================================================================
    #
    # An expected mutation that is covered but absent is a real conflict.
    #
    # An expected mutation outside the covered region is NOT a conflict.
    # ====================================================================

    if n_conflicts > 0:
        return "predefined_conflict, mismatching mutations"
    if n_missing > 0:
        return "predefined_conflict, missing expected mutations"

    # ====================================================================
    # TRUNCATED + OTHERWISE COMPATIBLE
    # ====================================================================

    is_truncated = truncated_left or truncated_right or coverage < min_coverage

    if is_truncated:

        if n_indel > 0:
            return "predefined_truncated_with_indel"

        if n_extra > 0:
            return "predefined_truncated_with_extra_snv"

        # This is the important case:
        #
        # The read is incomplete, but all observed mutations are
        # compatible with the predefined variant.
        #
        # Some expected mutations may simply be outside the read.
        if n_uncovered_expected > 0:
            return "predefined_truncated"

        return "predefined_truncated"

    # ====================================================================
    # FULL-LENGTH READ
    # ====================================================================

    if n_indel > 0:
        return "predefined_with_indel"

    if n_extra > 0:
        return "predefined_with_extra_snv"

    if n_missing == 0 and n_conflicts == 0:
        return "predefined_variant"

    return "predefined_conflict"


# # ============================================================================
# # ANALYZE AGGREGATED READS
# # ============================================================================


# def analyze_reads(
#     reads: dict[str, int],
#     wt: str,
#     predefined_sequences: dict[str, str],
#     aligner: PairwiseAligner | None = None,
#     min_coverage: float = 0.80,
# ) -> list[dict[str, Any]]:
#     """
#     Analyze unique read sequences while retaining their counts.

#     Parameters
#     ----------
#     reads:
#         {
#             sequence_1: count_1,
#             sequence_2: count_2,
#             ...
#         }

#     Returns
#     -------
#     List of result dictionaries.
#     """

#     if aligner is None:
#         aligner = create_aligner()

#     mutation_profiles = build_mutation_profiles(
#         predefined_sequences=predefined_sequences,
#         wt=wt,
#     )

#     results: list[dict[str, Any]] = []

#     for read, count in reads.items():

#         result = analyze_read(
#             read=read,
#             wt=wt,
#             aligner=aligner,
#             min_coverage=min_coverage,
#             predefined_sequences=predefined_sequences,
#             mutation_profiles=mutation_profiles,
#         )

#         result["count"] = count

#         results.append(result)

#     return results


# # ============================================================================
# # COMPATIBILITY HELPERS
# # ============================================================================


# def classify_against_predefined(
#     result: dict[str, Any],
#     mutation_lookup: dict[tuple[str, ...], str],
# ) -> dict[str, Any]:
#     """
#     Compatibility helper for the old API.

#     This performs the old exact-profile lookup, but now only on the
#     observed SNVs.

#     For the new workflow, prefer `analyze_read()` with
#     `mutation_profiles`.
#     """

#     mutation_profile = tuple(
#         sorted(
#             mutation["description"]
#             for mutation in result.get("mutations", [])
#             if mutation["type"] == "SNV"
#         )
#     )

#     matched = mutation_lookup.get(mutation_profile)

#     if matched is not None:

#         result["matched_sequence"] = matched

#         if result["classification"] == "single_snv_novel":
#             result["classification"] = "single_snv_predefined"

#     return result


def get_read_coordinates(
    wt_aligned: str,
    read_aligned: str,
) -> tuple[int, int] | None:
    """
    Determine the first and last 1-based WT positions covered by
    a gapped alignment string.

    This is retained for compatibility.

    For coverage calculations, however, use the alignment object
    directly via `get_alignment_coordinates()`.
    """

    if len(wt_aligned) != len(read_aligned):
        raise ValueError("wt_aligned and read_aligned must have equal length.")

    wt_pos = 0
    covered_positions: list[int] = []

    for wt_base, read_base in zip(
        wt_aligned,
        read_aligned,
    ):

        if wt_base != "-":
            wt_pos += 1

        if wt_base != "-" and read_base != "-":
            covered_positions.append(wt_pos)

    if not covered_positions:
        return None

    return (
        min(covered_positions),
        max(covered_positions),
    )


########################################################################################
# build profiles for our predefined sequences in the profile element
########################################################################################


def build_sequence_profile_lookup(
    predefined_sequences: dict[str, str],
    wt: str,
    aligner: PairwiseAligner | None = None,
) -> dict[str, str]:
    """
    Build:

        sequence -> complete mutation profile

    The profile contains SNVs, insertions and deletions.
    """

    prealigner = aligner or create_aligner(mode="global")

    lookup: dict[str, str] = {}

    for name, sequence in predefined_sequences.items():

        profile = get_mutation_profile(
            sequence=sequence,
            wt=wt,
            aligner=prealigner,
        )

        if sequence in lookup:
            raise ValueError(
                "Duplicate predefined sequence:\n"
                f"  sequence: {sequence}\n"
                f"  name: {name!r}"
            )

        lookup[sequence] = profile.profile_to_string()

    return lookup


################################################################################
# quick alignment to bam
################################################################################


def dataframe_to_igv_bam(
    df: pd.DataFrame,
    output_dir: str | Path = "igv_alignment",
    sample_name: str = "sample",
    reference_name: str = "WT",
    min_mapq: int = 60,
):
    """
    Generate an IGV-compatible BAM package from a read-analysis DataFrame.

    Expected columns
    ----------------
    wt
        Full WT/reference sequence.

    wt_aligned
        Reference sequence with '-' representing gaps introduced by the
        pairwise alignment.

    read_aligned
        Read sequence aligned against wt_aligned.

    count
        Number of original reads represented by this unique read sequence.

    Optional
    --------
    read
        Original read sequence.

    classification
        Classification of the read.

    Returns
    -------
    dict
        Paths to:
            - fasta
            - fai
            - bam
            - bai
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================
    # Validate input
    # =========================================================

    required = {
        "wt",
        "wt_aligned",
        "read_aligned",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if len(df) == 0:
        raise ValueError("DataFrame is empty.")

    # =========================================================
    # Reference
    # =========================================================

    wt_sequences = df["wt"].dropna().unique()

    if len(wt_sequences) != 1:
        raise ValueError("DataFrame contains multiple different WT sequences.")

    wt = wt_sequences[0].upper()

    reference_fasta = output_dir / f"{sample_name}.fasta"

    with open(reference_fasta, "w") as f:
        f.write(f">{reference_name}\n")

        # Wrap FASTA sequence at 80 bp
        for i in range(0, len(wt), 80):
            f.write(wt[i : i + 80] + "\n")

    # =========================================================
    # Index FASTA
    # =========================================================

    pysam.faidx(str(reference_fasta))

    reference_fai = Path(str(reference_fasta) + ".fai")

    # =========================================================
    # BAM header
    # =========================================================

    header = {
        "HD": {
            "VN": "1.6",
            "SO": "coordinate",
        },
        "SQ": [
            {
                "SN": reference_name,
                "LN": len(wt),
            }
        ],
        "PG": [
            {
                "ID": "dataframe_to_igv_bam",
                "PN": "dataframe_to_igv_bam",
                "VN": "1.0",
            }
        ],
    }

    bam_unsorted = output_dir / f"{sample_name}.unsorted.bam"
    bam = output_dir / f"{sample_name}.bam"
    bai = Path(str(bam) + ".bai")

    # =========================================================
    # Helpers
    # =========================================================

    def add_cigar_op(cigar, op, length=1):
        """Append/merge a CIGAR operation."""

        if length == 0:
            return

        if cigar and cigar[-1][0] == op:
            cigar[-1] = (
                op,
                cigar[-1][1] + length,
            )
        else:
            cigar.append((op, length))

    def alignment_to_cigar(
        wt_aligned,
        read_aligned,
    ):
        """
        Convert pairwise alignment strings into:

            CIGAR
            reference start
            reference end
            read sequence
        """

        if len(wt_aligned) != len(read_aligned):
            raise ValueError("wt_aligned and read_aligned have different lengths.")

        cigar = []

        read_seq = []

        ref_pos = 0

        # -----------------------------------------------------
        # Walk through aligned strings
        # -----------------------------------------------------

        for wt_base, read_base in zip(
            wt_aligned.upper(),
            read_aligned.upper(),
        ):

            # -------------------------------------------------
            # Reference gap -> insertion in read
            # -------------------------------------------------

            if wt_base == "-":

                if read_base == "-":
                    continue

                add_cigar_op(cigar, "I")
                read_seq.append(read_base)

                continue

            # -------------------------------------------------
            # Read gap -> deletion
            # -------------------------------------------------

            if read_base == "-":

                add_cigar_op(cigar, "D")

                ref_pos += 1

                continue

            # -------------------------------------------------
            # Normal base / mismatch
            # -------------------------------------------------

            add_cigar_op(cigar, "M")

            read_seq.append(read_base)

            ref_pos += 1

        # -----------------------------------------------------
        # Determine reference start from leading deletion
        # -----------------------------------------------------

        leading_ref_bases = 0

        for wt_base, read_base in zip(
            wt_aligned.upper(),
            read_aligned.upper(),
        ):

            if wt_base == "-":
                continue

            if read_base == "-":
                leading_ref_bases += 1
            else:
                break

        # -----------------------------------------------------
        # But more generally:
        # reference coordinate = number of WT bases before
        # first aligned read base
        # -----------------------------------------------------

        ref_before_first_read = 0

        for wt_base, read_base in zip(
            wt_aligned.upper(),
            read_aligned.upper(),
        ):

            if read_base != "-":
                break

            if wt_base != "-":
                ref_before_first_read += 1

        reference_start = ref_before_first_read

        # -----------------------------------------------------
        # Reference-consuming length
        # -----------------------------------------------------

        reference_consuming_ops = {
            "M",
            "=",
            "X",
            "D",
            "N",
        }

        reference_length = sum(
            length for op, length in cigar if op in reference_consuming_ops
        )

        reference_end = reference_start + reference_length

        return (
            cigar,
            reference_start,
            reference_end,
            "".join(read_seq),
        )

    # =========================================================
    # Write BAM
    # =========================================================

    with pysam.AlignmentFile(
        bam_unsorted,
        "wb",
        header=header,
    ) as bam_out:

        for idx, (_, row) in enumerate(df.iterrows()):

            wt_aligned = str(row["wt_aligned"]).upper()

            read_aligned = str(row["read_aligned"]).upper()

            (
                cigar,
                reference_start,
                reference_end,
                read_seq,
            ) = alignment_to_cigar(
                wt_aligned,
                read_aligned,
            )

            if not read_seq:
                continue

            # -------------------------------------------------
            # Read count
            # -------------------------------------------------

            count = int(row.get("count", 1))

            classification = str(
                row.get(
                    "classification",
                    "read",
                )
            )

            # -------------------------------------------------
            # Read name
            # -------------------------------------------------

            read_name = f"read_{idx + 1}" f"|count={count}" f"|class={classification}"

            # -------------------------------------------------
            # BAM record
            # -------------------------------------------------

            segment = pysam.AlignedSegment()

            segment.query_name = read_name

            segment.query_sequence = read_seq

            segment.flag = 0

            segment.reference_id = 0

            segment.reference_start = reference_start

            segment.mapping_quality = min_mapq

            segment.cigar = cigar

            segment.query_qualities = pysam.qualitystring_to_array("I" * len(read_seq))

            # -------------------------------------------------
            # Useful custom BAM tags
            # -------------------------------------------------

            segment.set_tag(
                "RC",
                count,
                value_type="i",
            )

            segment.set_tag(
                "CL",
                classification,
                value_type="Z",
            )

            # Coverage if available
            if "coverage" in row and pd.notna(row["coverage"]):
                segment.set_tag(
                    "CV",
                    float(row["coverage"]),
                    value_type="f",
                )

            # Mutation profile if available
            if "mutation_profile" in row:
                profile = row["mutation_profile"]

                if profile is not None:
                    profile_string = str(profile)

                    # BAM string tags cannot be unlimited;
                    # truncate to a sensible length.
                    profile_string = profile_string[:250]

                    segment.set_tag(
                        "MP",
                        profile_string,
                        value_type="Z",
                    )

            bam_out.write(segment)

    # =========================================================
    # Sort BAM
    # =========================================================

    pysam.sort(
        "-o",
        str(bam),
        str(bam_unsorted),
    )

    # Remove unsorted BAM
    bam_unsorted.unlink()

    # =========================================================
    # BAM index
    # =========================================================

    pysam.index(str(bam))

    # =========================================================
    # Return paths
    # =========================================================

    return {
        "fasta": reference_fasta,
        "fai": reference_fai,
        "bam": bam,
        "bai": bai,
    }


# def build_mutation_lookup_reverse(
#     predefined_sequences: dict[str, str],
#     wt: str,
# ) -> dict[str, tuple[str, ...]]:
#     """Map mutation profiles to predefined sequence names."""

#     lookup: dict[str, tuple[str, ...]] = {}

#     for name, sequence in predefined_sequences.items():

#         profile = get_snv_profile(sequence, wt)

#         if profile is None:
#             continue

#         profile = tuple(sorted(profile))
#         lookup[name] = profile
#     return lookup


# def plot_read_alignments(
#     analysis_df,
#     wt_length=None,
#     n_reads=50,
#     sort_by="count",
#     figsize=(14, 10),
# ):
#     """
#     Plot WT and read alignments relative to WT coordinates.

#     Parameters
#     ----------
#     analysis_df : pandas.DataFrame
#         Output of the read analysis.

#     wt_length : int, optional
#         Length of WT. If None, inferred from first row.

#     n_reads : int
#         Maximum number of reads to display.

#     sort_by : str
#         Column used for sorting, e.g. "count".

#     figsize : tuple
#         Matplotlib figure size.
#     """

#     df = analysis_df.copy()

#     if sort_by is not None:
#         df = df.sort_values(sort_by, ascending=False)

#     df = df.head(n_reads)

#     if wt_length is None:
#         wt_length = len(df.iloc[0]["wt"])

#     fig, ax = plt.subplots(figsize=figsize)

#     # ---------------------------------------------------------
#     # WT
#     # ---------------------------------------------------------

#     y_wt = len(df)

#     ax.plot(
#         [1, wt_length],
#         [y_wt, y_wt],
#         linewidth=5,
#         solid_capstyle="butt",
#     )

#     ax.text(
#         -5,
#         y_wt,
#         "WT",
#         ha="right",
#         va="center",
#         fontweight="bold",
#     )

#     # ---------------------------------------------------------
#     # Reads
#     # ---------------------------------------------------------

#     for i, (_, row) in enumerate(df.iterrows()):

#         y = len(df) - i - 1

#         coords = get_read_coordinates(
#             row["wt_aligned"],
#             row["read_aligned"],
#         )

#         if coords is None:
#             continue

#         start, end = coords

#         # Read coverage
#         ax.plot(
#             [start, end],
#             [y, y],
#             linewidth=6,
#             solid_capstyle="butt",
#         )

#         # Label
#         label = f"{row['classification']}  n={row['count']}"

#         ax.text(
#             -5,
#             y,
#             label,
#             ha="right",
#             va="center",
#             fontsize=8,
#         )

#         # -----------------------------------------------------
#         # Mutations
#         # -----------------------------------------------------

#         mutations = row["mutations"]

#         if isinstance(mutations, list):

#             for mutation in mutations:

#                 if mutation["type"] == "SNV":

#                     pos = mutation["position"]

#                     ax.scatter(
#                         pos,
#                         y,
#                         s=40,
#                         zorder=5,
#                     )

#                 elif mutation["type"] == "DEL":

#                     start = mutation["position"]
#                     end = mutation["end"]

#                     ax.plot(
#                         [start, end],
#                         [y, y],
#                         linewidth=8,
#                         alpha=0.5,
#                     )

#                 elif mutation["type"] == "INS":

#                     pos = mutation["position"]

#                     ax.scatter(
#                         pos,
#                         y,
#                         marker="^",
#                         s=50,
#                         zorder=5,
#                     )

#     # ---------------------------------------------------------
#     # Formatting
#     # ---------------------------------------------------------

#     ax.set_xlim(1, wt_length)

#     ax.set_ylim(-1, len(df) + 1)

#     ax.set_xlabel("Position in WT")

#     ax.set_yticks([])

#     ax.spines["left"].set_visible(False)
#     ax.spines["right"].set_visible(False)
#     ax.spines["top"].set_visible(False)

#     ax.grid(
#         axis="x",
#         alpha=0.2,
#     )

#     plt.tight_layout()

#     return fig, ax


########################################################################################


# def calculate_coverage(wt_aligned, read_aligned, wt_length):
#     """
#     Calculate coverage of the WT sequence by the read.
#     """

#     covered_wt_bases = sum(
#         1
#         for wt_base, read_base in zip(wt_aligned, read_aligned)
#         if wt_base != "-" and read_base != "-"
#     )

#     wt_coverage = covered_wt_bases / wt_length

#     return covered_wt_bases, wt_coverage


# def find_flank(sequence: str, flank: str, start: bool = True) -> int | None:
#     """
#     Find a flank sequence in a read.

#     Parameters
#     ----------
#     sequence : str
#         Sequencing read.
#     flank : str
#         Expected flank sequence.
#     start : bool
#         If True, search from the beginning.
#         If False, search from the end.

#     Returns
#     -------
#     int or None
#         Position of the flank, or None if not found.
#     """

#     sequence = sequence.upper()
#     flank = flank.upper()

#     if start:
#         pos = sequence.find(flank)
#     else:
#         pos = sequence.rfind(flank)

#     return pos if pos >= 0 else None


# def align_read_to_wt(
#     read: str, wt_sequence: str, flank_5: str | None = None, flank_3: str | None = None
# ) -> dict[str, int | str | list]:
#     """
#     Align one sequencing read against WT and annotate flanking regions.

#     Parameters
#     ----------
#     read : str
#         Sequencing read.

#     wt_sequence : str
#         WT/reference sequence.

#     flank_5 : str, optional
#         Expected 5' flanking sequence.

#     flank_3 : str, optional
#         Expected 3' flanking sequence.

#     Returns
#     -------
#     dict
#         Alignment and flank annotation.
#     """

#     read = read.upper()
#     wt_sequence = wt_sequence.upper()

#     # ---------------------------------------------------------
#     # Find flanks
#     # ---------------------------------------------------------

#     flank_5_start = None
#     flank_5_end = None

#     flank_3_start = None
#     flank_3_end = None

#     if flank_5 is not None:

#         flank_5 = flank_5.upper()

#         pos = read.find(flank_5)

#         if pos >= 0:
#             flank_5_start = pos
#             flank_5_end = pos + len(flank_5)

#     if flank_3 is not None:

#         flank_3 = flank_3.upper()

#         pos = read.rfind(flank_3)

#         if pos >= 0:
#             flank_3_start = pos
#             flank_3_end = pos + len(flank_3)

#     # ---------------------------------------------------------
#     # Align read to WT
#     # ---------------------------------------------------------

#     aligner = PairwiseAligner()

#     # Local alignment is useful because reads may be truncated
#     aligner.mode = "local"

#     aligner.match_score = 2
#     aligner.mismatch_score = -1
#     aligner.open_gap_score = -5
#     aligner.extend_gap_score = -1

#     alignment = aligner.align(wt_sequence, read)[0]

#     target_blocks = alignment.aligned[0]
#     query_blocks = alignment.aligned[1]

#     # ---------------------------------------------------------
#     # Create WT-relative representation
#     #
#     # "." = not covered by read
#     # "-" = deletion relative to WT
#     # base = observed base
#     # ---------------------------------------------------------

#     aligned_read = ["."] * len(wt_sequence)

#     insertions = []

#     for (wt_start, wt_end), (read_start, read_end) in zip(target_blocks, query_blocks):

#         for wt_pos, read_pos in zip(
#             range(wt_start, wt_end), range(read_start, read_end)
#         ):

#             aligned_read[wt_pos] = read[read_pos]

#     # ---------------------------------------------------------
#     # Detect gaps
#     # ---------------------------------------------------------

#     for i in range(len(target_blocks) - 1):

#         wt_end_previous = target_blocks[i][1]
#         wt_start_next = target_blocks[i + 1][0]

#         read_end_previous = query_blocks[i][1]
#         read_start_next = query_blocks[i + 1][0]

#         wt_gap = wt_start_next - wt_end_previous
#         read_gap = read_start_next - read_end_previous

#         # deletion relative to WT
#         if wt_gap > 0 and read_gap == 0:

#             for pos in range(wt_end_previous, wt_start_next):
#                 aligned_read[pos] = "-"

#         # insertion relative to WT
#         elif read_gap > 0 and wt_gap == 0:

#             inserted_sequence = read[read_end_previous:read_start_next]

#             insertions.append(
#                 {"wt_position": wt_end_previous, "sequence": inserted_sequence}
#             )

#     # ---------------------------------------------------------
#     # WT coordinates
#     # ---------------------------------------------------------

#     if len(target_blocks) > 0:

#         wt_start = target_blocks[0][0]
#         wt_end = target_blocks[-1][1]

#     else:

#         wt_start = None
#         wt_end = None

#     # ---------------------------------------------------------
#     # Flank status
#     # ---------------------------------------------------------

#     has_flank_5 = flank_5_start is not None
#     has_flank_3 = flank_3_start is not None

#     # ---------------------------------------------------------
#     # Build flank annotation string
#     # ---------------------------------------------------------

#     flank_annotation = []

#     if flank_5 is not None:

#         if has_flank_5:
#             flank_annotation.append("5prime_flank")
#         else:
#             flank_annotation.append("5prime_flank_missing")

#     if flank_3 is not None:

#         if has_flank_3:
#             flank_annotation.append("3prime_flank")
#         else:
#             flank_annotation.append("3prime_flank_missing")
#     return {
#         "read": read,
#         "aligned_read": "".join(aligned_read),
#         "wt_start": wt_start,
#         "wt_end": wt_end,
#         "aligned_length": (wt_end - wt_start if wt_start is not None else 0),
#         "alignment_score": alignment.score,
#         "insertions": insertions,
#         # flank information
#         "flank_5": flank_5,
#         "flank_5_start": flank_5_start,
#         "flank_5_end": flank_5_end,
#         "has_flank_5": has_flank_5,
#         "flank_3": flank_3,
#         "flank_3_start": flank_3_start,
#         "flank_3_end": flank_3_end,
#         "has_flank_3": has_flank_3,
#         "flank_annotation": ";".join(flank_annotation),
#     }


# def align_read_to_wt_no_flanks(read, wt_sequence):
#     """
#     Align one sequencing read against a WT sequence.

#     Parameters
#     ----------
#     read : str
#         Sequencing read.
#     wt_sequence : str
#         Reference / WT sequence.

#     Returns
#     -------
#     dict
#         Alignment information.
#     """

#     read = read.upper()
#     wt_sequence = wt_sequence.upper()

#     aligner = PairwiseAligner()
#     aligner.mode = "local"

#     # Scoring
#     aligner.match_score = 2
#     aligner.mismatch_score = -1
#     aligner.open_gap_score = -5
#     aligner.extend_gap_score = -1

#     alignment = aligner.align(wt_sequence, read)[0]

#     # Coordinates of aligned regions
#     target_blocks = alignment.aligned[0]
#     query_blocks = alignment.aligned[1]

#     # Alignment representation relative to WT
#     aligned_read = ["."] * len(wt_sequence)

#     insertions = []

#     for (wt_start, wt_end), (read_start, read_end) in zip(target_blocks, query_blocks):
#         for wt_pos, read_pos in zip(
#             range(wt_start, wt_end), range(read_start, read_end)
#         ):
#             aligned_read[wt_pos] = read[read_pos]

#     # Detect deletions between aligned blocks
#     for i in range(len(target_blocks) - 1):

#         wt_end_previous = target_blocks[i][1]
#         wt_start_next = target_blocks[i + 1][0]

#         read_end_previous = query_blocks[i][1]
#         read_start_next = query_blocks[i + 1][0]

#         wt_gap = wt_start_next - wt_end_previous
#         read_gap = read_start_next - read_end_previous

#         # Deletion relative to WT
#         if wt_gap > 0 and read_gap == 0:
#             for pos in range(wt_end_previous, wt_start_next):
#                 aligned_read[pos] = "-"

#         # Insertion relative to WT
#         elif read_gap > 0 and wt_gap == 0:
#             inserted_sequence = read[read_end_previous:read_start_next]

#             insertions.append(
#                 {"wt_position": wt_end_previous, "sequence": inserted_sequence}
#             )

#     # WT coordinates covered by the alignment
#     if len(target_blocks) > 0:
#         wt_start = target_blocks[0][0]
#         wt_end = target_blocks[-1][1]
#     else:
#         wt_start = None
#         wt_end = None

#     return {
#         "read": read,
#         "aligned_read": "".join(aligned_read),
#         "wt_start": wt_start,
#         "wt_end": wt_end,
#         "aligned_length": wt_end - wt_start if wt_start is not None else 0,
#         "insertions": insertions,
#         "alignment_score": alignment.score,
#     }


# # ============================================================================
# # PREDEFINED MATCHING
# # ============================================================================


# def parse_snv(
#     snv: str,
# ) -> tuple[int, str, str]:
#     """
#     Parse:

#         54G>C

#     into:

#         (54, "G", "C")
#     """

#     match = re.fullmatch(
#         r"(\d+)([ACGT])>([ACGT])",
#         snv,
#     )

#     if match is None:
#         raise ValueError(f"Invalid SNV string: {snv!r}")

#     return (
#         int(match.group(1)),
#         match.group(2),
#         match.group(3),
#     )
