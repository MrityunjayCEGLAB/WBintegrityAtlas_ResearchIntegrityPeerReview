#!/usr/bin/env python3
"""
WB-IntegrityAtlas: molecular-weight auditability analysis
=========================================================

Reads the output folder produced by WB-IntegrityAtlas and quantifies how
informative a protein's theoretical molecular weight is for protein identity.

Expected primary input:
    proteome_wb_integrity_atlas.csv

Optional input:
    mw_bin_summary.csv

Main outputs:
    TSV/
        protein_mw_auditability.tsv
        mw_bin_auditability.tsv
        mw_grid_information.tsv
        mw_auditability_summary.tsv
        top_ambiguous_proteins.tsv
        top_distinctive_proteins.tsv

    Figures/
        Figure_MW_distribution.*
        Figure_MW_density_and_effective_candidates.*
        Figure_protein_auditability_scatter.*
        Figure_ambiguity_distribution.*
        Figure_bin_rankings.*
        Figure_auditability_landscape.*
        Figure_nearest_neighbor_gap.*
        Figure_summary_multipanel.*

Usage:
    python wb_mw_auditability.py /path/to/human_wb_atlas

Optional:
    python wb_mw_auditability.py /path/to/human_wb_atlas \
        --bin-width 2 \
        --tolerance-kda 2 \
        --sigma-kda 1 \
        --kernel-cutoff-sigma 4 \
        --grid-step 0.25 \
        --query-mw 42.0,45.7

When --query-mw is supplied, all standard outputs are retained and additional
query-specific outputs are written under:

    mw_auditability_analysis/
        Query_MW/
            query_mw_summary_all.tsv
            42.0_kDa/
                query_mw_summary.tsv
                query_mw_candidates.tsv
                query_mw_local_landscape.tsv
                Figure_query_gaussian_compatibility.*
                Figure_query_top_candidates.*
                Figure_query_local_auditability.*
            45.7_kDa/
                ...

Notes
-----
1. "Ambiguity" is NOT a misconduct score. It quantifies the number/density of
   alternative proteome entries that are compatible with an apparent MW.
2. The current implementation is proteome-informed, not tissue-specific.
3. If all isoforms are present in the atlas, the results describe isoform-level
   ambiguity. If one canonical sequence per gene is present, they describe
   gene/canonical-protein-level ambiguity.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Iterable, Optional
from collections import Counter


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def first_existing_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    normalized = {norm_name(c): c for c in df.columns}
    for candidate in candidates:
        key = norm_name(candidate)
        if key in normalized:
            return normalized[key]
    return None


def find_mw_column(df: pd.DataFrame) -> str:
    candidates = [
        "molecular_weight_kda",
        "molecular_weight_kDa",
        "theoretical_mw_kda",
        "theoretical_molecular_weight_kda",
        "mw_kda",
        "molecular_weight",
        "protein_mw_kda",
        "mass_kda",
    ]
    col = first_existing_column(df, candidates)
    if col is not None:
        return col

    # Conservative heuristic fallback.
    for c in df.columns:
        n = norm_name(c)
        if ("molecularweight" in n or n.startswith("mw") or "masskda" in n) and "bin" not in n:
            vals = pd.to_numeric(df[c], errors="coerce")
            if vals.notna().sum() >= max(10, int(0.25 * len(df))):
                return c
    raise ValueError(
        "Could not identify a molecular-weight column. "
        "Please rename it to e.g. 'molecular_weight_kda' or use --mw-column."
    )


def numeric_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")
    # Extract first numeric token from strings such as "42.1 kDa"
    return pd.to_numeric(
        s.astype(str).str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", expand=False),
        errors="coerce",
    )


def save_figure(fig: plt.Figure, outbase: Path, dpi: int = 400) -> None:
    fig.savefig(outbase.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def empirical_percentile(values: np.ndarray) -> np.ndarray:
    """Percentile rank in [0,1], average-rank behavior for ties."""
    s = pd.Series(values)
    return s.rank(method="average", pct=True).to_numpy(dtype=float)


def robust_log10(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return np.log10(np.maximum(np.asarray(x, dtype=float), eps))


def parse_query_mw(value: Optional[str]) -> list[float]:
    """
    Parse --query-mw values.

    Accepted examples:
        --query-mw 42.0
        --query-mw 42.0,45.7,100
        --query-mw "42.0, 45.7, 100"

    Duplicate values are removed while preserving input order.
    """
    if value is None:
        return []

    raw = str(value).strip()
    if not raw:
        return []

    out = []
    seen = set()
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            mw = float(token)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"Invalid --query-mw value '{token}'. Use comma-separated numeric kDa values, "
                "for example --query-mw 42.0,45.7"
            ) from exc
        if not np.isfinite(mw) or mw <= 0:
            raise argparse.ArgumentTypeError(
                f"Query MW must be a finite positive value; got {token!r}."
            )
        if mw not in seen:
            out.append(mw)
            seen.add(mw)
    return out


def query_label(mw: float) -> str:
    """Filesystem-safe query label preserving useful decimal precision."""
    s = f"{mw:.6f}".rstrip("0").rstrip(".")
    if "." not in s:
        s += ".0"
    return f"{s}_kDa"


def percentile_of_value(value: float, reference: np.ndarray) -> float:
    """
    Empirical percentile of a scalar relative to a reference distribution.
    Uses midpoint treatment for exact ties:
        P(X < value) + 0.5 P(X = value)
    """
    ref = np.asarray(reference, dtype=float)
    ref = ref[np.isfinite(ref)]
    if ref.size == 0 or not np.isfinite(value):
        return np.nan
    less = np.count_nonzero(ref < value)
    equal = np.count_nonzero(np.isclose(ref, value, rtol=1e-12, atol=1e-12))
    return float((less + 0.5 * equal) / ref.size)



# ---------------------------------------------------------------------------
# Audit-relevant annotation context inherited from the full-report output
# ---------------------------------------------------------------------------

CASE_PROTEINS = [
    {"case": "p38delta / beta-actin neighborhood", "gene": "ACTB", "display": "ACTB / beta-actin", "role": "marker/source neighborhood"},
    {"case": "p38delta / beta-actin neighborhood", "gene": "MAPK13", "display": "MAPK13 / p38delta", "role": "claimed/target neighborhood"},
    {"case": "p38delta / beta-actin neighborhood", "gene": "MAPK14", "display": "MAPK14 / p38alpha", "role": "family context"},
    {"case": "MEK3 / p44 MAPK neighborhood", "gene": "MAP2K3", "display": "MAP2K3 / MEK3", "role": "claimed/target neighborhood"},
    {"case": "MEK3 / p44 MAPK neighborhood", "gene": "MAPK3", "display": "MAPK3 / ERK1-p44", "role": "source/family neighborhood"},
    {"case": "AP-1 / cell-cycle neighborhood", "gene": "JUN", "display": "JUN / AP-1", "role": "construct context"},
    {"case": "AP-1 / cell-cycle neighborhood", "gene": "FOSL1", "display": "FOSL1 / Fra-1", "role": "source/context neighborhood"},
    {"case": "AP-1 / cell-cycle neighborhood", "gene": "CCNA2", "display": "CCNA2 / Cyclin A2", "role": "source/context neighborhood"},
    {"case": "HSP25 / HSP27 neighborhood", "gene": "HSPB1", "display": "HSPB1 / HSP27", "role": "HSP relabeling context"},
    {"case": "Loading-control context", "gene": "GAPDH", "display": "GAPDH", "role": "common marker"},
    # -----------------------------------------------------------------------
    # Drosophila melanogaster documented protein-identity cases
    # -----------------------------------------------------------------------

    {"case": "Drosophila dPKB / Akt western-blot reuse", "gene": "AKT", "display": "dPKB / Akt", "role": "reported dPKB western-blot case"},

    {"case": "Drosophila dS6K western-blot reuse", "gene": "S6K", "display": "dS6K / S6K", "role": "reported dS6K western-blot case"},

    {"case": "Drosophila Mad phospho-western identity swap", "gene": "MAD", "display": "Mad / pMad", "role": "reported pMad blot-identity swap"},

    {"case": "Drosophila Flag-Mad-AVA blot swap", "gene": "MAD", "display": "Flag-Mad-AVA", "role": "reported Flag-Mad-AVA blot swap"},

    {"case": "Drosophila H3K9me3 western-blot orientation error", "gene": "H3", "display": "H3 / H3K9me3", "role": "reported H3K9me3 blot presentation error"},


    # -----------------------------------------------------------------------
    # Saccharomyces cerevisiae documented protein-identity cases
    # -----------------------------------------------------------------------

    {"case": "Yeast Leu1p western-blot lane assembly", "gene": "LEU1", "display": "Leu1p", "role": "reported cut-and-reassembled western-blot case"},

    {"case": "Yeast Rtn1-Cherry western-blot lane duplication", "gene": "RTN1", "display": "Rtn1-Cherry", "role": "reported duplicated Rtn1 western-blot lanes"},

    {"case": "Yeast Pgk1 loading-control lane mislabeling", "gene": "PGK1", "display": "Pgk1", "role": "reported loading-control lane mislabeling"},
]


def _truthy_series(s: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).to_numpy(dtype=float)
    vals = s.fillna("").astype(str).str.strip().str.lower()
    return vals.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=float)


def _nonempty_series(s: pd.Series) -> np.ndarray:
    vals = s.fillna("").astype(str).str.strip()
    return (~vals.isin({"", "nan", "none", "na", "n/a"})).to_numpy(dtype=float)


def _split_families(value) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "n/a"}:
        return []
    # Full-report output uses semicolon-separated family labels.
    return [x.strip() for x in re.split(r"\s*;\s*", text) if x.strip()]


def prepare_annotation_features(df: pd.DataFrame) -> tuple[dict[str, np.ndarray], list[list[str]], dict[str, int]]:
    """
    Recover annotation features already produced by wb_integrity_atlas_full_report.py.

    The auditability script is assumed to run after the full-report script, so
    re-annotation is intentionally avoided. This prevents divergence between two
    independently implemented marker/signaling/family matching pipelines.
    """
    n = len(df)

    if "is_marker" in df.columns:
        marker = _truthy_series(df["is_marker"])
    elif "marker_categories" in df.columns:
        marker = _nonempty_series(df["marker_categories"])
    else:
        marker = np.zeros(n, dtype=float)

    if "is_signaling_prefix" in df.columns:
        signaling = _truthy_series(df["is_signaling_prefix"])
    elif "signaling_categories" in df.columns:
        signaling = _nonempty_series(df["signaling_categories"])
    else:
        signaling = np.zeros(n, dtype=float)

    if "families" in df.columns:
        family_lists = [_split_families(v) for v in df["families"]]
    else:
        family_lists = [[] for _ in range(n)]

    family_memberships = np.array([len(set(x)) for x in family_lists], dtype=float)
    family_flag = (family_memberships > 0).astype(float)

    family_counts = Counter()
    for labels in family_lists:
        family_counts.update(set(labels))

    features = {
        "marker": marker,
        "signaling": signaling,
        # family_mass counts memberships, while family_protein_mass counts proteins
        # with >=1 configured family annotation.
        "family": family_memberships,
        "family_protein": family_flag,
    }
    return features, family_lists, dict(family_counts)


def gaussian_weighted_feature_for_centers(
    centers: np.ndarray,
    masses: np.ndarray,
    feature_values: np.ndarray,
    sigma: float,
    cutoff_sigma: float = 4.0,
) -> np.ndarray:
    """Gaussian-weighted local sum of an arbitrary per-protein feature."""
    order = np.argsort(masses)
    ms = np.asarray(masses, dtype=float)[order]
    fv = np.asarray(feature_values, dtype=float)[order]
    cutoff = cutoff_sigma * sigma
    out = np.zeros(len(centers), dtype=float)
    for idx, m in enumerate(np.asarray(centers, dtype=float)):
        lo = np.searchsorted(ms, m - cutoff, side="left")
        hi = np.searchsorted(ms, m + cutoff, side="right")
        if hi <= lo:
            continue
        d = ms[lo:hi] - m
        w = np.exp(-0.5 * (d / sigma) ** 2)
        out[idx] = float(np.dot(w, fv[lo:hi]))
    return out


def positive_context_percentile(values: np.ndarray) -> np.ndarray:
    """
    Rank positive annotation intensities in [0,1], with exact zeros fixed at 0.
    This avoids giving annotation-free regions an artificial ~0.5 percentile
    simply because many values are tied at zero.
    """
    arr = np.asarray(values, dtype=float)
    out = np.zeros_like(arr, dtype=float)
    mask = np.isfinite(arr) & (arr > 0)
    if mask.any():
        out[mask] = pd.Series(arr[mask]).rank(method="average", pct=True).to_numpy(dtype=float)
    out[~np.isfinite(arr)] = np.nan
    return out


def parse_context_weights(text: str) -> dict[str, float]:
    """Parse ambiguity,marker,signaling,family weights and normalize to sum 1."""
    try:
        vals = [float(x.strip()) for x in str(text).split(",")]
    except Exception as exc:
        raise argparse.ArgumentTypeError("--context-weights must contain four comma-separated numbers") from exc
    if len(vals) != 4 or any((not np.isfinite(v) or v < 0) for v in vals) or sum(vals) <= 0:
        raise argparse.ArgumentTypeError(
            "--context-weights must be four non-negative values with positive sum: ambiguity,marker,signaling,family"
        )
    vals = np.array(vals, dtype=float)
    vals = vals / vals.sum()
    return dict(zip(["ambiguity", "marker", "signaling", "family"], vals))


def augment_with_annotation_context(
    protein_df: pd.DataFrame,
    grid_df: pd.DataFrame,
    masses: np.ndarray,
    features: dict[str, np.ndarray],
    sigma: float,
    cutoff_sigma: float,
    context_weights: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Add Gaussian-smoothed marker, signaling and family context without changing
    the identity entropy itself. The identity entropy remains a property of the
    proteome MW candidate distribution; annotations are separate audit-context
    layers and are combined only after rank normalization.
    """
    protein_df = protein_df.copy()
    grid_df = grid_df.copy()

    for key, outname in [
        ("marker", "marker_gaussian_mass"),
        ("signaling", "signaling_gaussian_mass"),
        ("family", "family_membership_gaussian_mass"),
        ("family_protein", "family_protein_gaussian_mass"),
    ]:
        grid_df[outname] = gaussian_weighted_feature_for_centers(
            grid_df["mw_kda"].to_numpy(dtype=float), masses, features[key], sigma, cutoff_sigma
        )
        protein_df[outname] = gaussian_weighted_feature_for_centers(
            masses, masses, features[key], sigma, cutoff_sigma
        )

    # Grid-level scale-free context components.
    grid_df["marker_context_percentile"] = positive_context_percentile(grid_df["marker_gaussian_mass"].to_numpy())
    grid_df["signaling_context_percentile"] = positive_context_percentile(grid_df["signaling_gaussian_mass"].to_numpy())
    grid_df["family_context_percentile"] = positive_context_percentile(grid_df["family_membership_gaussian_mass"].to_numpy())

    # Existing ambiguity_percentile is already high when MW ambiguity is high.
    grid_df["contextual_audit_burden_index"] = (
        context_weights["ambiguity"] * grid_df["ambiguity_percentile"]
        + context_weights["marker"] * grid_df["marker_context_percentile"]
        + context_weights["signaling"] * grid_df["signaling_context_percentile"]
        + context_weights["family"] * grid_df["family_context_percentile"]
    )
    grid_df["contextual_auditability_index"] = 1.0 - grid_df["contextual_audit_burden_index"]
    grid_df["contextual_audit_burden_score_0_100"] = 100.0 * grid_df["contextual_audit_burden_index"]
    grid_df["contextual_auditability_score_0_100"] = 100.0 * grid_df["contextual_auditability_index"]

    # Rank-based priority categories chosen to mirror the interpretation of
    # high-documentation regions without reusing the old raw-score thresholds.
    grid_df["contextual_priority"] = np.where(
        grid_df["contextual_audit_burden_index"] >= grid_df["contextual_audit_burden_index"].quantile(0.90),
        "high",
        np.where(
            grid_df["contextual_audit_burden_index"] >= grid_df["contextual_audit_burden_index"].quantile(0.75),
            "elevated",
            "routine",
        ),
    )

    # Map protein-centered annotation masses to percentiles using the continuous
    # grid as the reference distribution, which keeps query and protein scores on
    # the same interpretation scale.
    for mass_col, pct_col in [
        ("marker_gaussian_mass", "marker_context_percentile"),
        ("signaling_gaussian_mass", "signaling_context_percentile"),
        ("family_membership_gaussian_mass", "family_context_percentile"),
    ]:
        ref = grid_df[mass_col].to_numpy(dtype=float)
        protein_df[pct_col] = [percentile_of_value(v, ref) if v > 0 else 0.0 for v in protein_df[mass_col]]

    # Protein ambiguity percentile already exists and is retained unchanged.
    protein_df["contextual_audit_burden_index"] = (
        context_weights["ambiguity"] * protein_df["mw_ambiguity_percentile"]
        + context_weights["marker"] * protein_df["marker_context_percentile"]
        + context_weights["signaling"] * protein_df["signaling_context_percentile"]
        + context_weights["family"] * protein_df["family_context_percentile"]
    )
    protein_df["contextual_auditability_index"] = 1.0 - protein_df["contextual_audit_burden_index"]
    protein_df["contextual_audit_burden_score_0_100"] = 100.0 * protein_df["contextual_audit_burden_index"]
    protein_df["contextual_auditability_score_0_100"] = 100.0 * protein_df["contextual_auditability_index"]

    return protein_df, grid_df


def build_family_profiles(
    grid_mw: np.ndarray,
    masses: np.ndarray,
    family_lists: list[list[str]],
    family_counts: dict[str, int],
    sigma: float,
    cutoff_sigma: float,
    min_members: int = 2,
) -> pd.DataFrame:
    rows = []
    labels = [k for k, v in sorted(family_counts.items(), key=lambda kv: (-kv[1], kv[0])) if v >= min_members]
    if not labels:
        return pd.DataFrame(columns=["family", "mw_kda", "gaussian_family_mass", "family_member_count"])
    sets = [set(x) for x in family_lists]
    for family in labels:
        feature = np.array([1.0 if family in s else 0.0 for s in sets], dtype=float)
        vals = gaussian_weighted_feature_for_centers(grid_mw, masses, feature, sigma, cutoff_sigma)
        part = pd.DataFrame({
            "family": family,
            "mw_kda": grid_mw,
            "gaussian_family_mass": vals,
            "family_member_count": family_counts[family],
        })
        rows.append(part)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def merge_legacy_bin_context(atlas_dir: Path, grid_df: pd.DataFrame) -> pd.DataFrame:
    """Create a comparison table with the original full-report 2 kDa bin metrics."""
    path = atlas_dir / "mw_bin_summary.csv"
    if not path.exists():
        return pd.DataFrame()
    legacy = pd.read_csv(path)
    required = {"bin_start_kda", "bin_end_kda"}
    if not required.issubset(legacy.columns):
        return pd.DataFrame()
    legacy = legacy.copy()
    legacy["bin_center_kda"] = (
        pd.to_numeric(legacy["bin_start_kda"], errors="coerce")
        + pd.to_numeric(legacy["bin_end_kda"], errors="coerce")
    ) / 2.0
    legacy = legacy.dropna(subset=["bin_center_kda"]).copy()
    x = grid_df["mw_kda"].to_numpy(dtype=float)
    for col in [
        "effective_candidate_number", "identity_entropy_bits", "ambiguity_percentile",
        "marker_gaussian_mass", "signaling_gaussian_mass", "family_membership_gaussian_mass",
        "contextual_audit_burden_index", "contextual_audit_burden_score_0_100",
    ]:
        if col in grid_df.columns:
            legacy[f"gaussian_{col}"] = np.interp(legacy["bin_center_kda"], x, grid_df[col].to_numpy(dtype=float))
    return legacy

# ---------------------------------------------------------------------------
# Core molecular-weight calculations
# ---------------------------------------------------------------------------

def hard_window_counts(masses: np.ndarray, tolerance: float) -> np.ndarray:
    """
    N_i(δ) = number of atlas proteins with |M_j - M_i| <= δ, including protein i.
    Implemented in O(N log N) using sorting + searchsorted.
    """
    order = np.argsort(masses)
    ms = masses[order]
    left = np.searchsorted(ms, ms - tolerance, side="left")
    right = np.searchsorted(ms, ms + tolerance, side="right")
    counts_sorted = right - left

    counts = np.empty_like(counts_sorted)
    counts[order] = counts_sorted
    return counts.astype(int)


def nearest_neighbor_gap(masses: np.ndarray) -> np.ndarray:
    """
    g_i = min_{j != i} |M_i - M_j|
    """
    order = np.argsort(masses)
    ms = masses[order]
    n = len(ms)

    gaps_sorted = np.full(n, np.nan, dtype=float)
    if n == 1:
        return gaps_sorted

    prev_gap = np.r_[np.inf, np.diff(ms)]
    next_gap = np.r_[np.diff(ms), np.inf]
    gaps_sorted = np.minimum(prev_gap, next_gap)

    gaps = np.empty(n, dtype=float)
    gaps[order] = gaps_sorted
    return gaps


def gaussian_metrics_for_centers(
    centers: np.ndarray,
    masses: np.ndarray,
    sigma: float,
    cutoff_sigma: float = 4.0,
):
    """
    For each center m:
        a(m) = sum_j exp(-(M_j-m)^2/(2 sigma^2))
        p_j(m) = w_j / sum_k w_k
        H(m) = -sum_j p_j log p_j
        N_eff(m) = exp(H(m))
        D(m) = 1 / N_eff(m)

    Only proteins within cutoff_sigma*sigma are evaluated. Tail omission at
    cutoff=4 sigma is negligible for practical visualization.
    """
    masses_sorted = np.sort(masses)
    cutoff = cutoff_sigma * sigma

    ambiguity = np.zeros(len(centers), dtype=float)
    entropy = np.zeros(len(centers), dtype=float)
    neff = np.zeros(len(centers), dtype=float)
    discriminability = np.zeros(len(centers), dtype=float)

    for idx, m in enumerate(centers):
        lo = np.searchsorted(masses_sorted, m - cutoff, side="left")
        hi = np.searchsorted(masses_sorted, m + cutoff, side="right")
        local = masses_sorted[lo:hi]

        if local.size == 0:
            ambiguity[idx] = 0.0
            entropy[idx] = np.nan
            neff[idx] = 0.0
            discriminability[idx] = np.nan
            continue

        d = local - m
        w = np.exp(-0.5 * (d / sigma) ** 2)
        a = float(w.sum())

        if a <= 0:
            ambiguity[idx] = 0.0
            entropy[idx] = np.nan
            neff[idx] = 0.0
            discriminability[idx] = np.nan
            continue

        p = w / a
        h = float(-(p * np.log(np.maximum(p, 1e-300))).sum())
        n_eff = float(np.exp(h))

        ambiguity[idx] = a
        entropy[idx] = h
        neff[idx] = n_eff
        discriminability[idx] = 1.0 / n_eff

    return ambiguity, entropy, neff, discriminability


def gaussian_self_excluded_ambiguity(
    masses: np.ndarray,
    sigma: float,
    cutoff_sigma: float = 4.0,
) -> np.ndarray:
    """
    A_i = sum_{j != i} exp(-(M_j-M_i)^2/(2 sigma^2))
    """
    a_total, _, _, _ = gaussian_metrics_for_centers(
        centers=masses,
        masses=masses,
        sigma=sigma,
        cutoff_sigma=cutoff_sigma,
    )
    # The target's self-weight is exactly exp(0)=1.
    return np.maximum(a_total - 1.0, 0.0)


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def build_bin_table(
    protein_df: pd.DataFrame,
    mw_col: str,
    bin_width: float,
    max_mw: float,
) -> pd.DataFrame:
    edges = np.arange(0, math.ceil(max_mw / bin_width) * bin_width + bin_width, bin_width)
    if len(edges) < 2:
        edges = np.array([0.0, bin_width])

    bins = pd.cut(
        protein_df[mw_col],
        bins=edges,
        right=False,
        include_lowest=True,
    )

    tmp = protein_df.copy()
    tmp["_mw_bin"] = bins

    agg = (
        tmp.groupby("_mw_bin", observed=False)
        .agg(
            protein_count=(mw_col, "count"),
            mean_mw_kda=(mw_col, "mean"),
            median_window_count=("mw_window_count", "median"),
            mean_window_count=("mw_window_count", "mean"),
            median_gaussian_ambiguity=("mw_gaussian_ambiguity", "median"),
            mean_gaussian_ambiguity=("mw_gaussian_ambiguity", "mean"),
            median_identity_entropy_nats=("mw_identity_entropy_nats", "median"),
            mean_identity_entropy_nats=("mw_identity_entropy_nats", "mean"),
            median_effective_candidates=("mw_effective_candidate_number", "median"),
            mean_effective_candidates=("mw_effective_candidate_number", "mean"),
            median_discriminability=("mw_discriminability", "median"),
            mean_discriminability=("mw_discriminability", "mean"),
            median_nearest_neighbor_gap_kda=("mw_nearest_neighbor_gap_kda", "median"),
        )
        .reset_index()
    )

    agg["bin_start_kda"] = agg["_mw_bin"].map(lambda x: float(x.left))
    agg["bin_end_kda"] = agg["_mw_bin"].map(lambda x: float(x.right))
    agg["mw_bin_label"] = agg.apply(
        lambda r: f"{r['bin_start_kda']:.1f}-{r['bin_end_kda']:.1f} kDa", axis=1
    )
    agg = agg.drop(columns=["_mw_bin"])

    agg["protein_density_per_kda"] = agg["protein_count"] / bin_width
    agg["ambiguity_percentile"] = empirical_percentile(
        agg["mean_effective_candidates"].fillna(0).to_numpy()
    )
    agg["auditability_percentile"] = 1.0 - agg["ambiguity_percentile"]

    first = [
        "mw_bin_label",
        "bin_start_kda",
        "bin_end_kda",
        "protein_count",
        "protein_density_per_kda",
        "mean_effective_candidates",
        "median_effective_candidates",
        "mean_identity_entropy_nats",
        "mean_discriminability",
        "ambiguity_percentile",
        "auditability_percentile",
    ]
    remaining = [c for c in agg.columns if c not in first]
    return agg[first + remaining]


def build_grid_table(
    masses: np.ndarray,
    grid_step: float,
    sigma: float,
    tolerance: float,
    cutoff_sigma: float,
    max_mw: float,
) -> pd.DataFrame:
    grid = np.arange(0, max_mw + grid_step, grid_step)

    a, h, n_eff, d = gaussian_metrics_for_centers(
        centers=grid,
        masses=masses,
        sigma=sigma,
        cutoff_sigma=cutoff_sigma,
    )

    sorted_masses = np.sort(masses)
    left = np.searchsorted(sorted_masses, grid - tolerance, side="left")
    right = np.searchsorted(sorted_masses, grid + tolerance, side="right")
    window_count = right - left

    out = pd.DataFrame(
        {
            "mw_kda": grid,
            "window_candidate_count": window_count.astype(int),
            "gaussian_weighted_candidate_mass": a,
            "identity_entropy_nats": h,
            "identity_entropy_bits": h / np.log(2.0),
            "effective_candidate_number": n_eff,
            "mw_discriminability": d,
        }
    )
    out["ambiguity_percentile"] = empirical_percentile(
        out["effective_candidate_number"].fillna(0).to_numpy()
    )
    out["auditability_percentile"] = 1.0 - out["ambiguity_percentile"]
    return out



# ---------------------------------------------------------------------------
# Query-MW analysis
# ---------------------------------------------------------------------------

def analyze_query_mw(
    query_mw: float,
    protein_df: pd.DataFrame,
    masses: np.ndarray,
    mw_col: str,
    tolerance: float,
    sigma: float,
    grid_df: pd.DataFrame,
    gene_col: Optional[str] = None,
    accession_col: Optional[str] = None,
    protein_name_col: Optional[str] = None,
):
    """
    Analyze an observed/claimed western-blot band centered at query_mw.

    Unlike protein-centered metrics, the query is an arbitrary observed MW and
    does not have to equal the theoretical MW of any atlas protein.

    Gaussian compatibility:
        w_j(q) = exp(-(M_j-q)^2 / (2 sigma^2))

    Normalized compatibility:
        p_j(q) = w_j(q) / sum_k w_k(q)

    Entropy:
        H(q) = -sum_j p_j(q) ln p_j(q)

    Effective candidate number:
        N_eff(q) = exp(H(q))

    MW discriminability:
        D(q) = 1 / N_eff(q)
    """
    if sigma <= 0:
        raise ValueError("--sigma-kda must be > 0.")
    if tolerance < 0:
        raise ValueError("--tolerance-kda must be >= 0.")

    delta = masses - query_mw
    abs_delta = np.abs(delta)

    # Every proteome entry is retained. Very distant proteins naturally receive
    # weights that underflow to zero, which is appropriate for the candidate TSV.
    gaussian_weight = np.exp(-0.5 * (delta / sigma) ** 2)
    gaussian_mass = float(gaussian_weight.sum())

    if gaussian_mass > 0:
        normalized = gaussian_weight / gaussian_mass
        nz = normalized > 0
        entropy_nats = float(-(normalized[nz] * np.log(normalized[nz])).sum())
        entropy_bits = float(entropy_nats / np.log(2.0))
        effective_candidates = float(np.exp(entropy_nats))
        discriminability = float(1.0 / effective_candidates)
    else:
        normalized = np.zeros_like(gaussian_weight)
        entropy_nats = np.nan
        entropy_bits = np.nan
        effective_candidates = 0.0
        discriminability = np.nan

    hard_mask = abs_delta <= tolerance
    hard_count = int(hard_mask.sum())

    # Auditability percentile is defined against the continuous MW landscape.
    # Higher N_eff means more ambiguity, so auditability reverses that percentile.
    ambiguity_percentile = percentile_of_value(
        effective_candidates,
        grid_df["effective_candidate_number"].to_numpy(dtype=float),
    )
    auditability_percentile = (
        1.0 - ambiguity_percentile if np.isfinite(ambiguity_percentile) else np.nan
    )
    auditability_score = (
        100.0 * auditability_percentile if np.isfinite(auditability_percentile) else np.nan
    )

    candidate = protein_df.copy()

    # Use paper-friendly names for query-specific fields while keeping all
    # original atlas annotations in the table.
    candidate.insert(0, "query_mw_kda", query_mw)
    candidate["mw_difference_kda"] = candidate[mw_col].to_numpy(dtype=float) - query_mw
    candidate["absolute_mw_difference_kda"] = np.abs(candidate["mw_difference_kda"])
    candidate["within_hard_window"] = hard_mask
    candidate["gaussian_weight"] = gaussian_weight
    candidate["normalized_compatibility"] = normalized

    # Rank 1 = strongest MW-compatible candidate.
    candidate["compatibility_rank"] = (
        pd.Series(normalized, index=candidate.index)
        .rank(method="first", ascending=False)
        .astype(int)
    )

    # Reorder the most useful fields to the front.
    preferred = [
        "query_mw_kda",
        gene_col,
        accession_col,
        protein_name_col,
        mw_col,
        "mw_difference_kda",
        "absolute_mw_difference_kda",
        "within_hard_window",
        "gaussian_weight",
        "normalized_compatibility",
        "compatibility_rank",
    ]
    preferred = [c for c in preferred if c is not None and c in candidate.columns]
    remaining = [c for c in candidate.columns if c not in preferred]
    candidate = candidate[preferred + remaining]

    candidate = candidate.sort_values(
        ["normalized_compatibility", "absolute_mw_difference_kda"],
        ascending=[False, True],
    ).reset_index(drop=True)

    summary = {
        "query_mw_kda": float(query_mw),
        "hard_window_tolerance_kda": float(tolerance),
        "hard_window_candidates": hard_count,
        "gaussian_sigma_kda": float(sigma),
        "gaussian_candidate_mass": gaussian_mass,
        "mw_entropy_nats": entropy_nats,
        "mw_entropy_bits": entropy_bits,
        "effective_candidate_number": effective_candidates,
        "mw_discriminability": discriminability,
        "mw_ambiguity_percentile": ambiguity_percentile,
        "mw_auditability_percentile": auditability_percentile,
        "mw_auditability_score_0_100": auditability_score,
        "closest_protein_distance_kda": float(abs_delta.min()) if len(abs_delta) else np.nan,
        "proteome_entries_evaluated": int(len(masses)),
    }
    return summary, candidate


def make_query_local_landscape(
    query_mw: float,
    masses: np.ndarray,
    sigma: float,
    tolerance: float,
    grid_df: pd.DataFrame,
    half_width_kda: float,
    step_kda: float,
) -> pd.DataFrame:
    """
    Build a dense local query-centered MW landscape for plotting/export.
    """
    lo = max(0.0, query_mw - half_width_kda)
    hi = query_mw + half_width_kda
    local_grid = np.arange(lo, hi + step_kda, step_kda)

    a, h, neff, discr = gaussian_metrics_for_centers(
        centers=local_grid,
        masses=masses,
        sigma=sigma,
        cutoff_sigma=8.0,
    )

    sorted_masses = np.sort(masses)
    left = np.searchsorted(sorted_masses, local_grid - tolerance, side="left")
    right = np.searchsorted(sorted_masses, local_grid + tolerance, side="right")
    hard_count = right - left

    ambiguity_ref = grid_df["effective_candidate_number"].to_numpy(dtype=float)
    ambiguity_pct = np.array(
        [percentile_of_value(v, ambiguity_ref) for v in neff],
        dtype=float,
    )

    return pd.DataFrame(
        {
            "mw_kda": local_grid,
            "distance_from_query_kda": local_grid - query_mw,
            "hard_window_candidate_count": hard_count.astype(int),
            "gaussian_candidate_mass": a,
            "mw_entropy_nats": h,
            "mw_entropy_bits": h / np.log(2.0),
            "effective_candidate_number": neff,
            "mw_discriminability": discr,
            "mw_ambiguity_percentile": ambiguity_pct,
            "mw_auditability_percentile": 1.0 - ambiguity_pct,
        }
    )


def make_query_figures(
    query_mw: float,
    candidate_df: pd.DataFrame,
    local_df: pd.DataFrame,
    out_dir: Path,
    mw_col: str,
    sigma: float,
    tolerance: float,
    dpi: int,
    top_n: int = 25,
):
    out_dir.mkdir(parents=True, exist_ok=True)

    # Figure 1: Gaussian compatibility around the claimed band.
    local_candidates = candidate_df[
        candidate_df["absolute_mw_difference_kda"] <= max(5.0 * sigma, tolerance)
    ].copy()

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    if len(local_candidates):
        ax.scatter(
            local_candidates[mw_col],
            local_candidates["gaussian_weight"],
            s=12,
            alpha=0.45,
            label="Proteome entries",
        )

    curve_x = np.linspace(
        max(0.0, query_mw - max(5.0 * sigma, tolerance)),
        query_mw + max(5.0 * sigma, tolerance),
        500,
    )
    curve_y = np.exp(-0.5 * ((curve_x - query_mw) / sigma) ** 2)
    ax.plot(curve_x, curve_y, linewidth=2, label="Gaussian MW compatibility")
    ax.axvline(query_mw, linestyle="--", linewidth=1.5, label=f"Query MW = {query_mw:g} kDa")
    ax.axvspan(
        max(0.0, query_mw - tolerance),
        query_mw + tolerance,
        alpha=0.10,
        label=f"Hard window ±{tolerance:g} kDa",
    )
    ax.set_xlabel("Theoretical molecular weight (kDa)")
    ax.set_ylabel("Gaussian compatibility weight")
    ax.set_title(f"MW compatibility centered on claimed band at {query_mw:g} kDa")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    save_figure(fig, out_dir / "Figure_query_gaussian_compatibility", dpi)

    # Figure 2: top normalized candidate proteins.
    top = candidate_df.head(top_n).copy()
    if len(top):
        label_col = None
        for c in ["gene_symbol", "gene", "primary_gene_symbol", "gene_name",
                  "uniprot_accession", "accession", "protein_name", "name"]:
            if c in top.columns:
                label_col = c
                break

        if label_col is None:
            top["_candidate_label"] = [
                f"Protein {i+1} ({mw:.2f} kDa)"
                for i, mw in enumerate(top[mw_col].to_numpy(dtype=float))
            ]
        else:
            top["_candidate_label"] = top[label_col].astype(str)
            top["_candidate_label"] = [
                f"{label} ({mw:.2f} kDa)"
                for label, mw in zip(
                    top["_candidate_label"],
                    top[mw_col].to_numpy(dtype=float),
                )
            ]

        top = top.iloc[::-1]
        fig, ax = plt.subplots(figsize=(9.5, max(6.0, 0.30 * len(top) + 1.8)))
        ax.barh(top["_candidate_label"], top["normalized_compatibility"])
        ax.set_xlabel("Normalized MW compatibility")
        ax.set_ylabel("Candidate protein")
        ax.set_title(
            f"Top {len(top)} MW-compatible proteome candidates for {query_mw:g} kDa"
        )
        ax.grid(axis="x", alpha=0.2)
        save_figure(fig, out_dir / "Figure_query_top_candidates", dpi)

    # Figure 3: local auditability landscape around query.
    fig, ax1 = plt.subplots(figsize=(9.5, 5.8))
    ax1.plot(
        local_df["mw_kda"],
        local_df["effective_candidate_number"],
        linewidth=1.8,
    )
    ax1.axvline(query_mw, linestyle="--", linewidth=1.5)
    ax1.set_xlabel("Apparent molecular weight (kDa)")
    ax1.set_ylabel("Effective candidate number")
    ax1.set_title(f"Local MW auditability landscape around {query_mw:g} kDa")
    ax1.grid(alpha=0.2)

    ax2 = ax1.twinx()
    ax2.plot(
        local_df["mw_kda"],
        local_df["mw_auditability_percentile"],
        alpha=0.7,
    )
    ax2.set_ylabel("MW auditability percentile")
    ax2.set_ylim(0, 1)
    fig.tight_layout()
    save_figure(fig, out_dir / "Figure_query_local_auditability", dpi)


def write_query_outputs(
    query_mws: list[float],
    protein_df: pd.DataFrame,
    masses: np.ndarray,
    mw_col: str,
    tolerance: float,
    sigma: float,
    grid_df: pd.DataFrame,
    query_root: Path,
    gene_col: Optional[str],
    accession_col: Optional[str],
    protein_name_col: Optional[str],
    dpi: int,
    local_half_width_kda: float,
    local_step_kda: float,
    query_top_n: int,
) -> pd.DataFrame:
    """
    Write combined and per-query outputs. Returns the combined summary table.
    """
    query_root.mkdir(parents=True, exist_ok=True)
    summaries = []

    for q in query_mws:
        label = query_label(q)
        qdir = query_root / label
        qfig = qdir / "Figures"
        qdir.mkdir(parents=True, exist_ok=True)
        qfig.mkdir(parents=True, exist_ok=True)

        summary, candidate = analyze_query_mw(
            query_mw=q,
            protein_df=protein_df,
            masses=masses,
            mw_col=mw_col,
            tolerance=tolerance,
            sigma=sigma,
            grid_df=grid_df,
            gene_col=gene_col,
            accession_col=accession_col,
            protein_name_col=protein_name_col,
        )
        summaries.append(summary)

        # Wide one-row summary for convenient downstream statistics.
        pd.DataFrame([summary]).to_csv(
            qdir / "query_mw_summary_wide.tsv",
            sep="\t",
            index=False,
        )

        # Human-readable Metric / Value table requested for each observed band.
        summary_long = pd.DataFrame(
            [
                ("Query MW", f"{summary['query_mw_kda']:.6g} kDa"),
                ("Hard-window candidates", summary["hard_window_candidates"]),
                ("Hard-window tolerance", f"±{summary['hard_window_tolerance_kda']:.6g} kDa"),
                ("Gaussian sigma", f"{summary['gaussian_sigma_kda']:.6g} kDa"),
                ("Gaussian candidate mass", f"{summary['gaussian_candidate_mass']:.8g}"),
                ("MW entropy", f"{summary['mw_entropy_nats']:.8g} nats"),
                ("MW entropy", f"{summary['mw_entropy_bits']:.8g} bits"),
                ("Effective candidates", f"{summary['effective_candidate_number']:.8g}"),
                ("MW discriminability", f"{summary['mw_discriminability']:.8g}"),
                ("MW ambiguity percentile", f"{summary['mw_ambiguity_percentile']:.8g}"),
                ("MW auditability percentile", f"{summary['mw_auditability_percentile']:.8g}"),
                ("MW auditability score", f"{summary['mw_auditability_score_0_100']:.4f} / 100"),
                ("Closest proteome MW distance", f"{summary['closest_protein_distance_kda']:.8g} kDa"),
                ("Proteome entries evaluated", summary["proteome_entries_evaluated"]),
            ],
            columns=["Metric", "Value"],
        )
        summary_long.to_csv(
            qdir / "query_mw_summary.tsv",
            sep="\t",
            index=False,
        )

        candidate.to_csv(
            qdir / "query_mw_candidates.tsv",
            sep="\t",
            index=False,
        )

        local_df = make_query_local_landscape(
            query_mw=q,
            masses=masses,
            sigma=sigma,
            tolerance=tolerance,
            grid_df=grid_df,
            half_width_kda=local_half_width_kda,
            step_kda=local_step_kda,
        )
        local_df.to_csv(
            qdir / "query_mw_local_landscape.tsv",
            sep="\t",
            index=False,
        )

        make_query_figures(
            query_mw=q,
            candidate_df=candidate,
            local_df=local_df,
            out_dir=qfig,
            mw_col=mw_col,
            sigma=sigma,
            tolerance=tolerance,
            dpi=dpi,
            top_n=query_top_n,
        )

    combined = pd.DataFrame(summaries)
    combined.to_csv(
        query_root / "query_mw_summary_all.tsv",
        sep="\t",
        index=False,
    )
    return combined




def _publication_axes(ax, xgrid=True, ygrid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if xgrid or ygrid:
        ax.grid(axis="both" if (xgrid and ygrid) else ("x" if xgrid else "y"), alpha=0.16, linewidth=0.6)
    ax.tick_params(labelsize=9)


def make_context_figures(
    protein_df: pd.DataFrame,
    grid_df: pd.DataFrame,
    family_profile_df: pd.DataFrame,
    legacy_context_df: pd.DataFrame,
    fig_dir: Path,
    mw_col: str,
    max_plot_mw: float,
    top_family_plots: int,
    dpi: int,
):
    """Generate annotation-aware Gaussian and composite publication figures."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    g = grid_df[grid_df["mw_kda"] <= max_plot_mw].copy()
    p = protein_df[protein_df[mw_col] <= max_plot_mw].copy()

    # 9. Annotation context on the continuous MW axis.
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 9.0), sharex=True)
    panels = [
        ("marker_gaussian_mass", "Gaussian-weighted marker mass", "A. Loading/compartment-marker context"),
        ("signaling_gaussian_mass", "Gaussian-weighted signaling mass", "B. Signaling-protein context"),
        ("family_membership_gaussian_mass", "Gaussian-weighted family memberships", "C. Protein-family context"),
    ]
    for ax, (col, ylabel, title) in zip(axes, panels):
        ax.fill_between(g["mw_kda"], 0, g[col], alpha=0.28)
        ax.plot(g["mw_kda"], g[col], linewidth=1.2)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
        _publication_axes(ax, xgrid=False, ygrid=True)
    axes[-1].set_xlabel("Apparent molecular weight (kDa)")
    fig.suptitle("Gaussian-smoothed audit-relevant annotation context across molecular weight", fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout()
    save_figure(fig, fig_dir / "Figure_annotation_context_gaussian", dpi)

    # 10. Contextual burden and identity-only ambiguity, deliberately shown together.
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    ax.plot(g["mw_kda"], g["ambiguity_percentile"], linewidth=1.35, label="MW ambiguity percentile")
    ax.plot(g["mw_kda"], g["contextual_audit_burden_index"], linewidth=1.65, label="Contextual audit-burden index")
    ax.plot(g["mw_kda"], g["contextual_auditability_index"], linewidth=1.15, alpha=0.75, label="Contextual auditability")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Apparent molecular weight (kDa)")
    ax.set_ylabel("Scale-free index")
    ax.set_title("Identity ambiguity and annotation-aware audit context")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    _publication_axes(ax)
    save_figure(fig, fig_dir / "Figure_contextual_audit_burden", dpi)

    # 11. How much contextual annotation changes prioritization at protein level.
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    sc = ax.scatter(
        p["mw_ambiguity_percentile"], p["contextual_audit_burden_index"],
        c=p[mw_col], cmap="viridis", s=9, alpha=0.30, linewidths=0,
    )
    lim = [0, 1]
    ax.plot(lim, lim, linestyle="--", linewidth=1, alpha=0.7)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("MW ambiguity percentile")
    ax.set_ylabel("Contextual audit-burden index")
    ax.set_title("Annotation context modifies MW-only prioritization")
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Theoretical MW (kDa)")
    _publication_axes(ax)
    save_figure(fig, fig_dir / "Figure_ambiguity_vs_contextual_burden", dpi)

    # 12. Family heatmap.
    if not family_profile_df.empty:
        fam_counts = (
            family_profile_df[["family", "family_member_count"]]
            .drop_duplicates()
            .sort_values("family_member_count", ascending=False)
            .head(top_family_plots)
        )
        fams = fam_counts["family"].tolist()
        fp = family_profile_df[
            family_profile_df["family"].isin(fams) & (family_profile_df["mw_kda"] <= max_plot_mw)
        ].copy()
        pivot = fp.pivot(index="family", columns="mw_kda", values="gaussian_family_mass").reindex(fams)
        # Normalize within family to emphasize location rather than family size.
        mat = pivot.to_numpy(dtype=float)
        denom = np.nanmax(mat, axis=1, keepdims=True)
        denom[~np.isfinite(denom) | (denom == 0)] = 1.0
        matn = mat / denom
        fig, ax = plt.subplots(figsize=(11.5, max(5.5, 0.42 * len(fams) + 2.0)))
        im = ax.imshow(matn, aspect="auto", interpolation="nearest", origin="upper", cmap="magma", extent=[pivot.columns.min(), pivot.columns.max(), len(fams)-0.5, -0.5])
        ax.set_yticks(np.arange(len(fams)))
        ax.set_yticklabels(fams, fontsize=8)
        ax.set_xlabel("Apparent molecular weight (kDa)")
        ax.set_ylabel("Configured protein family/context")
        ax.set_title("MW localization of configured protein-family contexts")
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label("Within-family normalized Gaussian intensity")
        fig.tight_layout()
        save_figure(fig, fig_dir / "Figure_family_context_heatmap", dpi)

    # 13. Legacy 2-kDa implementation versus Gaussian representation.
    if not legacy_context_df.empty and "protein_count" in legacy_context_df.columns:
        leg = legacy_context_df[legacy_context_df["bin_center_kda"] <= max_plot_mw].copy()
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), sharex=True)
        # Panel A shows discontinuous bin counts versus smooth Gaussian candidate mass.
        axes[0,0].step(leg["bin_center_kda"], leg["protein_count"], where="mid", linewidth=1.1, label="2 kDa protein count")
        axb = axes[0,0].twinx()
        axb.plot(g["mw_kda"], g["gaussian_weighted_candidate_mass"], linewidth=1.2, alpha=0.75, label="Gaussian candidate mass")
        axes[0,0].set_ylabel("Proteins / 2 kDa bin")
        axb.set_ylabel("Gaussian candidate mass")
        axes[0,0].set_title("A. Proteome density")
        _publication_axes(axes[0,0], xgrid=False, ygrid=True)

        for ax, old_col, new_col, title in [
            (axes[0,1], "marker_count", "marker_gaussian_mass", "B. Marker context"),
            (axes[1,0], "signaling_prefix_count", "signaling_gaussian_mass", "C. Signaling context"),
            (axes[1,1], "family_count", "family_membership_gaussian_mass", "D. Family context"),
        ]:
            if old_col in leg.columns:
                ax.step(leg["bin_center_kda"], leg[old_col], where="mid", linewidth=1.1, label="2 kDa count")
            ax2 = ax.twinx()
            ax2.plot(g["mw_kda"], g[new_col], linewidth=1.2, alpha=0.78, label="Gaussian context")
            ax.set_title(title)
            ax.set_ylabel("Legacy bin count")
            ax2.set_ylabel("Gaussian weighted mass")
            _publication_axes(ax, xgrid=False, ygrid=True)
        axes[1,0].set_xlabel("Molecular weight (kDa)")
        axes[1,1].set_xlabel("Molecular weight (kDa)")
        fig.suptitle("Fixed-bin and Gaussian representations of WB-IntegrityAtlas context", fontsize=13, fontweight="bold")
        fig.tight_layout(rect=[0,0,1,0.97])
        save_figure(fig, fig_dir / "Figure_legacy_2kda_vs_gaussian", dpi)

    # 14. Main-paper contextual multipanel figure.
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    axes[0,0].plot(g["mw_kda"], g["effective_candidate_number"], linewidth=1.25)
    axes[0,0].set_title("A. Effective MW candidate number")
    axes[0,0].set_ylabel("exp(H)")
    axes[0,0].set_xlabel("MW (kDa)")
    _publication_axes(axes[0,0])

    axes[0,1].plot(g["mw_kda"], g["marker_gaussian_mass"], label="Markers", linewidth=1.15)
    axes[0,1].plot(g["mw_kda"], g["signaling_gaussian_mass"], label="Signaling", linewidth=1.15)
    axes[0,1].plot(g["mw_kda"], g["family_membership_gaussian_mass"], label="Family memberships", linewidth=1.15)
    axes[0,1].set_title("B. Gaussian annotation context")
    axes[0,1].set_xlabel("MW (kDa)")
    axes[0,1].legend(frameon=False, fontsize=8)
    _publication_axes(axes[0,1])

    axes[1,0].plot(g["mw_kda"], g["ambiguity_percentile"], label="MW ambiguity", linewidth=1.15)
    axes[1,0].plot(g["mw_kda"], g["contextual_audit_burden_index"], label="Contextual burden", linewidth=1.5)
    axes[1,0].set_ylim(0,1)
    axes[1,0].set_title("C. Contextual prioritization")
    axes[1,0].set_xlabel("MW (kDa)")
    axes[1,0].legend(frameon=False, fontsize=8)
    _publication_axes(axes[1,0])

    axes[1,1].scatter(p["mw_effective_candidate_number"], p["contextual_audit_burden_index"], s=6, alpha=0.25)
    axes[1,1].set_xlabel("Effective candidate number")
    axes[1,1].set_ylabel("Contextual burden index")
    axes[1,1].set_title("D. Protein-level relationship")
    _publication_axes(axes[1,1])

    fig.suptitle("WB-IntegrityAtlas: Gaussian molecular-weight ambiguity and audit context", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.97])
    save_figure(fig, fig_dir / "Figure_contextual_summary_multipanel", dpi)


def generate_case_mapping_outputs(
    protein_df: pd.DataFrame,
    grid_df: pd.DataFrame,
    legacy_context_df: pd.DataFrame,
    out_dir: Path,
    mw_col: str,
    gene_col: Optional[str],
    dpi: int,
    mode: str,
) -> None:
    """Integrate the uploaded Figure-6 case map in both legacy-bin and Gaussian forms."""
    if mode == "off" or gene_col is None:
        return
    work = protein_df.copy()
    work["_gene_upper"] = work[gene_col].fillna("").astype(str).str.upper()
    rows = []
    for rec in CASE_PROTEINS:
        hits = work[work["_gene_upper"] == rec["gene"].upper()].copy()
        if hits.empty:
            continue
        # If multiple isoforms/entries exist, choose the one closest to the median MW
        # among matching entries to avoid selecting an extreme fragment.
        med = hits[mw_col].median()
        hit = hits.iloc[(hits[mw_col] - med).abs().argsort().iloc[0]]
        row = dict(rec)
        row.update({
            "protein_name": hit.get("protein_name", ""),
            "mw_kda": float(hit[mw_col]),
            "mw_effective_candidate_number": float(hit.get("mw_effective_candidate_number", np.nan)),
            "mw_auditability_score_0_100": float(hit.get("mw_auditability_score_0_100", np.nan)),
            "contextual_audit_burden_index": float(hit.get("contextual_audit_burden_index", np.nan)),
            "contextual_audit_burden_score_0_100": float(hit.get("contextual_audit_burden_score_0_100", np.nan)),
        })
        rows.append(row)
    cases = pd.DataFrame(rows)
    if cases.empty or (mode == "auto" and len(cases) < 3):
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    cases.sort_values(["case", "mw_kda"]).to_csv(out_dir / "Figure6_case_mapping_mapped_proteins.tsv", sep="\t", index=False)

    xmin = max(0.0, float(cases["mw_kda"].min()) - 8.0)
    xmax = float(cases["mw_kda"].max()) + 8.0
    g = grid_df[(grid_df["mw_kda"] >= xmin) & (grid_df["mw_kda"] <= xmax)].copy()

    # Recreate the original 2-kDa-background idea, but with labels placed using
    # a deterministic vertical ladder rather than fixed hard-coded pixel offsets.
    if not legacy_context_df.empty and "protein_count" in legacy_context_df.columns:
        leg = legacy_context_df[(legacy_context_df["bin_center_kda"] >= xmin) & (legacy_context_df["bin_center_kda"] <= xmax)].copy()
    else:
        edges = np.arange(math.floor(xmin/2)*2, math.ceil(xmax/2)*2 + 2, 2)
        cats = pd.cut(protein_df[mw_col], edges, right=False)
        counts = protein_df.groupby(cats, observed=False)[mw_col].count()
        leg = pd.DataFrame({"bin_center_kda": [c.mid for c in counts.index], "protein_count": counts.values})

    fig, ax = plt.subplots(figsize=(13.5, 7.5))
    ax.bar(leg["bin_center_kda"] - 1.0, leg["protein_count"], width=1.9, align="edge", alpha=0.24, label="Proteins per 2 kDa bin")
    ymax = max(float(leg["protein_count"].max()) if len(leg) else 1.0, 1.0)
    levels = np.linspace(ymax*0.58, ymax*0.96, max(len(cases), 2))
    for level, (_, row) in zip(levels, cases.sort_values("mw_kda").iterrows()):
        ax.scatter(row["mw_kda"], level, s=70, zorder=5)
        ax.plot([row["mw_kda"], row["mw_kda"]], [0, level*0.96], linewidth=0.55, alpha=0.35)
        ax.annotate(row["display"], (row["mw_kda"], level), xytext=(4,4), textcoords="offset points", fontsize=8)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Theoretical molecular weight (kDa)")
    ax.set_ylabel("Proteins per 2 kDa bin")
    ax.set_title("Documented protein-identity case examples mapped onto the legacy MW atlas")
    ax.legend(frameon=False)
    _publication_axes(ax, xgrid=False, ygrid=True)
    save_figure(fig, out_dir / "Figure6_case_mapping_legacy_2kda", dpi)

    # Gaussian case map: place each example on the continuous effective-candidate landscape.
    fig, ax = plt.subplots(figsize=(13.5, 7.5))
    ax.fill_between(g["mw_kda"], 0, g["effective_candidate_number"], alpha=0.22, label="Effective candidate landscape")
    ax.plot(g["mw_kda"], g["effective_candidate_number"], linewidth=1.4)
    for _, row in cases.iterrows():
        y = float(np.interp(row["mw_kda"], g["mw_kda"], g["effective_candidate_number"]))
        ax.scatter(row["mw_kda"], y, s=78, zorder=5)
        ax.annotate(row["display"], (row["mw_kda"], y), xytext=(5,6), textcoords="offset points", fontsize=8)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Theoretical molecular weight (kDa)")
    ax.set_ylabel("Effective candidate number, exp(H)")
    ax.set_title("Documented protein-identity case examples on the Gaussian MW-ambiguity landscape")
    ax.legend(frameon=False)
    _publication_axes(ax, xgrid=False, ygrid=True)
    save_figure(fig, out_dir / "Figure6_case_mapping_gaussian", dpi)

    # Side-by-side main manuscript version.
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.5), sharex=True)
    axes[0].bar(leg["bin_center_kda"] - 1.0, leg["protein_count"], width=1.9, align="edge", alpha=0.24)
    for _, row in cases.iterrows():
        bin_y = float(np.interp(row["mw_kda"], leg["bin_center_kda"], leg["protein_count"])) if len(leg) else 0
        axes[0].scatter(row["mw_kda"], bin_y, s=52)
        axes[0].annotate(row["gene"], (row["mw_kda"], bin_y), xytext=(3,5), textcoords="offset points", fontsize=7)
    axes[0].set_title("A. Original 2 kDa bin representation")
    axes[0].set_ylabel("Proteins per 2 kDa bin")
    axes[0].set_xlabel("MW (kDa)")
    _publication_axes(axes[0], xgrid=False, ygrid=True)

    axes[1].fill_between(g["mw_kda"], 0, g["effective_candidate_number"], alpha=0.20)
    axes[1].plot(g["mw_kda"], g["effective_candidate_number"], linewidth=1.25)
    for _, row in cases.iterrows():
        y = float(np.interp(row["mw_kda"], g["mw_kda"], g["effective_candidate_number"]))
        axes[1].scatter(row["mw_kda"], y, s=52)
        axes[1].annotate(row["gene"], (row["mw_kda"], y), xytext=(3,5), textcoords="offset points", fontsize=7)
    axes[1].set_title("B. Gaussian effective-candidate representation")
    axes[1].set_ylabel("Effective candidate number")
    axes[1].set_xlabel("MW (kDa)")
    _publication_axes(axes[1], xgrid=False, ygrid=True)
    fig.suptitle("Figure 6. Documented protein-identity examples mapped to WB-IntegrityAtlas", fontsize=13.5, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.96])
    save_figure(fig, out_dir / "Figure6_case_mapping_multipanel", dpi)


def add_query_annotation_outputs(
    query_mws: list[float],
    query_root: Path,
    protein_df: pd.DataFrame,
    grid_df: pd.DataFrame,
    family_profile_df: pd.DataFrame,
    mw_col: str,
    sigma: float,
    context_weights: dict[str, float],
    dpi: int,
):
    """Append annotation-aware metrics and figures to every --query-mw result."""
    if not query_mws:
        return
    masses = protein_df[mw_col].to_numpy(dtype=float)
    grid_x = grid_df["mw_kda"].to_numpy(dtype=float)
    combined = []
    for q in query_mws:
        qdir = query_root / query_label(q)
        if not qdir.exists():
            continue
        w = np.exp(-0.5 * ((masses - q) / sigma) ** 2)
        row = {
            "query_mw_kda": q,
            "marker_gaussian_mass": float(np.dot(w, protein_df["_audit_marker_value"].to_numpy(dtype=float))),
            "signaling_gaussian_mass": float(np.dot(w, protein_df["_audit_signaling_value"].to_numpy(dtype=float))),
            "family_membership_gaussian_mass": float(np.dot(w, protein_df["_audit_family_value"].to_numpy(dtype=float))),
        }
        for mass_col, pct_col in [
            ("marker_gaussian_mass", "marker_context_percentile"),
            ("signaling_gaussian_mass", "signaling_context_percentile"),
            ("family_membership_gaussian_mass", "family_context_percentile"),
        ]:
            row[pct_col] = percentile_of_value(row[mass_col], grid_df[mass_col].to_numpy(dtype=float)) if row[mass_col] > 0 else 0.0
        row["mw_ambiguity_percentile"] = float(np.interp(q, grid_x, grid_df["ambiguity_percentile"]))
        row["contextual_audit_burden_index"] = (
            context_weights["ambiguity"] * row["mw_ambiguity_percentile"]
            + context_weights["marker"] * row["marker_context_percentile"]
            + context_weights["signaling"] * row["signaling_context_percentile"]
            + context_weights["family"] * row["family_context_percentile"]
        )
        row["contextual_auditability_index"] = 1.0 - row["contextual_audit_burden_index"]
        row["contextual_audit_burden_score_0_100"] = 100.0 * row["contextual_audit_burden_index"]
        row["contextual_auditability_score_0_100"] = 100.0 * row["contextual_auditability_index"]
        combined.append(row)
        pd.DataFrame([row]).to_csv(qdir / "query_mw_annotation_context.tsv", sep="\t", index=False)

        # Extend the existing wide summary.
        wide_path = qdir / "query_mw_summary_wide.tsv"
        if wide_path.exists():
            wide = pd.read_csv(wide_path, sep="\t")
            for k, v in row.items():
                if k != "query_mw_kda": wide[k] = v
            wide.to_csv(wide_path, sep="\t", index=False)
        long_path = qdir / "query_mw_summary.tsv"
        if long_path.exists():
            long = pd.read_csv(long_path, sep="\t")
            extra = pd.DataFrame([
                ("Gaussian marker mass", f"{row['marker_gaussian_mass']:.8g}"),
                ("Gaussian signaling mass", f"{row['signaling_gaussian_mass']:.8g}"),
                ("Gaussian family-membership mass", f"{row['family_membership_gaussian_mass']:.8g}"),
                ("Marker context percentile", f"{row['marker_context_percentile']:.8g}"),
                ("Signaling context percentile", f"{row['signaling_context_percentile']:.8g}"),
                ("Family context percentile", f"{row['family_context_percentile']:.8g}"),
                ("Contextual audit-burden index", f"{row['contextual_audit_burden_index']:.8g}"),
                ("Contextual auditability index", f"{row['contextual_auditability_index']:.8g}"),
                ("Contextual audit-burden score", f"{row['contextual_audit_burden_score_0_100']:.4f} / 100"),
                ("Contextual auditability score", f"{row['contextual_auditability_score_0_100']:.4f} / 100"),
            ], columns=["Metric", "Value"])
            pd.concat([long, extra], ignore_index=True).to_csv(long_path, sep="\t", index=False)

        # Top families at the query MW.
        if not family_profile_df.empty:
            fam_rows = []
            for family, sub in family_profile_df.groupby("family", sort=False):
                val = float(np.interp(q, sub["mw_kda"], sub["gaussian_family_mass"]))
                fam_rows.append((family, val, int(sub["family_member_count"].iloc[0])))
            famq = pd.DataFrame(fam_rows, columns=["family", "gaussian_family_mass", "family_member_count"]).sort_values("gaussian_family_mass", ascending=False)
            famq.to_csv(qdir / "query_mw_family_context.tsv", sep="\t", index=False)

        # Compact query annotation figure.
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
        names = ["Marker", "Signaling", "Family"]
        vals = [row["marker_gaussian_mass"], row["signaling_gaussian_mass"], row["family_membership_gaussian_mass"]]
        axes[0].bar(names, vals, alpha=0.78)
        axes[0].set_ylabel("Gaussian-weighted annotation mass")
        axes[0].set_title(f"A. Annotation context at {q:g} kDa")
        _publication_axes(axes[0], xgrid=False, ygrid=True)
        pcts = [row["mw_ambiguity_percentile"], row["marker_context_percentile"], row["signaling_context_percentile"], row["family_context_percentile"], row["contextual_audit_burden_index"]]
        labels = ["MW ambiguity", "Marker", "Signaling", "Family", "Composite"]
        axes[1].barh(labels, pcts, alpha=0.78)
        axes[1].set_xlim(0,1)
        axes[1].set_xlabel("Percentile / scale-free index")
        axes[1].set_title("B. Context-normalized audit components")
        _publication_axes(axes[1], xgrid=True, ygrid=False)
        fig.suptitle(f"Query-centered audit context for an observed band at {q:g} kDa", fontsize=12.5, fontweight="bold")
        fig.tight_layout(rect=[0,0,1,0.95])
        save_figure(fig, qdir / "Figures" / "Figure_query_annotation_context", dpi)

    if combined:
        pd.DataFrame(combined).to_csv(query_root / "query_mw_annotation_context_all.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(
    protein_df: pd.DataFrame,
    bin_df: pd.DataFrame,
    grid_df: pd.DataFrame,
    fig_dir: Path,
    mw_col: str,
    bin_width: float,
    dpi: int,
    max_plot_mw: float,
):
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_proteins = protein_df[protein_df[mw_col] <= max_plot_mw].copy()
    plot_bins = bin_df[bin_df["bin_start_kda"] <= max_plot_mw].copy()
    plot_grid = grid_df[grid_df["mw_kda"] <= max_plot_mw].copy()

    # 1. MW distribution
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.hist(plot_proteins[mw_col].dropna(), bins=max(20, int(max_plot_mw / bin_width)))
    ax.set_xlabel("Theoretical molecular weight (kDa)")
    ax.set_ylabel("Protein count")
    ax.set_title("Proteome-wide molecular-weight distribution")
    ax.grid(alpha=0.2)
    save_figure(fig, fig_dir / "Figure_MW_distribution", dpi)

    # 2. Density + effective candidates
    fig, ax1 = plt.subplots(figsize=(9.5, 5.6))
    ax1.plot(plot_grid["mw_kda"], plot_grid["gaussian_weighted_candidate_mass"])
    ax1.set_xlabel("Apparent molecular weight (kDa)")
    ax1.set_ylabel("Gaussian-weighted candidate mass")
    ax1.set_title("Molecular-weight identity ambiguity across the proteome")
    ax1.grid(alpha=0.2)

    ax2 = ax1.twinx()
    ax2.plot(plot_grid["mw_kda"], plot_grid["effective_candidate_number"], alpha=0.7)
    ax2.set_ylabel("Effective candidate number, exp(H)")
    fig.tight_layout()
    save_figure(fig, fig_dir / "Figure_MW_density_and_effective_candidates", dpi)

    # 3. Protein-level auditability scatter
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.scatter(
        plot_proteins[mw_col],
        plot_proteins["mw_effective_candidate_number"],
        s=8,
        alpha=0.35,
    )
    ax.set_xlabel("Protein theoretical molecular weight (kDa)")
    ax.set_ylabel("Effective candidate number, exp(H)")
    ax.set_title("Protein-level molecular-weight ambiguity")
    ax.grid(alpha=0.2)
    save_figure(fig, fig_dir / "Figure_protein_auditability_scatter", dpi)

    # 4. Distribution of protein ambiguity
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    vals = protein_df["mw_effective_candidate_number"].replace([np.inf, -np.inf], np.nan).dropna()
    ax.hist(vals, bins=60)
    ax.set_xlabel("Effective candidate number, exp(H)")
    ax.set_ylabel("Protein count")
    ax.set_title("Distribution of MW-based protein-identity ambiguity")
    ax.grid(alpha=0.2)
    save_figure(fig, fig_dir / "Figure_ambiguity_distribution", dpi)

    # 5. Rank bins by ambiguity
    top = (
        plot_bins[plot_bins["protein_count"] > 0]
        .sort_values("mean_effective_candidates", ascending=False)
        .head(20)
        .sort_values("mean_effective_candidates", ascending=True)
    )
    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.barh(top["mw_bin_label"], top["mean_effective_candidates"])
    ax.set_xlabel("Mean effective candidate number")
    ax.set_ylabel("MW bin")
    ax.set_title("Molecular-weight bins with highest identity ambiguity")
    ax.grid(axis="x", alpha=0.2)
    save_figure(fig, fig_dir / "Figure_bin_rankings", dpi)

    # 6. Auditability landscape
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.plot(plot_grid["mw_kda"], plot_grid["auditability_percentile"])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Apparent molecular weight (kDa)")
    ax.set_ylabel("MW auditability percentile")
    ax.set_title("Proteome-wide western-blot MW auditability landscape")
    ax.grid(alpha=0.2)
    save_figure(fig, fig_dir / "Figure_auditability_landscape", dpi)

    # 7. Nearest-neighbor gap
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.scatter(
        plot_proteins[mw_col],
        plot_proteins["mw_nearest_neighbor_gap_kda"],
        s=8,
        alpha=0.35,
    )
    ax.set_xlabel("Protein theoretical molecular weight (kDa)")
    ax.set_ylabel("Nearest-neighbor MW gap (kDa)")
    ax.set_title("Local molecular-weight distinctiveness")
    ax.grid(alpha=0.2)
    save_figure(fig, fig_dir / "Figure_nearest_neighbor_gap", dpi)

    # 8. Multipanel summary
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    axes[0, 0].hist(
        plot_proteins[mw_col].dropna(),
        bins=max(20, int(max_plot_mw / bin_width)),
    )
    axes[0, 0].set_xlabel("Theoretical MW (kDa)")
    axes[0, 0].set_ylabel("Protein count")
    axes[0, 0].set_title("A. Proteome MW distribution")
    axes[0, 0].grid(alpha=0.2)

    axes[0, 1].plot(
        plot_grid["mw_kda"], plot_grid["effective_candidate_number"]
    )
    axes[0, 1].set_xlabel("Apparent MW (kDa)")
    axes[0, 1].set_ylabel("Effective candidate number")
    axes[0, 1].set_title("B. Identity ambiguity")
    axes[0, 1].grid(alpha=0.2)

    axes[1, 0].plot(
        plot_grid["mw_kda"], plot_grid["auditability_percentile"]
    )
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].set_xlabel("Apparent MW (kDa)")
    axes[1, 0].set_ylabel("Auditability percentile")
    axes[1, 0].set_title("C. MW auditability landscape")
    axes[1, 0].grid(alpha=0.2)

    axes[1, 1].scatter(
        plot_proteins[mw_col],
        plot_proteins["mw_nearest_neighbor_gap_kda"],
        s=6,
        alpha=0.3,
    )
    axes[1, 1].set_xlabel("Theoretical MW (kDa)")
    axes[1, 1].set_ylabel("Nearest-neighbor gap (kDa)")
    axes[1, 1].set_title("D. Local distinctiveness")
    axes[1, 1].grid(alpha=0.2)

    fig.suptitle("WB-IntegrityAtlas: molecular-weight auditability of the proteome")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_figure(fig, fig_dir / "Figure_summary_multipanel", dpi)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Quantify proteome-wide western-blot molecular-weight auditability."
    )
    p.add_argument("atlas_dir", type=Path, help="WB-IntegrityAtlas output folder")
    p.add_argument("--protein-file", default="proteome_wb_integrity_atlas.csv")
    p.add_argument("--bin-file", default="mw_bin_summary.csv")
    p.add_argument("--mw-column", default=None)
    p.add_argument("--bin-width", type=float, default=2.0, help="MW bin width in kDa")
    p.add_argument(
        "--tolerance-kda",
        type=float,
        default=2.0,
        help="Half-width δ for hard candidate window |M_j-M_i| <= δ",
    )
    p.add_argument(
        "--sigma-kda",
        type=float,
        default=1.0,
        help="Gaussian MW uncertainty sigma in kDa",
    )
    p.add_argument(
        "--kernel-cutoff-sigma",
        type=float,
        default=4.0,
        help="Gaussian evaluation cutoff in sigma units",
    )
    p.add_argument(
        "--grid-step",
        type=float,
        default=0.25,
        help="MW step for continuous auditability landscape",
    )
    p.add_argument(
        "--max-plot-mw",
        type=float,
        default=300.0,
        help="Maximum MW displayed in figures; calculations retain all proteins",
    )
    p.add_argument("--dpi", type=int, default=400)
    p.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of proteins in top ambiguous/distinctive TSVs",
    )
    p.add_argument(
        "--query-mw",
        default=None,
        help=(
            "Optional observed/claimed band MW or comma-separated list in kDa, "
            "e.g. --query-mw 42.0 or --query-mw 42.0,45.7,100"
        ),
    )
    p.add_argument(
        "--query-local-half-width-kda",
        type=float,
        default=8.0,
        help="Half-width in kDa for query-centered local landscape outputs/plots",
    )
    p.add_argument(
        "--query-local-step-kda",
        type=float,
        default=0.05,
        help="MW step in kDa for query-centered local landscape outputs/plots",
    )
    p.add_argument(
        "--query-top-n",
        type=int,
        default=25,
        help="Number of top MW-compatible proteins shown in each query figure",
    )
    p.add_argument(
        "--context-weights",
        default="1,1,1,1",
        help=(
            "Relative weights for ambiguity,marker,signaling,family in the contextual audit-burden index. "
            "Default equal weights: 1,1,1,1"
        ),
    )
    p.add_argument(
        "--family-min-members",
        type=int,
        default=2,
        help="Minimum annotated proteins required for an individual family Gaussian profile",
    )
    p.add_argument(
        "--top-family-plots",
        type=int,
        default=12,
        help="Number of most prevalent configured families shown in family heatmaps",
    )
    p.add_argument(
        "--case-mapping",
        choices=["auto", "on", "off"],
        default="auto",
        help="Generate Figure-6-style documented case mapping; auto generates when mapped genes are found",
    )
    return p.parse_args()


def main():
    args = parse_args()
    query_mws = parse_query_mw(args.query_mw)
    context_weights = parse_context_weights(args.context_weights)

    if args.query_local_half_width_kda <= 0:
        raise ValueError("--query-local-half-width-kda must be > 0.")
    if args.query_local_step_kda <= 0:
        raise ValueError("--query-local-step-kda must be > 0.")
    if args.query_top_n <= 0:
        raise ValueError("--query-top-n must be > 0.")

    atlas_dir = args.atlas_dir.expanduser().resolve()

    protein_path = atlas_dir / args.protein_file
    if not protein_path.exists():
        raise FileNotFoundError(f"Primary atlas file not found: {protein_path}")

    out_dir = atlas_dir / "mw_auditability_analysis"
    tsv_dir = out_dir / "TSV"
    fig_dir = out_dir / "Figures"
    tsv_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(protein_path)

    mw_col = args.mw_column or find_mw_column(df)
    df[mw_col] = numeric_series(df[mw_col])
    df = df[df[mw_col].notna() & (df[mw_col] > 0)].copy()
    df[mw_col] = df[mw_col].astype(float)

    masses = df[mw_col].to_numpy(dtype=float)
    if len(masses) < 2:
        raise ValueError("Need at least two proteins with valid positive molecular weights.")

    # Optional descriptive identifiers.
    gene_col = first_existing_column(
        df, ["gene_primary", "gene_symbol", "gene", "primary_gene_symbol", "gene_name"]
    )
    accession_col = first_existing_column(
        df, ["uniprot_accession", "accession", "protein_accession", "entry"]
    )
    protein_name_col = first_existing_column(
        df, ["protein_name", "name", "description"]
    )

    # Hard-window count.
    df["mw_window_count"] = hard_window_counts(masses, args.tolerance_kda)
    df["mw_window_alternative_count"] = np.maximum(df["mw_window_count"] - 1, 0)

    # Nearest-neighbor MW distance.
    df["mw_nearest_neighbor_gap_kda"] = nearest_neighbor_gap(masses)

    # Gaussian ambiguity excluding the target itself.
    df["mw_gaussian_ambiguity"] = gaussian_self_excluded_ambiguity(
        masses,
        sigma=args.sigma_kda,
        cutoff_sigma=args.kernel_cutoff_sigma,
    )

    # Information-theoretic identity metrics centered at each protein MW.
    total_a, entropy, neff, discr = gaussian_metrics_for_centers(
        centers=masses,
        masses=masses,
        sigma=args.sigma_kda,
        cutoff_sigma=args.kernel_cutoff_sigma,
    )
    df["mw_gaussian_candidate_mass_including_self"] = total_a
    df["mw_identity_entropy_nats"] = entropy
    df["mw_identity_entropy_bits"] = entropy / np.log(2.0)
    df["mw_effective_candidate_number"] = neff
    df["mw_discriminability"] = discr

    # Percentile interpretation.
    df["mw_ambiguity_percentile"] = empirical_percentile(
        df["mw_effective_candidate_number"].to_numpy()
    )
    df["mw_auditability_percentile"] = 1.0 - df["mw_ambiguity_percentile"]

    # A simple 0-100 paper/report score. High = greater MW distinctiveness.
    df["mw_auditability_score_0_100"] = 100.0 * df["mw_auditability_percentile"]

    # Preserve original columns first, metrics after.
    metric_cols = [
        "mw_window_count",
        "mw_window_alternative_count",
        "mw_nearest_neighbor_gap_kda",
        "mw_gaussian_ambiguity",
        "mw_gaussian_candidate_mass_including_self",
        "mw_identity_entropy_nats",
        "mw_identity_entropy_bits",
        "mw_effective_candidate_number",
        "mw_discriminability",
        "mw_ambiguity_percentile",
        "mw_auditability_percentile",
        "mw_auditability_score_0_100",
    ]

    # Bin-level analysis.
    max_mw = float(np.nanmax(masses))
    bin_df = build_bin_table(
        protein_df=df,
        mw_col=mw_col,
        bin_width=args.bin_width,
        max_mw=max_mw,
    )

    # Continuous MW information landscape.
    grid_df = build_grid_table(
        masses=masses,
        grid_step=args.grid_step,
        sigma=args.sigma_kda,
        tolerance=args.tolerance_kda,
        cutoff_sigma=args.kernel_cutoff_sigma,
        max_mw=max_mw,
    )

    # Recover marker, signaling and family annotations created by the preceding
    # wb_integrity_atlas_full_report.py run and smooth them on the same MW axis.
    annotation_features, family_lists, family_counts = prepare_annotation_features(df)
    # Internal numeric copies are retained in the auditability TSV because they make
    # downstream reproducibility explicit; original full-report annotation columns
    # are preserved unchanged.
    df["_audit_marker_value"] = annotation_features["marker"]
    df["_audit_signaling_value"] = annotation_features["signaling"]
    df["_audit_family_value"] = annotation_features["family"]
    df["_audit_family_protein_value"] = annotation_features["family_protein"]
    df, grid_df = augment_with_annotation_context(
        protein_df=df,
        grid_df=grid_df,
        masses=masses,
        features=annotation_features,
        sigma=args.sigma_kda,
        cutoff_sigma=args.kernel_cutoff_sigma,
        context_weights=context_weights,
    )

    family_profile_df = build_family_profiles(
        grid_mw=grid_df["mw_kda"].to_numpy(dtype=float),
        masses=masses,
        family_lists=family_lists,
        family_counts=family_counts,
        sigma=args.sigma_kda,
        cutoff_sigma=args.kernel_cutoff_sigma,
        min_members=args.family_min_members,
    )

    legacy_context_df = merge_legacy_bin_context(atlas_dir, grid_df)

    # Write TSVs.
    df.to_csv(tsv_dir / "protein_mw_auditability.tsv", sep="\t", index=False)
    bin_df.to_csv(tsv_dir / "mw_bin_auditability.tsv", sep="\t", index=False)
    grid_df.to_csv(tsv_dir / "mw_grid_information.tsv", sep="\t", index=False)
    if not family_profile_df.empty:
        family_profile_df.to_csv(tsv_dir / "family_gaussian_profiles.tsv", sep="\t", index=False)
    if not legacy_context_df.empty:
        legacy_context_df.to_csv(tsv_dir / "legacy_2kda_vs_gaussian_context.tsv", sep="\t", index=False)

    annotation_summary = pd.DataFrame([
        {
            "marker_annotated_proteins": int(annotation_features["marker"].sum()),
            "signaling_annotated_proteins": int(annotation_features["signaling"].sum()),
            "proteins_with_family_annotation": int(annotation_features["family_protein"].sum()),
            "total_family_memberships": int(annotation_features["family"].sum()),
            "distinct_configured_families_observed": int(len(family_counts)),
            "context_weight_ambiguity": context_weights["ambiguity"],
            "context_weight_marker": context_weights["marker"],
            "context_weight_signaling": context_weights["signaling"],
            "context_weight_family": context_weights["family"],
        }
    ])
    annotation_summary.to_csv(tsv_dir / "annotation_context_summary.tsv", sep="\t", index=False)

    id_cols = [c for c in [gene_col, accession_col, protein_name_col, mw_col] if c is not None]
    rank_cols = id_cols + [
        "mw_window_count",
        "mw_gaussian_ambiguity",
        "mw_identity_entropy_bits",
        "mw_effective_candidate_number",
        "mw_nearest_neighbor_gap_kda",
        "mw_auditability_score_0_100",
    ]

    df.sort_values(
        ["mw_effective_candidate_number", "mw_window_count"],
        ascending=[False, False],
    )[rank_cols].head(args.top_n).to_csv(
        tsv_dir / "top_ambiguous_proteins.tsv", sep="\t", index=False
    )

    df.sort_values(
        ["mw_effective_candidate_number", "mw_nearest_neighbor_gap_kda"],
        ascending=[True, False],
    )[rank_cols].head(args.top_n).to_csv(
        tsv_dir / "top_distinctive_proteins.tsv", sep="\t", index=False
    )

    # Summary statistics.
    summary_rows = [
        ("protein_count_analyzed", len(df)),
        ("mw_column", mw_col),
        ("minimum_mw_kda", float(df[mw_col].min())),
        ("median_mw_kda", float(df[mw_col].median())),
        ("mean_mw_kda", float(df[mw_col].mean())),
        ("maximum_mw_kda", float(df[mw_col].max())),
        ("hard_window_tolerance_kda", args.tolerance_kda),
        ("gaussian_sigma_kda", args.sigma_kda),
        ("bin_width_kda", args.bin_width),
        ("median_window_candidate_count", float(df["mw_window_count"].median())),
        ("mean_window_candidate_count", float(df["mw_window_count"].mean())),
        ("median_effective_candidate_number", float(df["mw_effective_candidate_number"].median())),
        ("mean_effective_candidate_number", float(df["mw_effective_candidate_number"].mean())),
        ("median_identity_entropy_bits", float(df["mw_identity_entropy_bits"].median())),
        ("median_nearest_neighbor_gap_kda", float(df["mw_nearest_neighbor_gap_kda"].median())),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(
        tsv_dir / "mw_auditability_summary.tsv", sep="\t", index=False
    )

    # Publication figures.
    make_figures(
        protein_df=df,
        bin_df=bin_df,
        grid_df=grid_df,
        fig_dir=fig_dir,
        mw_col=mw_col,
        bin_width=args.bin_width,
        dpi=args.dpi,
        max_plot_mw=args.max_plot_mw,
    )
    make_context_figures(
        protein_df=df,
        grid_df=grid_df,
        family_profile_df=family_profile_df,
        legacy_context_df=legacy_context_df,
        fig_dir=fig_dir,
        mw_col=mw_col,
        max_plot_mw=args.max_plot_mw,
        top_family_plots=args.top_family_plots,
        dpi=args.dpi,
    )

    case_dir = out_dir / "Case_Mapping"
    generate_case_mapping_outputs(
        protein_df=df,
        grid_df=grid_df,
        legacy_context_df=legacy_context_df,
        out_dir=case_dir,
        mw_col=mw_col,
        gene_col=gene_col,
        dpi=args.dpi,
        mode=args.case_mapping,
    )

    # Optional observed/claimed band analysis.
    if query_mws:
        query_root = out_dir / "Query_MW"
        write_query_outputs(
            query_mws=query_mws,
            protein_df=df,
            masses=masses,
            mw_col=mw_col,
            tolerance=args.tolerance_kda,
            sigma=args.sigma_kda,
            grid_df=grid_df,
            query_root=query_root,
            gene_col=gene_col,
            accession_col=accession_col,
            protein_name_col=protein_name_col,
            dpi=args.dpi,
            local_half_width_kda=args.query_local_half_width_kda,
            local_step_kda=args.query_local_step_kda,
            query_top_n=args.query_top_n,
        )
        add_query_annotation_outputs(
            query_mws=query_mws,
            query_root=query_root,
            protein_df=df,
            grid_df=grid_df,
            family_profile_df=family_profile_df,
            mw_col=mw_col,
            sigma=args.sigma_kda,
            context_weights=context_weights,
            dpi=args.dpi,
        )

    # Parameter record for reproducibility.
    parameter_text = f"""WB-IntegrityAtlas molecular-weight auditability analysis
protein_file={protein_path}
mw_column={mw_col}
n_proteins={len(df)}
bin_width_kda={args.bin_width}
hard_window_tolerance_kda={args.tolerance_kda}
gaussian_sigma_kda={args.sigma_kda}
kernel_cutoff_sigma={args.kernel_cutoff_sigma}
grid_step_kda={args.grid_step}
max_plot_mw_kda={args.max_plot_mw}
query_mw={",".join(f"{q:g}" for q in query_mws) if query_mws else "none"}
query_local_half_width_kda={args.query_local_half_width_kda}
query_local_step_kda={args.query_local_step_kda}
query_top_n={args.query_top_n}
context_weights_ambiguity_marker_signaling_family={args.context_weights}
family_min_members={args.family_min_members}
top_family_plots={args.top_family_plots}
case_mapping={args.case_mapping}

Interpretation:
- Higher effective_candidate_number = greater MW identity ambiguity.
- Higher identity_entropy = greater uncertainty among proteins compatible with the MW.
- Higher mw_discriminability = stronger MW-based identity specificity.
- Higher mw_auditability_score_0_100 = more distinctive MW neighborhood.
- These values measure auditability/identity ambiguity, NOT evidence of misconduct.
"""
    (out_dir / "analysis_parameters.txt").write_text(parameter_text)

    print(f"Done. Results written to: {out_dir}")
    print(f"Proteins analyzed: {len(df)}")
    print(f"MW column used: {mw_col}")
    if query_mws:
        print("Query MW values: " + ", ".join(f"{q:g} kDa" for q in query_mws))
        print(f"Query-specific outputs: {out_dir / 'Query_MW'}")


if __name__ == "__main__":
    main()
