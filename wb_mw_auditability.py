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
        --grid-step 0.25

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
    return p.parse_args()


def main():
    args = parse_args()
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
        df, ["gene_symbol", "gene", "primary_gene_symbol", "gene_name"]
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

    # Write TSVs.
    df.to_csv(tsv_dir / "protein_mw_auditability.tsv", sep="\t", index=False)
    bin_df.to_csv(tsv_dir / "mw_bin_auditability.tsv", sep="\t", index=False)
    grid_df.to_csv(tsv_dir / "mw_grid_information.tsv", sep="\t", index=False)

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


if __name__ == "__main__":
    main()
