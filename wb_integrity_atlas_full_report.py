#!/usr/bin/env python3
"""
Western Blot Integrity Atlas: complete report generator
=======================================================

DEFENSIVE USE ONLY.

This is the combined report runner for the Western Blot Integrity Atlas package.
It preserves the command-line flags of wb_integrity_atlas_externalized.py and,
from the same FASTA/configuration inputs, generates:

1. proteome_wb_integrity_atlas.csv
2. mw_bin_summary.csv
3. reported_malpractice_bins.csv, when --malpractice-file is supplied
4. four standalone plot PNG files in plots/
5. western_blot_integrity_atlas_full_report_with_plots.pdf

No pre-generated CSV is required. In particular, there is no --bin_csv flag.
All tables and figures are generated internally from the FASTA and TSV inputs.

Dependencies:
    pip install pandas matplotlib reportlab
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
from pathlib import Path
from typing import Iterable, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image,
    LongTable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from wb_integrity_atlas_externalized import (
    annotate,
    attach_malpractice_bins,
    html_escape,
    load_family_pattern_file,
    load_malpractice_file,
    load_marker_file,
    load_signaling_prefix_file,
    proteome_dataframe_from_fasta,
    safe_str,
    write_outputs,
)


REPORT_FILENAME = "western_blot_integrity_atlas_full_report_with_plots.pdf"
COMPATIBILITY_REPORT_FILENAME = "western_blot_integrity_atlas_report.pdf"
PLOT_DIRNAME = "plots"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Generate the complete defensive western-blot molecular-weight integrity atlas, "
            "including CSVs, plots, and the consolidated PDF report."
        )
    )
    parser.add_argument(
        "--proteome-fasta",
        required=True,
        help="Proteome FASTA file, preferably a UniProt reference-proteome FASTA",
    )
    parser.add_argument(
        "--marker-file",
        required=True,
        help="TSV with columns gene_symbol, category, note",
    )
    parser.add_argument(
        "--signaling-prefix-file",
        required=True,
        help="TSV with columns prefix, category, note",
    )
    parser.add_argument(
        "--family-pattern-file",
        required=True,
        help="TSV with columns family, regex, category, note",
    )
    parser.add_argument(
        "--malpractice-file",
        default=None,
        help="Optional TSV of reported malpractice cases",
    )
    parser.add_argument(
        "--species-name",
        default="Unspecified species",
        help="Species/report label",
    )
    parser.add_argument(
        "--outdir",
        default="wb_integrity_atlas_out",
        help="Output directory",
    )
    parser.add_argument(
        "--bin-kda",
        type=float,
        default=2.0,
        help="Molecular-weight bin width in kDa",
    )
    return parser.parse_args(argv)


def _numeric(series: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(default)


def normalize_report_dataframes(
    df: pd.DataFrame,
    bins: pd.DataFrame,
    mal_df: Optional[pd.DataFrame] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Normalize data types so report generation also works on reloaded CSVs."""
    df = df.copy()
    bins = bins.copy()
    mal_df = pd.DataFrame() if mal_df is None else mal_df.copy()

    for col in ["mw_kda", "mass_da", "length_aa", "protein_count", "marker_count", "signaling_prefix_count", "family_count", "audit_burden_score"]:
        if col in df.columns:
            df[col] = _numeric(df[col])
    for col in ["bin_start_kda", "bin_end_kda", "protein_count", "marker_count", "signaling_prefix_count", "family_count", "audit_burden_score"]:
        if col in bins.columns:
            bins[col] = _numeric(bins[col])

    for col in ["is_marker", "is_signaling_prefix"]:
        if col in df.columns and df[col].dtype != bool:
            df[col] = (
                df[col]
                .astype(str)
                .str.strip()
                .str.lower()
                .isin({"true", "1", "yes"})
            )

    text_cols_df = [
        "accession", "entry_name", "gene_primary", "protein_name", "families",
        "mw_bin_label", "defensive_note", "marker_categories", "signaling_categories",
    ]
    for col in text_cols_df:
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str)

    text_cols_bins = [
        "mw_bin_label", "families", "marker_genes", "signaling_genes",
        "examples_audit_families", "audit_priority",
    ]
    for col in text_cols_bins:
        if col in bins.columns:
            bins[col] = bins[col].fillna("").astype(str)

    if not mal_df.empty:
        for col in mal_df.columns:
            if col != "approx_mw_kda" and mal_df[col].dtype == object:
                mal_df[col] = mal_df[col].fillna("").astype(str)
        if "approx_mw_kda" in mal_df.columns:
            mal_df["approx_mw_kda"] = _numeric(mal_df["approx_mw_kda"], default=math.nan)

    return df, bins, mal_df


def _save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def generate_plots(bins: pd.DataFrame, outdir: str) -> dict[str, str]:
    """Generate the four figures used in the consolidated PDF."""
    plot_dir = Path(outdir) / PLOT_DIRNAME
    plot_dir.mkdir(parents=True, exist_ok=True)

    bins = bins.sort_values(["bin_start_kda", "bin_end_kda"]).copy()
    plot_bins = bins[bins["bin_start_kda"] <= 250].copy()
    if plot_bins.empty:
        plot_bins = bins.copy()

    paths: dict[str, str] = {}

    # Figure 2: distribution by molecular-weight bin.
    fig, ax = plt.subplots(figsize=(10.8, 4.2))
    ax.bar(plot_bins["bin_start_kda"], plot_bins["protein_count"], width=(plot_bins["bin_end_kda"] - plot_bins["bin_start_kda"]).median() * 0.88)
    ax.set_title("Reference proteome distribution across molecular-weight bins")
    ax.set_xlabel("Theoretical molecular-weight bin (kDa)")
    ax.set_ylabel("Number of proteins")
    ax.grid(axis="y", alpha=0.22)
    top = plot_bins.nlargest(5, "protein_count")
    for _, r in top.iterrows():
        ax.text(r["bin_start_kda"], r["protein_count"] + max(plot_bins["protein_count"].max() * 0.015, 1), f"{int(r['protein_count'])}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    path = plot_dir / "figure_2_molecular_weight_distribution.png"
    _save_figure(fig, path)
    paths["figure_2"] = str(path)

    # Figure 3: crowded neighborhoods.
    crowded = bins.nlargest(20, ["protein_count", "audit_burden_score"]).sort_values("protein_count")
    fig, ax = plt.subplots(figsize=(10.8, 6.5))
    bars = ax.barh(crowded["mw_bin_label"], crowded["protein_count"])
    ax.set_title("Top 20 most crowded molecular-weight bins")
    ax.set_xlabel("Number of proteins")
    ax.set_ylabel("2 kDa molecular-weight bin")
    ax.grid(axis="x", alpha=0.22)
    for bar, value in zip(bars, crowded["protein_count"]):
        ax.text(bar.get_width() + max(crowded["protein_count"].max() * 0.01, 1), bar.get_y() + bar.get_height() / 2, f"{int(value)}", va="center", fontsize=7)
    fig.tight_layout()
    path = plot_dir / "figure_3_most_crowded_bins.png"
    _save_figure(fig, path)
    paths["figure_3"] = str(path)

    # Figure 4: marker and signaling distributions over total protein density.
    x = plot_bins["bin_start_kda"].to_numpy()
    density = plot_bins["protein_count"].to_numpy()
    marker = plot_bins["marker_count"].to_numpy()
    signal = plot_bins["signaling_prefix_count"].to_numpy()

    fig = plt.figure(figsize=(10.8, 6.6))
    ax1 = fig.add_axes([0.08, 0.56, 0.84, 0.34])
    ax1b = ax1.twinx()
    ax1b.bar(x, density, width=1.7, alpha=0.18)
    ax1.plot(x, marker, marker="o", markersize=2.5, linewidth=1.0)
    ax1.set_ylabel("Marker proteins per bin")
    ax1b.set_ylabel("Total proteins per bin")
    ax1.set_title("A. Common loading and compartment-marker proteins", loc="left", fontsize=9)
    ax1.grid(axis="y", alpha=0.18)
    ax1.tick_params(labelbottom=False)

    ax2 = fig.add_axes([0.08, 0.13, 0.84, 0.34])
    ax2b = ax2.twinx()
    ax2b.bar(x, density, width=1.7, alpha=0.18)
    ax2.plot(x, signal, marker="o", markersize=2.5, linewidth=1.0)
    ax2.set_ylabel("Signaling proteins per bin")
    ax2b.set_ylabel("Total proteins per bin")
    ax2.set_xlabel("Theoretical molecular-weight bin (kDa)")
    ax2.set_title("B. Common signaling proteins", loc="left", fontsize=9)
    ax2.grid(axis="y", alpha=0.18)
    fig.suptitle("Audit-relevant annotations across molecular-weight bins", y=0.98, fontsize=11)
    path = plot_dir / "figure_4_marker_signaling_distribution.png"
    _save_figure(fig, path)
    paths["figure_4"] = str(path)

    # Figure 5: audit burden, retaining the priority colors used in the prior report.
    priority_colors = {"routine": "#1f77b4", "elevated": "#ff7f0e", "high": "#2ca02c"}
    colors_for_bins = [priority_colors.get(str(v).lower(), "#7f7f7f") for v in plot_bins["audit_priority"]]
    fig, ax = plt.subplots(figsize=(10.8, 4.5))
    ax.bar(plot_bins["bin_start_kda"], plot_bins["audit_burden_score"], width=1.7, color=colors_for_bins)
    ax.set_title("Audit-burden score across molecular-weight bins")
    ax.set_xlabel("Theoretical molecular-weight bin (kDa)")
    ax.set_ylabel("Audit-burden score")
    ax.grid(axis="y", alpha=0.22)
    handles = [plt.Rectangle((0, 0), 1, 1, color=priority_colors[k]) for k in ["routine", "elevated", "high"]]
    ax.legend(handles, ["Routine", "Elevated", "High"], loc="upper right", frameon=False, ncol=3)
    for _, r in plot_bins.nlargest(8, "audit_burden_score").iterrows():
        ax.text(r["bin_start_kda"], r["audit_burden_score"] + max(plot_bins["audit_burden_score"].max() * 0.015, 0.3), r["mw_bin_label"].replace(" kDa", ""), ha="center", va="bottom", fontsize=6, rotation=0)
    fig.tight_layout()
    path = plot_dir / "figure_5_audit_burden.png"
    _save_figure(fig, path)
    paths["figure_5"] = str(path)

    return paths


def _page_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(colors.HexColor("#777777"))
    canvas.drawString(1.0 * cm, 0.65 * cm, "Western Blot Integrity Atlas - defensive audit report only")
    canvas.drawRightString(A4[0] - 1.0 * cm, 0.65 * cm, f"Page {doc.page}")
    canvas.restoreState()


def _styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="AtlasTitle", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=18, leading=21, alignment=TA_CENTER, spaceAfter=5))
    styles.add(ParagraphStyle(name="AtlasSubtitle", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=12, leading=14, alignment=TA_LEFT, spaceAfter=5))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=12, leading=14, spaceBefore=2, spaceAfter=5))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=8.2, leading=10.2, spaceAfter=4))
    styles.add(ParagraphStyle(name="Warning", parent=styles["BodyText"], fontSize=8.2, leading=10.0, textColor=colors.HexColor("#7a2e00"), spaceAfter=5))
    styles.add(ParagraphStyle(name="Cell", parent=styles["BodyText"], fontSize=5.35, leading=6.45, spaceAfter=0))
    styles.add(ParagraphStyle(name="CellSmall", parent=styles["BodyText"], fontSize=4.85, leading=5.75, spaceAfter=0))
    styles.add(ParagraphStyle(name="CellTiny", parent=styles["BodyText"], fontSize=4.35, leading=5.15, spaceAfter=0))
    styles.add(ParagraphStyle(name="Head", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=5.2, leading=6.0, textColor=colors.white, spaceAfter=0))
    styles.add(ParagraphStyle(name="HeadTiny", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=4.6, leading=5.3, textColor=colors.white, spaceAfter=0))
    styles.add(ParagraphStyle(name="Caption", parent=styles["BodyText"], fontSize=8.2, leading=10.2, spaceAfter=6))
    return styles


def _pcell(value: object, style: ParagraphStyle) -> Paragraph:
    return Paragraph(html_escape(value), style)


def _make_table(
    rows,
    col_widths,
    header_bg="#17324d",
    font_size=5.3,
    repeat_rows=1,
    long=True,
    row_backgrounds=True,
):
    cls = LongTable if long else Table
    table = cls(rows, colWidths=col_widths, repeatRows=repeat_rows, splitByRow=True, hAlign="LEFT")
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("GRID", (0, 0), (-1, -1), 0.22, colors.HexColor("#d0d4d8")),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.2),
        ("TOPPADDING", (0, 0), (-1, -1), 2.1),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.1),
    ]
    if row_backgrounds:
        style_cmds.append(("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7f9fb")]))
    table.setStyle(TableStyle(style_cmds))
    return table


def _image_for_pdf(path: str, max_width: float, max_height: float) -> Image:
    img = Image(path)
    scale = min(max_width / img.imageWidth, max_height / img.imageHeight)
    img.drawWidth = img.imageWidth * scale
    img.drawHeight = img.imageHeight * scale
    img.hAlign = "CENTER"
    return img


def _metric_rows(
    df: pd.DataFrame,
    bins: pd.DataFrame,
    fasta_path: str,
    bin_kda: float,
) -> list[list[str]]:
    total = max(len(df), 1)
    le50 = int((df["mw_kda"] <= 50).sum())
    mid = int(((df["mw_kda"] > 50) & (df["mw_kda"] <= 100)).sum())
    gt100 = int((df["mw_kda"] > 100).sum())
    priorities = bins["audit_priority"].astype(str).str.lower().value_counts()
    return [
        ["Metric", "Value"],
        ["Proteins parsed", f"{len(df):,}"],
        ["Molecular-weight bins", f"{len(bins):,}"],
        ["Bin width", f"{bin_kda:g} kDa"],
        ["MW range", f"{df['mw_kda'].min():.2f} to {df['mw_kda'].max():.2f} kDa"],
        ["Median MW", f"{df['mw_kda'].median():.2f} kDa"],
        ["Mean MW", f"{df['mw_kda'].mean():.2f} kDa"],
        ["Proteins <= 50 kDa", f"{le50:,} ({100.0 * le50 / total:.1f}%)"],
        ["Proteins 50-100 kDa", f"{mid:,} ({100.0 * mid / total:.1f}%)"],
        ["Proteins > 100 kDa", f"{gt100:,} ({100.0 * gt100 / total:.1f}%)"],
        ["High-priority bins", f"{int(priorities.get('high', 0)):,}"],
        ["Elevated-priority bins", f"{int(priorities.get('elevated', 0)):,}"],
        ["Routine-priority bins", f"{int(priorities.get('routine', 0)):,}"],
        ["Source FASTA", os.path.basename(fasta_path)],
    ]


def build_full_pdf(
    df: pd.DataFrame,
    bins: pd.DataFrame,
    mal_df: pd.DataFrame,
    outdir: str,
    species_name: str,
    fasta_path: str,
    marker_file: str,
    signaling_file: str,
    family_file: str,
    bin_kda: float,
    plot_paths: dict[str, str],
) -> str:
    """Build the consolidated portrait-A4 report with figures and all tables."""
    df, bins, mal_df = normalize_report_dataframes(df, bins, mal_df)
    out_path = os.path.join(outdir, REPORT_FILENAME)
    doc = SimpleDocTemplate(
        out_path,
        pagesize=A4,
        rightMargin=0.8 * cm,
        leftMargin=0.8 * cm,
        topMargin=0.75 * cm,
        bottomMargin=1.0 * cm,
        title="Western Blot Integrity Atlas",
        author="Western Blot Integrity Atlas generator",
    )
    styles = _styles()
    story = []

    # Cover/summary page.
    story.append(Paragraph("Western Blot Integrity Atlas", styles["AtlasTitle"]))
    story.append(Paragraph("Full reference-proteome molecular-weight audit report with embedded figures and complete tables", styles["AtlasSubtitle"]))
    story.append(
        Paragraph(
            "Defensive use only. This report helps evaluate whether a claimed western blot result is plausible and what controls are needed. It must not be used to plan, enable, or disguise blot substitution or fabrication.",
            styles["Warning"],
        )
    )
    metric_rows = [[_pcell(c, styles["Head"]) for c in _metric_rows(df, bins, fasta_path, bin_kda)[0]]]
    for r in _metric_rows(df, bins, fasta_path, bin_kda)[1:]:
        metric_rows.append([_pcell(r[0], styles["BodySmall"]), _pcell(r[1], styles["BodySmall"])])
    story.append(_make_table(metric_rows, [6.6 * cm, 12.0 * cm], header_bg="#17324d", font_size=7.2, long=False))
    story.append(Spacer(1, 0.18 * cm))
    story.append(Paragraph("Method in one breath", styles["Section"]))
    story.append(
        Paragraph(
            "For each protein sequence, the script computes theoretical average molecular weight from residue masses, assigns a molecular-weight bin, applies the external marker, signaling-prefix, and family-pattern annotations, maps supplied public malpractice cases to bins, and calculates an audit-burden priority. Similar apparent molecular weight is never treated as proof of identity; it only increases the need for full-blot provenance and orthogonal controls.",
            styles["BodySmall"],
        )
    )

    # Full figures.
    figure_specs = [
        (
            "Molecular-weight distribution of the reference proteome.",
            "Proteins were grouped into theoretical molecular-weight bins. Dense gel regions are especially relevant for western blot review because many unrelated proteins occupy nearby molecular-weight neighborhoods.",
            "figure_2",
            17.8 * cm,
            16.5 * cm,
        ),
        (
            "Most crowded molecular-weight neighborhoods.",
            "The top bins are ranked by protein count. These bins are not suspicious by themselves; they indicate where molecular weight alone is a weak identity criterion and where validation documentation becomes more important.",
            "figure_3",
            17.8 * cm,
            19.5 * cm,
        ),
        (
            "Distribution of common western blot markers and signaling proteins across molecular-weight bins.",
            "The pale background bars show total protein density, while the overlaid lines show audit-relevant annotations. Panel A shows configured loading or compartment markers; Panel B shows configured signaling-prefix matches.",
            "figure_4",
            17.8 * cm,
            18.5 * cm,
        ),
        (
            "Audit-burden score across molecular-weight bins.",
            "The score integrates protein density, marker count, signaling-prefix count, and configured family/context information. High-scoring bins define high-documentation zones for western blot identity claims.",
            "figure_5",
            17.8 * cm,
            17.5 * cm,
        ),
    ]
    for title, caption, key, w, h in figure_specs:
        story.append(PageBreak())
        story.append(Paragraph(title, styles["Section"]))
        story.append(Paragraph(caption, styles["Caption"]))
        story.append(_image_for_pdf(plot_paths[key], w, h))

    # Compact score table retained from the plot report.
    story.append(PageBreak())
    story.append(Paragraph("Highest-audit molecular-weight neighborhoods: score summary", styles["Section"]))
    top_score = bins.sort_values(["audit_burden_score", "protein_count"], ascending=False).head(25)
    headers = ["MW bin", "N", "Markers", "Signal", "Priority", "Score"]
    rows = [[_pcell(x, styles["Head"]) for x in headers]]
    for _, r in top_score.iterrows():
        rows.append([
            _pcell(r["mw_bin_label"], styles["Cell"]),
            _pcell(int(r["protein_count"]), styles["Cell"]),
            _pcell(int(r["marker_count"]), styles["Cell"]),
            _pcell(int(r["signaling_prefix_count"]), styles["Cell"]),
            _pcell(r["audit_priority"], styles["Cell"]),
            _pcell(f"{r['audit_burden_score']:.1f}", styles["Cell"]),
        ])
    story.append(_make_table(rows, [3.1 * cm, 2.0 * cm, 2.3 * cm, 2.3 * cm, 3.2 * cm, 2.2 * cm], header_bg="#17324d", font_size=6.6))

    # Detailed highest-audit table from the original full tabular report.
    story.append(PageBreak())
    story.append(Paragraph("Highest-audit molecular-weight neighborhoods", styles["Section"]))
    story.append(
        Paragraph(
            "These bins combine protein density with configured markers, signaling-prefix matches, and family-rich regions. They are not substitution recommendations; they are places where reviewers should be most insistent about source files and controls.",
            styles["BodySmall"],
        )
    )
    top_bins = bins.sort_values(["audit_burden_score", "protein_count"], ascending=False).head(30)
    headers = ["MW bin", "N", "Markers", "Signal", "Families", "Priority", "Examples / audit families"]
    rows = [[_pcell(x, styles["HeadTiny"]) for x in headers]]
    for _, r in top_bins.iterrows():
        rows.append([
            _pcell(r["mw_bin_label"], styles["CellTiny"]),
            _pcell(int(r["protein_count"]), styles["CellTiny"]),
            _pcell(int(r["marker_count"]), styles["CellTiny"]),
            _pcell(int(r["signaling_prefix_count"]), styles["CellTiny"]),
            _pcell(r["families"], styles["CellTiny"]),
            _pcell(r["audit_priority"], styles["CellTiny"]),
            _pcell(r["examples_audit_families"], styles["CellTiny"]),
        ])
    story.append(_make_table(rows, [2.1 * cm, 0.65 * cm, 0.8 * cm, 0.75 * cm, 4.6 * cm, 1.2 * cm, 8.4 * cm], header_bg="#17324d", font_size=4.5))

    # Common marker table retained from the earlier full report.
    story.append(PageBreak())
    story.append(Paragraph("Common loading/compartment markers present in this proteome file", styles["Section"]))
    markers = df[df.get("is_marker", False) == True].sort_values(["mw_kda", "gene_primary", "accession"]).copy()  # noqa: E712
    headers = ["MW kDa", "Gene", "Accession", "Entry", "Protein", "MW bin"]
    rows = [[_pcell(x, styles["Head"]) for x in headers]]
    for _, r in markers.iterrows():
        rows.append([
            _pcell(f"{r['mw_kda']:.2f}", styles["CellSmall"]),
            _pcell(r["gene_primary"], styles["CellSmall"]),
            _pcell(r["accession"], styles["CellSmall"]),
            _pcell(r["entry_name"], styles["CellSmall"]),
            _pcell(r["protein_name"], styles["CellSmall"]),
            _pcell(r["mw_bin_label"], styles["CellSmall"]),
        ])
    if len(rows) == 1:
        rows.append([_pcell("No configured marker genes were present in the parsed proteome.", styles["Cell"])] + [""] * 5)
    story.append(_make_table(rows, [1.3 * cm, 1.8 * cm, 2.1 * cm, 2.4 * cm, 7.9 * cm, 2.7 * cm], header_bg="#17324d", font_size=4.9))

    # Complete bin table, all rows.
    story.append(PageBreak())
    story.append(Paragraph("Complete molecular-weight bin summary, low to high", styles["Section"]))
    story.append(Paragraph("All molecular-weight bins are shown. Long Families and Examples / audit families cells wrap and continue across pages.", styles["BodySmall"]))
    headers = ["MW bin", "N", "Markers", "Signal", "Families", "Priority", "Examples / audit families"]
    rows = [[_pcell(x, styles["HeadTiny"]) for x in headers]]
    for _, r in bins.sort_values(["bin_start_kda", "bin_end_kda"]).iterrows():
        rows.append([
            _pcell(r["mw_bin_label"], styles["CellTiny"]),
            _pcell(int(r["protein_count"]), styles["CellTiny"]),
            _pcell(int(r["marker_count"]), styles["CellTiny"]),
            _pcell(int(r["signaling_prefix_count"]), styles["CellTiny"]),
            _pcell(r["families"], styles["CellTiny"]),
            _pcell(r["audit_priority"], styles["CellTiny"]),
            _pcell(r["examples_audit_families"], styles["CellTiny"]),
        ])
    story.append(_make_table(rows, [2.05 * cm, 0.62 * cm, 0.78 * cm, 0.72 * cm, 4.65 * cm, 1.15 * cm, 8.35 * cm], header_bg="#244c2a", font_size=4.4))

    # Protein table from original report, first 120 records.
    story.append(PageBreak())
    story.append(Paragraph("Proteins sorted from lowest to highest molecular weight", styles["Section"]))
    story.append(Paragraph("The PDF shows the first 120 proteins. The complete table is in proteome_wb_integrity_atlas.csv.", styles["BodySmall"]))
    headers = ["MW kDa", "Gene", "Accession", "Protein", "MW bin", "Families", "Audit note"]
    rows = [[_pcell(x, styles["HeadTiny"]) for x in headers]]
    for _, r in df.sort_values(["mw_kda", "gene_primary", "accession"]).head(120).iterrows():
        rows.append([
            _pcell(f"{r['mw_kda']:.2f}", styles["CellTiny"]),
            _pcell(r["gene_primary"], styles["CellTiny"]),
            _pcell(r["accession"], styles["CellTiny"]),
            _pcell(r["protein_name"], styles["CellTiny"]),
            _pcell(r["mw_bin_label"], styles["CellTiny"]),
            _pcell(r["families"], styles["CellTiny"]),
            _pcell(r["defensive_note"], styles["CellTiny"]),
        ])
    story.append(_make_table(rows, [1.2 * cm, 1.65 * cm, 1.85 * cm, 5.2 * cm, 2.15 * cm, 2.8 * cm, 4.1 * cm], header_bg="#4b2c62", font_size=4.35))

    # Malpractice table from original report.
    if not mal_df.empty:
        story.append(PageBreak())
        story.append(Paragraph("Reported western-blot malpractice cases mapped to MW bins", styles["Section"]))
        story.append(
            Paragraph(
                "These are reported public cases supplied in the malpractice input TSV. They are evidence examples for integrity-audit results sections, not substitution recommendations.",
                styles["Warning"],
            )
        )
        headers = ["Case", "Claim/protein", "MW bin", "Reported issue", "Results sentence"]
        rows = [[_pcell(x, styles["HeadTiny"]) for x in headers]]
        for _, r in mal_df.iterrows():
            rows.append([
                _pcell(r.get("case_id", ""), styles["CellSmall"]),
                _pcell(r.get("reported_protein_claim", ""), styles["CellSmall"]),
                _pcell(r.get("mw_bin_label", ""), styles["CellSmall"]),
                _pcell(r.get("reported_issue", ""), styles["CellSmall"]),
                _pcell(r.get("results_sentence", ""), styles["CellSmall"]),
            ])
        story.append(_make_table(rows, [2.8 * cm, 3.0 * cm, 2.1 * cm, 5.3 * cm, 5.7 * cm], header_bg="#7a2e00", font_size=4.8))

    # Final checklist and generated-file inventory.
    story.append(PageBreak())
    story.append(Paragraph("Integrity-audit checklist", styles["Section"]))
    checklist = [
        "Ask for the full, uncropped blot with ladder and lane labels. Cropped panels are not enough.",
        "Confirm expected molecular weight and apparent molecular weight. Do not treat molecular-weight agreement as proof of identity.",
        "For configured marker proteins, request source-imager metadata and original antibody-validation records.",
        "Check whether a marker band is near common markers such as GAPDH, actin, tubulin, HSPs, lamins, VDAC, PCNA, or other configured controls.",
        "For signaling proteins, verify phospho/total controls, controls for activation state, and an independent identity control when central to the claim.",
        "In crowded or high-audit bins, prioritize raw-data provenance, full blot context, and orthogonal validation.",
        "Interpret mapped malpractice cases only as documented audit precedents. They do not imply that proteins in the same bin are interchangeable.",
    ]
    for item in checklist:
        story.append(Paragraph("- " + html_escape(item), styles["BodySmall"]))

    story.append(Spacer(1, 0.2 * cm))
    story.append(Paragraph("Generated files", styles["Section"]))
    generated = [
        "proteome_wb_integrity_atlas.csv - complete per-protein atlas sorted by theoretical MW.",
        "mw_bin_summary.csv - complete molecular-weight-bin audit summary.",
        "reported_malpractice_bins.csv - generated when a malpractice TSV is supplied.",
        f"{PLOT_DIRNAME}/ - four standalone PNG figures embedded in this report.",
        f"{REPORT_FILENAME} - this consolidated report.",
        f"{COMPATIBILITY_REPORT_FILENAME} - byte-identical compatibility copy of the consolidated report.",
    ]
    for item in generated:
        story.append(Paragraph("- " + html_escape(item), styles["BodySmall"]))

    story.append(Spacer(1, 0.2 * cm))
    story.append(Paragraph("Input configuration", styles["Section"]))
    config_rows = [
        ["Input", "File"],
        ["Proteome FASTA", os.path.basename(fasta_path)],
        ["Marker file", os.path.basename(marker_file)],
        ["Signaling-prefix file", os.path.basename(signaling_file)],
        ["Family-pattern file", os.path.basename(family_file)],
        ["Species/report label", species_name],
    ]
    rows = [[_pcell(x, styles["Head"]) for x in config_rows[0]]]
    rows.extend([[_pcell(r[0], styles["BodySmall"]), _pcell(r[1], styles["BodySmall"])] for r in config_rows[1:]])
    story.append(_make_table(rows, [5.2 * cm, 12.8 * cm], header_bg="#555555", font_size=7.0, long=False))

    doc.build(story, onFirstPage=_page_footer, onLaterPages=_page_footer)

    compatibility_path = os.path.join(outdir, COMPATIBILITY_REPORT_FILENAME)
    shutil.copyfile(out_path, compatibility_path)
    return out_path


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.bin_kda <= 0:
        raise SystemExit("--bin-kda must be greater than zero")

    os.makedirs(args.outdir, exist_ok=True)

    marker_df = load_marker_file(args.marker_file)
    signaling_df = load_signaling_prefix_file(args.signaling_prefix_file)
    family_df = load_family_pattern_file(args.family_pattern_file)
    mal_df = load_malpractice_file(args.malpractice_file)

    df = proteome_dataframe_from_fasta(args.proteome_fasta)
    df, bins = annotate(df, marker_df, signaling_df, family_df, args.bin_kda)
    mal_df = attach_malpractice_bins(mal_df, args.bin_kda)

    paths = write_outputs(df, bins, mal_df, args.outdir)
    plot_paths = generate_plots(bins, args.outdir)
    report_path = build_full_pdf(
        df=df,
        bins=bins,
        mal_df=mal_df,
        outdir=args.outdir,
        species_name=args.species_name,
        fasta_path=args.proteome_fasta,
        marker_file=args.marker_file,
        signaling_file=args.signaling_prefix_file,
        family_file=args.family_pattern_file,
        bin_kda=args.bin_kda,
        plot_paths=plot_paths,
    )

    paths.update(plot_paths)
    paths["full_pdf"] = report_path
    paths["compatibility_pdf"] = os.path.join(args.outdir, COMPATIBILITY_REPORT_FILENAME)

    print("Generated:")
    for key, value in paths.items():
        print(f"  {key}: {value}")
    print(f"Proteins parsed: {len(df):,}; MW bins: {len(bins):,}")
    print("No --bin_csv input was used; all tables and figures were generated internally.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
