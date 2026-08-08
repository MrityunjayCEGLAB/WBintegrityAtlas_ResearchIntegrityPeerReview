#!/usr/bin/env python3
"""
Western Blot Integrity Atlas Generator
======================================

DEFENSIVE USE ONLY.
This script clusters proteins by theoretical molecular weight to help reviewers,
editors, research-integrity officers, and lab members audit western blot figures.
It is NOT a guide for substituting one protein blot for another.

Requested updates implemented:
1. Marker genes, signalling prefixes, and family patterns are no longer hard-coded.
   They are supplied as separate TSV files at runtime.
2. The same script can be used for human, cow, mouse, chicken, or any proteome FASTA.
3. PDF tables wrap long text in the "Examples / audit families" and "Families" columns.
   Dynamic row heights are used through ReportLab Paragraph objects.
4. Bin-summary tables use wider columns and landscape A4 pages to avoid overlap.
5. Optional malpractice case files can be supplied to map reported cases onto MW bins.

Typical command:
    python wb_integrity_atlas_externalized.py \
      --proteome-fasta human_UP000005640.fasta \
      --marker-file configs/human_marker_genes.tsv \
      --signaling-prefix-file configs/human_signaling_prefixes.tsv \
      --family-pattern-file configs/human_family_patterns.tsv \
      --malpractice-file malpractice/human_reported_malpractice_cases.tsv \
      --species-name "Homo sapiens" \
      --outdir human_wb_atlas

Dependencies:
    pip install pandas reportlab
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
import textwrap
from dataclasses import dataclass
from typing import Iterable, Optional

try:
    import pandas as pd
except Exception as exc:
    raise SystemExit("This script requires pandas: pip install pandas") from exc

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import (
        SimpleDocTemplate,
        Paragraph,
        Spacer,
        Table,
        TableStyle,
        PageBreak,
    )
except Exception as exc:
    raise SystemExit("This script requires reportlab: pip install reportlab") from exc


# Average amino-acid residue masses in Da, residue form inside a peptide chain.
# Add water mass once per complete protein.
AA_AVG_MASS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
    "U": 150.0388,
    "O": 237.3018,
}
WATER_MASS = 18.01528


@dataclass
class FastaProtein:
    accession: str
    entry_name: str
    gene_primary: str
    protein_name: str
    sequence: str


def safe_str(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def html_escape(text: object) -> str:
    s = safe_str(text)
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("\n", "<br/>")
    )


def short(text: object, n: int = 72) -> str:
    s = safe_str(text)
    if len(s) <= n:
        return s
    return s[: max(0, n - 1)] + "…"


def protein_mass_da(seq: str) -> float:
    seq = re.sub(r"[^A-Z]", "", seq.upper())
    total = WATER_MASS
    for aa in seq:
        total += AA_AVG_MASS.get(aa, 0.0)
    return total


def parse_uniprot_header(header: str) -> tuple[str, str, str, str]:
    """Parse a UniProt-style FASTA header.

    Handles headers such as:
    >sp|P04406|G3P_HUMAN Glyceraldehyde-3-phosphate dehydrogenase OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3
    """
    h = header.lstrip(">").strip()
    accession = ""
    entry_name = ""
    rest = h

    parts = h.split("|", 2)
    if len(parts) == 3 and parts[0] in {"sp", "tr", "up"}:
        accession = parts[1].strip()
        rest = parts[2].strip()
        toks = rest.split(None, 1)
        entry_name = toks[0].strip() if toks else ""
        desc = toks[1].strip() if len(toks) > 1 else ""
    else:
        toks = h.split(None, 1)
        accession = toks[0].strip()
        desc = toks[1].strip() if len(toks) > 1 else ""

    gene = ""
    m = re.search(r"\bGN=([^\s]+)", h)
    if m:
        gene = m.group(1).strip()

    # Protein name is usually description before OS=.
    protein_name = re.split(r"\sOS=", desc)[0].strip() if desc else ""
    if not protein_name:
        protein_name = entry_name or accession
    return accession, entry_name, gene, protein_name


def read_fasta(path: str) -> list[FastaProtein]:
    proteins: list[FastaProtein] = []
    header: Optional[str] = None
    seq_parts: list[str] = []

    def flush():
        nonlocal header, seq_parts
        if header is None:
            return
        accession, entry_name, gene, protein_name = parse_uniprot_header(header)
        sequence = "".join(seq_parts).strip().upper()
        if sequence:
            proteins.append(FastaProtein(accession, entry_name, gene, protein_name, sequence))
        header = None
        seq_parts = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()
    return proteins


def proteome_dataframe_from_fasta(path: str) -> pd.DataFrame:
    proteins = read_fasta(path)
    rows = []
    for p in proteins:
        mass = protein_mass_da(p.sequence)
        rows.append(
            {
                "accession": p.accession,
                "entry_name": p.entry_name,
                "gene_primary": p.gene_primary,
                "protein_name": p.protein_name,
                "length_aa": len(p.sequence),
                "mass_da": mass,
                "mw_kda": mass / 1000.0,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"No protein sequences were parsed from FASTA: {path}")
    df = df.sort_values(["mw_kda", "gene_primary", "accession"]).reset_index(drop=True)
    return df


def read_tsv_required(path: str, required_columns: list[str]) -> pd.DataFrame:
    if not path:
        raise ValueError("A required TSV file path was not supplied.")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    return df


def load_marker_file(path: str) -> pd.DataFrame:
    return read_tsv_required(path, ["gene_symbol", "category", "note"])


def load_signaling_prefix_file(path: str) -> pd.DataFrame:
    return read_tsv_required(path, ["prefix", "category", "note"])


def load_family_pattern_file(path: str) -> pd.DataFrame:
    return read_tsv_required(path, ["family", "regex", "category", "note"])


def load_malpractice_file(path: Optional[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    df = read_tsv_required(
        path,
        [
            "case_id",
            "paper_or_case_title",
            "reported_protein_claim",
            "reported_issue",
            "approx_mw_kda",
            "source_url",
            "results_sentence",
        ],
    )
    df["approx_mw_kda"] = pd.to_numeric(df["approx_mw_kda"], errors="coerce")
    return df


def compile_regex(pattern: str):
    return re.compile(pattern, flags=re.IGNORECASE)


def annotate(
    df: pd.DataFrame,
    marker_df: pd.DataFrame,
    signaling_df: pd.DataFrame,
    family_df: pd.DataFrame,
    bin_kda: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df["gene_upper"] = df["gene_primary"].fillna("").astype(str).str.upper()
    marker_df = marker_df.copy()
    marker_df["gene_upper"] = marker_df["gene_symbol"].astype(str).str.upper()

    marker_map = marker_df.groupby("gene_upper").agg(
        marker_categories=("category", lambda x: "; ".join(sorted(set([v for v in x if v])))),
        marker_notes=("note", lambda x: "; ".join(sorted(set([v for v in x if v])))),
    )
    df = df.merge(marker_map, left_on="gene_upper", right_index=True, how="left")
    df["is_marker"] = df["marker_categories"].notna()
    df["marker_categories"] = df["marker_categories"].fillna("")
    df["marker_notes"] = df["marker_notes"].fillna("")

    prefixes = []
    for _, row in signaling_df.iterrows():
        prefix = safe_str(row["prefix"]).strip()
        if prefix:
            prefixes.append((prefix.upper(), safe_str(row["category"]), safe_str(row["note"])))

    def match_prefixes(gene: str) -> tuple[str, str]:
        g = safe_str(gene).upper()
        cats = []
        notes = []
        for prefix, cat, note in prefixes:
            if g.startswith(prefix):
                if cat:
                    cats.append(cat)
                if note:
                    notes.append(note)
        return "; ".join(sorted(set(cats))), "; ".join(sorted(set(notes)))

    prefix_matches = df["gene_primary"].apply(match_prefixes)
    df["signaling_categories"] = [x[0] for x in prefix_matches]
    df["signaling_notes"] = [x[1] for x in prefix_matches]
    df["is_signaling_prefix"] = df["signaling_categories"].astype(bool)

    regexes = []
    for _, row in family_df.iterrows():
        fam = safe_str(row["family"]).strip()
        regex = safe_str(row["regex"]).strip()
        cat = safe_str(row["category"]).strip()
        note = safe_str(row["note"]).strip()
        if fam and regex:
            regexes.append((fam, compile_regex(regex), cat, note))

    def match_families(row: pd.Series) -> tuple[str, str, str]:
        gene = safe_str(row["gene_primary"])
        protein = safe_str(row["protein_name"])
        text = f"{gene} {protein}"
        fams = []
        cats = []
        notes = []
        for fam, rgx, cat, note in regexes:
            if rgx.search(text):
                fams.append(fam)
                if cat:
                    cats.append(cat)
                if note:
                    notes.append(note)
        return (
            "; ".join(sorted(set(fams))),
            "; ".join(sorted(set(cats))),
            "; ".join(sorted(set(notes))),
        )

    family_matches = df.apply(match_families, axis=1)
    df["families"] = [x[0] for x in family_matches]
    df["family_categories"] = [x[1] for x in family_matches]
    df["family_notes"] = [x[2] for x in family_matches]

    df["mw_bin_start_kda"] = (df["mw_kda"].apply(lambda x: math.floor(x / bin_kda) * bin_kda)).round(3)
    df["mw_bin_end_kda"] = (df["mw_bin_start_kda"] + bin_kda).round(3)
    df["mw_bin_label"] = df.apply(lambda r: f"{r.mw_bin_start_kda:.1f}-{r.mw_bin_end_kda:.1f} kDa", axis=1)

    bin_rows = []
    for label, g in df.groupby("mw_bin_label", sort=False):
        examples = []
        for _, r in g.head(18).iterrows():
            fam = safe_str(r["families"])
            tag = f" [{fam}]" if fam else ""
            examples.append(f"{safe_str(r['gene_primary'])}{tag}")
        families = sorted(set("; ".join(g["families"].dropna().astype(str)).split("; ")) - {""})
        marker_genes = sorted(g.loc[g["is_marker"], "gene_primary"].dropna().astype(str).unique().tolist())
        signal_genes = sorted(g.loc[g["is_signaling_prefix"], "gene_primary"].dropna().astype(str).unique().tolist())
        protein_count = len(g)
        marker_count = int(g["is_marker"].sum())
        signaling_count = int(g["is_signaling_prefix"].sum())
        family_count = len(families)
        audit_burden_score = round(
            min(protein_count, 400) * 0.05
            + marker_count * 5.0
            + signaling_count * 1.5
            + family_count * 0.75,
            2,
        )
        bin_rows.append(
            {
                "mw_bin_label": label,
                "bin_start_kda": g["mw_bin_start_kda"].iloc[0],
                "bin_end_kda": g["mw_bin_end_kda"].iloc[0],
                "protein_count": protein_count,
                "marker_count": marker_count,
                "signaling_prefix_count": signaling_count,
                "family_count": family_count,
                "families": "; ".join(families),
                "marker_genes": ", ".join(marker_genes),
                "signaling_genes": ", ".join(signal_genes),
                "examples_audit_families": ", ".join(examples),
                "audit_burden_score": audit_burden_score,
            }
        )
    bins = pd.DataFrame(bin_rows).sort_values(["bin_start_kda", "bin_end_kda"]).reset_index(drop=True)
    q75 = bins["audit_burden_score"].quantile(0.75) if len(bins) else 0
    q90 = bins["audit_burden_score"].quantile(0.90) if len(bins) else 0

    def priority(score: float) -> str:
        if score >= q90:
            return "high"
        if score >= q75:
            return "elevated"
        return "routine"

    bins["audit_priority"] = bins["audit_burden_score"].apply(priority)
    df = df.merge(
        bins[
            [
                "mw_bin_label",
                "protein_count",
                "marker_count",
                "signaling_prefix_count",
                "family_count",
                "audit_burden_score",
                "audit_priority",
            ]
        ],
        on="mw_bin_label",
        how="left",
    )
    df["defensive_note"] = df.apply(make_defensive_note, axis=1)
    return df, bins


def make_defensive_note(row: pd.Series) -> str:
    notes = []
    if row.get("is_marker", False):
        notes.append("configured marker/control; verify source-file provenance if used near target bands")
    if row.get("is_signaling_prefix", False):
        notes.append("configured signaling gene; check phospho/total controls and antibody validation")
    if safe_str(row.get("families")):
        notes.append("configured audit family match: " + safe_str(row.get("families")))
    if row.get("audit_priority") == "high":
        notes.append("crowded MW neighborhood; require full blot, ladder, raw metadata, and controls")
    return "; ".join(notes) if notes else "standard full-blot and antibody-validation checks"


def attach_malpractice_bins(mal_df: pd.DataFrame, bin_kda: float) -> pd.DataFrame:
    if mal_df.empty:
        return mal_df
    out = mal_df.copy()
    out["mw_bin_start_kda"] = out["approx_mw_kda"].apply(lambda x: math.floor(x / bin_kda) * bin_kda if pd.notna(x) else math.nan)
    out["mw_bin_end_kda"] = out["mw_bin_start_kda"] + bin_kda
    out["mw_bin_label"] = out.apply(
        lambda r: f"{r.mw_bin_start_kda:.1f}-{r.mw_bin_end_kda:.1f} kDa" if pd.notna(r["mw_bin_start_kda"]) else "NA",
        axis=1,
    )
    return out


def write_outputs(df: pd.DataFrame, bins: pd.DataFrame, mal_df: pd.DataFrame, outdir: str) -> dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    paths = {
        "protein_csv": os.path.join(outdir, "proteome_wb_integrity_atlas.csv"),
        "bin_csv": os.path.join(outdir, "mw_bin_summary.csv"),
    }
    df.to_csv(paths["protein_csv"], index=False)
    bins.to_csv(paths["bin_csv"], index=False)
    if not mal_df.empty:
        paths["malpractice_csv"] = os.path.join(outdir, "reported_malpractice_bins.csv")
        mal_df.to_csv(paths["malpractice_csv"], index=False)
    return paths


def page_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.grey)
    canvas.drawString(1.2 * cm, 0.8 * cm, "Western Blot Integrity Atlas - defensive audit report only")
    canvas.drawRightString(landscape(A4)[0] - 1.2 * cm, 0.8 * cm, f"Page {doc.page}")
    canvas.restoreState()


def pcell(text: object, style: ParagraphStyle) -> Paragraph:
    return Paragraph(html_escape(text), style)


def make_table(rows, col_widths, header_bg="#17324d", font_size=7):
    table = Table(rows, colWidths=col_widths, repeatRows=1, splitByRow=True)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                ("FONTSIZE", (0, 0), (-1, -1), font_size),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7f9fb")]),
            ]
        )
    )
    return table


def build_pdf(
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
) -> str:
    os.makedirs(outdir, exist_ok=True)
    pdf_path = os.path.join(outdir, "western_blot_integrity_atlas_report.pdf")
    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=landscape(A4),
        rightMargin=1.0 * cm,
        leftMargin=1.0 * cm,
        topMargin=1.0 * cm,
        bottomMargin=1.2 * cm,
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=8, leading=10))
    styles.add(ParagraphStyle(name="Tiny", parent=styles["BodyText"], fontSize=6.5, leading=8))
    styles.add(ParagraphStyle(name="Warn", parent=styles["BodyText"], fontSize=9, leading=11, textColor=colors.HexColor("#7a2e00")))
    styles.add(ParagraphStyle(name="TableCell", parent=styles["BodyText"], fontSize=6.5, leading=8))
    styles.add(ParagraphStyle(name="TableHead", parent=styles["BodyText"], fontSize=7, leading=8, textColor=colors.white))

    story = []
    story.append(Paragraph("Western Blot Integrity Atlas", styles["Title"]))
    story.append(Paragraph(species_name, styles["Heading2"]))
    story.append(
        Paragraph(
            "This report is a defensive audit aid. It groups proteins by theoretical molecular weight to identify crowded gel regions where full blots, molecular-weight ladders, antibody validation, knockdown/knockout controls, and raw-image provenance are especially important. It must not be used to plan blot substitution or fabrication.",
            styles["Warn"],
        )
    )
    story.append(Spacer(1, 0.25 * cm))

    meta_rows = [
        ["Species/report label", species_name],
        ["Proteome FASTA", os.path.basename(fasta_path)],
        ["Proteins parsed", f"{len(df):,}"],
        ["MW bin width", f"{bin_kda:g} kDa"],
        ["Lowest MW", f"{df['mw_kda'].min():.2f} kDa"],
        ["Highest MW", f"{df['mw_kda'].max():.2f} kDa"],
        ["Marker file", os.path.basename(marker_file)],
        ["Signaling-prefix file", os.path.basename(signaling_file)],
        ["Family-pattern file", os.path.basename(family_file)],
    ]
    table = make_table(meta_rows, [5.0 * cm, 20.0 * cm], header_bg="#555555")
    story.append(table)
    story.append(PageBreak())

    story.append(Paragraph("Highest-audit molecular-weight neighborhoods", styles["Heading2"]))
    story.append(
        Paragraph(
            "This table now wraps long text in the Examples / audit families column. Row heights expand dynamically so no gene-family text is clipped.",
            styles["Small"],
        )
    )
    top_bins = bins.sort_values(["audit_burden_score", "protein_count"], ascending=False).head(30)
    rows = [[pcell(x, styles["TableHead"]) for x in ["MW bin", "N", "Markers", "Signal", "Families", "Priority", "Examples / audit families"]]]
    for _, r in top_bins.iterrows():
        rows.append(
            [
                pcell(r["mw_bin_label"], styles["TableCell"]),
                pcell(int(r["protein_count"]), styles["TableCell"]),
                pcell(int(r["marker_count"]), styles["TableCell"]),
                pcell(int(r["signaling_prefix_count"]), styles["TableCell"]),
                pcell(r["families"], styles["TableCell"]),
                pcell(r["audit_priority"], styles["TableCell"]),
                pcell(r["examples_audit_families"], styles["TableCell"]),
            ]
        )
    story.append(make_table(rows, [2.6 * cm, 1.2 * cm, 1.6 * cm, 1.6 * cm, 4.4 * cm, 2.0 * cm, 14.0 * cm]))
    story.append(PageBreak())

    story.append(Paragraph("Complete molecular-weight bin summary, low to high", styles["Heading2"]))
    story.append(
        Paragraph(
            "The Families and Examples / audit families columns are Paragraph-wrapped with wider columns and dynamic row heights to avoid overlap.",
            styles["Small"],
        )
    )
    rows = [[pcell(x, styles["TableHead"]) for x in ["MW bin", "N", "Markers", "Signal", "Families", "Priority", "Examples / audit families"]]]
    for _, r in bins.sort_values(["bin_start_kda", "bin_end_kda"]).iterrows():
        rows.append(
            [
                pcell(r["mw_bin_label"], styles["TableCell"]),
                pcell(int(r["protein_count"]), styles["TableCell"]),
                pcell(int(r["marker_count"]), styles["TableCell"]),
                pcell(int(r["signaling_prefix_count"]), styles["TableCell"]),
                pcell(r["families"], styles["TableCell"]),
                pcell(r["audit_priority"], styles["TableCell"]),
                pcell(r["examples_audit_families"], styles["TableCell"]),
            ]
        )
    story.append(make_table(rows, [2.5 * cm, 1.1 * cm, 1.5 * cm, 1.5 * cm, 4.7 * cm, 1.8 * cm, 14.0 * cm], header_bg="#244c2a"))
    story.append(PageBreak())

    story.append(Paragraph("Proteins sorted from lowest to highest molecular weight", styles["Heading2"]))
    story.append(Paragraph("The PDF shows the first 120 proteins. The complete table is in proteome_wb_integrity_atlas.csv.", styles["Small"]))
    rows = [[pcell(x, styles["TableHead"]) for x in ["MW kDa", "Gene", "Accession", "Protein", "MW bin", "Families", "Audit note"]]]
    for _, r in df.sort_values("mw_kda").head(120).iterrows():
        rows.append(
            [
                pcell(f"{r['mw_kda']:.2f}", styles["TableCell"]),
                pcell(r["gene_primary"], styles["TableCell"]),
                pcell(r["accession"], styles["TableCell"]),
                pcell(r["protein_name"], styles["TableCell"]),
                pcell(r["mw_bin_label"], styles["TableCell"]),
                pcell(r["families"], styles["TableCell"]),
                pcell(r["defensive_note"], styles["TableCell"]),
            ]
        )
    story.append(make_table(rows, [1.7 * cm, 2.0 * cm, 2.5 * cm, 6.5 * cm, 2.4 * cm, 4.5 * cm, 8.3 * cm], header_bg="#4b2c62"))

    if not mal_df.empty:
        story.append(PageBreak())
        story.append(Paragraph("Reported western-blot malpractice cases mapped to MW bins", styles["Heading2"]))
        story.append(
            Paragraph(
                "These are reported public cases supplied in the malpractice input TSV. They are evidence examples for integrity-audit results sections, not substitution recommendations.",
                styles["Warn"],
            )
        )
        rows = [[pcell(x, styles["TableHead"]) for x in ["Case", "Claim/protein", "MW bin", "Reported issue", "Results sentence"]]]
        for _, r in mal_df.iterrows():
            rows.append(
                [
                    pcell(r["case_id"], styles["TableCell"]),
                    pcell(r["reported_protein_claim"], styles["TableCell"]),
                    pcell(r["mw_bin_label"], styles["TableCell"]),
                    pcell(r["reported_issue"], styles["TableCell"]),
                    pcell(r["results_sentence"], styles["TableCell"]),
                ]
            )
        story.append(make_table(rows, [3.0 * cm, 4.0 * cm, 2.5 * cm, 7.5 * cm, 11.0 * cm], header_bg="#7a2e00"))

    story.append(PageBreak())
    story.append(Paragraph("Audit checklist", styles["Heading2"]))
    checklist = [
        "Check the claimed protein's expected molecular weight against the full blot and marker ladder, not just the cropped panel.",
        "For crowded MW neighborhoods, request raw imager files, exposure settings, sample maps, and antibody validation records.",
        "Treat similar molecular weight as an audit-risk signal only. It is not evidence of protein identity or substitution by itself.",
        "Use knockdown/knockout, overexpression, peptide blocking, or an independent antibody when protein identity is central to the claim.",
        "Pay special attention to configured marker genes, signaling prefixes, and family-pattern matches supplied in the external TSV files.",
        "When malpractice case bins are supplied, use them to motivate stricter review of those MW neighborhoods without implying that unrelated proteins are interchangeable.",
    ]
    for item in checklist:
        story.append(Paragraph("• " + html_escape(item), styles["BodyText"]))

    doc.build(story, onFirstPage=page_footer, onLaterPages=page_footer)
    return pdf_path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Generate a defensive western-blot molecular-weight integrity atlas.")
    parser.add_argument("--proteome-fasta", required=True, help="Proteome FASTA file, preferably UniProt reference proteome FASTA")
    parser.add_argument("--marker-file", required=True, help="TSV with columns gene_symbol, category, note")
    parser.add_argument("--signaling-prefix-file", required=True, help="TSV with columns prefix, category, note")
    parser.add_argument("--family-pattern-file", required=True, help="TSV with columns family, regex, category, note")
    parser.add_argument("--malpractice-file", default=None, help="Optional TSV of reported malpractice cases")
    parser.add_argument("--species-name", default="Unspecified species", help="Species/report label")
    parser.add_argument("--outdir", default="wb_integrity_atlas_out", help="Output directory")
    parser.add_argument("--bin-kda", type=float, default=2.0, help="Molecular-weight bin width in kDa")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    marker_df = load_marker_file(args.marker_file)
    signaling_df = load_signaling_prefix_file(args.signaling_prefix_file)
    family_df = load_family_pattern_file(args.family_pattern_file)
    mal_df = load_malpractice_file(args.malpractice_file)

    df = proteome_dataframe_from_fasta(args.proteome_fasta)
    df, bins = annotate(df, marker_df, signaling_df, family_df, args.bin_kda)
    mal_df = attach_malpractice_bins(mal_df, args.bin_kda)

    paths = write_outputs(df, bins, mal_df, args.outdir)
    pdf_path = build_pdf(
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
    )
    paths["pdf"] = pdf_path

    print("Generated:")
    for key, value in paths.items():
        print(f"  {key}: {value}")
    print(f"Proteins parsed: {len(df):,}; MW bins: {len(bins):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
