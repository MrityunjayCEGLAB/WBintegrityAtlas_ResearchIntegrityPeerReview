# WB-IntegrityAtlas

## A Biology-Informed Framework for Evaluating Western Blot Protein-Identity Claims Using Proteome-Scale Molecular-Weight Neighbourhoods

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

---

## Overview

WB-IntegrityAtlas is a proteome-scale computational framework designed to support biology-informed evaluation of western blot protein-identity claims.

Unlike conventional western blot image-forensics tools, which primarily detect image duplication, manipulation, inappropriate splicing, or other visual irregularities, WB-IntegrityAtlas evaluates the **biological context** surrounding protein identity.

The software systematically characterises molecular-weight neighbourhoods across an entire reference proteome using protein density, common loading controls, signalling proteins, and audit-relevant protein families to identify molecular-weight regions where protein-identity claims may require stronger supporting biological evidence.

The framework is intended as a **defensive research-integrity resource** to support:

- Journal editors
- Peer reviewers
- Research-integrity professionals
- Institutional investigation committees
- Biomedical researchers

during the evaluation of western blot evidence.

---

## Why WB-IntegrityAtlas?

Western blot integrity assessment has traditionally focused on detecting:

- duplicated panels
- inappropriate image manipulation
- image splicing
- image reuse
- AI-generated scientific images

Although these approaches are extremely valuable, they do not determine whether a western blot band actually represents the reported protein.

Protein identity depends upon biological evidence including

- theoretical molecular weight
- neighbouring proteins with similar migration
- loading controls
- signalling proteins
- protein-family relationships
- antibody specificity
- orthogonal validation

WB-IntegrityAtlas addresses this complementary problem by providing proteome-wide biological context for western blot interpretation.

---

# Key Features

✓ Proteome-wide molecular-weight atlas

✓ Automatic theoretical molecular-weight calculation

✓ Protein-density analysis

✓ Common loading-control annotation

✓ Compartment-marker annotation

✓ Signalling-protein annotation

✓ Audit-relevant protein-family annotation

✓ Molecular-weight neighbourhood construction

✓ Audit-priority scoring

✓ Comprehensive PDF report generation

✓ Molecular-weight auditability analysis

✓ Species-independent design

---

# Workflow

```
Reference proteome FASTA
            │
            ▼
Sequence parsing
            │
            ▼
Theoretical molecular-weight calculation
            │
            ▼
Protein assignment into 2-kDa neighbourhoods
            │
            ▼
Biological annotation

    • Protein density
    • Common markers
    • Signalling proteins
    • Protein families

            │
            ▼
Neighbourhood summarisation
            │
            ▼
Audit-burden scoring
            │
            ▼
Audit-priority classification
            │
            ▼
Generated outputs

• Protein atlas
• Bin summary
• Figures
• PDF report
• Auditability analysis
```

---

# Repository Structure

```
WB-IntegrityAtlas/

configs/
    Species-specific annotation files

malpractice/
    Publicly documented protein-substitution cases

demo_full_output/
    Example output generated from demo FASTA

reference_full_output/
    Reference atlas output

wb_integrity_atlas_full_report.py
    Main atlas-generation pipeline

wb_integrity_atlas_externalized.py
    Legacy implementation

wb_mw_auditability.py
    Molecular-weight auditability analysis

requirements.txt

README.md
```

---

# Installation

Clone the repository

```bash
git clone https://github.com/USERNAME/WB-IntegrityAtlas.git

cd WB-IntegrityAtlas
```

Install dependencies

```bash
python -m pip install -r requirements.txt
```

or

```bash
pip install pandas matplotlib reportlab numpy
```

---

# Required Inputs

The main pipeline requires:

### 1. Reference proteome FASTA

Example

```
human_UP000005640.fasta
```

---

### 2. Marker annotation

```
configs/human_marker_genes.tsv
```

Contains

- gene symbols
- marker categories
- notes

---

### 3. Signalling prefixes

```
configs/human_signaling_prefixes.tsv
```

Contains

- signalling prefixes
- pathway annotations

---

### 4. Protein-family annotations

```
configs/human_family_patterns.tsv
```

Contains

- family names
- regular-expression patterns
- biological categories

---

### 5. Malpractice dataset (optional)

```
malpractice/human_reported_malpractice_cases.tsv
```

Used for mapping documented protein-substitution cases onto the atlas.

---

# Running WB-IntegrityAtlas

Example:

```bash
python wb_integrity_atlas_full_report.py \
--proteome-fasta human_UP000005640.fasta \
--marker-file configs/human_marker_genes.tsv \
--signaling-prefix-file configs/human_signaling_prefixes.tsv \
--family-pattern-file configs/human_family_patterns.tsv \
--malpractice-file malpractice/human_reported_malpractice_cases.tsv \
--species-name "Homo sapiens" \
--outdir human_wb_atlas
```

---

# Outputs

Each analysis generates the following outputs.

---

## Protein Atlas

```
proteome_wb_integrity_atlas.csv
```

Contains one row per protein including

- accession
- protein name
- gene symbol
- sequence length
- theoretical molecular weight
- molecular-weight neighbourhood
- marker annotation
- signalling annotation
- protein-family annotation
- audit-burden score
- audit-priority category
- defensive review recommendation

---

## Molecular-Weight Bin Summary

```
mw_bin_summary.csv
```

Contains

- protein count
- marker count
- signalling count
- family count
- representative proteins
- dominant families
- audit-burden score
- audit priority

---

## Figures

High-resolution publication-quality figures including

- Molecular-weight distribution

- Protein-density distribution

- Marker distribution

- Signalling-protein distribution

- Audit-burden distribution

---

## PDF Report

```
western_blot_integrity_atlas_full_report_with_plots.pdf
```

Contains

- summary statistics

- methods

- publication-quality figures

- molecular-weight neighbourhood tables

- highest-priority regions

- integrity recommendations

---

# Molecular-Weight Auditability Analysis

Following atlas generation, molecular-weight auditability can be analysed using

```bash
python wb_mw_auditability.py human_wb_atlas
```

This module quantifies how informative a protein's theoretical molecular weight is for distinguishing it from other proteins within the reference proteome.

Generated outputs include

- protein auditability tables

- neighbourhood auditability

- ambiguity rankings

- information-theoretic metrics

- publication-quality figures

---

# Species Support

The framework is fully configurable.

Changing

- reference proteome

- marker annotations

- signalling annotations

- family annotations

allows straightforward adaptation to additional organisms without modifying the underlying workflow.

---

# Intended Applications

WB-IntegrityAtlas may assist

- journal editors

- peer reviewers

- research-integrity investigators

- laboratory researchers

- educators

during

- western blot review

- manuscript assessment

- editorial screening

- biological interpretation of protein identity

- research-integrity training

---

# Important Interpretation

WB-IntegrityAtlas is **not** an image-forensics tool.

The framework **does not**

- detect image manipulation

- identify scientific misconduct

- determine protein identity

- replace expert scientific judgement

Instead, it identifies molecular-weight neighbourhoods where additional supporting evidence may reasonably be requested.

A high audit-priority score **does not imply**

- wrongdoing

- incorrect protein identity

- fabricated data

It indicates that the surrounding molecular-weight neighbourhood contains numerous biologically plausible alternative proteins and therefore warrants more careful evaluation.

---

# Citation

If you use WB-IntegrityAtlas in your research, please cite

Dwivedi M., Vijay N.

**Beyond Image Duplication: A Proteome-Scale Framework for Auditing Western Blot Protein-Identity Claims**

(Manuscript submitted.)

---

# License

MIT License

---

# Contact

**Mrityunjay Dwivedi**

Computational Evolutionary Genomics Laboratory

Department of Biological Sciences

Indian Institute of Science Education and Research (IISER) Bhopal

India

---

# Acknowledgements

The authors thank members of the Computational Evolutionary Genomics Laboratory, IISER Bhopal, for discussions and feedback during the development of WB-IntegrityAtlas.

---

## Disclaimer

WB-IntegrityAtlas is intended exclusively as a defensive research-integrity and educational resource.

It should be used to support transparent, evidence-informed evaluation of western blot protein-identity claims alongside experimental validation and expert scientific judgement. The framework is not intended to diagnose scientific misconduct or replace established editorial or institutional investigation procedures.
