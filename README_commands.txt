Western Blot Integrity Atlas - Complete Report Package
======================================================

Purpose
-------
This package generates the complete Western Blot Integrity Atlas from a
proteome FASTA and the external marker/signaling/family TSV configuration
files. The updated runner combines the previous table-rich report with all
four molecular-weight plots in one consolidated PDF.

The runner is defensive-audit software. It is intended to support review of
western blot identity claims, source-file provenance, antibody validation,
and appropriate controls. It must not be used to plan, enable, or disguise
blot substitution or fabrication.

Primary runner
--------------
Use:

  wb_integrity_atlas_full_report.py

The original runner is retained as:

  wb_integrity_atlas_externalized.py

The updated runner keeps the same command-line inputs as the original runner.
It does NOT require a pre-generated mw_bin_summary CSV, and there is no
--bin_csv flag. The protein atlas, bin summary, malpractice-bin table, plot
PNGs, and complete PDF are all generated internally in one run.

Install dependencies
--------------------
From the wb_atlas_update directory:

  python -m pip install -r requirements.txt

Equivalent direct installation:

  python -m pip install pandas matplotlib reportlab

Human example
-------------
Place the human reference-proteome FASTA in this directory, or give its full
path, then run:

  python wb_integrity_atlas_full_report.py \
    --proteome-fasta human_UP000005640.fasta \
    --marker-file configs/human_marker_genes.tsv \
    --signaling-prefix-file configs/human_signaling_prefixes.tsv \
    --family-pattern-file configs/human_family_patterns.tsv \
    --malpractice-file malpractice/human_reported_malpractice_cases.tsv \
    --species-name "Homo sapiens" \
    --outdir human_wb_atlas

Cow example
-----------
  python wb_integrity_atlas_full_report.py \
    --proteome-fasta cow_UP000009136.fasta \
    --marker-file configs/cow_marker_genes.tsv \
    --signaling-prefix-file configs/cow_signaling_prefixes.tsv \
    --family-pattern-file configs/cow_family_patterns.tsv \
    --malpractice-file malpractice/cow_reported_malpractice_cases.tsv \
    --species-name "Bos taurus" \
    --outdir cow_wb_atlas

Mouse example
-------------
  python wb_integrity_atlas_full_report.py \
    --proteome-fasta mouse_UP000000589.fasta \
    --marker-file configs/mouse_marker_genes.tsv \
    --signaling-prefix-file configs/mouse_signaling_prefixes.tsv \
    --family-pattern-file configs/mouse_family_patterns.tsv \
    --malpractice-file malpractice/mouse_reported_malpractice_cases.tsv \
    --species-name "Mus musculus" \
    --outdir mouse_wb_atlas

Chicken example
---------------
  python wb_integrity_atlas_full_report.py \
    --proteome-fasta chicken_UP000000539.fasta \
    --marker-file configs/chicken_marker_genes.tsv \
    --signaling-prefix-file configs/chicken_signaling_prefixes.tsv \
    --family-pattern-file configs/chicken_family_patterns.tsv \
    --malpractice-file malpractice/chicken_reported_malpractice_cases.tsv \
    --species-name "Gallus gallus" \
    --outdir chicken_wb_atlas

Demo run included in the package
--------------------------------
  python wb_integrity_atlas_full_report.py \
    --proteome-fasta demo_human_subset.fasta \
    --marker-file configs/human_marker_genes.tsv \
    --signaling-prefix-file configs/human_signaling_prefixes.tsv \
    --family-pattern-file configs/human_family_patterns.tsv \
    --malpractice-file malpractice/human_reported_malpractice_cases.tsv \
    --species-name "Homo sapiens demo subset" \
    --outdir demo_full_output

Command-line flags
------------------
--proteome-fasta
  Required. Proteome FASTA, preferably a UniProt reference-proteome FASTA.

--marker-file
  Required. TSV containing gene_symbol, category, and note columns.

--signaling-prefix-file
  Required. TSV containing prefix, category, and note columns.

--family-pattern-file
  Required. TSV containing family, regex, category, and note columns.

--malpractice-file
  Optional. TSV of reported public malpractice cases. When supplied, the
  script generates reported_malpractice_bins.csv and embeds its table in the
  PDF.

--species-name
  Optional report label. Default: "Unspecified species".

--outdir
  Optional output directory. Default: wb_integrity_atlas_out.

--bin-kda
  Optional molecular-weight bin width. Default: 2.0 kDa.

Outputs
-------
Each run creates the following inside --outdir:

  proteome_wb_integrity_atlas.csv
    Complete per-protein atlas with calculated mass, annotations, bin metrics,
    audit priority, and defensive audit note.

  mw_bin_summary.csv
    Complete molecular-weight-bin summary with protein counts, marker counts,
    signaling-prefix counts, family annotations, examples, burden score, and
    priority.

  reported_malpractice_bins.csv
    Generated only when --malpractice-file is supplied.

  plots/figure_2_molecular_weight_distribution.png
  plots/figure_3_most_crowded_bins.png
  plots/figure_4_marker_signaling_distribution.png
  plots/figure_5_audit_burden.png
    Standalone high-resolution PNGs embedded in the PDF.

  western_blot_integrity_atlas_full_report_with_plots.pdf
    The complete consolidated report. It includes:
      - summary metrics and method
      - all four updated plots
      - highest-audit score summary
      - detailed highest-audit table
      - common loading/compartment marker table
      - complete molecular-weight-bin table
      - first 120 proteins sorted by molecular weight
      - reported malpractice case table, when supplied
      - integrity-audit checklist and generated-file inventory

  western_blot_integrity_atlas_report.pdf
    Byte-identical compatibility copy of the consolidated report, preserving
    the previous report filename for downstream workflows.

Package layout
--------------
configs/
  Species-specific marker, signaling-prefix, and family-pattern TSV files.

malpractice/
  Public reported-case TSV files and source tables.

wb_integrity_atlas_externalized.py
  Original externalized table-report runner retained unchanged.

wb_integrity_atlas_full_report.py
  Updated one-pass runner for CSVs, all plots, and the consolidated PDF.

requirements.txt
  Python dependencies.

demo_human_subset.fasta
  Small input for a quick installation test.

demo_full_output/
  Output generated from the included demo FASTA by the updated runner.

reference_full_output/
  Consolidated report and plots generated from the full uploaded atlas CSVs.
  These files demonstrate the intended final layout. The normal production
  workflow still starts from FASTA/configuration inputs and does not use CSV
  inputs.

Important interpretation notes
------------------------------
- Theoretical molecular weight is an audit context, not proof of identity.
- High-priority bins indicate where documentation and orthogonal controls are
  especially important. They do not indicate wrongdoing.
- Similar molecular weight never makes proteins interchangeable.
- Cow, mouse, and chicken case files may include orthology/protein-family audit
  projections. Check species_interpretation before describing a case as direct
  evidence for a species.
- Edit the species TSV configuration files when the FASTA uses alternative gene
  symbols or when the audit categories need to be updated.

Troubleshooting
---------------
1. "No protein sequences were parsed"
   Confirm that --proteome-fasta points to a readable FASTA with records that
   begin with >.

2. Missing required TSV columns
   Keep the column names exactly as documented above. TSV files must be
   tab-delimited.

3. Matplotlib or ReportLab import error
   Run: python -m pip install -r requirements.txt

4. Empty marker table
   Confirm that the FASTA headers include GN= gene symbols matching the marker
   TSV. Alternative symbols may require edits to the species config file.

5. Very large PDF
   This is expected for the full reference proteome because the consolidated
   report contains all bin rows, multiple tables, and four high-resolution
   figures. The complete per-protein table remains in CSV to keep the PDF
   finite; the PDF displays the first 120 proteins exactly as in the previous
   table-rich report.
