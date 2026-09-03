# `derived_data/reference_tables/` -- input files, **not** the supplementary tables

Checked and decided on **2026-08-24** (task 7.3 of session prompt V).

## What this is

Eighteen CSV files numbered **S1-S18**. That is the numbering of an
**earlier** version of the paper and has nothing to do with today's
supplementary tables **TS1-TS14**. The files are **inputs** of the running
chain, not deliverables:

| File | read by |
|---|---|
| `S5_konvergente_gene.csv` | `00_shared/_module.py` (the 173 module genes), `06_orthogonal_layers/61_gene_set_enrichment.R`, `08_disease_gene_orthogonality/50_disease_genes_two_layers.py` |
| `S1_sichtung_alle_datensaetze.csv` | `09_figures/11_panel_data_supplement.py`, `08_disease_gene_orthogonality/21_diagnosis_shear.py` |
| `S7_kohorte_18_datensaetze.csv` | `09_figures/11_panel_data_supplement.py` |
| `S6_gensets.csv` | cited by name in the methods part of the manuscript |

**The folder therefore stays where it is.** Moving it to `_archive/` would
pull the module list out from under the chain -- the same caution that the
session prompt names for `derived_data/manuscript/` (`PDAT`) applies here.

## Why the numbering does not matter anyway

The supplementary tables of the manuscript are called **TS1-TS14** and live
in `figures/data/`; they are built by `09_figures/10_panel_data_main.py` and
`09_figures/11_panel_data_supplement.py`, delivered as
`submission/Supplementary_Tables.xlsx` (one sheet per table, English column
headers). There is no mapping between the two number ranges, and the
manuscript refers only to the TS range.

Whoever wants to tidy up here must touch the three reading scripts at the
same time -- and some of those are scripts of the analysis chain that this
session did not touch.
