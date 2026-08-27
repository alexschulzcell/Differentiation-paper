# `derived_data/reference_tables/` -- input files, **not** the supplementary tables

Checked and decided on **2026-08-24** (task 7.3 of session prompt V).

## What this is

Eighteen CSV files numbered **S1-S18**, plus the compact public `dWT` matrix
used to verify the module. The S-numbering is that of an
**earlier** version of the paper and has nothing to do with today's
supplementary tables **TS1-TS14**. The files are **inputs** of the running
chain, not deliverables:

| File | read by |
|---|---|
| `S5_konvergente_gene.csv` | `reference_implementations/_module.py` (the 173 module genes), `code/24_gene_sets_v2.R`, `reference_implementations/followup/ws1_zwei_schichten.py` |
| `20d_dWT_matrix.csv.gz` | `reference_implementations/manuscript/methods/20f_convergence_dwt.py` (the public provenance check for S5) |
| `S1_sichtung_alle_datensaetze.csv` | `code/51_supplement_data.py`, `reference_implementations/53_diagnosis_screening.py` |
| `S7_kohorte_18_datensaetze.csv` | `code/51_supplement_data.py` |
| `S6_gensets.csv` | cited by name in the methods part of the manuscript |

**The folder therefore stays where it is.** Moving it to `_archive/` would
pull the module list out from under the chain -- the same caution that the
session prompt names for `derived_data/manuscript/` (`PDAT`) applies here.

## Why the numbering does not matter anyway

The supplementary tables of the manuscript are called **TS1-TS14** and live
in `figures/data/`; they are built by `code/50_panel_data.py` and
`code/51_supplement_data.py`, delivered as
`submission/Supplementary_Tables.xlsx` (one sheet per table, English column
headers). There is no mapping between the two number ranges, and the
manuscript refers only to the TS range.

Whoever wants to tidy up here must touch the three reading scripts at the
same time -- and some of those are scripts of the analysis chain that this
session did not touch.

## Provenance of `S5_konvergente_gene.csv`

`S5_konvergente_gene.csv` is the frozen output of an explicitly exploratory
analysis, not an externally defined or preregistered gene set. The underlying
per-gene `dWT` values are supplied in `20d_dWT_matrix.csv.gz`. The public
`20f_convergence_dwt.py` script reproduces the archived selection rule: first
retain genes measurable in at least 14 of 18 data sets, then median-centre each
data-set column within those 12,563 genes, and retain genes for which at least
90% of the non-zero centred signs agree. It verifies the resulting 173 genes,
their `n` and `v` values, and their directions `ri` against S5. The list is
then used unchanged by downstream analyses.
