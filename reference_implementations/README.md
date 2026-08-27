# `reference_implementations/` — one implementation per metric

These scripts compute the orthogonal measurement levels and the human-genetics
anchor. They are the upstream half of the pipeline: they read raw data from
`data_raw/` and write their results into `derived_data/`. Everything
downstream — the panel data under `figures/data/`, the figures and the
submission package — is built by `code/` from those results and recomputes
nothing.

The rule that governs this folder: **no metric is implemented twice.** The
module concordance, the two-set contrast and the detection limit exist once
each, in `_module.py`, and every script calls them.

A note on names: the file names and the short column names of the outputs are
the internal vocabulary of the pipeline. `column_glossary.csv` in the
repository root gives the English meaning of every column name that appears in
`derived_data/` and `results/`.

## Shared foundation

| file | content |
|---|---|
| `_module.py` | the 173 genes identified in the exploratory formation analysis, frozen thereafter with direction `ri`; **the only** implementation of the directional concordance (`konkordanz`), of the two-set contrast (`kontrast`) and of the baseline-stratified null |
| `_marker.py` | the canonical lineage marker sets (osteogenic, adipogenic, myogenic, chondrogenic, undifferentiated). Textbook knowledge, fixed before any computation and checked for disjointness. They serve **only** to calibrate a measurement level |
| `_enrichment.R` | the one implementation of the enrichment test used for the gene-set analyses |

## Formation of the 173-gene module

The 173-gene module is an internal, data-derived result. It was identified in
the explicitly exploratory `20_Exploration` analysis of the 18 perturbation
data sets, not supplied by an external catalogue or defined before the paper.
`20d_genes.R` produces the per-gene `dWT` matrix. The selection rule is
preserved in section (3) of
`manuscript/methods/20i_dexamethasone.R` and is independently re-run by:

```bash
python reference_implementations/manuscript/methods/20f_convergence_dwt.py
```

The rule first retains genes measurable in at least 14 of 18 data sets, then
median-centres each data-set column within that 12,563-gene universe and keeps
genes with at least 90% concordant non-zero signs. It yields 173 genes, with
129 `ri = +1` and 44 `ri = -1`, exactly matching
`derived_data/reference_tables/S5_konvergente_gene.csv`. The module is frozen
only after this exploratory formation step; downstream analyses never
reselect it. The 15-of-18 filter in `01_internal_gene_map.py` belongs to that
script's separate continuous gene-map analysis and does not form S5.

---

## Order

```bash
python reference_implementations/01_internal_gene_map.py
```
The internal RNA level as a **continuous** quantity per gene (11 581 genes
with at least 15 evaluable data sets).
→ `derived_data/R_intern/`

```bash
python reference_implementations/10_dnam_GSE33896.py
python reference_implementations/11_dnam_sensitivity.py
python reference_implementations/12_dnam450_GSE129266.py
```
The promoter methylomes. `10_` the 27K data set (osteogenic and myogenic,
3 paired donors), `11_` its sensitivities and positive control, `12_` the 450K
data set (chondrogenic, 2 paired donors).
→ `derived_data/A_dnam/`

```bash
Rscript reference_implementations/22_atac_window_calibration.R
python  reference_implementations/24_atac_calibration.py
```
Chromatin accessibility in GSE332758: four measurement windows, the
calibration per axis, the module test.
→ `derived_data/B_atac/`

```bash
python reference_implementations/26_atac_GSE224251_strict.py
```
Chromatin, GSE224251 (n = 3 per arm), eight preparations crossed. All eight are
reported.
→ `derived_data/B_atac/`

```bash
SPUR=H3K27ac Rscript reference_implementations/27_h3k27ac_GSE129031.R
SPUR=input   Rscript reference_implementations/27_h3k27ac_GSE129031.R
python       reference_implementations/28_h3k27ac_module_test.py
```
H3K27ac in the **same** in vitro chondrogenesis as the 450K data set
GSE129266 (the same laboratory and protocol, though not demonstrably the same
donors). Two donor lines. This level makes the comparison of chromatin against
promoter methylome **within one differentiation axis** possible.
→ `derived_data/B_atac/`

```bash
python reference_implementations/30_integration.py
python reference_implementations/31_module_vs_markers.py
python reference_implementations/32_dose_response.py
python reference_implementations/40_figure7_data.py
```
The integration: internal convergence as a dose (`30_`, with a free and a
baseline-stratified null), the module against the canonical marker set on the
same level (`31_`, with a size control), and the dose-response curve across the
convergence octiles (`32_`). `40_` collects the numbers.
→ `derived_data/Z_integration/`, `derived_data/manuscript/`

Total runtime about 40 minutes; `24_` and `30_` account for most of it
(baseline-stratified permutation nulls). The seed is 20260821 throughout.

---

## The patient and human-genetics phase (`50_`–`53_`)

Preregistered in `preregistrations/PRAEREG_M_B.md` and `PRAEREG_M_A.md`, both
dated **before the respective first download**.

| script | what it does | output |
|---|---|---|
| `50_cohort_search.py` | ten search axes against GEO; **nothing is downloaded**, only accessions and metadata | `derived_data/M_patienten/` |
| `50b_screening.py` | mechanical pre-screen with the preregistered exclusion codes; in case of doubt a series stays in | `sichtung_mechanisch.csv`, `kandidaten.csv` |
| `50c_check_candidates.py` | GSM metadata per candidate -- **only this** makes the inclusion criteria decidable; a series title never suffices | `data_raw/_meta/` |
| `51a_fetch.py` | downloads the seven cohorts that survive the hand screen | `data_raw/<GSE>/` |
| `51_patient_variability.py` | **calibration first**, the main computation only where the calibration passes | `streuung.csv`, `eichung.csv`, `synthese.csv` |
| `52a_go_sets.R` | exports the two secretion sets **verbatim** from the reference implementation of the metric | `derived_data/M_humangenetik/` |
| `52b_fetch_panels.py` | seven reference panels with their retrieval date (HPO including ontology propagation, GWAS catalogue, PanelApp) | `panels.csv` |
| `52_human_genetics_anchor.py` | the positive control first, then the enrichment against the expression- and length-matched null | `anker.csv`, `eichung_A.csv` |
| `53_diagnosis_screening.py` | the audit of all 127 hand-checked candidates; the verdicts stand by name in the script | `derived_data/M_diagnosen/` |

**Three things that matter when rebuilding this:**

1. **No second implementation.** `51_` and `52_` call `_module.konkordanz` and
   `_module.kontrast`; the marker sets come unchanged from `_marker.py` and the
   module unchanged from its stored definition.
2. **The calibration runs before the main computation**, not after it. Reverse
   the order and the assignment of the marker sets can no longer credibly be
   called fixed in advance.
3. **The exclusions belong in the result.** The screening files contain
   **every** screened series with its verdict, not only the included ones -- in
   the diagnosis screen the screening is the whole finding.

Runtime about 25 minutes. The seed is 20260821 throughout.

---

## The donor-resolved phase (`54_`, `56_`)

Preregistered in `preregistrations/PRAEREG_M_D.md`, dated **before the first
download**, with an addendum dated **before the first statistic**.

| script | what it does | output |
|---|---|---|
| `54_donor_search.py` | three search axes against GEO, aimed at the study design rather than at entities; **nothing is downloaded** | `derived_data/M_donoren/` |
| `54e_screening_by_hand.py` | the hand screen of every hit and of the candidates named in advance, each with a code and a reason | `sichtung_hand.csv` |
| `54a_self_test.py` | **before the first real number**: is the null of the new statistics calibrated, and do they fire at 0.35 z? | `selbsttest.csv` |
| `54b_cells.py` | builds the 14 cells (donor x axis x study) and runs the **built-in calibration** | `eichung.csv`, `gene_level.csv.gz` |
| `54c_ladder.py` | the preregistered statistic ladder S1 to S3 on the **calibrated** cells, the study synthesis against the donor-flip null, and the complete leave-one-out | `statistik.csv`, `je_zelle.csv`, `auslassung.csv` |
| `54d_circularity.py` | **added afterwards, not preregistered**: the same ladder on the cells that did **not** help define the module | `zirkularitaet.csv` |
| `56_calibration_eighteen.py` | the per-data-set calibration of the 18 perturbation data sets, on unfiltered marker sets | `derived_data/M_kalibrierung/` |
| `56b_calibration_sensitivity.py` | the same calibration on the subset of markers reachable in vitro | `derived_data/M_kalibrierung/` |

**Four things that matter when rebuilding this:**

1. **No second implementation.** The statistic ladder lives **in**
   `_module.py` and uses the same baseline-stratified null as the concordance.
   `_marker.py` is unchanged.
2. **The calibration runs before the main computation.** `54c_` computes only
   on cells whose calibration passed; cells that failed carry no number. Seven
   of fourteen cells pass.
3. **Clones are not donors.** `54b_` assigns the donor identifier per cell
   explicitly; the replicate experiments of GSE221128 and the clones of the
   urine-derived stem cell series are **one** donor each. The null permutes
   cells of the same donor together.
4. **Two degeneracies are named in the code rather than hidden:** the flip null
   is exactly degenerate for S2 (the principal-component share is invariant
   under row signs), and studies with only one donor carry no between-donor
   statistic. In both cases `NaN` is reported rather than an apparent number.

Runtime about 50 minutes, almost all of it in the baseline-stratified nulls of
`54c_` and `54d_` (20 000 draws each) and in the self-test. The seed is
20260822 throughout. **Nothing is downloaded** — all six study units were
already present.

---

## `manuscript/`

`reference_implementations/manuscript/methods/03_metric.R` is the reference
implementation of the metric for Figures 1, 2 and S1 to S5; the other scripts
in that folder produce the derived inputs under `derived_data/manuscript/`.
`reference_implementations/followup/` holds the follow-up analyses whose
outputs live in `derived_data/followup/`.
