# 00_setup — environment, paths and data acquisition

## 1 · Paths

No path is hard-coded. Every script derives the repository root from its own
location, so a clone works wherever it is put. Two environment variables
override that, for the case where the code and the data live apart:

| variable | meaning | default |
|---|---|---|
| `PAPER_V2_ROOT` | the repository root | the parent directory of the script |
| `LIMB_ATLAS` | path to the human fetal limb atlas `.h5ad` | inside `data_raw/` |

```bash
export PAPER_V2_ROOT="/path/to/repo"
```

## 2 · Software

R 4.4.3 — `ggplot2`, `patchwork`, `ragg`, `systemfonts`, `DESeq2`,
`matrixStats`, `rtracklayer`, `GenomicRanges`.
Python 3.12 — `numpy`, `pandas`, `scipy`.

Every script writes `sessionInfo()` or the package versions into its log; the
figure scripts write a session-information file into `results/`.

## 3 · Raw data

**The 96 GB of raw data are not in this repository.** They are public. Every
series used is listed with its accession in Supplementary Table 1
(`figures/data/TS1_eighteen_datasets.csv`); every screened series, including
the excluded ones with their exclusion code, is in Supplementary Table 2.

Reference data with retrieval dates:

| source | what | retrieved |
|---|---|---|
| GENCODE v46 (hg38) | gene annotation, union exon lengths | — |
| gnomAD | LOEUF constraint | — |
| Genomics England PanelApp | skeletal dysplasia panels (309 green; 1 471 broad) | — |
| Nosology of Genetic Skeletal Disorders | core and broad panels | — |
| NCBI `gene2pubmed` and `gene2ensembl` | publications per gene, tax_id 9606 | **2026-08-22** |
| Human fetal limb atlas | post-conception weeks 5.1–9.3, 136 311 cells | — |

The repository itself is the archive for the code and the derived
tables: **https://github.com/alexschulzcell/Differentiation-paper**.

## 4 · Order of execution

```bash
python 03_lineage_calibration/10_calibration_18_datasets.py  # the per-data-set calibration, 2 of 18
python 03_lineage_calibration/12_calibration_gene_space.py                     # which gene space belongs under it
python 07_in_vivo_growth_plate/12_fetal_donor_trend_test.py                         # the donor-stratified in vivo test
Rscript 06_orthogonal_layers/60_gene_sets_build.R                        # fix the broad gene sets (writes once)
Rscript 06_orthogonal_layers/61_gene_set_enrichment.R                               # Figure 2D against broad external sets
python 04_programme_decomposition/10_decomposition_18_datasets.py                     # the decomposition on all 18 data sets
python 07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py                         # atlas pseudobulk, cached (writes once)
python 07_in_vivo_growth_plate/13_fetal_gene_decomposition.py                 # is the in vivo trend broadly carried?
python 07_in_vivo_growth_plate/20_postnatal_growth_plate_annotation.py          # postnatal zones (needs the raw archive)
python 07_in_vivo_growth_plate/21_postnatal_growth_plate_test.py                # postnatal anchor: not calibratable
python 07_in_vivo_growth_plate/14_hypertrophic_zone_sensitivity.py              # does the fetal anchor hang on the terminal zone?
python 05_programme_validation/10_heldout_and_robustness.py                          # held-out validation and robustness (Figure 2G,H)
python 05_programme_validation/11_external_differentiation_systems.py                        # external validation on independent data (Figure 2I)
python 09_figures/10_panel_data_main.py                                 # one CSV per main-figure panel
python 09_figures/11_panel_data_supplement.py                            # supplement panels and Tables S1-S14
python 09_figures/12_panel_data_second_cohort.py                                    # the panels of supplementary Figure S9
Rscript 09_figures/20_figures_main.R                               # F1 to F6, PDF and PNG
Rscript 09_figures/21_figures_supplement.R                         # S1 to S9, PDF and PNG
python 09_figures/30_graphical_abstract.py                         # graphical abstract, PDF and PNG
python 10_manuscript_checks/10_check_numbers.py                              # every number against its panel file
```

Submission packaging (`10_manuscript_checks/20_key_resources_table.py`,
`21_build_submission.py`) and the manuscript-side checks
(`11_check_references.py`, `12_check_language.py`) live in
`10_manuscript_checks/` and complete the reproducible chain.

The panel and figure steps read only files that already exist under
`derived_data/` and `results/`; they recompute nothing and take seconds. Steps
`20_`, `24_`, `25_`, `27_`, `35_` and `36_` recompute one analysis each from stored
per-sample or per-gene values and take seconds to a few minutes. `35_` reads
only the per-dataset dWT matrix (`derived_data/reference_tables/20d_dWT_matrix.csv.gz`),
the dataset-to-study map (`dataset_study_map.csv`) and the frozen gene table,
all of them in the repository. `36_` reads only the archived external
log-fold-changes (`derived_data/reference_tables/external_differentiation_logfc.csv.gz`,
built once from the public GEO processed matrices of GSE37558, GSE283759 and
GSE214987) and the frozen gene table.

Two steps are **exceptions to "reads only stored files"** and are marked as
such:

| step | needs | note |
|---|---|---|
| `60_gene_sets_build.R` | `msigdbr`, `org.Hs.eg.db`, `GO.db` | writes the frozen gene sets **once** and refuses to overwrite them; set `GENSAETZE_V2_NEU=1` to rebuild |
| `11_fetal_atlas_pseudobulk_store.py` | the 7.6 GB limb atlas `.h5ad` and about 3 GB of memory | writes `results/invivo_pseudobulk.csv.gz` **once**; set `INVIVO_PSEUDOBULK_NEU=1` to rebuild. Path override: `LIMB_ATLAS` |

`07_in_vivo_growth_plate/20_postnatal_growth_plate_annotation.py` additionally needs the unpacked
GEO archive of the postnatal growth-plate series, plus scanpy with leidenalg.
(The literature-derived zone-marker sets are frozen in the repository at
`derived_data/reference_tables/growth_plate_zone_markers.csv`, so the postnatal
*test* in step `21` reproduces without raw data.) Its decision rule stands in
its own header, before any number. The outcome is fixed by guards in
`10_check_numbers.py`: that level is **not calibratable**, and should a re-run
pass its positive control, those guards fire and the manuscript text has to be
revisited.

Both cached outputs are checked against their own source on every use: `25_`
asserts that its own-axis contrast reproduces the stored calibration
character for character, and `27_` asserts the same against the stored
per-sample module values. If either assertion fails, the script stops.

Re-running the upstream analyses — the ones that touch raw data — needs the raw
data; those scripts are the numbered upstream stages (`01`–`03`, `06`–`08`) and
`data_acquisition/`, each with its own header and each marked "needs raw data"
in its stage `README.md`.

## 5 · Deposit

The repository is archived on Zenodo on submission. After deposition the DOI
is entered here and in the *Data and code availability* section of the
manuscript:

- **DOI 10.5281/zenodo.XXXXXXX (assigned on deposition)**

The deposit contains no raw data; the accessions are in Supplementary
Table 1.
