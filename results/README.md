# `results/` — result files, and what is known about each

Every file here is 1:1 with a figure panel or a supplementary table. The
authoritative panel-by-panel mapping is in `figures/data/` (one CSV per panel,
named after the panel) and in `results/panel_data_log.txt` /
`results/supplement_data_log.txt`, which record how many rows and columns each
panel file has and which source it came from.

| file | produced by | used in |
|---|---|---|
| `invivo_spendertest.csv` | `07_in_vivo_growth_plate/12_fetal_donor_trend_test.py` | Fig. 3C, panel `F3C` |
| `eichung_achtzehn.csv` | `03_lineage_calibration/10_calibration_18_datasets.py` | Fig. 1D, Table S3 |
| `eichung.csv` | `08_disease_gene_orthogonality/31_donor_cells_build_calibrate.py` | Fig. 1E, Table S3 |
| `statistik.csv`, `auslassung.csv`, `zirkularitaet.csv` | `08_disease_gene_orthogonality/32_donor_statistic_ladder.py`, `33_donor_circularity_control.py` | Fig. S5A, S5B |
| `B_atac_modultest_final.csv`, `B_atac_eichung_je_achse.csv` | `06_orthogonal_layers/21_atac_calibration.py` | Fig. 2E/2F, Fig. S6A |
| `figures_*_sessionInfo.txt` | the two figure scripts | reproducibility record |
| `panel_data_log.txt`, `supplement_data_log.txt` | `09_figures/10_panel_data_main.py`, `11_panel_data_supplement.py` | reproducibility record |
| `gensaetze_v2_*.csv` | `06_orthogonal_layers/61_gene_set_enrichment.R` | Fig. 2D, Tables S9/S9b |
| `zerlegung_achtzehn*.csv/.txt` | `04_programme_decomposition/10_decomposition_18_datasets.py` | Fig. 2F, Table S10 |
| `invivo_pseudobulk.csv.gz` | `07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py` | input to `13_fetal_gene_decomposition.py`; cached atlas pseudobulk |
| `invivo_genzerlegung*.csv/.txt` | `07_in_vivo_growth_plate/13_fetal_gene_decomposition.py` | Fig. 3C legend, Tables S11/S11b |
| `derived_data/reference_tables/geo_primary_publications*.csv` | frozen bibliographic metadata (no fetch script in the repo) | Table S12 |
| `data_raw/_referenz/gensaetze_v2/` | `06_orthogonal_layers/60_gene_sets_build.R` | the fixed broad gene sets, with source/version/date |
| `numbers_check.txt`, `reference_check.txt`, `language_check.txt` | `10_manuscript_checks/10_`, `11_`, `12_` | the three self-tests |

## Files that need a warning

**These four are flagged deliberately.** Marking them is cheaper than
discovering them again.

| file | status |
|---|---|
| `derived_data/followup/ws8_atac_linienunabhaengigkeit.csv` | **orphan, provenance resolved.** No script in the project produces it. It is a condensation of `B_atac_modultest_final.csv` (rows with the background null) and `B_atac_eichung_je_achse.csv`, verified character-for-character by `09_figures/10_panel_data_main.py`, which now regenerates it as `figures/data/F2E_atac_per_axis.csv`. **Do not cite the orphan file as a source.** |
| `derived_data/followup/ws5_bilanz.csv` | **hand-made**, the companion table of the WS5 audit. Not a computation and not reproducible by script. Marked as such in the audit itself. |
| `derived_data/M_diagnosen/punkte.csv` | **header row only, no data.** That is not an omission — it is the finding: the diagnosis screen produced no new points. Kept for the record, used nowhere. |
| `derived_data/manuscript/f6_forest.csv` | no unambiguous producing script. Its content is contained in `f6_s12_fixed173_by_study.csv`, which **was** reproduced bit-identically on 2026-08-23. Use that file instead. |

## A file that was corrected, and the one it replaced

**`derived_data/M_kalibrierung/eichung_achtzehn.csv` — resolved 2026-08-24.**
The stored file did not reproduce under its own script. The cause was not the
annotation (both use Gencode v46) but a **filter**: the calibration had been
run on the internal gene map, which keeps only genes measurable in ≥ 15 of the
18 datasets. That filter removes the **terminal** differentiation markers
preferentially while leaving NAIV complete, so it acts asymmetrically against
the half of the calibration that measures arrival at the lineage. The
calibration now runs on unfiltered marker sets and an unfiltered background,
and `bestanden` follows the preregistered z ≥ 2 rule. Head count **2 of 18**
instead of 3. The full account, with the decomposition that isolates the
cause, is kept in the project's internal audit records (working archive, not
published); the rule itself is stated in the header of
`03_lineage_calibration/12_calibration_gene_space.py`.

| file | status |
|---|---|
| `eichung_achtzehn_MIN_N15_ALT.csv` | **superseded.** The pre-2026-08-24 run under the MIN_N ≥ 15 filter. Kept for the record; **not a source of numbers.** |
| `eichung_genraum.csv` | the four-variant decomposition behind that decision (`03_lineage_calibration/12_calibration_gene_space.py`) |

## Files that must not be used as a source of numbers

**`derived_data/manuscript/f7_ebenen.csv` is blocked.** Its `modul_z` column holds a
**rank-based** statistic that has no detection limit, and it sits next to MDE80
values belonging to the **concordance** statistic. The two have been quoted
side by side as if they belonged together, which is where the H3K27ac
"z +4.88" came from (the concordance value is z +3.29 to +3.63). For
GSE33896 osteogenic the two statistics even differ in sign (+0.54 rank versus
−0.46 concordance). Use the per-level module-test files instead, all of which
carry their own MDE80. See `Neu/KONSISTENZ_PROTOKOLL.md` §9a.

**`derived_data/manuscript/m6_nachweisgrenzen.csv`, row `external`, column `z`** holds the
value 1.21, which is not a z value but the observed statistic of the
study-level synthesis. Retired; see `KONSISTENZ_PROTOKOLL.md` §1.

## Numbers that changed on 2026-08-23

Four numbers in earlier internal narratives were not supported by their source
files and have been corrected. Old value, new value, reason and source file for
each are in **`Neu/KONSISTENZ_PROTOKOLL.md`**, which is the file to read before
quoting anything from a document dated before 2026-08-23.
