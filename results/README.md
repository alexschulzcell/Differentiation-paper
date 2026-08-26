# `results/` — result files, and what is known about each

Every file here is 1:1 with a figure panel or a supplementary table. The
authoritative panel-by-panel mapping is in `figures/data/` (one CSV per panel,
named after the panel) and in `results/panel_data_log.txt` /
`results/supplement_data_log.txt`, which record how many rows and columns each
panel file has and which source it came from.

| file | produced by | used in |
|---|---|---|
| `invivo_spendertest.csv` | `code/20_in_vivo_donor_test.py` | Fig. 3C, panel `F3C` |
| `eichung_achtzehn.csv` | `reference_implementations/56_calibration_eighteen.py` | Fig. 1D, Table S3 |
| `eichung.csv` | `reference_implementations/54b_cells.py` | Fig. 1E, Table S3 |
| `statistik.csv`, `auslassung.csv`, `zirkularitaet.csv` | `reference_implementations/54c_ladder.py`, `54d_circularity.py` | Fig. S5A, S5B |
| `B_atac_modultest_final.csv`, `B_atac_eichung_je_achse.csv` | `reference_implementations/24_atac_calibration.py` | Fig. 2E/2F, Fig. S6A |
| `figures_*_sessionInfo.txt` | the two figure scripts | reproducibility record |
| `panel_data_log.txt`, `supplement_data_log.txt` | `code/50_`, `code/51_` | reproducibility record |
| `gensaetze_v2_*.csv` | `code/24_gene_sets_v2.R` | Fig. 2D, Tables S9/S9b |
| `zerlegung_achtzehn*.csv/.txt` | `code/25_decomposition_eighteen.py` | Fig. 2F, Table S10 |
| `invivo_pseudobulk.csv.gz` | `code/26_in_vivo_pseudobulk.py` | input to `27_`; cached atlas pseudobulk |
| `invivo_genzerlegung*.csv/.txt` | `code/27_in_vivo_gene_decomposition.py` | Fig. 3C legend, Tables S11/S11b |
| `data_raw/_referenz/geo_primaerarbeiten.csv` | `code/28_geo_primary_publications.py` | Table S12 |
| `data_raw/_referenz/gensaetze_v2/` | `code/24a_gene_sets_v2_build.R` | the fixed broad gene sets, with source/version/date |
| `numbers_check.txt`, `reference_check.txt`, `language_check.txt` | `code/70_`, `71_`, `72_` | the three self-tests |

## Files that need a warning

**These four are flagged deliberately.** Marking them is cheaper than
discovering them again.

| file | status |
|---|---|
| `derived_data/followup/ws8_atac_linienunabhaengigkeit.csv` | **orphan, provenance resolved.** No script in the project produces it. It is a condensation of `B_atac_modultest_final.csv` (rows with the background null) and `B_atac_eichung_je_achse.csv`, verified character-for-character by `code/50_panel_data.py`, which now regenerates it as `figures/data/F2E_atac_per_axis.csv`. **Do not cite the orphan file as a source.** |
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
`code/29_calibration_gene_space.py`.

| file | status |
|---|---|
| `eichung_achtzehn_MIN_N15_ALT.csv` | **superseded.** The pre-2026-08-24 run under the MIN_N ≥ 15 filter. Kept for the record; **not a source of numbers.** |
| `eichung_genraum.csv` | the four-variant decomposition behind that decision (`code/29_calibration_gene_space.py`) |

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
