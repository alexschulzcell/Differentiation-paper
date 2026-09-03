# 06_orthogonal_layers — Orthogonal molecular layers: methylation, accessibility, H3K27ac

Tests the frozen programme on measurement layers other than RNA — promoter methylation, chromatin accessibility and H3K27ac — each calibrated against its own positive control before any programme verdict is read, plus the convergence dose-response and the lineage contrast.

## Scripts

- `10_methylation_osteogenic_27k.py` — _needs raw data_
- `11_methylation_sensitivity.py` — _needs raw data_
- `12_methylation_chondrogenic_450k.py` — _needs raw data_
- `20_atac_window_calibration.R` — _needs raw data_
- `21_atac_calibration.py` — _needs raw data_
- `22_atac_accessibility_osteogenic.py` — _needs raw data_
- `23_h3k27ac_chondrogenic_build.R` — _needs raw data_
- `24_h3k27ac_chondrogenic_test.py` — _needs raw data_
- `30_convergence_dose_integration.py` — _needs raw data_
- `31_programme_vs_markers.py` — _needs raw data_
- `32_dose_response.py` — _needs raw data_
- `40_collect_layer_numbers.py` — _needs raw data_
- `50_second_chromatin_cohort_windows.R` — _needs raw data_
- `51_lineage_contrast.py` — _needs raw data_
- `60_gene_sets_build.R` — _needs raw data_
- `61_gene_set_enrichment.R` — _repo-runnable_

_“repo-runnable” = reproduces from files in this repository alone. “needs raw data” = archival; requires the ~98 GB public raw data under `data_raw/` (see `00_setup.md`)._
