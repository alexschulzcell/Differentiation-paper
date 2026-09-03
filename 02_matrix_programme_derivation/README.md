# 02_matrix_programme_derivation — Derive and freeze the 173-gene matrix programme

Loads the 18 data sets in one common format, establishes the metric's detection limit by an empirical noise scan, and derives the convergent gene programme at gene level. The resulting 173-gene list and its per-gene directions are then FROZEN (`derived_data/reference_tables/S5_konvergente_gene.csv`) and never refitted downstream.

## Scripts

- `10_load_reference_metric.R` — _needs raw data_
- `11_load_18_datasets.R` — _needs raw data_
- `12_metric_reference.R` — _needs raw data_
- `13_geo_matrices_to_metric_format.R` — _needs raw data_
- `14_geo_matrices_s5_format.R` — _needs raw data_
- `20_detection_limit_scan.R` — _needs raw data_
- `21_detection_limit_analysis.py` — _needs raw data_
- `30_gene_level_convergence_build.R` — _needs raw data_
- `31_derive_matrix_programme.py` — _needs raw data_
- `32_dexamethasone_confounder.R` — _needs raw data_
- `33_aggregate_direction_replication.py` — _needs raw data_

_“repo-runnable” = reproduces from files in this repository alone. “needs raw data” = archival; requires the ~98 GB public raw data under `data_raw/` (see `00_setup.md`)._
