# 00_shared — Shared foundation — one implementation per metric

Every downstream script imports the matrix programme, the marker sets and the enrichment test from here. Nothing statistical is defined anywhere else, so a result cannot silently diverge between analyses.

- `_module.py` — the fixed 173 convergent genes with their frozen directions, and the only implementations of the directional concordance test, the two-set contrast and the baseline-stratified null.
- `_marker.py` — the canonical, textbook lineage marker sets (osteogenic, adipogenic, myogenic, chondrogenic, undifferentiated). They exist only to calibrate a measurement level.
- `_enrichment.R` — the single enrichment (Fisher) test used for all gene-set analyses.

## Scripts


_“repo-runnable” = reproduces from files in this repository alone. “needs raw data” = archival; requires the ~98 GB public raw data under `data_raw/` (see `00_setup.md`)._
