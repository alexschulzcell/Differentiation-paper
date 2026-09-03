# data_acquisition — GEO cohort search and screening — metadata only

How the cohorts were found and screened. These scripts query GEO for accessions and sample metadata and apply the preregistered exclusion codes; they perform NO analysis and (except the explicit fetch step) download nothing. Kept for provenance, clearly separated from the analysis.

## Scripts

- `10_cohort_search.py` — _needs raw data_
- `11_cohort_prescreen.py` — _needs raw data_
- `12_cohort_sample_metadata.py` — _needs raw data_
- `13_cohort_fetch.py` — _needs raw data_
- `14_reference_panels_fetch.py` — _needs raw data_
- `20_donor_search.py` — _needs raw data_
- `21_donor_manual_screen.py` — _needs raw data_

_“repo-runnable” = reproduces from files in this repository alone. “needs raw data” = archival; requires the ~98 GB public raw data under `data_raw/` (see `00_setup.md`)._
