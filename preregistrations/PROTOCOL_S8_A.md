> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S8 step A — the release of the orthogonal single-cell test

Dated **2026-08-19**, before any new S8 number and before any new GEO download.
The preregistration in force: `preregistrations/PREREG_S8.md`.

## The state of knowledge before this step

- `GSE196652`: 317 cells passing quality control, S7-AB2 binding; **not**
  attempted again as a replication.
- `GSE337700`: downloaded, eight Seurat objects analysed as context; 98-99 %
  `transitioning` under the six C8 signatures; no undifferentiated arm, no perturbation.
- `s7n1_geo_summaries.csv`: 945 hits of the single broad E-Search; the R1
  candidates (GSE166824, GSE324998, GSE255646) and the R2 sensitivity contexts
  (GSE241505, GSE150768) fixed in advance with their roles (`PREREG_S8.md` §1).
- No R1, R2 or R3 number has been computed so far; no R1 or sensitivity data set
  has been downloaded so far.

## What this step fixes (by reference)

- Roles, data sets, the 173-gene set as confirmatory, the cell-type signatures
  and the assignment rule, the quality control, replicate rule and author
  exclusions, the primary and sensitivity analysis per role, the MDE80 and
  confidence-interval reasoning, the noise expectation, the step E outcomes and
  the stopping criteria: `PREREG_S8.md` §1-§8.
- The base data-set size and the SHA-256 rule for the three R1 and the two R2
  sensitivity downloads: `PREREG_S8.md` addendum 1.

## Release

S8 may now begin step B: first the temporary directories and the download of
the R1 data sets `GSE166824`, `GSE324998` and `GSE255646` (and the R2
sensitivity contexts `GSE241505`, `GSE150768`, as soon as their sources are
documented), each with its size and SHA-256 logged, **before** any number or
any cell-type assignment is looked at. The primary R2 analysis in `GSE337700`
begins only after the download and SHA step of the R1 data, so that every
number arises afterwards without result distortion.
