> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S9 step A — the release of the R2 impact scan

Dated **2026-08-19**, before the new search query and before any new download.

## The state of knowledge before this step

- S8 is complete (R1 only one preparation, giving S8-AB1; R2 and R3 in
  `GSE337700` not decidable because of the cell type, giving S8-AB2).
- Already screened and known **not** to be R2 impact-carrying candidates (from
  the S7 screening table and from S8): the non-union cohort (GSE337700), whose
  annotation collapses; the OPLL cohort (GSE241505), with no non-union against
  healed axis and technical wells; a murine AHR model (GSE150768); and the three
  R1 cohorts (GSE166824, GSE324998, GSE255646), which have no lesion against
  control axis.
- No new S9 search query has been run so far and no new data set has been
  assessed or downloaded.

## What is fixed

The frozen search query, the inclusion rules and the stopping criteria are given
in `PREREG_S9.md` §3-§6. The scan is run as the S9 GEO scan script; the hit
list, the search URL, the date and the return code are stored with the S9
results. Only after a complete screen and design verification is anything
downloaded, with size and SHA-256 recorded before the analysis.

## Release

S9 may now run the search query.
