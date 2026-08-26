> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S10 step A — the release of the R2 computation with a pinned reference

Dated **2026-08-19**, before any new S10 number. In force:
`preregistrations/PREREG_S10.md`.

## The state of knowledge

- S9: `GSE255460` (8 osteoarthritis against 3 control donors, 135 896 cells) is
  the only design-valid R2 candidate with a sufficient control arm. Under C8 the
  assignment collapsed to 93.4 % `transitioning` (a cross-reference stability
  against HCA of 0.148).
- The count matrix and the metadata are already loaded; the marker rule from the
  S9 marker export is implemented and unchanged.
- No S10 number has been computed yet.

## Release

S10 may now export the reference markers (HCA and BlueprintEncode), stream the
matrix again and run the R2 analysis following `PREREG_S10.md` §6-§7.
