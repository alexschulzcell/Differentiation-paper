> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S7 step A — the release before the computation

Written on **2026-08-19**, before the first S7 quality-control, aggregation or
metric number.

## In force

- `PREREG_whole_study.md` addendum 2 opens the pseudobulk route only for
  `GSE196652`; `A9` stays in force for all other single-cell data sets.
- `PREREG_S7.md` fixes in advance the quality control, the replicate rule, the
  external cell-type definition, the thresholds for step C, the detection limit
  and every outcome.
- The confirmatory axis question from `PREREG_S6.md` §1 stays closed.

## Individual decisions fixed in advance

- Quality control: at least 1 000 detected genes per cell, at most 20 %
  mitochondrial reads.
- Symbol to Ensembl: `org.Hs.eg.db` 3.20.0; unique assignments only, with
  ambiguities and unmappable symbols logged.
- Pseudobulk: one sample per source x exact time point x matched patient, at
  least 10 cells passing quality control, no artificial splitting and no
  merging of time points.
- Cell populations: six MSigDB C8 signatures named in advance, version
  `2026.1.Hs`, loaded with `msigdbr` 26.1.0 on 2026-08-19; a score of at least
  10 % and a distance of at least 2 percentage points to the second-best score.
- The existing reference implementation is called after the aggregation; the
  aggregation step is not a second metric.
- Every reported number carries its MDE80 and its noise expectation. No new
  axis scan and no new gene-set scan is admissible.

## The state of knowledge

At the time of dating, the GEO description, the two raw matrices and their
column structure were known. No S7 quality control, no pseudobulk sample, no
cell-type figure and no new `iv` or `dWT` value had been computed or looked at.

Step A was thereby complete and step B was allowed to begin. The result of
step B is recorded in the session protocol.
