> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S7 addendum 2 — the context analysis of `GSE337700`

Dated **2026-08-19**, after the download and before the Seurat objects were
read.

## The question and its limit

`GSE337700` is a new, independent clinical single-cell cohort with three
fracture non-union and three healed control samples. An undifferentiated arm is missing;
the data set can therefore deliver neither `dWT` nor the interaction term `iv`
of the reference implementation. It is **not a nineteenth point**, not a
replication of Figure 4 and not a confirmation of the bulk convergence
statement.

It may be used only as an orthogonal, exploratory context for the cell-type and
lesion question: do cell-type shares differ between non-union and healed
fracture, and are the cell-type markers named in advance present in the
respective populations? A positive or negative result is not written back into
the main metric or the existing story.

## The extent of the data

All eight patient-level Seurat objects contained in the GEO TAR are taken into
account. The two non-union samples the original team described as
batch-disturbed (`NUBF01`, `NUBF03`) are not quietly removed; they are marked
as an author exclusion. The primary descriptive table uses the final six-sample
set defined by the original team (three non-union, three controls) and carries
a sensitivity row with the two excluded samples.

## The analysis fixed in advance

- no recomputation of `kern()`, `mk_zieh()`, `kontrast_f()` or `einzel_f()`;
- no new axis, no self-chosen cluster and no result filter;
- cell types are described solely through the external or original annotation
  already present in the objects or, where it is missing, through the six
  MSigDB C8 signatures frozen in `PREREG_S7.md` §3;
- per patient the share of each population is computed; the comparison is made
  at sample level between non-union and healed control, with median, range and
  a Wilson 95 % interval. An MDE80 for the cell shares is not reported as a
  biological power figure, because only three biological samples per group
  exist; the Wilson interval describes only the counting uncertainty within one
  sample;
- all six marker populations are reported, even where one of them is flat or
  cannot be found.

The results are filed in a context protocol of their own. Without a further
preregistration, neither a confirmatory statement nor a new main figure may
arise from them.
