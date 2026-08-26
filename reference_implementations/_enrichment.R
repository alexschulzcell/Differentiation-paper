# =============================================================================
# _enrichment.R -- the ONE implementation of the enrichment test.
#
# Project rule 3 ("no second implementation of a statistic"): the Fisher test
# of a gene category against the background measurable on this layer stands
# exactly here. `reference_implementations/followup/ws6_p3_annotation.R` (narrow
# GO sets) and `code/24_gene_sets_v2.R` (broad, independently curated sets)
# call the same function; they differ only in the gene sets they pass in.
#
# Extracted on 2026-08-24 from ws6_p3_annotation.R, unchanged in its
# computation. Equivalence is proven: after extraction ws6_p3_annotation.R
# produces `ws6_p3_go_annotation.csv` sign-identical to the retired version
# (see Neu/KONSISTENZ_PROTOKOLL.md §11).
# =============================================================================

ft <- function(satz, fg, hg) {
  satz <- intersect(satz, hg)
  fg <- intersect(fg, hg)
  a <- length(intersect(satz, fg))
  b <- length(fg) - a
  c <- length(satz) - a
  d <- length(hg) - length(fg) - c
  f <- fisher.test(matrix(c(a, b, c, d), 2))
  list(k = a, n_fg = length(fg), n_satz = length(satz), n_hg = length(hg),
       OR = unname(f$estimate), lo = f$conf.int[1], hi = f$conf.int[2],
       p = f$p.value)
}
