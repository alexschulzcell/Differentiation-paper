# =============================================================================
# 20d_genes.R -- EXPLORATORY. Changing the level: not gene sets but GENES.
# =============================================================================
# *** EXPRESSLY EXPLORATORY. *** No confirmatory claim.
#
# WHY. The clustering step showed that on these data the 105 gene sets have
# effectively only about ELEVEN degrees of freedom (101 sets fall into 11
# groups at r >= 0.90, one of them holding 80 sets), and that the one dominant
# dimension is the shared offset of the reference side. A scan over gene sets
# can therefore find almost nothing here structurally -- not because nothing is
# there but because the resolution is missing. Everything in this work since S1
# has been computed at set level.
#
# This script produces the GENE-WISE matrix for the first time: the interaction
# term iv per gene and per data set. The question "do the data sets converge"
# can thereby be put at the level where the resolution sits -- about 10 000
# genes instead of 11 effective sets. The enrichment (skeletal dysplasia, short
# stature, matrisome, GO) comes AFTERWARDS, in 20e, and then against the right
# null.
#
# Only kern() is computed per data set -- NO bootstrap, no NB. That is cheap;
# the expensive part of this work was always the null per gene set, and it is
# not needed yet here.
#
# Written per point: gene, iv, dWT, basis, im_pool. Everything is thereby
# reconstructible in 20e, including the control on the baseline expression.
#
# DECISIONS:
#  (1) ALL genes with a valid iv are written out, not only the pool. Pool
#      membership (dWT >= 0.5) stands beside it as a column -- 20e can then
#      compute both, and the pool rule stays visible instead of baked in.
#  (2) No filtering by effect size here. Filtering is analysis.
#  (3) One partial file per point; existing ones are skipped.
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
PUNKTE <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:18

# The session tree is not part of the public archive; set PAPER_V2_SESSIONS.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
NUR_PUNKTE <- PUNKTE
set.seed(20260819)
source(file.path(SESSIONS, "13_Konvergenzachsen", "reference_implementations",
                 "13_load18.R"))
OUT20 <- file.path(SESSIONS, "20_Exploration", "derived_data")  # AFTER the source()
dir.create(OUT20, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("20d_genes.R -- gene-wise matrix, points %s\n",
            paste(PUNKTE, collapse = ",")))

for (ii in intersect(PUNKTE, seq_along(LABELS_18))) {
  lab  <- LABELS_18[ii]
  teil <- file.path(OUT20, sprintf("20d_gene_%02d.csv", ii))
  if (file.exists(teil)) { cat(sprintf("[%d] %s already present.\n", ii, lab)); next }
  if (is.null(DATEN[[lab]])) { cat(sprintf("[%d] not loaded.\n", ii)); next }
  d <- DATEN[[lab]]
  K0 <- kern(d$Z, d$meta)
  g <- names(K0$iv)
  ok <- is.finite(K0$iv) & is.finite(K0$dWT)
  T <- data.frame(punkt = ii, gen = g[ok], iv = K0$iv[ok], dWT = K0$dWT[ok],
                  basis = K0$basis[ok],
                  im_pool = g[ok] %in% K0$pool, stringsAsFactors = FALSE)
  write.csv(T, teil, row.names = FALSE)
  cat(sprintf("[%2d] %-30s %6d genes, of which %5d in pool\n", ii, lab, nrow(T),
              sum(T$im_pool)))
  DATEN[[lab]] <- NULL; rm(d, K0); invisible(gc())
}
cat("20d_genes.R done.\n")
