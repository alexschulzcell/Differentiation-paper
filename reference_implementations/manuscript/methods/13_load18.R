# =============================================================================
# 13_load18.R -- Loads all EIGHTEEN points: the eleven from 04_load.R and the
#                seven from step B, in exactly the same format.
# =============================================================================
# As in 04_load.R: do not rewrite it, call it. This script does not implement
# the metric -- it appends the seven new data sets to DATEN in the same
# preparation that vorb2() in 03_metric.R applies to the eleven old ones (rlog
# for raw counts, otherwise a variance filter; then zmat()).
#
# DATEN[[i]] holds: Z, meta, expr, label, arm, rolle, klasse -- unchanged.
#
# The calling convention follows the distal-axis script: through the variable
# NUR_PUNKTE (set before the source()) the script can be told to discard
# everything else after loading. That keeps memory under control -- 18 data
# sets in one process exceed the 15.4 GB when several processes run at the
# same time (section 12 of the S5 brief).
# =============================================================================
suppressMessages({ library(DESeq2); library(matrixStats) })
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- The path parameter: the session tree ---------------------------------
# This script used to point at a directory that survived only as a Windows
# directory junction onto a backup -- a junction cannot be cloned or
# versioned, and its target was documented nowhere. SESSIONS points at the
# tree of raw analysis sessions, which holds the reference implementation of
# the metric and the loading block. That tree is not part of the public
# archive; set PAPER_V2_SESSIONS to point at it.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
stopifnot(dir.exists(SESSIONS))

# 04_load.R is sourced from the repaired copy IN THE REPOSITORY (identical in
# content to the session version; only the path block is parameterised).
.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
METHODEN <- Sys.getenv("PAPER_V2_METHODS")
if (!nzchar(METHODEN)) METHODEN <- .skriptordner()
source(file.path(METHODEN, "04_load.R"))
stopifnot(length(DATEN) == 11)

source(file.path(SESSIONS, "03_Metrik_Elf_Punkte", "reference_implementations",
       "03b_data_s5.R"))

# Word for word from 03_metric.R (zmat / vorb2), only without the loading part
zmat18 <- function(X) { Z <- t(scale(t(X))); Z[!is.na(rowSums(Z)), , drop = FALSE] }
vorb18 <- function(X, roh) {
  X <- X[!duplicated(rownames(X)), , drop = FALSE]
  if (roh) {
    X <- X[rowSums(X >= 5) >= max(3, floor(ncol(X) / 4)), , drop = FALSE]
    d <- DESeqDataSetFromMatrix(round(X), data.frame(x = rep(1, ncol(X))), ~ 1)
    X <- assay(rlog(d, blind = TRUE))
  } else X <- X[rowVars(as.matrix(X)) > 0, , drop = FALSE]
  list(Z = zmat18(as.matrix(X)), expr = rowMeans(as.matrix(X)))
}

NEUE <- names(NEU_S5)
IDX_NEU <- 11L + seq_along(NEUE)
BRAUCHE <- if (exists("NUR_PUNKTE")) NUR_PUNKTE else seq_len(11L + length(NEUE))

for (k in seq_along(NEUE)) {
  n <- NEUE[k]; ii <- IDX_NEU[k]
  if (!(ii %in% BRAUCHE)) { cat(sprintf("  skipped: %s (point %d)\n", n, ii)); next }
  D <- NEU_S5[[n]]()
  V <- vorb18(D$X, isTRUE(D$roh))
  M <- NEUMETA_S5[[n]]
  DATEN[[M$lab]] <- list(Z = V$Z, meta = D$meta, expr = V$expr, label = M$lab,
                         arm = M$arm, rolle = "Pruefling", klasse = M$kl)
  cat(sprintf("  loaded: %-32s %-11s class %s  (%d genes x %d samples)\n",
              M$lab, M$arm, M$kl, nrow(V$Z), ncol(V$Z)))
  rm(D, V); invisible(gc())
}
# The numbering of the points is FIXED and does not depend on what this
# process has loaded -- otherwise splitting the work over processes would
# permute the points. LABELS_18[i] is point i, always.
LABELS_18 <- c(names(DATEN)[1:11],
               vapply(NEUE, function(n) NEUMETA_S5[[n]]$lab, character(1)))
names(LABELS_18) <- NULL
stopifnot(length(LABELS_18) == 18, !any(duplicated(LABELS_18)))

# Discard the old points this process does not need
for (ii in seq_len(11)) if (!(ii %in% BRAUCHE)) DATEN[[LABELS_18[ii]]] <- NULL
invisible(gc())

cat(sprintf("13_load18.R: %d data sets in memory, points %s.\n",
            length(DATEN), paste(BRAUCHE, collapse = ",")))
