# =============================================================================
# 40_noise_floor_contrasts.R -- phase M-E, test 1: single-sample contrasts per dataset.
# =============================================================================
# Preregistration: `preregistrations/PRAEREG_M_E.md` §4, dated before the
# first statistic.
#
# Per dataset and per draw, THREE quantities are built from **single samples**,
# i.e. with identical sample count and identical noise structure:
#
#   iv_1x1        = [KO_diff(s1) - KO_undiff(s2)] - [WT_diff(s3) - WT_undiff(s4)]
#   pseudo_iv_1x1 = [WT_diff(a1) - WT_undiff(a2)] - [WT_diff(b1) - WT_undiff(b2)]
#   dWT_1x1       =  WT_diff(s3) - WT_undiff(s4)
#
# `pseudo_iv_1x1` has the same algebra, the same order and the same degrees
# of freedom as `iv_1x1` -- it just crosses NO genotype boundary. It is the
# noise floor against which the lesion response should be measured.
#
# This script DOES NOT IMPLEMENT THE METRIC ANEW: loading uses the reference
# loader `11_load_18_datasets.R`, computation runs on its `Z` (per-gene z per dataset)
# with the same algebra as `kern()` in `12_metric_reference.R`.
#
# Only the SIGN (-1/0/+1) is written, because the convergence rule of the
# main part is sign-based. Counting happens in `41_noise_floor_tests.py`.
#
# Output: derived_data/M_kalibrierung/kontraste/<groesse>_<punkt>.csv.gz
# =============================================================================
# --- path parameters -------------------------------------------------------
# Two separate things, kept in separate variables: where the output goes (the
# repository, found from this file's location or from PAPER_V2_ROOT) and where
# the reference loader comes from (the session tree, PAPER_V2_SESSIONS), which
# is not part of this repository.
.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL   <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
METHODEN <- Sys.getenv("PAPER_V2_METHODEN",
  file.path(WURZEL, "02_matrix_programme_derivation"))

args   <- commandArgs(trailingOnly = TRUE)
NUR_PUNKTE <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:18

source(file.path(METHODEN, "11_load_18_datasets.R"))

# WARNING: the reference loader occupies short names (among others `B`, `D`,
# `Z`). All own constants therefore come AFTER the source() and carry a
# suffix.
NZIEH_E <- 200L
SEED_E  <- 20260823L
AUS_E   <- file.path(WURZEL, "derived_data", "M_kalibrierung", "kontraste")
dir.create(AUS_E, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n40_noise_floor_contrasts.R -- %d Datensaetze, B = %d, Seed %d\n",
            length(DATEN), NZIEH_E, SEED_E))

schreib_E <- function(M, gene, groesse, punkt) {
  d <- data.frame(gen = gene, M, check.names = FALSE)
  con <- gzfile(file.path(AUS_E, sprintf("%s_%02d.csv.gz", groesse, punkt)), "wt")
  write.csv(d, con, row.names = FALSE, quote = FALSE)
  close(con)
}

bilanz <- list()

for (i in seq_along(DATEN)) {
  D <- DATEN[[i]]
  Z  <- D$Z
  me <- D$meta
  if (!is.null(me$sample) && all(colnames(Z) %in% me$sample)) {
    me <- me[match(colnames(Z), me$sample), ]
  }
  stopifnot(nrow(me) == ncol(Z))
  gt <- as.character(me$genotype)
  cd <- as.character(me$condition)
  cd <- ifelse(cd == "naiv", "undiff", cd)

  KOd <- which(gt == "KO" & cd == "diff");   KOu <- which(gt == "KO" & cd == "undiff")
  WTd <- which(gt == "WT" & cd == "diff");   WTu <- which(gt == "WT" & cd == "undiff")

  ok <- length(KOd) >= 1 && length(KOu) >= 1 &&
        length(WTd) >= 2 && length(WTu) >= 2
  cat(sprintf("  %2d %-22s Gene %6d | KO %d/%d  WT %d/%d | %s\n",
              i, D$label, nrow(Z), length(KOd), length(KOu),
              length(WTd), length(WTu),
              if (ok) "ok" else "EXCLUDED -- < 2 WT per condition"))
  bilanz[[length(bilanz) + 1]] <- data.frame(
    lauf_index = i, label = D$label, n_gene = nrow(Z),
    n_KO_diff = length(KOd), n_KO_undiff = length(KOu),
    n_WT_diff = length(WTd), n_WT_undiff = length(WTu),
    verwendbar = ok)
  if (!ok) next

  set.seed(SEED_E + i)
  n <- nrow(Z)
  IV <- matrix(0L, n, NZIEH_E); PIV <- matrix(0L, n, NZIEH_E); DW <- matrix(0L, n, NZIEH_E)
  for (b in seq_len(NZIEH_E)) {
    s1 <- KOd[sample.int(length(KOd), 1)]; s2 <- KOu[sample.int(length(KOu), 1)]
    s3 <- WTd[sample.int(length(WTd), 1)]; s4 <- WTu[sample.int(length(WTu), 1)]
    a  <- WTd[sample.int(length(WTd), 2)]; cc <- WTu[sample.int(length(WTu), 2)]
    dw <- Z[, s3] - Z[, s4]
    IV[,  b] <- as.integer(sign((Z[, s1] - Z[, s2]) - dw))
    PIV[, b] <- as.integer(sign((Z[, a[1]] - Z[, cc[1]]) -
                                (Z[, a[2]] - Z[, cc[2]])))
    DW[,  b] <- as.integer(sign(dw))
  }
  colnames(IV) <- colnames(PIV) <- colnames(DW) <- sprintf("b%03d", seq_len(NZIEH_E))
  schreib_E(IV,  rownames(Z), "iv",     i)
  schreib_E(PIV, rownames(Z), "pseudo", i)
  schreib_E(DW,  rownames(Z), "dwt",    i)
  rm(Z, IV, PIV, DW); gc(verbose = FALSE)
}

B2 <- do.call(rbind, bilanz)
write.csv(B2, file.path(WURZEL, "derived_data", "M_kalibrierung",
                        "test1_datensaetze.csv"), row.names = FALSE)
cat(sprintf("\nWritten: %d datasets usable of %d -> %s\n",
            sum(B2$verwendbar), nrow(B2), AUS_E))
