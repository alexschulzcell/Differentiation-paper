# =============================================================================
# 27_h3k27ac_GSE129031.R -- layer B3: H3K27ac during CHONDROGENESIS.
# =============================================================================
# Why this dataset closes the gap
# -------------------------------
# So far the chromatin statement stands on the OSTEOGENIC axis (GSE224251)
# and the methylation statement on the CHONDROGENIC one (GSE129266). The
# comparison "accessibility yes, promoter methylation no" thus runs across
# two different differentiations and is attackable: perhaps the difference
# is due to the axis, not to the measurement layer.
#
# GSE129031 closes exactly this gap. It is H3K27ac (the active regulatory
# mark) in **the same** in-vitro chondrogenesis, from **the same lab** and
# with **the same protocol** as GSE129266 (Figure 5, methylome). This makes
# the layer comparison feasible within one axis.
#
# IMPORTANT, and it stands this way in the report: the donor labels of the
# two series do not match (`8A`/`2454e` vs `Donor 1..4`). It is the same
# system and the same lab, but these are NOT verifiably the same donors.
# The wording remains "same axis, same protocol", never "same donors".
#
# Design: 2 donor lines (8A, 2454e) x 2 states (MSC naive, CHON
# differentiated), H3K27ac, hg38. Two biological units -> the layer is
# DESCRIPTIVE at cohort level; only the gene level is inferential, against
# the background draw, as for GSE332758.
#
# Output: derived_data/B_atac/B3_GSE129031_matrix_<Fenster>.csv
# =============================================================================
suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

# --- path parameters (2026-08-23) ------------------------------------------
# Previously a hard-coded path pointed to ".../DFG Antrag/
# Scherenpaper_Folgeprojekt" -- a directory that only survived as a Windows
# junction onto a backup and is therefore not clonable. Moreover, the
# `data_raw/` read here and the `derived_data/` written here never actually
# lived there, but under `Paper v2` itself.
# WURZEL is from now on the project root; overridable via the environment
# variable PAPER_V2_ROOT so the script runs from within a clone.
.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
DAT <- file.path(WURZEL, "data_raw")
AUS <- file.path(WURZEL, "derived_data", "B_atac")
dir.create(AUS, recursive = TRUE, showWarnings = FALSE)

cat("Layer B3 -- H3K27ac, GSE129031 (chondrogenesis, hg38)\n")
cat(strrep("=", 70), "\n")

gtf <- import(file.path(DAT, "_referenz", "gencode.v46.annotation.hg38.gtf.gz"),
              format = "gtf", feature.type = "gene")
gtf <- gtf[gtf$gene_type == "protein_coding"]
gtf <- gtf[!duplicated(gtf$gene_name)]
neg <- as.character(strand(gtf)) == "-"
tss <- ifelse(neg, end(gtf), start(gtf))

mk <- function(a, b) {
  s <- ifelse(neg, tss - b, tss - a)
  e <- ifelse(neg, tss + a, tss + b)
  g <- GRanges(seqnames(gtf), IRanges(pmax(1L, s), e))
  g$symbol <- gtf$gene_name
  g
}
# H3K27ac sits at promoters AND at distal enhancers. Hence the same windows
# as for chromatin accessibility, so both layers are compared with the same
# geometry.
FENSTER <- list(P = mk(2000L, 500L), T10 = mk(10000L, 10000L),
                T50 = mk(50000L, 50000L))

# Two tracks: the mark itself and the input control. The input measures
# mappability, copy number and fragmentation propensity -- without it an
# H3K27ac difference can be a pure coverage difference. Selection via the
# environment variable SPUR (default: H3K27ac).
SPUR <- Sys.getenv("SPUR"); if (!nzchar(SPUR)) SPUR <- "H3K27ac"
bws <- sort(list.files(file.path(DAT, "GSE129031"),
                       pattern = paste0(SPUR, ".*\\.bigwig$"), full.names = TRUE))
stopifnot(length(bws) == 4)
namen <- sub("^GSM\\d+_", "",
             sub(paste0("_", SPUR), "", sub("\\.bigwig$", "", basename(bws))))
cat("track:", SPUR, "| samples:", paste(namen, collapse = ", "), "\n")

for (fn in names(FENSTER)) {
  ziel <- file.path(AUS, sprintf("B3_GSE129031_%s_matrix_%s.csv", SPUR, fn))
  if (file.exists(ziel)) { cat("exists:", basename(ziel), "\n"); next }
  gr <- FENSTER[[fn]]
  M <- matrix(NA_real_, length(gr), length(bws), dimnames = list(gr$symbol, namen))
  for (i in seq_along(bws)) {
    bf <- BigWigFile(bws[i])
    ok <- as.character(seqnames(gr)) %in% seqnames(seqinfo(bf))
    g2 <- gr[ok]
    sl <- seqlengths(seqinfo(bf))[as.character(seqnames(g2))]
    end(g2) <- pmin(end(g2), as.integer(sl))
    v <- unlist(summary(bf, g2, type = "mean", defaultValue = 0))
    M[g2$symbol, i] <- v$score
    cat(sprintf("  [%s] %-14s Median %.4f\n", fn, namen[i],
                median(v$score, na.rm = TRUE)))
    rm(v); invisible(gc())
  }
  M <- M[rowSums(is.finite(M)) == ncol(M), , drop = FALSE]
  med <- apply(M, 2, median, na.rm = TRUE)
  M <- sweep(M, 2, med / mean(med), "/")
  write.csv(data.frame(symbol = rownames(M), M, check.names = FALSE),
            ziel, row.names = FALSE)
  cat(sprintf("window %s written: %d genes\n", fn, nrow(M)))
}
cat("done.\n")
