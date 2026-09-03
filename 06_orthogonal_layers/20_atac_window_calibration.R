# =============================================================================
# 20_atac_window_calibration.R -- layer B, step 2: calibrate the window BEFORE
#                             the module is tested.
# =============================================================================
# 21_atac_moduletest.py finds no module direction at the promoter window -- but
# the positive control of the layer fails there as well: canonical osteogenic
# markers (RUNX2, SP7, ALPL, ...) do not separate from adipogenic ones at the
# promoter. A null finding on a layer that fails its own positive control says
# nothing about the module, only about the window.
#
# This script measures the same ATAC signal in four windows:
#     P    promoter         TSS -2000 / +500
#     T10  TSS +- 10 kb
#     T50  TSS +- 50 kb     (regulatory neighborhood, enhancers)
#     GB   gene body        start to end of the gene
# and passes all four to the calibration in 23_. Only a window that passes the
# positive control may carry the module test. The order is binding and stands
# before any module result.
#
# Technically: `summary(BigWigFile, ...)` reads across zoom levels and is the
# only workable route for wide windows.
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

cat("Layer B -- window calibration GSE332758\n"); cat(strrep("=", 70), "\n")

gtf <- import(file.path(DAT, "_referenz", "gencode.v46.annotation.hg38.gtf.gz"),
              format = "gtf", feature.type = "gene")
gtf <- gtf[gtf$gene_type == "protein_coding"]
gtf <- gtf[!duplicated(gtf$gene_name)]
neg <- as.character(strand(gtf)) == "-"
tss <- ifelse(neg, end(gtf), start(gtf))

mk <- function(a, b) {
  # a bp upstream, b bp downstream of the TSS, strand-aware
  s <- ifelse(neg, tss - b, tss - a)
  e <- ifelse(neg, tss + a, tss + b)
  g <- GRanges(seqnames(gtf), IRanges(pmax(1L, s), e))
  g$symbol <- gtf$gene_name
  g
}

FENSTER <- list(
  P   = mk(2000L, 500L),
  T10 = mk(10000L, 10000L),
  T50 = mk(50000L, 50000L)
)
gb <- GRanges(seqnames(gtf), IRanges(start(gtf), end(gtf)))
gb$symbol <- gtf$gene_name
FENSTER$GB <- gb

for (n in names(FENSTER)) {
  cat(sprintf("window %-4s: %d regions, median width %.0f bp\n",
              n, length(FENSTER[[n]]), median(width(FENSTER[[n]]))))
}

bws <- sort(list.files(file.path(DAT, "GSE332758", "bw"), pattern = "\\.bw$",
                       full.names = TRUE))
namen <- sub("^GSM\\d+_", "", sub("\\.bw$", "", basename(bws)))

for (fn in names(FENSTER)) {
  ziel <- file.path(AUS, sprintf("B_atac_matrix_%s.csv", fn))
  if (file.exists(ziel)) { cat("exists:", basename(ziel), "\n"); next }
  gr <- FENSTER[[fn]]
  M <- matrix(NA_real_, length(gr), length(bws), dimnames = list(gr$symbol, namen))
  for (i in seq_along(bws)) {
    bf <- BigWigFile(bws[i])
    ok <- as.character(seqnames(gr)) %in% seqnames(seqinfo(bf))
    g2 <- gr[ok]
    # clip to contig length, otherwise summary() refuses
    sl <- seqlengths(seqinfo(bf))[as.character(seqnames(g2))]
    end(g2) <- pmin(end(g2), as.integer(sl))
    v <- unlist(summary(bf, g2, type = "mean", defaultValue = 0))
    M[g2$symbol, i] <- v$score
    cat(sprintf("  [%s] %-8s Median %.3f\n", fn, namen[i], median(v$score, na.rm = TRUE)))
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
