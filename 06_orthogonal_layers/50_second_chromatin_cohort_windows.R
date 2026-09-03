# =============================================================================
# 50_second_chromatin_cohort_windows.R -- window matrices of the SECOND chromatin cohort
# =============================================================================
# Purpose  GSE332758 cannot carry the lineage contrast: its osteogenic axis
#          fails its own calibration in every window, and on the difference
#          axis the osteogenic markers do not move at all (z -0.47 to +1.66)
#          while the adipogenic ones move correctly (z -3.67 in the promoter
#          window). A null result on such an axis means "not measurable".
#
#          This script fetches the level from a second, independent source:
#          GSE151311 (ATAC) and GSE151315 (H3K27ac), both from one study of
#          chromatin rearrangement in adipogenesis and osteogenesis:
#
#            GSE151311  ATAC-seq   AC rep1/2, OB rep1/2    -- no undifferentiated arm
#            GSE151315  H3K27ac    hMSC, AC, OB, rep1/2    -- WITH an undifferentiated arm
#
#          GSE151315 is the real gain: its undifferentiated arm makes both
#          differentiation axes calibratable (calibration D), and the lineage
#          contrast OB minus AC gets its own calibration L.
#
# Genome   Both cohorts are in **hg19**. The windows are built from GENCODE
#          v46 (hg38) exactly as in the reference implementation
#          `06_orthogonal_layers/20_atac_window_calibration.R` and then
#          lifted to hg19 point by point. Transcription start sites and gene
#          body boundaries are lifted separately; a gene is dropped when its
#          parts land on different contigs or when the gene body grows by
#          more than threefold.
#
# Inputs   data_raw/GSE151311/bw/*.bw, data_raw/GSE151315/bw/*.bw
#          data_raw/_referenz/gencode.v46.annotation.hg38.gtf.gz
#          data_raw/_referenz/hg38ToHg19.over.chain.gz
# Outputs  derived_data/B_atac2/<GSE>_matrix_<window>.csv
#          derived_data/B_atac2/22_zweitkohorte_log.txt
# Runtime  about 15 minutes (10 BigWigs x 4 windows)
# =============================================================================
suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
DAT <- file.path(WURZEL, "data_raw")
AUS <- file.path(WURZEL, "derived_data", "B_atac2")
dir.create(AUS, recursive = TRUE, showWarnings = FALSE)

LOG <- character(0)
log <- function(...) {
  s <- sprintf(...)
  cat(s, "\n")
  LOG <<- c(LOG, s)
}

log("=== Second cohort: window matrices (hg19) ===")

# ------------------------------------------------------- Annotation (hg38)
gtf <- import(file.path(DAT, "_referenz", "gencode.v46.annotation.hg38.gtf.gz"),
              format = "gtf", feature.type = "gene")
gtf <- gtf[gtf$gene_type == "protein_coding"]
gtf <- gtf[!duplicated(gtf$gene_name)]
neg <- as.character(strand(gtf)) == "-"
tss38 <- ifelse(neg, end(gtf), start(gtf))
log("GENCODE v46: %d protein-coding genes", length(gtf))

# ------------------------------------------------------- point-wise to hg19
# import.chain() wants a file path, not a connection, so the chain file is
# unpacked into results/ once.
chain_gz <- file.path(DAT, "_referenz", "hg38ToHg19.over.chain.gz")
chain_txt <- file.path(WURZEL, "results", "hg38ToHg19.over.chain")
if (!file.exists(chain_txt)) {
  dir.create(dirname(chain_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(readLines(gzfile(chain_gz)), chain_txt)
}
chain <- import.chain(chain_txt)

hebe <- function(chrom, pos, sym) {
  # Lift one point per gene; only unique hits are kept.
  g <- GRanges(chrom, IRanges(pos, pos))
  g$symbol <- sym
  l <- liftOver(g, chain)
  n <- lengths(l)
  ok <- n == 1
  out <- unlist(l[ok])
  data.frame(symbol = sym[ok], chr = as.character(seqnames(out)),
             pos = start(out), stringsAsFactors = FALSE)
}

ch <- as.character(seqnames(gtf))
t19 <- hebe(ch, tss38, gtf$gene_name)
s19 <- hebe(ch, start(gtf), gtf$gene_name)
e19 <- hebe(ch, end(gtf), gtf$gene_name)
log("lifted to hg19: TSS %d, gene start %d, gene end %d",
    nrow(t19), nrow(s19), nrow(e19))

strang <- setNames(ifelse(neg, "-", "+"), gtf$gene_name)
laenge38 <- setNames(width(gtf), gtf$gene_name)

# windows anchored on the transcription start site
tt <- t19
tt$neg <- strang[tt$symbol] == "-"
mk <- function(a, b) {
  s <- ifelse(tt$neg, tt$pos - b, tt$pos - a)
  e <- ifelse(tt$neg, tt$pos + a, tt$pos + b)
  g <- GRanges(tt$chr, IRanges(pmax(1L, as.integer(s)), as.integer(e)))
  g$symbol <- tt$symbol
  g
}

# Gene body: start and end must land on the same contig, and the lifted body
# must not grow by more than threefold.
gb <- merge(s19, e19, by = "symbol", suffixes = c("_s", "_e"))
gb <- gb[gb$chr_s == gb$chr_e, ]
lo <- pmin(gb$pos_s, gb$pos_e); hi <- pmax(gb$pos_s, gb$pos_e)
plaus <- (hi - lo + 1) <= 3 * laenge38[gb$symbol]
gb <- gb[plaus, ]; lo <- lo[plaus]; hi <- hi[plaus]
grb <- GRanges(gb$chr_s, IRanges(lo, hi)); grb$symbol <- gb$symbol
log("gene bodies hg19 plausible: %d of %d", length(grb), nrow(s19))

FENSTER <- list(P = mk(2000L, 500L), T10 = mk(10000L, 10000L),
                T50 = mk(50000L, 50000L), GB = grb)
for (n in names(FENSTER))
  log("window %-4s: %d regions, median width %.0f bp",
      n, length(FENSTER[[n]]), median(width(FENSTER[[n]])))

# ------------------------------------------------------------------ read-out
# Word for word as in the reference implementation: mean per region across
# the zoom levels, genes without a complete row are dropped, column-median
# normalisation. Nothing here is invented for this cohort.
kohorten <- list(GSE151311 = file.path(DAT, "GSE151311", "bw"),
                 GSE151315 = file.path(DAT, "GSE151315", "bw"))

for (ko in names(kohorten)) {
  bws <- sort(list.files(kohorten[[ko]], pattern = "\\.bw$", full.names = TRUE))
  namen <- sub("\\.bw$", "", basename(bws))
  log("\n--- %s: %d BigWigs (%s)", ko, length(bws), paste(namen, collapse = ", "))
  if (!length(bws)) { log("  no files -- skipped"); next }
  for (fn in names(FENSTER)) {
    ziel <- file.path(AUS, sprintf("%s_matrix_%s.csv", ko, fn))
    if (file.exists(ziel)) { log("  present: %s", basename(ziel)); next }
    gr <- FENSTER[[fn]]
    M <- matrix(NA_real_, length(gr), length(bws),
                dimnames = list(gr$symbol, namen))
    for (i in seq_along(bws)) {
      bf <- BigWigFile(bws[i])
      ok <- as.character(seqnames(gr)) %in% seqnames(seqinfo(bf))
      g2 <- gr[ok]
      sl <- seqlengths(seqinfo(bf))[as.character(seqnames(g2))]
      end(g2) <- pmin(end(g2), as.integer(sl))
      g2 <- g2[start(g2) <= end(g2)]
      v <- unlist(summary(bf, g2, type = "mean", defaultValue = 0))
      M[g2$symbol, i] <- v$score
      log("  [%s] %-16s Median %.4f", fn, namen[i],
          median(v$score, na.rm = TRUE))
      rm(v); invisible(gc())
    }
    M <- M[rowSums(is.finite(M)) == ncol(M), , drop = FALSE]
    med <- apply(M, 2, median, na.rm = TRUE)
    M <- sweep(M, 2, med / mean(med), "/")
    write.csv(data.frame(symbol = rownames(M), M, check.names = FALSE),
              ziel, row.names = FALSE)
    log("  -> %s: %d genes", basename(ziel), nrow(M))
  }
}

writeLines(LOG, file.path(AUS, "22_zweitkohorte_log.txt"))
capture.output(sessionInfo(),
               file = file.path(WURZEL, "results",
                                "zweitkohorte_sessionInfo.txt"))
cat("done ->", AUS, "\n")
