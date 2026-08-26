# =============================================================================
# 03a_data.R -- Brings the four newly loaded GEO matrices into the format the
#               reference implementation expects (a matrix plus meta).
# =============================================================================
# Preregistration: PREREG_whole_study.md; the screen: the screening record.
# NOTHING is filtered that does not stand in the screening record. The
# assignment of genotype and condition follows the screening entry, not the
# result.
#
# THE PROJECT CONVENTION, taken over from the data sets already included: if
# several differentiation time points exist, ALL of them enter as "diff" and
# the undifferentiated time point as "undiff". GSE245585 (day 0/7/14/21) and
# GSE227512 (day 0/7/14/21) are thereby handled without a selection -- no time
# point is chosen.
#
# A DEVIATION FROM section 2.1 of the preregistration, recorded in advance in
# the screening record: for GSE218101 (CPM), GSE221128 (TPM), GSE205432 and
# GSE245585 (normalised counts) neither a raw count matrix nor a
# variance-stabilised matrix is publicly available. log2(x+1) is used.
# =============================================================================
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi)
                   library(matrixStats) })

# The raw GEO matrices live in the tree of raw analysis sessions, which is not
# part of the public archive; set PAPER_V2_SESSIONS to point at it.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
DAT <- file.path(SESSIONS, "03_Metrik_Elf_Punkte", "data_raw")

ent <- function(x) sub("\\..*$", "", x)
entdup <- function(X) X[!duplicated(rownames(X)) & rownames(X) != "" &
                        !is.na(rownames(X)), , drop = FALSE]

# --------------------------------------------------------------- GSE218101
# ARSB (MPS VI), chondrogenic. 4 patient lines, each empty vector (EV =
# ARSB-deficient = KO) against gene-corrected (GE = WT). iPS = undiff,
# D14 = diff. The two "Healthy control" columns do not belong to the isogenic
# 2x2 and exist only at day 0 -- they are not used.
gse218101 <- function() {
  D <- read.csv(gzfile(file.path(DAT, "GSE218101_CPM.csv.gz")),
                check.names = FALSE, stringsAsFactors = FALSE)
  # Symbols occur more than once; the first entry wins, as in vorb()
  D <- D[!duplicated(D[[1]]), , drop = FALSE]
  X <- D[, -1, drop = FALSE]; rownames(X) <- D[[1]]
  X <- X[, !grepl("^Healthy", colnames(X)), drop = FALSE]
  # Symbols -> Ensembl (the matrix is deposited at symbol level)
  s2e <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(X),
           keytype = "SYMBOL", columns = "ENSEMBL")
  s2e <- s2e[!is.na(s2e$ENSEMBL) & !duplicated(s2e$SYMBOL), ]
  X <- as.matrix(X[s2e$SYMBOL, , drop = FALSE])
  rownames(X) <- s2e$ENSEMBL          # a matrix, so duplicates are allowed
  X <- entdup(X)                      # the first entry wins, as in vorb()
  sp <- colnames(X)
  meta <- data.frame(sample = sp,
    genotype  = ifelse(grepl("_EV", sp), "KO", "WT"),
    condition = ifelse(grepl("_iPS$", sp), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}

# --------------------------------------------------------------- GSE221128
# ACVR1 (FOP), chondrogenic. FOP = KO, resFOP (corrected) = WT.
# Day 0 = undiff (3/3), day 6 = diff (FOP 2, resFOP 3).
gse221128 <- function() {
  D <- read.delim(gzfile(file.path(DAT, "GSE221128_TPM_allsamples.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(D[, -(1:2), drop = FALSE])
  rownames(X) <- ent(D$gene_id); X <- entdup(X)
  sp <- colnames(X)
  meta <- data.frame(sample = sp,
    genotype  = ifelse(grepl("_FOP_", sp), "KO", "WT"),
    condition = ifelse(grepl("Day0", sp), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}

# --------------------------------------------------------------- GSE205432
# RNF4 siRNA against scrambled, osteogenic. untreated = undiff, 7d = diff.
gse205432 <- function() {
  D <- read.csv(file.path(DAT, "GSE205432_norm.csv"), row.names = 1,
                check.names = FALSE)
  X <- entdup(as.matrix(D)); rownames(X) <- ent(rownames(X)); X <- entdup(X)
  sp <- colnames(X)
  meta <- data.frame(sample = sp,
    genotype  = ifelse(grepl("^siRNAaRNF4", sp), "KO", "WT"),
    condition = ifelse(grepl("untreated", sp), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}

# --------------------------------------------------------------- GSE245585
# RB1 (RB+/-), osteogenic. One file per sample, joined here into one matrix.
# The assignment comes from the series matrix: file prefix U = WT1,
# B = MT1 (RB+/-MUT). Day 0 = undiff (2/2), day 7/14/21 = diff (6/6).
gse245585 <- function() {
  fs <- sort(list.files(file.path(DAT, "GSE245585"), full.names = TRUE))
  L <- lapply(fs, function(f) {
    d <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    setNames(d[[2]], ent(d[[1]])) })
  nm <- sub("^GSM[0-9]+_GT_SO_10449_", "", sub("\\.txt\\.gz$", "", basename(fs)))
  gemein <- Reduce(intersect, lapply(L, names))
  X <- do.call(cbind, lapply(L, function(v) v[gemein]))
  rownames(X) <- gemein; colnames(X) <- nm
  X <- entdup(X)
  meta <- data.frame(sample = nm,
    genotype  = ifelse(grepl("^B_", nm), "KO", "WT"),
    condition = ifelse(grepl("_D0_", nm), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}

NEU <- list(GSE218101 = gse218101, GSE221128 = gse221128,
            GSE205432 = gse205432, GSE245585 = gse245585)
