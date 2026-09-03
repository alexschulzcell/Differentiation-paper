# =============================================================================
# 14_geo_matrices_s5_format.R -- Brings the GEO matrices newly included in S5 step B into
#                  the format of the reference implementation.
# =============================================================================
# The screen was completed BEFORE any computation on these data sets.
# Preregistration: PREREG_S5.md section 4, with the criteria of
# PREREG_whole_study.md 3.1 and 3.2 unchanged.
#
# The same PROJECT CONVENTION applies as in 13_geo_matrices_to_metric_format.R: if several
# differentiation time points exist, ALL of them enter as "diff" and the
# undifferentiated time point as "undiff". No time point is chosen.
#
# A DEVIATION FROM section 2.1 of the preregistration, recorded in advance here
# as in 13_geo_matrices_to_metric_format.R: for GSE102732 (RPKM, log2), GSE145235 (FPKM), GSE137035
# (log-CPM) and GSE226565 (counts, partly non-integer) neither a raw count
# matrix in integer form NOR a variance-stabilised matrix is publicly
# available. log2(x+1) is used, or the already logarithmised matrix unchanged.
# GSE190542 and GSE247528 deliver real counts and run through rlog.
#
# DECISIONS of this script, all made before the computation and recorded in the
# protocol:
#  (1) GSE226565: only the classes "Nonunion" (KO) and "Healthy" (WT). The
#      third group "Callus" is a state of healing, not of failure, and as a
#      second point would share the same control arm; a second, dependent point
#      would violate the independence assumption of the sign test. Chosen by
#      substance, not by result.
#  (2) GSE190542: only "Expansion MEM-a media" (undiff) against "Osteogenic
#      media" (diff), and only cell type MSC. Adipogenic, granulocytic and
#      erythroid are other lineages and not the differentiation of this study.
#  (3) GSE102732: only the arms WT and LFS. The arms shSFRP2 and sFRP2OE have
#      no non-targeting control arm of their own; setting them against the
#      parental WT line would be exactly the clone confounding documented in
#      the record of fallen hypotheses. Exclusion A2 for these sub-arms, noted
#      in the protocol.
#  (4) GSE145235: two points, RBM against cRBM and RBD against cRBD. Both are
#      ISOGENIC (a gene-corrected line as the control) -- the shared WT arm is
#      therefore not needed and not used. Two mutations of the same GSE count
#      as ONE study under section 5.2.
#  (5) GSE247528: transcript-wise kallisto abundances. est_counts are summed
#      per gene (the standard transcript-to-gene route). That is preparation,
#      not a second implementation of the metric.
# =============================================================================
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi)
                   library(matrixStats) })

# The session tree is not part of the public archive; set PAPER_V2_SESSIONS.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
DAT_S5 <- file.path(SESSIONS, "03_Metrik_Elf_Punkte", "data_raw")

ent5    <- function(x) sub("\\..*$", "", x)
entdup5 <- function(X) X[!duplicated(rownames(X)) & rownames(X) != "" &
                         !is.na(rownames(X)), , drop = FALSE]

# Symbol -> Ensembl, word-for-word procedure as in 13_geo_matrices_to_metric_format.R (gse218101)
sym2ens5 <- function(X) {
  s2e <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
           keys = rownames(X), keytype = "SYMBOL", columns = "ENSEMBL"))
  s2e <- s2e[!is.na(s2e$ENSEMBL) & !duplicated(s2e$SYMBOL), ]
  X <- as.matrix(X[s2e$SYMBOL, , drop = FALSE])
  rownames(X) <- s2e$ENSEMBL
  entdup5(X)
}

# Read the SOFT header: returns title and description per GSM
soft5 <- function(datei) {
  z <- readLines(file.path(DAT_S5, datei), warn = FALSE)
  tit <- sub("^!Sample_title = ", "", grep("^!Sample_title = ", z, value = TRUE))
  acc <- sub("^!Sample_geo_accession = ", "",
             grep("^!Sample_geo_accession = ", z, value = TRUE))
  # collect the descriptions per sample (there can be several per sample)
  idx <- grep("^!Sample_title = ", z)
  ende <- c(idx[-1] - 1, length(z))
  des <- vapply(seq_along(idx), function(i) {
    blk <- z[idx[i]:ende[i]]
    d <- sub("^!Sample_description = ", "",
             grep("^!Sample_description = ", blk, value = TRUE))
    if (!length(d)) NA_character_ else d[length(d)]   # last one = column key
  }, character(1))
  data.frame(gsm = acc, titel = tit, beschreibung = des, stringsAsFactors = FALSE)
}

NEU_S5 <- list()

# --------------------------------------------------------------- GSE226565
# Osteoprogenitors, nonunion (KO) against healthy (WT). Day 0 = undiff,
# day 7/14/21/28 = diff. The matrix column names stand in Sample_description.
NEU_S5$GSE226565 <- function() {
  D <- read.delim(gzfile(file.path(DAT_S5, "GSE226565_raw_counts.tsv.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  S <- soft5("GSE226565_soft.txt")
  S$typ  <- ifelse(grepl("_NU_", S$titel), "KO",
             ifelse(grepl("_HB_", S$titel), "WT", NA_character_))
  S$tag  <- as.integer(sub(".*Day ", "", S$titel))
  S$dono <- sub("_.*", "", S$titel)
  S <- S[!is.na(S$typ) & S$beschreibung %in% colnames(D), ]
  X <- as.matrix(D[, S$beschreibung, drop = FALSE])
  colnames(X) <- S$titel
  X <- sym2ens5(X)
  meta <- data.frame(sample = S$titel, genotype = S$typ,
                     condition = ifelse(S$tag == 0, "undiff", "diff"),
                     block = S$dono, stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}

# --------------------------------------------------------------- GSE247528
# SERPINA3-KD, osteogenic. Transcript-wise kallisto abundances per sample;
# est_counts summed per gene. Day 0 = undiff, day 3/7 = diff.
NEU_S5$GSE247528 <- function() {
  fs <- list.files(file.path(DAT_S5, "GSE247528"), pattern = "abundances.txt.gz$",
                   full.names = TRUE)
  stopifnot(length(fs) == 18)
  nm <- sub("^GSM[0-9]+_", "", sub("\\.abundances\\.txt\\.gz$", "", basename(fs)))
  L <- lapply(fs, function(f) {
    d <- read.delim(gzfile(f), stringsAsFactors = FALSE)
    setNames(d$est_counts, ent5(d$target_id)) })
  tx <- Reduce(union, lapply(L, names))
  M <- vapply(L, function(v) { o <- setNames(rep(0, length(tx)), tx)
                               o[names(v)] <- v; o }, numeric(length(tx)))
  colnames(M) <- nm
  t2g <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
           keys = rownames(M), keytype = "ENSEMBLTRANS", columns = "ENSEMBL"))
  t2g <- t2g[!is.na(t2g$ENSEMBL) & !duplicated(t2g$ENSEMBLTRANS), ]
  M <- M[t2g$ENSEMBLTRANS, , drop = FALSE]
  X <- rowsum(M, group = t2g$ENSEMBL)          # transcript -> gene
  meta <- data.frame(sample = nm,
    genotype  = ifelse(grepl("_KD_", nm), "KO", "WT"),
    condition = ifelse(grepl("^Day0_", nm), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = X, meta = meta, roh = TRUE)
}

# --------------------------------------------------------------- GSE190542
# ERCC6L2 knockdown against scrambled, MSC, osteogenic.
# Expansion MEM-a = undiff, Osteogenic media = diff. n = 2 per cell.
NEU_S5$GSE190542 <- function() {
  fs <- list.files(file.path(DAT_S5, "GSE190542"), pattern = "Counts.txt.gz$",
                   full.names = TRUE)
  gsm <- sub("_.*", "", basename(fs))
  S <- soft5("GSE190542_soft.txt")
  S <- S[match(gsm, S$gsm), ]
  keep <- grepl("^Expansion MEM-a media|^Osteogenic media", S$titel) &
          grepl("Scramble|Knockdown", S$titel)
  fs <- fs[keep]; S <- S[keep, ]
  stopifnot(nrow(S) == 8)
  L <- lapply(fs, function(f) { d <- read.delim(gzfile(f), header = FALSE,
                                  stringsAsFactors = FALSE)
                                setNames(d$V2, ent5(d$V1)) })
  g <- Reduce(intersect, lapply(L, names))
  X <- vapply(L, function(v) v[g], numeric(length(g)))
  rownames(X) <- g; colnames(X) <- S$titel
  meta <- data.frame(sample = S$titel,
    genotype  = ifelse(grepl("Knockdown", S$titel), "KO", "WT"),
    condition = ifelse(grepl("^Expansion", S$titel), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = X, meta = meta, roh = TRUE)
}

# --------------------------------------------------------------- GSE102732
# Li-Fraumeni (TP53 G245D) against WT, osteogenic. D0 MSCs = undiff,
# D7/D14/D17 pre-osteoblasts = diff. Only the WT and LFS arms (see (3)).
NEU_S5$GSE102732 <- function() {
  D <- read.delim(gzfile(file.path(DAT_S5, "GSE102732_RPKM_log2.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  S <- soft5("GSE102732_soft.txt")
  S <- S[grepl("^(WT|LFS) D", S$titel), ]
  S <- S[S$beschreibung %in% colnames(D), ]
  stopifnot(nrow(S) == 16)          # 2 arms x 4 time points x 2 replicates
  X <- as.matrix(D[, S$beschreibung, drop = FALSE])
  colnames(X) <- S$titel
  X <- sym2ens5(X)
  meta <- data.frame(sample = S$titel,
    genotype  = ifelse(grepl("^LFS", S$titel), "KO", "WT"),
    condition = ifelse(grepl(" D0 ", S$titel), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = X, meta = meta)          # already log2(RPKM)
}

# --------------------------------------------------------------- GSE137035
# YAP/TAZ siRNA against scrambled siRNA in hFOB. 0 nM BMP2 = undiff (not
# induced), 5 nM BMP2 = diff. n = 3 per cell. The matrix is log-CPM with
# comma as decimal separator.
NEU_S5$GSE137035 <- function() {
  D <- read.delim(gzfile(file.path(DAT_S5,
         "GSE137035_DESEQ_cpm_values_w_pvals.txt.gz")),
         check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
         dec = ",")
  S <- soft5("GSE137035_soft.txt")
  z <- readLines(file.path(DAT_S5, "GSE137035_soft.txt"), warn = FALSE)
  idx <- grep("^!Sample_title = ", z); ende <- c(idx[-1] - 1, length(z))
  merk <- function(schl) vapply(seq_along(idx), function(i) {
    blk <- z[idx[i]:ende[i]]
    v <- grep(paste0("^!Sample_characteristics_ch1 = ", schl, ": "), blk,
              value = TRUE)
    if (!length(v)) NA_character_ else sub(".*: ", "", v[1]) }, character(1))
  S$sirna <- merk("sirna"); S$bmp2 <- merk("treatment")
  S <- S[S$titel %in% colnames(D) & !is.na(S$sirna), ]
  stopifnot(nrow(S) == 12)
  X <- as.matrix(D[, S$titel, drop = FALSE])
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  X <- sym2ens5(X)
  meta <- data.frame(sample = S$titel,
    genotype  = ifelse(grepl("YAP", S$sirna), "KO", "WT"),
    condition = ifelse(grepl("^0nM", S$bmp2), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = X, meta = meta)          # already log-CPM
}

# ------------------------------------------------- GSE145235, two points
# RB1-mutant or RB1-deleted iPSC-MSC against the respectively GENE-CORRECTED
# line -- isogenic. D0 = undiff, D15/D24 = diff. n = 3 per cell.
lade145235 <- function(mut, ktrl) function() {
  D <- read.csv(gzfile(file.path(DAT_S5, "GSE145235_TJ_RBiPSC_RNAseq_fpkm.csv.gz")),
                check.names = FALSE, stringsAsFactors = FALSE)
  names(D)[1] <- "gene_id"
  X <- as.matrix(D[, -(1:2), drop = FALSE])
  rownames(X) <- D$gene_id_alias                    # symbols
  X <- entdup5(X)
  sp <- colnames(X)
  sel <- grepl(paste0("^", mut, "[0-9]_"), sp) | grepl(paste0("^", ktrl, "[0-9]_"), sp)
  X <- X[, sel, drop = FALSE]; sp <- colnames(X)
  X <- sym2ens5(X)
  meta <- data.frame(sample = sp,
    genotype  = ifelse(grepl(paste0("^", ktrl), sp), "WT", "KO"),
    condition = ifelse(grepl("_D0$", sp), "undiff", "diff"),
    stringsAsFactors = FALSE)
  list(X = log2(X + 1), meta = meta)
}
# Mind the order: "cRBM" contains "RBM" -- therefore the control arm is
# tested first and the mutant arm determined as the remainder.
NEU_S5$GSE145235_RBM <- lade145235("RBM", "cRBM")
NEU_S5$GSE145235_RBD <- lade145235("RBD", "cRBD")

# =============================================================================
# Metadata of the seven new points (class from 12_klassen.R, external)
# =============================================================================
NEUMETA_S5 <- list(
  GSE226565     = list(lab = "Nonunion (GSE226565)",     arm = "osteogen", kl = "N"),
  GSE247528     = list(lab = "SERPINA3-KD osteo (GSE247528)", arm = "osteogen", kl = "N"),
  GSE190542     = list(lab = "ERCC6L2-KD (GSE190542)",   arm = "osteogen", kl = "N"),
  GSE102732     = list(lab = "TP53 LFS (GSE102732)",     arm = "osteogen", kl = "M"),
  GSE137035     = list(lab = "YAP/TAZ-KD (GSE137035)",   arm = "osteogen", kl = "N"),
  GSE145235_RBM = list(lab = "RB1-mut isogen (GSE145235)", arm = "osteogen", kl = "N"),
  GSE145235_RBD = list(lab = "RB1-del isogen (GSE145235)", arm = "osteogen", kl = "N"))
