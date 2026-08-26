# =============================================================================
# 20i_dexamethasone.R -- EXPLORATORY. The reservation raised in the
#                        exploration protocol: is the convergence only the
#                        effect of the medium?
# =============================================================================
# Osteogenic and chondrogenic differentiation media both contain
# dexamethasone, and among the 129 convergently up-running genes stand TSC22D3
# (GILZ), ZBTB16 (PLZF), FOXO1/3, IRS2, GADD45B and SAA1 -- canonical
# glucocorticoid targets. If the convergence is only that, it carries no paper.
#
# THE TEST RUNS IN THREE STEPS:
#  (1) How strongly are the convergent genes enriched for an EXTERNAL
#      glucocorticoid signature at all?
#  (2) THE DECISIVE QUESTION: does the matrix enrichment (OR 3.40, p 0.0038)
#      hold if all glucocorticoid-responsive genes are removed BEFOREHAND from
#      the foreground AND the background?
#  (3) Does the convergence itself survive the removal of the GC genes -- that
#      is, is the count of 173 against 7.9 carried by them?
#
# THE SOURCE, external and version-fixed: MSigDB 2026.1.Hs through msigdbr
# 26.1.0, retrieved 2026-08-19. No fresh GEO download -- a curated, versioned
# signature is more traceable than a data set analysed by ourselves, and it
# sidesteps the question of which cell line and which dose would be "right".
#   RHEIN_ALL_GLUCOCORTICOID_THERAPY_UP   88 genes  (in vivo, patient therapy)
#   RHEIN_ALL_GLUCOCORTICOID_THERAPY_DN  423 genes
#   WP_GLUCOCORTICOID_RECEPTOR_PATHWAY    70 genes
#   GO:0051384 response to glucocorticoid           (org.Hs.eg.db, as before)
# The union of these four is the exclusion set in (2) and (3) -- deliberately
# cut GENEROUSLY, because too small an exclusion set would not dispel the
# reservation.
# =============================================================================
suppressMessages({ library(msigdbr); library(org.Hs.eg.db) })
# The session tree is not part of the public archive; set PAPER_V2_SESSIONS.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
ERG <- file.path(SESSIONS, "20_Exploration", "derived_data")

g2e <- as.list(org.Hs.egGO2ALLEGS); e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) { eg <- g2e[[id]]; if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE))) }

MS  <- msigdbr(species = "Homo sapiens", collection = "C2")
MSV <- paste(unique(MS$db_version), collapse = "/")
ms  <- function(nm) unique(na.omit(MS$ensembl_gene[MS$gs_name == nm]))
GC_UP <- ms("RHEIN_ALL_GLUCOCORTICOID_THERAPY_UP")
GC_DN <- ms("RHEIN_ALL_GLUCOCORTICOID_THERAPY_DN")
GC_WP <- ms("WP_GLUCOCORTICOID_RECEPTOR_PATHWAY")
GC_GO <- GOr("GO:0051384")
GC    <- unique(c(GC_UP, GC_DN, GC_WP, GC_GO))
ECM   <- GOr("GO:0030198")

K  <- read.csv(file.path(ERG, "20f_konvergente_dWT.csv"), stringsAsFactors = FALSE)
HG <- read.csv(file.path(ERG, "20f_hintergrund.csv"), stringsAsFactors = FALSE)$gen
UP <- K$gen[K$ri > 0]; DN <- K$gen[K$ri < 0]

cat("################################################################\n")
cat(sprintf("# 20i_dexamethasone.R | MSigDB %s | msigdbr %s | retrieved 2026-08-19\n",
            MSV, as.character(packageVersion("msigdbr"))))
cat(sprintf("# GC exclusion set: UP %d + DN %d + WP %d + GO %d -> union %d\n",
            length(GC_UP), length(GC_DN), length(GC_WP), length(GC_GO), length(GC)))
cat(sprintf("# of these in the background (%d genes): %d\n", length(HG),
            length(intersect(GC, HG))))
cat("################################################################\n\n")

ft <- function(menge, ziel, hg, nm) {
  menge <- intersect(menge, hg); ziel <- intersect(ziel, hg)
  a <- length(intersect(menge, ziel)); b <- length(menge) - a
  c <- length(ziel) - a; d <- length(hg) - length(menge) - c
  f <- fisher.test(matrix(c(a, b, c, d), 2))
  cat(sprintf("  %-44s %3d/%4d (%5.2f%%)  HG %5.2f%%  OR %5.2f [%.2f-%.2f]  p %.4g\n",
      nm, a, length(menge), 100*a/max(length(menge),1),
      100*length(ziel)/length(hg), f$estimate, f$conf.int[1], f$conf.int[2],
      f$p.value))
  invisible(c(OR = unname(f$estimate), p = f$p.value, k = a, n = length(menge)))
}

cat("--- (1) Are the convergent genes glucocorticoid-responsive? -----------\n")
ft(UP, GC_UP, HG, "convergent UP    x  RHEIN GC therapy UP")
ft(UP, GC,    HG, "convergent UP    x  GC overall")
ft(DN, GC,    HG, "convergent DOWN  x  GC overall")
cat("\n  Which of the 129 they are, by name:\n")
sy <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys = intersect(UP, GC), keytype = "ENSEMBL", columns = "SYMBOL"))
cat("   ", paste(sort(unique(na.omit(sy$SYMBOL))), collapse = ", "), "\n")

cat("\n--- (2) THE DECISIVE TEST ---------------------------------------------\n")
cat("Matrix enrichment BEFORE and AFTER removing all GC-responsive genes\n")
cat("(from the foreground AND the background, or the test would be biased):\n\n")
a1 <- ft(UP, ECM, HG, "ECM organisation, all genes")
HG2 <- setdiff(HG, GC); UP2 <- setdiff(UP, GC)
cat(sprintf("\n  Background %d -> %d genes | convergent UP %d -> %d genes\n\n",
            length(HG), length(HG2), length(UP), length(UP2)))
a2 <- ft(UP2, ECM, HG2, "ECM organisation, WITHOUT GC genes")

cat("\n--- (3) Does GC carry the convergence itself? -------------------------\n")
G <- do.call(rbind, lapply(list.files(ERG, "^20d_gene_.*csv$", full.names = TRUE),
                           read.csv, stringsAsFactors = FALSE))
mat <- function(sp) {
  M <- reshape(G[, c("gen", "punkt", sp)], idvar = "gen", timevar = "punkt",
               direction = "wide")
  rownames(M) <- M$gen; M <- as.matrix(M[, -1]); M[, order(colnames(M))]
}
DW <- mat("dWT")
DW <- DW[rowSums(!is.na(DW)) >= 14, , drop = FALSE]
zent <- sweep(DW, 2, apply(DW, 2, median, na.rm = TRUE), "-")
S <- sign(zent); S[is.na(S)] <- 0
kons <- function(rows) {
  s <- S[rows, , drop = FALSE]
  p <- rowSums(s > 0); n <- rowSums(s != 0); v <- pmax(p, n - p)
  sum(v / n >= 0.90) }
alle <- rownames(DW); ohne <- setdiff(alle, GC)
set.seed(20260819)
rausch <- replicate(500, {
  Sp <- S %*% diag(sample(c(-1, 1), ncol(S), TRUE))
  p <- rowSums(Sp > 0); n <- rowSums(Sp != 0); v <- pmax(p, n - p)
  sum((v / n >= 0.90)[match(ohne, alle)]) })
cat(sprintf("  convergent (>= 90%%) with all %d genes    : %d\n",
            length(alle), kons(alle)))
cat(sprintf("  convergent (>= 90%%) WITHOUT the %d GC genes: %d  (noise %.1f)\n",
            length(alle) - length(ohne), kons(ohne), mean(rausch)))

cat("\n################################################################\n")
cat("VERDICT: the reservation is dispelled if (2) still carries after\n")
cat("the removal and (3) loses only a few genes.\n")
