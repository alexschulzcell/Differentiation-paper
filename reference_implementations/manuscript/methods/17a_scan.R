# =============================================================================
# 17a_scan.R -- S6 step C: the scan in the DISCOVERY HALF, its empirical noise
#               arm, and the detection limit of the metric.
# =============================================================================
# Preregistration: preregistrations/PREREG_S6.md, dated 2026-08-19, before the
# first metric of this script. This script fixes NOTHING anew.
#
# Three parts, per point in this order, three separate CSV files:
#   (A) SCAN         4.1/4.3  K = 105 candidates against REF = S_MASCHINE, NB 2000
#   (B) NOISE        4.4 (b)  the same scan over NPERM = 10 permutations of the
#                             genotype labels within the condition, NB 500
#   (C) SENSITIVITY  4.7      neutral pairs with a built-in offset delta, NB 2000
#
# It touches only the nine discovery points. The validation half is not loaded
# (NUR_PUNKTE).
#
# DECISIONS of this script (section 12 of the S6 brief -- decide, give the
# reason, carry on):
#  (1) The candidate list is NOT rebuilt but read from
#      17_setcheck_kandidaten.csv (column testbar). K = 105 is thereby exactly
#      the same K that stands in PREREG_S6.md 4.3 -- a second construction
#      could deviate from it.
#  (2) The confounder checks do NOT run in the scan. Under 5.1 they are part of
#      step D and belong to the ONE validated hypothesis; computing them 105
#      times costs twenty times as much and carries no verdict.
#  (3) The noise arm permutes genotype WITHIN the condition. The induction part
#      dWT is thereby structurally preserved and only the interaction term iv is
#      destroyed -- that is the null the scan is tested against. The identity
#      permutation is rejected (up to 50 attempts); if that does not succeed,
#      the round is marked with ident = TRUE.
#  (4) NB is lowered to 500 ONLY in the noise arm (PREREG_S6.md 4.2). The value
#      is set before and after the arm and written into every row.
#  (5) The result is written PER POINT as its own CSV; existing partial files
#      are skipped. The work can be split through the first argument.
#
# The gene-set source: GO Biological Process through org.Hs.eg.db 3.20.0
# (GO2ALLEGS), retrieved 2026-08-19, through the same function GOr() as the
# reference implementation. No gene sets by hand.
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
PUNKTE <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else NULL

# The session tree is not part of the public archive; set PAPER_V2_SESSIONS.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
AUF <- read.csv(file.path(SESSIONS, "16_Pseudobulk", "derived_data",
                          "16_aufteilung.csv"), stringsAsFactors = FALSE)
ENT <- AUF$punkt[AUF$haelfte == "Entdeckung"]
stopifnot(length(ENT) == 9)
if (is.null(PUNKTE)) PUNKTE <- ENT
PUNKTE <- intersect(PUNKTE, ENT)          # never the validation half
stopifnot(length(PUNKTE) >= 1)
NUR_PUNKTE <- PUNKTE

set.seed(20260819)
HIER17 <- file.path(SESSIONS, "17_Scan")
source(file.path(SESSIONS, "13_Konvergenzachsen", "reference_implementations",
                 "13_load18.R"))
OUT17 <- file.path(HIER17, "derived_data")   # AFTER the source() -- a path trap
dir.create(OUT17, showWarnings = FALSE, recursive = TRUE)

OUTLOG <- character(0)
log <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); OUTLOG <<- c(OUTLOG, s) }
schreibe_log <- function() writeLines(OUTLOG,
  file.path(OUT17, sprintf("17a_log_%s.txt", paste(range(PUNKTE), collapse = "-"))))

KH      <- as.integer(KHAUPT)
NB_HAUPT <- 2000L
NB_PERM  <- 500L
NPERM    <- 10L
DELTA    <- c(0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00)
GROESSEN <- c(15L, 30L, 60L)
NRUND    <- 5L

# ------------------------------------------------------ the candidate list (1)
KAND <- read.csv(file.path(OUT17, "17_setcheck_kandidaten.csv"),
                 stringsAsFactors = FALSE)
KAND <- KAND[KAND$testbar %in% c(TRUE, "TRUE"), ]
K <- nrow(KAND)
stopifnot(K == 105)
REF <- S_MASCHINE
GENE_A <- lapply(KAND$go, function(id) setdiff(GOr(id), REF))
GENE_B <- lapply(KAND$go, function(id) setdiff(REF, GOr(id)))
names(GENE_A) <- names(GENE_B) <- KAND$go

# Jaccard against the sets already used (4.6 (c) 2) -- once, independent of point
S_STRESS_roh <- unique(c(S_UPR, S_ERAD))
S_STRESS     <- setdiff(S_STRESS_roh, S_MASCHINE)
S_ABBAU_roh  <- unique(c(GOr("GO:0006914"), GOr("GO:0005764"), GOr("GO:0000502")))
S_ABBAU      <- setdiff(S_ABBAU_roh, c(S_MASCHINE, S_FRACHT0))
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
ALT <- list(S_STRESS = S_STRESS, S_ABBAU = S_ABBAU, S_DISTAL = S_DISTAL,
            S_BIOSYN = S_BIOSYN, S_FRACHT0 = S_FRACHT0)
KAND$jaccard_max <- vapply(GENE_A, function(g)
  max(vapply(ALT, function(s) jac(g, s), numeric(1))), numeric(1))
write.csv(KAND[, c("go", "name", "groesse", "jaccard_max")],
          file.path(OUT17, "17a_kandidaten_jaccard.csv"), row.names = FALSE)

log("################################################################")
log("# 17a_scan.R -- S6 step C | PREREG_S6.md section 4, dated 2026-08-19")
log("# Seed 20260819 | NB scan %d | NB noise arm %d | NVIF %d | main null %s",
    NB_HAUPT, NB_PERM, NVIF, KHAUPT)
log("# K = %d candidates | REF = S_MASCHINE (%d genes) | GO BP org.Hs.eg.db %s",
    K, length(REF), as.character(packageVersion("org.Hs.eg.db")))
log("# NPERM %d | delta %s | sizes %s | rounds %d", NPERM,
    paste(DELTA, collapse = "/"), paste(GROESSEN, collapse = "/"), NRUND)
log("# discovery half: %s | to be computed: %s",
    paste(ENT, collapse = ","), paste(PUNKTE, collapse = ","))
log("################################################################")
schreibe_log()

t0 <- Sys.time()
for (ii in PUNKTE) {
  lab <- LABELS_18[ii]
  fA <- file.path(OUT17, sprintf("17a_scan_%02d.csv", ii))
  fB <- file.path(OUT17, sprintf("17b_perm_%02d.csv", ii))
  fC <- file.path(OUT17, sprintf("17c_empfind_%02d.csv", ii))
  if (is.null(DATEN[[lab]])) { log("\n[%d] %s: not loaded.", ii, lab); next }
  d <- DATEN[[lab]]
  K0 <- kern(d$Z, d$meta)
  iv <- K0$iv; pool <- K0$pool; rq <- K0$rho_quer; bas <- K0$basis
  Zi <- mk_zieh(K0$dWT, pool, KH)
  # d$expr is missing at some points (the convergence script checks the same).
  # d_niveau then stays NA; d_basis is untouched by that and carries 4.5.
  HAT_EXPR <- !is.null(d$expr) && length(intersect(names(d$expr), pool)) > 100
  sd_expr  <- if (HAT_EXPR) sd(d$expr[intersect(names(d$expr), pool)]) else NA_real_
  mexpr <- function(g) if (!HAT_EXPR) NA_real_ else
    mean(d$expr[intersect(g, names(d$expr))], na.rm = TRUE)
  log("\n============================================================")
  log("--- [%d] %s (%s, class %s) | pool %d", ii, lab, d$arm, d$klasse, length(pool))

  MR <- function(...) data.frame(punkt = ii, datensatz = lab, arm = d$arm,
                                 klasse = d$klasse, ..., stringsAsFactors = FALSE)

  # =========================================================== (A) the scan
  if (file.exists(fA)) log("  (A) scan: partial file present, skipped.") else {
    NB <<- NB_HAUPT
    RA <- vector("list", K); ta <- Sys.time()
    for (k in seq_len(K)) {
      gA <- intersect(GENE_A[[k]], pool); gB <- intersect(GENE_B[[k]], pool)
      set.seed(20260819L + 100L * ii + k)
      z <- kontrast_f(iv, Zi, gA, gB, pool, rq,
             MR(go = KAND$go[k], name = KAND$name[k], teil = "scan", nb = NB))
      if (is.null(z)) next
      # 4.5 the control on the baseline expression, per candidate and point
      z$d_basis  <- mean(bas[gA], na.rm = TRUE) - mean(bas[gB], na.rm = TRUE)
      z$d_niveau <- (mexpr(gA) - mexpr(gB)) / sd_expr
      RA[[k]] <- z
    }
    TA <- do.call(rbind, RA); write.csv(TA, fA, row.names = FALSE)
    log("  (A) scan: %d rows, %.1f min.", nrow(TA),
        as.numeric(difftime(Sys.time(), ta, units = "mins")))
    schreibe_log()
  }

  # ================================================ (B) the empirical noise arm
  # The noise arm writes its own file PER ROUND (17b_perm_<point>_<r>.csv) and
  # skips existing ones. The reason, for the record: the first run was
  # interrupted on 2026-08-19 at about 08:06 because the machine was closed,
  # and lost an almost complete series of ten in the process. Seeds are keyed
  # to (point, round), not to the loop order -- the result is therefore the
  # same whether a round falls in the first or in the second run.
  {
    NB <<- NB_PERM
    tb <- Sys.time()
    for (r in seq_len(NPERM)) {
      fBr <- file.path(OUT17, sprintf("17b_perm_%02d_%02d.csv", ii, r))
      if (file.exists(fBr)) { log("    perm %d/%d: present.", r, NPERM); next }
      RB <- list()
      set.seed(20260819L + 3000L * r + 100L * ii)
      mp <- d$meta; ident <- TRUE
      for (v in 1:50) {
        for (cc in unique(mp$condition)) {
          s <- which(mp$condition == cc)
          mp$genotype[s] <- sample(d$meta$genotype[s])
        }
        if (!identical(mp$genotype, d$meta$genotype)) { ident <- FALSE; break }
      }
      KP <- kern(d$Z, mp)
      if (length(KP$pool) < 1000) { log("    perm %d: pool %d < 1000, discarded.",
                                        r, length(KP$pool)); next }
      ZiP <- mk_zieh(KP$dWT, KP$pool, KH)
      for (k in seq_len(K)) {
        gA <- intersect(GENE_A[[k]], KP$pool); gB <- intersect(GENE_B[[k]], KP$pool)
        set.seed(20260819L + 7000L * r + 100L * ii + k)
        z <- kontrast_f(KP$iv, ZiP, gA, gB, KP$pool, KP$rho_quer,
               MR(go = KAND$go[k], name = KAND$name[k], teil = "perm", nb = NB,
                  runde = r, ident = ident))
        if (!is.null(z)) RB[[length(RB) + 1]] <- z
      }
      write.csv(do.call(rbind, RB), fBr, row.names = FALSE)
      log("    perm %d/%d done and written (%.1f min).", r, NPERM,
          as.numeric(difftime(Sys.time(), tb, units = "mins")))
      rm(KP, ZiP, RB); invisible(gc()); schreibe_log()
    }
    log("  (B) noise arm: %d rounds, %.1f min.", NPERM,
        as.numeric(difftime(Sys.time(), tb, units = "mins")))
    NB <<- NB_HAUPT
    schreibe_log()
  }

  # ============================================= (C) the detection limit (4.7)
  if (file.exists(fC)) log("  (C) sensitivity: partial file present, skipped.") else {
    NB <<- NB_HAUPT
    npool <- intersect(S_NEUTRAL, pool)
    gr <- unique(pmin(GROESSEN, floor(length(npool) / 2)))
    gr <- gr[gr >= 8]
    RC <- list(); tc <- Sys.time()
    if (!length(gr)) log("  (C) S_NEUTRAL in pool %d -- no size >= 8, not computable.",
                         length(npool))
    for (g in gr) for (dl in DELTA) for (r in seq_len(NRUND)) {
      set.seed(20260819L + 11000L * r + 100L * ii + 7L * g + round(1000 * dl))
      p <- sample(npool, 2 * g); gA <- p[1:g]; gB <- p[(g + 1):(2 * g)]
      iv2 <- iv; iv2[gA] <- iv2[gA] + dl
      z <- kontrast_f(iv2, Zi, gA, gB, pool, rq,
             MR(teil = "empfindlichkeit", nb = NB, groesse_seite = g,
                delta = dl, runde = r, n_neutral_pool = length(npool)))
      if (!is.null(z)) RC[[length(RC) + 1]] <- z
    }
    if (length(RC)) {
      TC <- do.call(rbind, RC); write.csv(TC, fC, row.names = FALSE)
      log("  (C) sensitivity: %d rows, sizes %s, %.1f min.", nrow(TC),
          paste(gr, collapse = "/"),
          as.numeric(difftime(Sys.time(), tc, units = "mins")))
    }
    schreibe_log()
  }

  DATEN[[lab]] <- NULL; rm(d, K0, Zi); invisible(gc())
  log("  point %d done (%.1f min in total).", ii,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
  schreibe_log()
}

log("\nRun finished after %.1f min.", as.numeric(difftime(Sys.time(), t0, units = "mins")))
schreibe_log()
