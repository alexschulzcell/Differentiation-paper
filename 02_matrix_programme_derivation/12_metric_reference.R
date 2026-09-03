# =============================================================================
# 12_metric_reference.R -- The metric per PRAEREG_Gesamtstudie.md §2 across all
#                    eleven points, plus the confounder checks of §7.
# =============================================================================
# Preregistration: preregistrations/PRAEREG_Gesamtstudie.md, in force since
# 2026-08-18, including addendum 1 (A9). Screening: SICHTUNG.md, completed.
#
# EXTENSION of the frozen reference implementation
# 01_FN1_Induktionsnull/reference_implementations/01_induktionsnull.R. The computational core -- pool,
# interaction, induction-matched null, VIF correction, MDE80 -- is
# adopted word for word. What is new is exclusively:
#
#   (E1) einzel() now also runs with ZIEH[["20"]]. In 01_ einzel() was
#        called only with ZIEH[["4"]] and with the T2 residualisation. PRAEREG
#        §2.1 defines apparatus-z and cargo-z as the z_korr of the single sets
#        against the 20-class null (§2.6). Those are the axes of the level in
#        §5 -- that is why this row here is the governing one. The old
#        4-class values keep running unchanged and are reported alongside.
#   (E2) five further single sets against the same null: distal, biosynthetic,
#        neutral (negative control), UPR and ERAD separately with sign (§7).
#   (E3) the confounder checks from §7 as a repetition of the two
#        statistics from §2.1 -- without the 500 most frequent genes, without
#        the cell-cycle genes, against an additionally length-matched null.
#   (E4) one positive control per data set: did the WT arm differentiate at
#        all? Tested is the induction term dWT of the arm-specific
#        differentiation set against a null matched to the baseline
#        expression -- the same procedure as in S2.
#   (E5) all eleven data sets instead of six. The six old points are
#        recomputed with the same implementation; in particular LAMA5
#        chondrogenic no longer sits on the 4-class null.
#
# NOT changed: gene sets, pool rule (§2.4, A6), threshold THR, NB, NVIF,
# seed, VIF machinery, MDE80 formula.
# =============================================================================
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi)
                   library(DESeq2); library(matrixStats)
                   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
                   library(GenomicFeatures) })
set.seed(20260818)

# --- path parameters ---------------------------------------------------------
# As in 10_load_reference_metric.R, the tree of raw analysis sessions is an
# explicit, overridable parameter. It holds 03_Metrik_Elf_Punkte (HIER) and
# the older sessions read through ALT below. That tree is not part of this
# repository; set PAPER_V2_SESSIONS to point at it.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
stopifnot(dir.exists(SESSIONS))
ALT  <- SESSIONS
HIER <- file.path(SESSIONS, "03_Metrik_Elf_Punkte")
OUT  <- file.path(HIER, "derived_data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

THR <- 0.5; NB <- 2000; NVIF <- 100
KLASSEN <- c(4, 20, 50)
KHAUPT  <- "20"          # PRAEREG §2.6: the 20-class null is the governing one
NTOP    <- 500           # §7 library composition

OUTLOG <- character(0)
log <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); OUTLOG <<- c(OUTLOG, s) }

log("################################################################")
log("# 12_metric_reference.R -- PRAEREG §2 ueber alle elf Punkte")
log("# Seed 20260818 | NB %d | NVIF %d | Klassenleiter %s | Hauptnull %s Klassen",
    NB, NVIF, paste(KLASSEN, collapse = "/"), KHAUPT)
log("################################################################")

# ------------------------------------------------- Sets, identical to 01_/02_
g2e <- as.list(org.Hs.egGO2ALLEGS); e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) { eg <- g2e[[id]]; if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE))) }
S_MASCHINE <- unique(unlist(lapply(c("GO:0006396","GO:0042254","GO:0006412",
  "GO:0005730","GO:0006888","GO:0048193","GO:0006486"), GOr)))
S_BIOSYN <- unique(unlist(lapply(c("GO:0006396","GO:0042254","GO:0006412",
                                   "GO:0005730"), GOr)))
S_DISTAL <- unique(unlist(lapply(c("GO:0006888","GO:0048193","GO:0006486",
                                   "GO:0030968"), GOr)))
S_FRACHT <- GOr("GO:0005615")
S_NEUTRAL <- GOr("GO:0007268")                       # synaptic transmission
ueb <- intersect(S_MASCHINE, S_FRACHT)
S_MASCHINE <- setdiff(S_MASCHINE, ueb); S_FRACHT0 <- setdiff(S_FRACHT, ueb)
S_DISTAL <- setdiff(S_DISTAL, intersect(S_DISTAL, S_FRACHT))
S_BIOSYN <- setdiff(S_BIOSYN, intersect(S_BIOSYN, S_FRACHT))
S_BIOSYN <- setdiff(S_BIOSYN, S_DISTAL)
S_NEUTRAL <- setdiff(S_NEUTRAL, c(S_MASCHINE, S_FRACHT, S_DISTAL, S_BIOSYN))
# §7: stress response separated and with direction
S_UPR  <- GOr("GO:0030968")
S_ERAD <- GOr("GO:0036503")
# §7: Proliferation
S_ZYKLUS <- GOr("GO:0007049")
# (E4) arm-specific differentiation sets, only as the data set's positive control
S_OSSIF  <- GOr("GO:0001503")                        # ossification
S_KNORPEL <- GOr("GO:0051216")                       # cartilage development
log("\nSets: Maschine %d | biosyn %d | distal %d | Fracht %d | Neutral %d",
    length(S_MASCHINE), length(S_BIOSYN), length(S_DISTAL),
    length(S_FRACHT0), length(S_NEUTRAL))
log("      UPR %d | ERAD %d | Zellzyklus %d | Ossifikation %d | Knorpel %d",
    length(S_UPR), length(S_ERAD), length(S_ZYKLUS), length(S_OSSIF),
    length(S_KNORPEL))

# ------------------------------------------- Protein length (§7 gene size)
# Word for word from 9_Leithypothese_Kapazitaet/reference_implementations/41_frachtlast.R.
cdsl <- sum(width(cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "tx")))
tx2g <- suppressMessages(
  AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene,
                        keys = names(cdsl), keytype = "TXID",
                        columns = c("TXID", "GENEID")))
tx2g <- tx2g[!is.na(tx2g$GENEID), ]
tx2g$len <- cdsl[as.character(tx2g$TXID)]
per_eg <- tapply(tx2g$len, tx2g$GENEID, max)
ens <- e2ens[names(per_eg)]
AS <- setNames(rep(per_eg / 3, lengths(ens)), unlist(ens, use.names = FALSE))
AS <- AS[!is.na(names(AS))]; AS <- tapply(AS, names(AS), max)
log("Proteinlaengen (AS, laengste codierende Isoform) fuer %d Gene", length(AS))

# =============================================================================
# The computational core. Taken from 01_induktionsnull.R; Z, meta as there.
# Returns the quantities all tests need -- so that the confounder variants
# from §7 can call the same core with a modified Z.
# =============================================================================
kern <- function(Z, meta) {
  meta <- meta[match(colnames(Z), meta$sample), ]
  mZ <- function(s) rowMeans(Z[, s, drop = FALSE])
  s_of <- function(g, c) meta$sample[meta$genotype == g & meta$condition == c]
  iv  <- (mZ(s_of("KO","diff")) - mZ(s_of("KO","undiff"))) -
         (mZ(s_of("WT","diff")) - mZ(s_of("WT","undiff")))
  dWT <-  mZ(s_of("WT","diff")) - mZ(s_of("WT","undiff"))

  mm <- model.matrix(~ genotype * condition, data = meta)
  H <- mm %*% solve(crossprod(mm)) %*% t(mm)
  RES <- Z %*% (diag(ncol(Z)) - H); sdr <- sqrt(rowSums(RES^2)); ok <- sdr > 1e-8
  RESn <- RES[ok, ] / sdr[ok]
  rho_quer <- function(g) { g <- intersect(g, rownames(RESn)); m <- length(g)
    if (m < 3) return(NA_real_)
    s <- colSums(RESn[g, , drop = FALSE]); (sum(s^2) - m) / (m * (m - 1)) }

  pool <- names(dWT)[dWT >= THR & !is.na(dWT)]
  list(iv = iv, dWT = dWT, pool = pool, rho_quer = rho_quer,
       basis = mZ(s_of("WT","undiff")))
}

# an induction-matched draw with K classes
mk_zieh <- function(dWT, pool, K) {
  brk <- unique(quantile(dWT[pool], seq(0, 1, length.out = K + 1)))
  bin <- cut(dWT[pool], breaks = brk, include.lowest = TRUE, labels = FALSE)
  names(bin) <- pool; pb <- split(pool, bin)
  function(g) { cnt <- table(bin[g]); cnt <- cnt[cnt > 0]
    unlist(lapply(names(cnt), function(b)
      sample(pb[[b]], cnt[[b]], replace = FALSE))) } }

# additionally length-matched: cell = induction bin x length quartile
# (procedure from 41_frachtlast.R, there 4 x 4; here the same)
mk_zieh_L <- function(dWT, pool) {
  brk <- unique(quantile(dWT[pool], seq(0, 1, 0.25)))
  bin <- cut(dWT[pool], breaks = brk, include.lowest = TRUE, labels = FALSE)
  names(bin) <- pool
  lp <- AS[intersect(pool, names(AS))]
  lbrk <- unique(quantile(lp, seq(0, 1, 0.25)))
  lbin <- cut(lp, breaks = lbrk, include.lowest = TRUE, labels = FALSE)
  names(lbin) <- names(lp)
  key <- paste(bin[names(lbin)], lbin); plb <- split(names(lbin), key)
  list(gene = names(lbin),
       zieh = function(g) { g <- intersect(g, names(lbin))
         k <- paste(bin[g], lbin[g]); cnt <- table(k)
         unlist(lapply(names(cnt), function(kk) {
           p <- plb[[kk]]; n <- min(cnt[[kk]], length(p))
           if (n == 0) character(0) else sample(p, n, replace = FALSE) })) }) }

# ------------------------------------------------------ one single-set test
# Word for word from 01_, extended by MDE80 (in 01_ it was NA in einzel(); PRAEREG
# §2.7 requires it for EVERY number, and the formula already stands there in kontrast()).
einzel_f <- function(v, zieh, g, pool, rho_quer, meta_row) {
  g <- intersect(g, pool); if (length(g) < 8) return(NULL)
  obs <- mean(v[g]); nl <- replicate(NB, mean(v[zieh(g)]))
  z <- (obs - mean(nl)) / sd(nl)
  r_s <- rho_quer(g)
  r_n <- mean(replicate(NVIF, rho_quer(zieh(g))), na.rm = TRUE)
  f <- sqrt((1 + (length(g)-1)*r_s) / (1 + (length(g)-1)*r_n))
  mde <- mean(nl) + (1.96 + 0.84) * sd(nl) * f
  cbind(meta_row, data.frame(n1 = length(g), n2 = NA_integer_, diff = obs,
        z_roh = z, faktor = f, z_korr = z / f,
        p_korr = 2 * pnorm(-abs(z / f)), mde80 = mde,
        stringsAsFactors = FALSE)) }

kontrast_f <- function(v, zieh, g1, g2, pool, rho_quer, meta_row) {
  g1 <- intersect(g1, pool); g2 <- intersect(g2, pool)
  if (length(g1) < 8 || length(g2) < 8) return(NULL)
  obs <- mean(v[g1]) - mean(v[g2])
  nl  <- replicate(NB, mean(v[zieh(g1)]) - mean(v[zieh(g2)]))
  z   <- (obs - mean(nl)) / sd(nl)
  ff <- function(g) { r_s <- rho_quer(g)
    r_n <- mean(replicate(NVIF, rho_quer(zieh(g))), na.rm = TRUE)
    sqrt((1 + (length(g)-1)*r_s) / (1 + (length(g)-1)*r_n)) }
  f <- sqrt((ff(g1)^2 + ff(g2)^2) / 2)
  mde <- mean(nl) + (1.96 + 0.84) * sd(nl) * f
  cbind(meta_row, data.frame(n1 = length(g1), n2 = length(g2), diff = obs,
        z_roh = z, faktor = f, z_korr = z / f,
        p_korr = 2 * pnorm(-abs(z / f)), mde80 = mde,
        stringsAsFactors = FALSE)) }

# =============================================================================
# One data set, complete.
#   Z     z-standardised matrix of the contrast (project standard)
#   expr  mean expression per gene, only for ranking the 500 most frequent (§7)
# =============================================================================
auswerten <- function(Z, meta, expr, label, arm, rolle, klasse) {
  K0 <- kern(Z, meta)
  if (length(K0$pool) < 1000) {
    log("\n--- %s --- POOL %d < 1000 => AUSSCHLUSS A6 (PRAEREG §2.4/§3.2)",
        label, length(K0$pool)); return(NULL) }

  iv <- K0$iv; dWT <- K0$dWT; pool <- K0$pool; rq <- K0$rho_quer
  abstumpfung <- mean(iv[pool]); sd_pool <- sd(iv[pool])

  log("\n============================================================")
  log("--- %s (%s, Klasse %s, %s) ---", label, arm, klasse, rolle)
  log("  Pool %d Gene | Abstumpfung %+.3f | sd(iv) im Pool %.3f",
      length(pool), abstumpfung, sd_pool)

  # T2 residualisation (§2.2, for EVERY data set)
  lo <- lowess(dWT[pool], iv[pool], f = 0.3)
  fit <- approx(lo$x, lo$y, xout = dWT[pool], rule = 2)$y
  iv_res <- setNames(iv[pool] - fit, pool)
  log("  T2 Loess iv~dWT: erklaerte Varianz %.1f %% | Rest-Abstumpfung %+.4f",
      100 * (1 - var(iv_res) / var(iv[pool])), mean(iv_res))

  ZIEH <- lapply(KLASSEN, function(K) mk_zieh(dWT, pool, K))
  names(ZIEH) <- as.character(KLASSEN)

  MR <- function(test, null) data.frame(datensatz = label, arm = arm,
    klasse = klasse, rolle = rolle, test = test, null = null,
    stringsAsFactors = FALSE)
  R <- list()
  add <- function(x) if (!is.null(x)) R[[length(R)+1]] <<- x

  # --- unchanged from 01_: contrasts over the class ladder -----------------
  for (K in as.character(KLASSEN)) {
    nn <- paste0(K, "_Klassen")
    add(kontrast_f(iv, ZIEH[[K]], S_FRACHT0, S_MASCHINE, pool, rq,
                   MR("Fracht-Maschine", nn)))
    add(kontrast_f(iv, ZIEH[[K]], S_FRACHT0, S_DISTAL, pool, rq,
                   MR("Fracht-distal", nn))) }
  add(kontrast_f(iv_res, ZIEH[["4"]], S_FRACHT0, S_MASCHINE, pool, rq,
                 MR("Fracht-Maschine", "T2_residualisiert")))
  add(kontrast_f(iv_res, ZIEH[["4"]], S_FRACHT0, S_DISTAL, pool, rq,
                 MR("Fracht-distal", "T2_residualisiert")))

  # --- unchanged from 01_: single sets against the 4-class null -----------
  add(einzel_f(iv, ZIEH[["4"]], S_FRACHT0,  pool, rq, MR("nur_Fracht",   "4_Klassen")))
  add(einzel_f(iv, ZIEH[["4"]], S_MASCHINE, pool, rq, MR("nur_Maschine", "4_Klassen")))
  add(einzel_f(iv_res, ZIEH[["4"]], S_FRACHT0,  pool, rq, MR("nur_Fracht",   "T2_residualisiert")))
  add(einzel_f(iv_res, ZIEH[["4"]], S_MASCHINE, pool, rq, MR("nur_Maschine", "T2_residualisiert")))

  # --- (E1) THE GAP: apparatus-z and cargo-z per §2.1 against the ----------
  # ---      20-class null. Those are the axes of the level from §5.  --------
  nnH <- paste0(KHAUPT, "_Klassen")
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_MASCHINE, pool, rq, MR("Apparat_z",  nnH)))
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_FRACHT0,  pool, rq, MR("Fracht_z",   nnH)))
  # --- (E2) accompanying single sets against the same null (§2.2, §7) -------
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_DISTAL,  pool, rq, MR("nur_distal",  nnH)))
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_BIOSYN,  pool, rq, MR("nur_biosyn",  nnH)))
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_NEUTRAL, pool, rq, MR("NEG_neutral", nnH)))
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_UPR,     pool, rq, MR("UPR",         nnH)))
  add(einzel_f(iv, ZIEH[[KHAUPT]], S_ERAD,    pool, rq, MR("ERAD",        nnH)))
  # decoupling as a contrast against the same null (§2.1, supplement; §5.4)
  add(kontrast_f(iv, ZIEH[[KHAUPT]], S_FRACHT0, S_MASCHINE, pool, rq,
                 MR("Entkopplung", nnH)))

  # --- (E3) confounders per §7 ---------------------------------------------
  # (a) length-matched null -- same target quantity, different null
  ZL <- mk_zieh_L(dWT, pool)
  add(einzel_f(iv, ZL$zieh, intersect(S_MASCHINE, ZL$gene), pool, rq,
               MR("Apparat_z", "S7_laengenangeglichen")))
  add(einzel_f(iv, ZL$zieh, intersect(S_FRACHT0, ZL$gene), pool, rq,
               MR("Fracht_z",  "S7_laengenangeglichen")))
  # (b) without the 500 most frequent genes: they are removed from the matrix,
  #     every sample is recentred on its median over the remaining genes
  #     (that is the composition effect that is to be tested), and the core
  #     runs completely anew.
  if (!is.null(expr)) {
    top <- names(sort(expr[intersect(names(expr), rownames(Z))],
                      decreasing = TRUE))[seq_len(min(NTOP, length(expr)))]
    Z2 <- Z[setdiff(rownames(Z), top), , drop = FALSE]
    Z2 <- sweep(Z2, 2, colMedians(Z2), "-")
    K2 <- kern(Z2, meta)
    if (length(K2$pool) >= 1000) {
      Zi2 <- mk_zieh(K2$dWT, K2$pool, as.integer(KHAUPT))
      add(einzel_f(K2$iv, Zi2, S_MASCHINE, K2$pool, K2$rho_quer,
                   MR("Apparat_z", "S7_ohne_Top500")))
      add(einzel_f(K2$iv, Zi2, S_FRACHT0, K2$pool, K2$rho_quer,
                   MR("Fracht_z",  "S7_ohne_Top500")))
    } else log("  S7 ohne Top500: Pool %d < 1000, nicht auswertbar", length(K2$pool))
  }
  # (c) without the cell-cycle genes
  Z3 <- Z[setdiff(rownames(Z), S_ZYKLUS), , drop = FALSE]
  K3 <- kern(Z3, meta)
  if (length(K3$pool) >= 1000) {
    Zi3 <- mk_zieh(K3$dWT, K3$pool, as.integer(KHAUPT))
    add(einzel_f(K3$iv, Zi3, S_MASCHINE, K3$pool, K3$rho_quer,
                 MR("Apparat_z", "S7_ohne_Zellzyklus")))
    add(einzel_f(K3$iv, Zi3, S_FRACHT0, K3$pool, K3$rho_quer,
                 MR("Fracht_z",  "S7_ohne_Zellzyklus")))
  } else log("  S7 ohne Zellzyklus: Pool %d < 1000, nicht auswertbar", length(K3$pool))

  # --- (E4) positive control: did the WT arm differentiate at all? ----------
  # Procedure as in S2: induction term dWT against a null matched to the
  # baseline expression (10 classes). Arm-specific set, never across arms (§4).
  SDIFF <- if (arm == "osteogen") S_OSSIF else S_KNORPEL
  alle <- names(dWT)[!is.na(dWT) & !is.na(K0$basis)]
  zieh_b <- mk_zieh(K0$basis, alle, 10)
  add(einzel_f(dWT, zieh_b, SDIFF, alle, rq,
               MR(if (arm == "osteogen") "POS_ossifikation" else "POS_knorpel",
                  "Ausgangsexpression_10")))
  add(einzel_f(dWT, zieh_b, S_NEUTRAL, alle, rq,
               MR("NEG_neutral_dWT", "Ausgangsexpression_10")))

  tab <- do.call(rbind, R)
  tab$diff_norm <- tab$diff / sd_pool
  tab$abstumpfung <- abstumpfung; tab$sd_pool <- sd_pool; tab$pool <- length(pool)

  for (i in seq_len(nrow(tab)))
    log("  %-17s %-24s n=%9s | diff %+.3f | z_korr %+7.2f | p %9.3g | MDE80 %+.3f",
        tab$test[i], tab$null[i],
        ifelse(is.na(tab$n2[i]), as.character(tab$n1[i]),
               paste0(tab$n1[i], "/", tab$n2[i])),
        tab$diff[i], tab$z_korr[i], tab$p_korr[i], tab$mde80[i])
  tab
}

# =============================================================================
# Die elf Datensaetze
# =============================================================================
zmat <- function(X) { Z <- t(scale(t(X))); Z[!is.na(rowSums(Z)), , drop = FALSE] }
lade <- function(pfad, nm) {
  cs <- read.csv(file.path(pfad, paste0(nm, "_counts.csv")), row.names = 1,
                 check.names = FALSE)
  mt <- read.csv(file.path(pfad, paste0(nm, "_meta.csv")), stringsAsFactors = FALSE)
  list(X = as.matrix(cs), meta = mt) }
# vorb() unchanged from 01_, additionally returns the pre-z matrix so that
# §7 (b) can rank the 500 most frequent genes.
vorb2 <- function(X, roh) {
  rownames(X) <- sub("\\..*$", "", rownames(X))
  X <- X[!duplicated(rownames(X)), , drop = FALSE]
  if (roh) {
    X <- X[rowSums(X >= 5) >= max(3, floor(ncol(X) / 4)), , drop = FALSE]
    d <- DESeqDataSetFromMatrix(round(X), data.frame(x = rep(1, ncol(X))), ~ 1)
    X <- assay(rlog(d, blind = TRUE))
  } else X <- X[rowVars(X) > 0, , drop = FALSE]
  list(Z = zmat(X), expr = rowMeans(X)) }

erg <- list()

# ---- 1/2: LAMA5-KO, own data, osteogenic and chondrogenic -----------------
B <- readRDS(file.path(ALT, "8_Leithypothese_Kontaktschritt", "derived_data",
                       "bulk_results.rds"))
dat <- B$dat
me <- data.frame(sample = dat$sample, genotype = as.character(dat$genotype),
                 condition = ifelse(as.character(dat$condition) == "naiv", "undiff",
                                    as.character(dat$condition)),
                 stringsAsFactors = FALSE)
# Rank of the most frequent genes for §7 (b): baseMean of the same bulk data set
bg <- read.csv(file.path(ALT, "12_Konvergenz_Kontaktschritt", "derived_data",
                         "A0_hintergrund.csv"), stringsAsFactors = FALSE)
expr_lama5 <- setNames(bg$baseMean, bg$ensembl)
for (arm in c("chondro", "osteo")) {
  sel <- me[me$condition %in% c("undiff", arm), ]
  sel$condition <- ifelse(sel$condition == "undiff", "undiff", "diff")
  erg[[paste0("LAMA5_", arm)]] <-
    auswerten(B$Z[, sel$sample, drop = FALSE], sel, expr_lama5,
              paste0("LAMA5-KO (", substr(arm, 1, 6), ")"),
              if (arm == "osteo") "osteogen" else "chondrogen",
              "Positivkontrolle", "M")
}

# ---- 3/4: FN1 C123R and C231W (GSE251698) --------------------------------
D20 <- file.path(ALT, "20_Generalisierung_Matrixdefekte", "data_raw")
for (nm in c("FN1_FNC123R", "FN1_FNC231W")) {
  L <- lade(D20, nm); V <- vorb2(L$X, FALSE)
  erg[[nm]] <- auswerten(V$Z, L$meta, V$expr, sub("FN1_", "FN1 ", nm),
                         "chondrogen", "Pruefling", "M") }

# ---- 5/6: the two negative controls of the predecessor project ------------
D18 <- file.path(ALT, "18_Spezifitaet_Schere", "data_raw")
L <- lade(D18, "D2_GSE247491"); V <- vorb2(L$X, TRUE)
erg[["SERPINA3"]] <- auswerten(V$Z, L$meta, V$expr, "SERPINA3-KD",
                               "chondrogen", "Negativkontrolle", "N")
L <- lade(D18, "D1_GSE184087"); V <- vorb2(L$X, FALSE)
erg[["MIR181A1HG"]] <- auswerten(V$Z, L$meta, V$expr, "MIR181A1HG-KD",
                                 "osteogen", "Negativkontrolle", "N")

# ---- 7: GSE227512, LINC01638-KD, osteogenic, raw counts -------------------
L <- lade(D18, "D3_GSE227512"); V <- vorb2(L$X, TRUE)
erg[["GSE227512"]] <- auswerten(V$Z, L$meta, V$expr, "LINC01638-KD (GSE227512)",
                                "osteogen", "Pruefling", "N")

# ---- 8-11: the four newly loaded GEO matrices -----------------------------
source(file.path(WURZEL, "02_matrix_programme_derivation", "13_geo_matrices_to_metric_format.R"))
NEUMETA <- list(
  GSE218101 = list(lab = "ARSB / MPS VI (GSE218101)", arm = "chondrogen", kl = "M"),
  GSE221128 = list(lab = "ACVR1 / FOP (GSE221128)",   arm = "chondrogen", kl = "M"),
  GSE205432 = list(lab = "RNF4-KD (GSE205432)",       arm = "osteogen",   kl = "N"),
  GSE245585 = list(lab = "RB1 +/- (GSE245585)",       arm = "osteogen",   kl = "N"))
for (n in names(NEUMETA)) {
  L <- NEU[[n]]()
  M <- L$X[rowVars(L$X) > 0, , drop = FALSE]
  erg[[n]] <- auswerten(zmat(M), L$meta, rowMeans(M), NEUMETA[[n]]$lab,
                        NEUMETA[[n]]$arm, "Pruefling", NEUMETA[[n]]$kl) }

TAB <- do.call(rbind, erg)
write.csv(TAB, file.path(OUT, "03_metrik_elf.csv"), row.names = FALSE)

# =============================================================================
# Overview: the two statistics per §2.1 for each data set
# =============================================================================
hol <- function(d, t, nl) TAB[TAB$datensatz == d & TAB$test == t & TAB$null == nl, ]
DS <- unique(TAB$datensatz)
nnH <- paste0(KHAUPT, "_Klassen")
log("\n\n################ §2.1 Apparat-z und Fracht-z (%s) ################", nnH)
log("%-28s %-11s %-2s %8s %8s %10s %10s", "Datensatz", "Arm", "Kl",
    "Apparat", "Fracht", "Entkopplg", "MDE80(App)")
for (d in DS) {
  a <- hol(d, "Apparat_z", nnH); f <- hol(d, "Fracht_z", nnH)
  e <- hol(d, "Entkopplung", nnH)
  if (!nrow(a) || !nrow(f)) next
  log("%-28s %-11s %-2s %+8.2f %+8.2f %+10.2f %+10.3f", d, a$arm, a$klasse,
      a$z_korr, f$z_korr, if (nrow(e)) e$z_korr else NA_real_, a$mde80) }

log("\n################ Entkopplung Fracht-Maschine ueber alle Nullen ################")
log("%-28s %-24s %8s %8s %9s %10s", "Datensatz", "Null", "diff", "z_korr", "p", "MDE80")
for (d in DS) for (nl in c(paste0(KLASSEN, "_Klassen"), "T2_residualisiert")) {
  r <- hol(d, "Fracht-Maschine", nl); if (!nrow(r)) next
  log("%-28s %-24s %+8.3f %+8.2f %9.2g %+10.3f", d, nl, r$diff, r$z_korr,
      r$p_korr, r$mde80) }

log("\n################ §7 Stoergroessen: kippt die Zuweisung? ################")
feld <- function(a, f) {
  if (is.na(a) || is.na(f)) return("?")
  A <- if (a > 2) "hoch" else if (a < -2) "tief" else "flach"
  Fr <- if (f > 2) "hoch" else if (f < -2) "tief" else "flach"
  if (A == "tief" && Fr == "tief") return("Globalausfall")
  if (A == "tief") return("Kapazitaetsversagen")
  if (A == "hoch" && Fr != "hoch") return("Programmversagen")
  if (A == "flach" && Fr == "tief") return("Programmversagen")
  if (Fr == "hoch") return("intakt/ausbauend")
  "keine Laesion" }
NULLEN <- c(nnH, "S7_laengenangeglichen", "S7_ohne_Top500", "S7_ohne_Zellzyklus",
            "T2_residualisiert")
log("%-28s %-24s %8s %8s  %s", "Datensatz", "Pruefung", "Apparat", "Fracht", "Feld")
STAB <- list()
for (d in DS) {
  ref <- NA_character_
  for (nl in NULLEN) {
    # T2 exists only as a 4-class single set (nur_Maschine / nur_Fracht)
    tA <- if (nl == "T2_residualisiert") "nur_Maschine" else "Apparat_z"
    tF <- if (nl == "T2_residualisiert") "nur_Fracht"   else "Fracht_z"
    a <- hol(d, tA, nl); f <- hol(d, tF, nl); if (!nrow(a) || !nrow(f)) next
    q <- feld(a$z_korr, f$z_korr); if (is.na(ref)) ref <- q
    log("%-28s %-24s %+8.2f %+8.2f  %s%s", d, nl, a$z_korr, f$z_korr, q,
        if (nl != nnH && q != ref) "   <-- KIPPT" else "")
    STAB[[length(STAB)+1]] <- data.frame(datensatz = d, pruefung = nl,
      apparat_z = a$z_korr, fracht_z = f$z_korr, feld = q,
      kippt = (nl != nnH && q != ref), stringsAsFactors = FALSE) } }
STAB <- do.call(rbind, STAB)
write.csv(STAB, file.path(OUT, "03_stabilitaet.csv"), row.names = FALSE)

log("\n################ Kontrollen je Datensatz (§7 letzte Zeile) ################")
log("%-28s %10s %10s %10s", "Datensatz", "POS(dWT)", "NEG(dWT)", "NEG(iv)")
for (d in DS) {
  p <- rbind(hol(d, "POS_ossifikation", "Ausgangsexpression_10"),
             hol(d, "POS_knorpel", "Ausgangsexpression_10"))
  n1 <- hol(d, "NEG_neutral_dWT", "Ausgangsexpression_10")
  n2 <- hol(d, "NEG_neutral", nnH)
  log("%-28s %+10.2f %+10.2f %+10.2f", d,
      if (nrow(p)) p$z_korr else NA_real_,
      if (nrow(n1)) n1$z_korr else NA_real_,
      if (nrow(n2)) n2$z_korr else NA_real_) }

log("\n################ UPR / ERAD, mit Richtung (§7) ################")
log("%-28s %10s %10s", "Datensatz", "UPR", "ERAD")
for (d in DS) {
  u <- hol(d, "UPR", nnH); e <- hol(d, "ERAD", nnH)
  log("%-28s %+10.2f %+10.2f", d, if (nrow(u)) u$z_korr else NA_real_,
      if (nrow(e)) e$z_korr else NA_real_) }

writeLines(OUTLOG, file.path(OUT, "03_log.txt"))
cat("\ngeschrieben:", file.path(OUT, "03_metrik_elf.csv"), "\n")
