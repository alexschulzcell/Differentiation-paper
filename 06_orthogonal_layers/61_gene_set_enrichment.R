# =============================================================================
# 61_gene_set_enrichment.R -- put Figure 2D on gene sets that can carry it.
#
# The gene sets are fixed and stored BEFORE this run
# (`06_orthogonal_layers/60_gene_sets_build.R`), each with its source, version and
# retrieval date. This script READS them and recomputes nothing. The
# enrichment test itself is the one implementation of the project
# (`00_shared/_enrichment.R`), the same one that produced the
# narrow run.
#
# -----------------------------------------------------------------------------
# THE DECISION RULE -- written here before the first number was computed, and
# not changed afterwards.
# -----------------------------------------------------------------------------
# For each of the six categories the enrichment is recomputed with the broad
# set, using the same statistic (Fisher exact, two-sided), the same null (the
# 11 581 genes of the internal gene map) and the same detection limit as
# before.
#
#   (a) CONFIRMED: same direction, and the 95 % CI excludes 1 (for the
#       secretory machinery, excludes 1 from above). The broad set then
#       replaces the narrow one in Figure 2D, and the legend names both.
#   (b) WEAKENED: same direction, CI includes 1. The direction stays in the
#       text and THE NUMBER LEAVES THE RUNNING TEXT -- it appears only in the
#       supplement, with both sets side by side.
#   (c) REVERSED OR NULL: the category leaves Figure 2D and the discussion is
#       rewritten. That is the case which would hit the orthogonality
#       argument, and it would be reported as such.
#
# "Same direction" means OR > 1 in the narrow run and OR > 1 in the broad one
# (or < 1 and < 1). Categories that were already null in the narrow run can
# only give (a) or, in the sense of "was already null", (b)/(c); they are
# recorded as "was already null".
#
# Bonferroni over the number of PRIMARY categories = 6, alpha = 0.05 / 6 =
# 0.00833. The Reactome variants of the two matrix axes run as a sensitivity
# analysis; they do not enter the Bonferroni factor and replace nothing.
#
# The narrow set is reported in EVERY case -- nothing is replaced without both
# staying visible.
#
# n_module is carried per run. In the narrow GO run it is 147 of the 173
# module genes, namely those that lie in the GO background. The background is
# the same for all sets here -- the 11 581 genes of the internal gene map --
# so n_module is identical for narrow and broad sets. That is checked and
# written out, because the narrow and the broad number must never be mixed.
# =============================================================================

.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
source(file.path(WURZEL, "00_shared", "_enrichment.R"))

ABLAGE <- file.path(WURZEL, "data_raw", "_referenz", "gensaetze_v2")
ERG    <- file.path(WURZEL, "results")
dir.create(ERG, recursive = TRUE, showWarnings = FALSE)

GK  <- read.csv(file.path(WURZEL, "derived_data", "R_intern",
                          "R_interne_genkarte.csv"), stringsAsFactors = FALSE)
MOD <- read.csv(file.path(WURZEL, "derived_data", "reference_tables",
                          "S5_konvergente_gene.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(MOD) == 173)

HG <- unique(GK$ensembl)
FG <- unique(MOD$gen)
stopifnot(length(FG) == 173, length(HG) == 11581)
cat(sprintf("background %d genes | module %d, of which in the background %d\n",
            length(HG), length(FG), length(intersect(FG, HG))))

if (!file.exists(file.path(ABLAGE, "mitglieder.csv")))
  stop("The frozen MSigDB gene sets are not in this repository. Run ",
       "06_orthogonal_layers/60_gene_sets_build.R once (needs MSigDB and ",
       "org.Hs.eg.db, see 00_setup.md). The outputs of this script are ",
       "committed under results/, and the figure steps read those.",
       call. = FALSE)

MIT <- read.csv(file.path(ABLAGE, "mitglieder.csv"), stringsAsFactors = FALSE)
UEB <- read.csv(file.path(ABLAGE, "uebersicht.csv"), stringsAsFactors = FALSE)

satz <- function(kat, var) MIT$ensembl[MIT$kategorie == kat & MIT$variante == var]

KATEGORIEN <- sort(unique(UEB$kategorie))
zeilen <- list()
for (kat in KATEGORIEN) {
  for (var in c("eng", "breit", "empfindlichkeit")) {
    s <- satz(kat, var)
    if (!length(s)) next
    r <- ft(s, FG, HG)
    zeilen[[length(zeilen) + 1]] <- data.frame(
      kategorie = kat, variante = var,
      quelle = UEB$quelle[UEB$kategorie == kat & UEB$variante == var],
      bestandteile = UEB$bestandteile[UEB$kategorie == kat & UEB$variante == var],
      k = r$k, n_module = r$n_fg, n_satz_gesamt = length(s),
      n_satz_im_hg = r$n_satz, n_hintergrund = r$n_hg,
      OR = r$OR, lo = r$lo, hi = r$hi, p = r$p, stringsAsFactors = FALSE)
  }
}
T <- do.call(rbind, zeilen)

# Bonferroni over the six PRIMARY categories (broad) only.
N_PRIMAER <- length(KATEGORIEN)
ALPHA <- 0.05 / N_PRIMAER
T$p_bonferroni <- ifelse(T$variante == "empfindlichkeit", NA,
                         pmin(1, T$p * N_PRIMAER))
T$ki_schliesst_1_aus <- (T$lo > 1) | (T$hi < 1)

# n_module must not move between the variants.
stopifnot(length(unique(T$n_module)) == 1)
N_MODUL <- unique(T$n_module)
cat(sprintf("n_module (in the background of the internal gene map): %d of 173\n",
            N_MODUL))

# ---- apply the decision rule ------------------------------------------------
eng   <- T[T$variante == "eng", ]
breit <- T[T$variante == "breit", ]
rownames(eng) <- eng$kategorie; rownames(breit) <- breit$kategorie

urteil <- data.frame(kategorie = KATEGORIEN, stringsAsFactors = FALSE)
urteil$OR_eng        <- eng[KATEGORIEN, "OR"]
urteil$lo_eng        <- eng[KATEGORIEN, "lo"]
urteil$hi_eng        <- eng[KATEGORIEN, "hi"]
urteil$k_eng         <- eng[KATEGORIEN, "k"]
urteil$n_satz_eng    <- eng[KATEGORIEN, "n_satz_im_hg"]
urteil$OR_breit      <- breit[KATEGORIEN, "OR"]
urteil$lo_breit      <- breit[KATEGORIEN, "lo"]
urteil$hi_breit      <- breit[KATEGORIEN, "hi"]
urteil$k_breit       <- breit[KATEGORIEN, "k"]
urteil$n_satz_breit  <- breit[KATEGORIEN, "n_satz_im_hg"]
urteil$p_breit       <- breit[KATEGORIEN, "p"]
urteil$p_bonf_breit  <- breit[KATEGORIEN, "p_bonferroni"]

richtung <- function(or) ifelse(or > 1, "auf", ifelse(or < 1, "ab", "null"))
urteil$richtung_eng   <- richtung(urteil$OR_eng)
urteil$richtung_breit <- richtung(urteil$OR_breit)
urteil$ki_eng_aus     <- (urteil$lo_eng > 1) | (urteil$hi_eng < 1)
urteil$ki_breit_aus   <- (urteil$lo_breit > 1) | (urteil$hi_breit < 1)

urteil$fall <- with(urteil, ifelse(
  !ki_eng_aus, "war schon null",
  ifelse(richtung_breit != richtung_eng, "(c) umgekehrt oder null",
  ifelse(ki_breit_aus, "(a) bestaetigt", "(b) abgeschwaecht"))))

cat("\n== narrow vs broad set ==\n")
print(urteil[, c("kategorie", "k_eng", "n_satz_eng", "OR_eng",
                 "k_breit", "n_satz_breit", "OR_breit", "lo_breit",
                 "hi_breit", "p_breit", "p_bonf_breit", "fall")],
      digits = 3)

cat(sprintf("\nBonferroni: %d primary tests, alpha = %.5f\n", N_PRIMAER, ALPHA))
cat("\n== Sensitivity (Reactome, second independent source) ==\n")
e <- T[T$variante == "empfindlichkeit", ]
if (nrow(e)) print(e[, c("kategorie", "k", "n_satz_im_hg", "OR", "lo", "hi", "p")],
                   digits = 3)

write.csv(T, file.path(ERG, "gensaetze_v2_anreicherung.csv"), row.names = FALSE)
write.csv(urteil, file.path(ERG, "gensaetze_v2_urteil.csv"), row.names = FALSE)

# The panel file F2D is shaped -- like every panel file in this repository --
# by `09_figures/10_panel_data_main.py` from `results/gensaetze_v2_anreicherung.csv`, not
# here. One panel with two writers would be exactly the second implementation
# that the project rules forbid.

# ---- which module genes carry which category? ------------------------------
sym <- setNames(MOD$symbol, MOD$gen)
gz <- list()
for (kat in KATEGORIEN) for (var in c("eng", "breit", "empfindlichkeit")) {
  s <- satz(kat, var); if (!length(s)) next
  g <- sort(intersect(intersect(FG, HG), s))
  if (!length(g)) next
  gz[[length(gz) + 1]] <- data.frame(kategorie = kat, variante = var,
    ensembl = g, symbol = unname(sym[g]), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, gz), file.path(ERG, "gensaetze_v2_modulgene.csv"),
          row.names = FALSE)

cat(sprintf("\n-> %s\n-> %s\n-> %s\n",
            file.path(ERG, "gensaetze_v2_anreicherung.csv"),
            file.path(ERG, "gensaetze_v2_urteil.csv"),
            file.path(ERG, "gensaetze_v2_modulgene.csv")))
writeLines(capture.output(sessionInfo()),
           file.path(ERG, "gensaetze_v2_sessionInfo.txt"))
