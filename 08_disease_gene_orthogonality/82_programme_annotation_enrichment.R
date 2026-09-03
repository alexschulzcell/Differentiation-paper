# =============================================================================
# 82_programme_annotation_enrichment.R -- WS6 check P3: what is the 173-gene module, substantively?
#
# Annotates the 173 module genes against EXTERNAL, curated gene sets
# (GO.db / org.Hs.eg.db / reactome.db -- no hand-built set), six predefined,
# non-overlapping substantive categories:
#   1 matrix components        GO:0005201 ECM structural constituent
#   2 matrix remodeling        GO:0022617 ECM disassembly (MMP/TIMP/LOX route)
#   3 secretion machinery      GO:0016192 vesicle-mediated transport
#                               + GO:0006888 ER-to-Golgi vesicle transport
#   4 TGFb/BMP signaling       GO:0007179 TGFbeta receptor signaling
#   5 hypoxia/stress           GO:0036293 response to decreased oxygen levels
#   6 cell-cycle exit          GO:0008285 negative reg. of cell proliferation
#                               + GO:0090398 cellular senescence
#
# Background: the 11581 genes of the internal gene map (derived_data/R_intern/
# R_interne_genkarte.csv) -- the space measurable on this layer at all.
# Test: Fisher exact, two-sided, against this background; Bonferroni over the
# EXPLICIT 6 tests (not over all GO terms -- no fishing).
# =============================================================================
suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(GO.db)
})

WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  here <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
  WURZEL <- normalizePath(file.path(here, "..", ".."))
}
AUS <- file.path(WURZEL, "derived_data", "followup")
dir.create(AUS, recursive = TRUE, showWarnings = FALSE)

GK <- read.csv(file.path(WURZEL, "derived_data", "R_intern", "R_interne_genkarte.csv"),
                stringsAsFactors = FALSE)
MOD <- read.csv(file.path(WURZEL, "manuscript", "Tabellen", "S5_konvergente_gene.csv"),
                 stringsAsFactors = FALSE)
stopifnot(nrow(MOD) == 173)

HG <- unique(GK$ensembl)
FG <- unique(MOD$gen)
stopifnot(length(FG) == 173)
cat(sprintf("background (gene map): %d genes | module: %d genes, of which in background: %d\n",
            length(HG), length(FG), length(intersect(FG, HG))))

g2e <- as.list(org.Hs.egGO2ALLEGS)
e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) {
  eg <- g2e[[id]]
  if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE)))
}

KATEGORIEN <- list(
  "1_Matrixbestandteile"   = GOr("GO:0005201"),
  "2_Matrixremodellierung" = GOr("GO:0022617"),
  "3_Sekretionsmaschine"   = union(GOr("GO:0016192"), GOr("GO:0006888")),
  "4_TGFb_BMP"             = GOr("GO:0007179"),
  "5_Hypoxie_Stress"       = GOr("GO:0036293"),
  "6_Zellzyklusausstieg"   = union(GOr("GO:0008285"), GOr("GO:0090398"))
)

# The enrichment test has lived since 2026-08-24 in 00_shared/_enrichment.R --
# one implementation for the narrow (here) and the broad sets
# (06_orthogonal_layers/61_gene_set_enrichment.R). Computation unchanged.
source(file.path(WURZEL, "00_shared", "_enrichment.R"))

zeilen <- list()
for (nm in names(KATEGORIEN)) {
  r <- ft(KATEGORIEN[[nm]], FG, HG)
  zeilen[[nm]] <- data.frame(kategorie = nm, k = r$k, n_module = r$n_fg,
                             n_satz_gesamt = r$n_satz, n_satz_im_hg = length(intersect(KATEGORIEN[[nm]], HG)),
                             n_hintergrund = r$n_hg, OR = r$OR, lo = r$lo, hi = r$hi, p = r$p)
}
T <- do.call(rbind, zeilen)
n_tests <- nrow(T)
T$p_bonferroni <- pmin(1, T$p * n_tests)
T$mde80_OR_ansatz <- NA  # detection limit separately via the power note below

cat(sprintf("\n%d explicit tests, Bonferroni factor %d\n", n_tests, n_tests))
print(T[, c("kategorie", "k", "n_module", "n_satz_im_hg", "OR", "lo", "hi", "p", "p_bonferroni")])

write.csv(T, file.path(AUS, "ws6_p3_go_annotation.csv"), row.names = FALSE)

# Which module genes fall into none of the six categories? (remainder)
alle_kat <- unique(unlist(KATEGORIEN))
rest <- setdiff(FG, alle_kat)
cat(sprintf("\nmodule genes in NONE of the six categories: %d of %d (%.0f %%)\n",
            length(rest), length(FG), 100 * length(rest) / length(FG)))

# write out symbols per category (for the report)
karte_sym <- setNames(MOD$symbol, MOD$gen)
je_kat <- lapply(KATEGORIEN, function(s) {
  g <- intersect(FG, s)
  sort(karte_sym[g])
})
sink(file.path(AUS, "ws6_p3_gene_je_kategorie.txt"))
for (nm in names(je_kat)) {
  cat(sprintf("\n== %s (n=%d) ==\n", nm, length(je_kat[[nm]])))
  cat(paste(je_kat[[nm]], collapse = ", "), "\n")
}
cat(sprintf("\n== no category (n=%d) ==\n", length(rest)))
cat(paste(sort(karte_sym[rest]), collapse = ", "), "\n")
sink()

cat(sprintf("\n-> %s\n-> %s\n", file.path(AUS, "ws6_p3_go_annotation.csv"),
            file.path(AUS, "ws6_p3_gene_je_kategorie.txt")))
