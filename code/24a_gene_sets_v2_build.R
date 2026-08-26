# =============================================================================
# 24a_gene_sets_v2_build.R -- FIX the broad gene sets and store them.
#
# Separate from the analysis (`code/24_gene_sets_v2.R`) and run before it: the
# sets are written once, with source, version and retrieval date, and are not
# touched afterwards. A set that is changed after a result has been seen is no
# longer an external set.
#
# The script does NOT overwrite an existing store. Rebuilding it requires
# GENSAETZE_V2_NEU=1, and that then appears in the log.
#
# Why these sets (the reasoning before the run, not after)
# --------------------------------------------------------
# The narrow run uses single, very narrow GO terms. GO:0022617 (ECM
# disassembly) has 33 genes in the background, GO:0005201 (ECM structural
# constituent) has 104. With 147 module genes, k cannot become large there,
# and the confidence intervals are correspondingly wide (1.53-25.9). More
# propagation does not help -- `org.Hs.egGO2ALLEGS` already propagates to all
# child terms. What helps are BROADER and INDEPENDENTLY curated sets.
#
#   matrix components     -> the matrisome, core division (collagens, ECM
#                            glycoproteins, proteoglycans). Naba et al.,
#                            MatrisomeDB, through MSigDB. Curated
#                            independently of GO and built for this question.
#   matrix remodelling    -> the matrisome, ECM regulators (MMP, TIMP, LOX,
#                            ADAMTS).
#   secretory machinery   -> STAYS on GO. The background there already holds
#                            1 127 genes; there is nothing to broaden.
#   TGFb/BMP              -> Reactome, the TGFb family and BMP together.
#   hypoxia/stress        -> Reactome, cellular response to hypoxia, in
#                            addition to the existing GO set.
#   cell-cycle exit       -> the existing GO terms PLUS Reactome cellular
#                            senescence.
#
# As a SECOND source for the two matrix axes, independent of the matrisome,
# Reactome sets are stored as well (ECM organization, degradation of the ECM,
# collagen formation, ECM proteoglycans). They run as a sensitivity analysis,
# not as a primary analysis -- this is stated in the header of 24_ and does
# not change the Bonferroni factor of the primary analysis.
#
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION is deliberately NOT used: that
# set mixes matrix, cytoskeleton and signalling pathways and carries a cancer
# biology history that has no place here.
#
# Gene bridge: for EVERY set -- narrow and broad alike -- the members are
# translated to Ensembl through the Entrez identifier and `org.Hs.egENSEMBL`,
# that is, over exactly the same bridge as in the narrow run. The narrow and
# the broad run therefore differ in set content and in nothing else.
# =============================================================================
suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(GO.db); library(msigdbr)
})

.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
ABLAGE <- file.path(WURZEL, "data_raw", "_referenz", "gensaetze_v2")

if (dir.exists(ABLAGE) && !nzchar(Sys.getenv("GENSAETZE_V2_NEU"))) {
  cat("Store already exists -- nothing written.\n",
      "To rebuild: GENSAETZE_V2_NEU=1\n", sep = "")
  quit(save = "no", status = 0)
}
dir.create(ABLAGE, recursive = TRUE, showWarnings = FALSE)

ABRUF <- format(Sys.Date())

# ---- bridge Entrez -> Ensembl (the same as in the narrow run) --------------
e2ens <- as.list(org.Hs.egENSEMBL)
nach_ensembl <- function(entrez) {
  unique(na.omit(unlist(e2ens[unique(as.character(entrez))], use.names = FALSE)))
}

g2e <- as.list(org.Hs.egGO2ALLEGS)
GOr <- function(id) {
  eg <- g2e[[id]]
  if (is.null(eg)) return(character(0))
  nach_ensembl(unname(eg))
}

MSIG <- rbind(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
MSIG_VERSION <- unique(MSIG$db_version)
stopifnot(length(MSIG_VERSION) == 1)
MS <- function(name) {
  z <- MSIG[MSIG$gs_name == name, ]
  stopifnot(nrow(z) > 0)
  nach_ensembl(z$ncbi_gene)
}

# ---- the sets ---------------------------------------------------------------
# One row each: category, variant ("eng" | "breit" | "empfindlichkeit"),
# source, version, members.
DEF <- list(
  list("1_Matrixbestandteile", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = "GO:0005201")),
  list("1_Matrixbestandteile", "breit", "MSigDB / Matrisom (Naba et al.)",
       MSIG_VERSION,
       list(MSIGDB = c("NABA_COLLAGENS", "NABA_ECM_GLYCOPROTEINS",
                       "NABA_PROTEOGLYCANS"))),
  list("1_Matrixbestandteile", "empfindlichkeit", "MSigDB / Reactome",
       MSIG_VERSION,
       list(MSIGDB = c("REACTOME_COLLAGEN_FORMATION",
                       "REACTOME_ECM_PROTEOGLYCANS"))),

  list("2_Matrixremodellierung", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = "GO:0022617")),
  list("2_Matrixremodellierung", "breit", "MSigDB / Matrisom (Naba et al.)",
       MSIG_VERSION, list(MSIGDB = "NABA_ECM_REGULATORS")),
  list("2_Matrixremodellierung", "empfindlichkeit", "MSigDB / Reactome",
       MSIG_VERSION,
       list(MSIGDB = c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
                       "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"))),

  list("3_Sekretionsmaschine", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = c("GO:0016192", "GO:0006888"))),
  list("3_Sekretionsmaschine", "breit", "GO via org.Hs.eg.db (unveraendert)",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = c("GO:0016192", "GO:0006888"))),

  list("4_TGFb_BMP", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = "GO:0007179")),
  list("4_TGFb_BMP", "breit", "MSigDB / Reactome", MSIG_VERSION,
       list(MSIGDB = c("REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
                       "REACTOME_SIGNALING_BY_BMP"))),

  list("5_Hypoxie_Stress", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = "GO:0036293")),
  list("5_Hypoxie_Stress", "breit", "GO + MSigDB / Reactome",
       paste(as.character(packageVersion("org.Hs.eg.db")), MSIG_VERSION),
       list(GO = "GO:0036293",
            MSIGDB = c("REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
                       "REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR"))),

  list("6_Zellzyklusausstieg", "eng", "GO via org.Hs.eg.db",
       as.character(packageVersion("org.Hs.eg.db")),
       list(GO = c("GO:0008285", "GO:0090398"))),
  list("6_Zellzyklusausstieg", "breit", "GO + MSigDB / Reactome",
       paste(as.character(packageVersion("org.Hs.eg.db")), MSIG_VERSION),
       list(GO = c("GO:0008285", "GO:0090398"),
            MSIGDB = "REACTOME_CELLULAR_SENESCENCE"))
)

zeilen <- list()
uebersicht <- list()
for (d in DEF) {
  kat <- d[[1]]; var <- d[[2]]; quelle <- d[[3]]; version <- d[[4]]
  teile <- d[[5]]
  g <- character(0)
  for (typ in names(teile)) {
    for (id in teile[[typ]]) {
      g <- union(g, if (typ == "GO") GOr(id) else MS(id))
    }
  }
  g <- sort(unique(g))
  bestandteile <- paste(unlist(teile), collapse = " | ")
  zeilen[[length(zeilen) + 1]] <- data.frame(
    kategorie = kat, variante = var, ensembl = g, stringsAsFactors = FALSE)
  uebersicht[[length(uebersicht) + 1]] <- data.frame(
    kategorie = kat, variante = var, quelle = quelle, version = version,
    bestandteile = bestandteile, n_gene = length(g), abrufdatum = ABRUF,
    stringsAsFactors = FALSE)
  cat(sprintf("%-24s %-16s n=%5d  %s\n", kat, var, length(g), bestandteile))
}

MIT <- do.call(rbind, zeilen)
UEB <- do.call(rbind, uebersicht)
write.csv(MIT, file.path(ABLAGE, "mitglieder.csv"), row.names = FALSE)
write.csv(UEB, file.path(ABLAGE, "uebersicht.csv"), row.names = FALSE)

writeLines(c(
  "# Gene sets v2 -- sources, versions, retrieval date",
  "",
  sprintf("Created %s by `code/24a_gene_sets_v2_build.R`. **Not changed after",
          ABRUF),
  "creation.** The analysis (`code/24_gene_sets_v2.R`) reads",
  "`mitglieder.csv` only and recomputes nothing.",
  "",
  "| Source | Version | Retrieved |",
  "|---|---|---|",
  sprintf("| GO via `org.Hs.eg.db` / `org.Hs.egGO2ALLEGS` | %s | %s |",
          as.character(packageVersion("org.Hs.eg.db")), ABRUF),
  sprintf("| MSigDB via `msigdbr` (Matrisome NABA_*, Reactome) | %s (msigdbr %s) | %s |",
          MSIG_VERSION, as.character(packageVersion("msigdbr")), ABRUF),
  "",
  "Matrisome: Naba A. *et al.*, PMID 22159717, <http://matrisome.org>.",
  "The gene bridge is the same for every set: Entrez ID ->",
  "`org.Hs.egENSEMBL`. The narrow and the broad run therefore differ in",
  "set content and in nothing else.",
  "",
  "## Sets",
  "",
  "| Category | Variant | Source | Components | n genes |",
  "|---|---|---|---|---|",
  sprintf("| %s | %s | %s | `%s` | %d |", UEB$kategorie, UEB$variante,
          UEB$quelle, UEB$bestandteile, UEB$n_gene)
), file.path(ABLAGE, "QUELLEN.md"))

cat(sprintf("\n-> %s\n", file.path(ABLAGE, "mitglieder.csv")))
