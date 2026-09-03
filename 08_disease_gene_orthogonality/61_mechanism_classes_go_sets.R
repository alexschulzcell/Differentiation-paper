# =============================================================================
# 61_mechanism_classes_go_sets.R -- work package WS2: mechanism classes of the disease genes.
#
# Exports seven external, versioned GO gene sets (one mechanism axis per set).
# Same acquisition route as `08_disease_gene_orthogonality/10_secretion_go_sets.R`
# (org.Hs.egGO2ALLEGS -> Ensembl, GO.db version as installed in the project).
# Multiple membership is explicitly allowed -- the sets are NOT cleaned
# against each other (unlike S_DISTAL/S_BIOSYN in 52a, which form a fixed
# positive control and remain unchanged).
#
# This is NOT a search for a new convergence axis: the sets classify only the
# DISEASE-GENE side (panels), not the differentiation program. The program
# (173/147 genes) and its definition remain untouched.
# =============================================================================
suppressPackageStartupMessages({library(org.Hs.eg.db); library(GO.db)})

g2e <- as.list(org.Hs.egGO2ALLEGS); e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) { eg <- g2e[[id]]; if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE))) }

klassen <- list(
  ECM_STRUKTUR   = c("GO:0005201",   # extracellular matrix structural constituent (MF)
                      "GO:0030198"), # extracellular matrix organization (BP)
  GLYKO_LINKER   = c("GO:0015012",   # heparan sulfate proteoglycan biosynthetic process
                      "GO:0050650",  # chondroitin sulfate proteoglycan biosynthetic process
                      "GO:0030166"), # proteoglycan biosynthetic process
  TF_DNABINDEND  = c("GO:0003700"),  # DNA-binding transcription factor activity
  CILIUM         = c("GO:0060271",   # cilium assembly
                      "GO:0005929"), # cilium (CC)
  SIGNAL_FWBI    = c("GO:0008543",   # fibroblast growth factor receptor signaling pathway
                      "GO:0016055",  # Wnt signaling pathway
                      "GO:0030509",  # BMP signaling pathway
                      "GO:0007224"), # smoothened signaling pathway (Hedgehog/IHH)
  VESIKEL_SEKRET = c("GO:0016192",   # vesicle-mediated transport
                      "GO:0006887"), # exocytosis
  LYSOSOM        = c("GO:0005764")   # lysosome (CC)
)

aus <- lapply(names(klassen), function(kn) {
  ens <- unique(unlist(lapply(klassen[[kn]], GOr)))
  data.frame(klasse = kn, ensembl = ens)
})
D <- do.call(rbind, aus)

dir.create("derived_data/followup", recursive = TRUE, showWarnings = FALSE)
write.csv(D, "derived_data/followup/ws2_mechanismusklassen_go.csv", row.names = FALSE)

cat("mechanism classes (genes per set, GO2ALLEGS -> Ensembl):\n")
for (kn in names(klassen)) {
  cat(sprintf("  %-14s %5d genes  (GO terms: %s)\n", kn,
              sum(D$klasse == kn), paste(klassen[[kn]], collapse = ", ")))
}

# overlap matrix (multiple membership is allowed and is reported)
kn <- names(klassen)
M <- matrix(0L, length(kn), length(kn), dimnames = list(kn, kn))
sets <- split(D$ensembl, D$klasse)
for (i in kn) for (j in kn) M[i, j] <- length(intersect(sets[[i]], sets[[j]]))
cat("\noverlap matrix (number of shared genes, diagonal = set size):\n")
print(M)
write.csv(as.data.frame(M), "derived_data/followup/ws2_ueberlappungsmatrix.csv")
