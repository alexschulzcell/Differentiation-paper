# =============================================================================
# ws2_go_klassen.R -- WS2: seven mechanism classes of the skeletal-dysplasia
# disease genes, defined via external, versioned GO term sets.
#
# Follows LITERALLY the same pattern as reference_implementations/52a_go_sets.R (the same
# packages org.Hs.eg.db + GO.db, the same GOr() mechanism: GO term -> all
# EntrezGene via GO2ALLEGS [incl. all child terms] -> Ensembl). Nothing is
# changed about the program definition (173/147 genes) -- what is classified
# is exclusively the disease-gene side (the external panels).
#
# Seven axes, with the GO terms they carry:
#   structure/ECM              GO:0005201 (ECM structural constituent),
#                              GO:0031012 (extracellular matrix)
#   glycosylation/PG-linker    GO:0006486 (protein glycosylation),
#                              GO:0030203 (glycosaminoglycan metabolic process),
#                              GO:0015012 (heparan sulfate PG biosynthesis),
#                              GO:0050650 (chondroitin sulfate PG biosynthesis)
#   TF/DNA-binding             GO:0003700 (DNA-binding TF activity)
#   ciliopathy/cilium          GO:0005929 (cilium), GO:0060271 (cilium assembly)
#   signaling FGF/WNT/BMP/IHH  GO:0008543 (FGFR signaling),
#                              GO:0016055 (Wnt signaling),
#                              GO:0030509 (BMP signaling),
#                              GO:0007224 (smoothened/Hedgehog signaling)
#   vesicle transport/secretion GO:0016192 (vesicle-mediated transport),
#                              GO:0006887 (exocytosis)
#   lysosome                   GO:0005764 (lysosome)
# =============================================================================
suppressPackageStartupMessages({library(org.Hs.eg.db); library(GO.db)})

g2e <- as.list(org.Hs.egGO2ALLEGS); e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) { eg <- g2e[[id]]; if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE))) }

klassen <- list(
  Struktur_ECM   = c("GO:0005201", "GO:0031012"),
  Glykosylierung = c("GO:0006486", "GO:0030203", "GO:0015012", "GO:0050650"),
  TF_DNAbindend  = c("GO:0003700"),
  Ziliopathie    = c("GO:0005929", "GO:0060271"),
  Signal_FGF_WNT_BMP_IHH = c("GO:0008543", "GO:0016055", "GO:0030509", "GO:0007224"),
  Vesikel_Sekretion = c("GO:0016192", "GO:0006887"),
  Lysosom        = c("GO:0005764")
)

aus <- do.call(rbind, lapply(names(klassen), function(kl) {
  ens <- unique(unlist(lapply(klassen[[kl]], GOr)))
  data.frame(klasse = kl, ensembl = ens, stringsAsFactors = FALSE)
}))

dir.create("derived_data/followup", recursive = TRUE, showWarnings = FALSE)
write.csv(aus, "derived_data/followup/ws2_go_klassen.csv", row.names = FALSE)

cat("mechanism classes (GO sets, all genes in the genome):\n")
for (kl in names(klassen)) {
  n <- sum(aus$klasse == kl)
  cat(sprintf("  %-26s %5d genes | GO terms: %s\n", kl, n,
              paste(klassen[[kl]], collapse = ", ")))
}
