# =============================================================================
# 10_secretion_go_sets.R -- exports the GO sets of the positive control M-A.
#
# The definition is LITERALLY the one from
# `02_matrix_programme_derivation/12_metric_reference.R` and is not
# changed here -- otherwise the anchor (OR 4.97) would not be the same anchor.
# =============================================================================
suppressPackageStartupMessages({library(org.Hs.eg.db); library(GO.db)})

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
ueb <- intersect(S_MASCHINE, S_FRACHT)
S_DISTAL <- setdiff(S_DISTAL, intersect(S_DISTAL, S_FRACHT))
S_BIOSYN <- setdiff(S_BIOSYN, intersect(S_BIOSYN, S_FRACHT))
S_BIOSYN <- setdiff(S_BIOSYN, S_DISTAL)

dir.create("derived_data/M_humangenetik", recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(satz = c(rep("S_DISTAL", length(S_DISTAL)),
                              rep("S_BIOSYN", length(S_BIOSYN))),
                     ensembl = c(S_DISTAL, S_BIOSYN)),
          "derived_data/M_humangenetik/go_saetze.csv", row.names = FALSE)
cat(sprintf("S_DISTAL %d | S_BIOSYN %d genes\n",
            length(S_DISTAL), length(S_BIOSYN)))
