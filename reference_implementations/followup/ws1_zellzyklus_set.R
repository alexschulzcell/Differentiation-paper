suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
g2e <- as.list(org.Hs.egGO2ALLEGS); e2ens <- as.list(org.Hs.egENSEMBL)
GOr <- function(id) { eg <- g2e[[id]]; if (is.null(eg)) return(character(0))
  unique(na.omit(unlist(e2ens[unique(unname(eg))], use.names = FALSE))) }
S_ZYKLUS <- GOr("GO:0007049")
cat("cell cycle (GO:0007049):", length(S_ZYKLUS), "genes\n")
write.csv(data.frame(ensembl = S_ZYKLUS, satz = "S_ZELLZYKLUS"),
          "derived_data/followup/ws1_zellzyklus.csv", row.names = FALSE)
