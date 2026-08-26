NUR_PUNKTE <- 1:2
source("reference_implementations/manuscript/methods/13_load18.R")
cat(head(names(DATEN[[1]]$expr)),"\n")
cat(head(rownames(DATEN[[1]]$Z)),"\n")
cat(class(DATEN[[1]]$klasse), DATEN[[1]]$klasse, "\n")
