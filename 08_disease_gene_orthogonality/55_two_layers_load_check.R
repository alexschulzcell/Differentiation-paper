NUR_PUNKTE <- 1:2
source("02_matrix_programme_derivation/11_load_18_datasets.R")
for (i in seq_along(DATEN)) {
  cat("== point", i, LABELS_18[i], "==\n")
  cat("arm:", DATEN[[i]]$arm, " klasse:", DATEN[[i]]$klasse, "\n")
  cat("meta cols:", paste(colnames(DATEN[[i]]$meta), collapse=","), "\n")
  print(head(DATEN[[i]]$meta,3))
  cat("expr len:", length(DATEN[[i]]$expr), " Z dim:", paste(dim(DATEN[[i]]$Z),collapse="x"), "\n\n")
}
