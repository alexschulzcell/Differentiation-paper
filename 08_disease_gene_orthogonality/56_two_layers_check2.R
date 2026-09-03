NUR_PUNKTE <- 1:4
source("02_matrix_programme_derivation/11_load_18_datasets.R")
for (i in seq_along(DATEN)) {
  cat("== point", i, LABELS_18[i], "arm=",DATEN[[i]]$arm,"==\n")
  print(table(DATEN[[i]]$meta$condition))
}
