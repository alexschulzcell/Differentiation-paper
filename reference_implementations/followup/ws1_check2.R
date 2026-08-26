NUR_PUNKTE <- 1:4
source("reference_implementations/manuscript/methods/13_load18.R")
for (i in seq_along(DATEN)) {
  cat("== point", i, LABELS_18[i], "arm=",DATEN[[i]]$arm,"==\n")
  print(table(DATEN[[i]]$meta$condition))
}
