# install_r_packages.R -- install the R packages the figure chain needs.
#
#   Rscript install_r_packages.R
#
# r_packages.txt is a record of the versions actually used, in the
# `name_version` form that sessionInfo() prints. It is therefore not a list of
# installable names, and `install.packages(readLines("r_packages.txt"))` does
# not work: it would try to install "ggplot2_4.0.3" and the comment lines.
# This script reads that record, strips the version suffix, and routes each
# package to CRAN or to Bioconductor.
#
# Only four of the recorded packages are needed for `python reproduce.py`
# (ggplot2, patchwork, ragg, systemfonts). The rest belong to the stages
# that read the raw data or query MSigDB; r_packages.txt says which is
# which, and this script installs whatever is missing from either group.
#
# Versions are not pinned here. The figure scripts do not depend on the exact
# versions; r_packages.txt states the ones the published figures were drawn
# with, should an exact reconstruction be needed.

CRAN_MIRROR <- "https://cloud.r-project.org"

# Bioconductor packages; install.packages() cannot find them. matrixStats and
# msigdbr look like Bioconductor packages but are both on CRAN.
BIOC <- c("DESeq2", "rtracklayer", "GenomicRanges", "GenomicFeatures",
          "AnnotationDbi", "org.Hs.eg.db", "GO.db",
          "TxDb.Hsapiens.UCSC.hg38.knownGene")

root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
record <- file.path(root, "r_packages.txt")
if (!file.exists(record)) stop("r_packages.txt not found next to this script")

lines <- readLines(record, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
# "ggplot2_4.0.3" -> "ggplot2"; a name without a version is left alone.
wanted <- unique(sub("_.*$", "", lines))
if (!length(wanted)) stop("no package names found in r_packages.txt")

cat("Packages recorded in r_packages.txt:\n  ", paste(wanted, collapse = ", "),
    "\n\n", sep = "")

missing <- wanted[!vapply(wanted, requireNamespace, logical(1),
                          quietly = TRUE)]
if (!length(missing)) {
  cat("All recorded packages are already installed. Nothing to do.\n")
  quit(save = "no", status = 0)
}
cat("Missing, will be installed:\n  ", paste(missing, collapse = ", "),
    "\n\n", sep = "")

from_cran <- setdiff(missing, BIOC)
from_bioc <- intersect(missing, BIOC)

if (length(from_cran)) {
  cat("--- CRAN:", paste(from_cran, collapse = ", "), "---\n")
  install.packages(from_cran, repos = CRAN_MIRROR)
}

if (length(from_bioc)) {
  cat("--- Bioconductor:", paste(from_bioc, collapse = ", "), "---\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = CRAN_MIRROR)
  }
  BiocManager::install(from_bioc, update = FALSE, ask = FALSE)
}

# Report rather than assume: a compile failure on Linux is usually a missing
# system library, and the reviewer needs to see which package it was.
still <- wanted[!vapply(wanted, requireNamespace, logical(1), quietly = TRUE)]
cat("\n")
if (length(still)) {
  cat("STILL MISSING:\n  ", paste(still, collapse = ", "), "\n", sep = "")
  cat("These usually fail on a missing system library. ragg and systemfonts\n")
  cat("need freetype, harfbuzz and fribidi; on Debian/Ubuntu:\n")
  cat("  sudo apt-get install libfreetype6-dev libharfbuzz-dev ",
      "libfribidi-dev libpng-dev\n", sep = "")
  quit(save = "no", status = 1)
}
cat("All recorded R packages are installed.\n")
