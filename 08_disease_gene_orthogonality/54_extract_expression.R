# =============================================================================
# 54_extract_expression.R -- WS1: absolute expression per gene from the 18 datasets,
# blockwise (memory), incl. rank per dataset (scale-free, no mean over
# cohorts on the raw scale needed). Writes one CSV per point:
#   ensembl, expr, expr_rank, arm, klasse, punkt, label, n_undiff, n_diff
# to derived_data/followup/ws1_expr/punkt_<i>.csv
# =============================================================================
AUS <- "derived_data/followup/ws1_expr"
dir.create(AUS, showWarnings = FALSE, recursive = TRUE)

bloecke <- list(1:5, 6:11, 12:18)
for (blk in bloecke) {
  NUR_PUNKTE <- blk
  source("02_matrix_programme_derivation/11_load_18_datasets.R")
  for (i in seq_along(DATEN)) {
    d <- DATEN[[i]]
    idx <- blk[i]  # since DATEN in this session only contains those of the block, in order
    ex <- d$expr
    r <- rank(ex, na.last = "keep") / sum(!is.na(ex))
    cond <- if (!is.null(d$meta$condition)) table(d$meta$condition) else NA
    n_undiff <- if (!is.null(d$meta$condition)) sum(d$meta$condition == "undiff") else NA
    n_diff   <- if (!is.null(d$meta$condition)) sum(d$meta$condition != "undiff") else NA
    df <- data.frame(ensembl = names(ex), expr = as.numeric(ex),
                      expr_rank = as.numeric(r[names(ex)]),
                      arm = d$arm, klasse = d$klasse,
                      punkt = idx, label = d$label,
                      n_undiff = n_undiff, n_diff = n_diff)
    fn <- file.path(AUS, sprintf("punkt_%02d.csv", idx))
    write.csv(df, fn, row.names = FALSE)
    cat(sprintf("written: %s (%d genes)\n", fn, nrow(df)))
  }
  rm(DATEN); invisible(gc())
}
cat("done.\n")
