# =============================================================================
# 20_figures_main.R -- the six main figures
# =============================================================================
# Purpose  Draws F1 to F6 from the panel CSV files in figures/data/. The style
#          comes unchanged from figure_style/publication_style.R; PUB_DIR is
#          overwritten after the source() call, so the style script itself is
#          never touched.
#
# Inputs   figures/data/*.csv  (written by 09_figures/10_panel_data_main.py)
# Outputs  figures/F1.pdf to F6.pdf and the same as PNG at 600 dpi
#          results/abbildungen_haupt_sessionInfo.txt
# Runtime  about a minute
#
# The rules from figure_style/FIGURE_RULES.md that apply here:
#   * no title and no prose inside the image -- everything explanatory goes
#     into the legend
#   * panel letters are placed on the sheet with grid.text, not via plot.tag
#   * no scale_*_continuous(limits=); the view is set with coord_cartesian
#   * no colour vectors passed to a geom; always aes() plus scale_*_manual()
#   * the key is a direct label in the field, not a ggplot legend
# =============================================================================
suppressMessages({
  library(ggplot2); library(patchwork); library(grid)
})

.skriptordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
WURZEL <- Sys.getenv("PAPER_V2_ROOT")
if (!nzchar(WURZEL)) WURZEL <- dirname(.skriptordner())
source(file.path(WURZEL, "figure_style", "publication_style.R"))

# PUB_DIR is set after the source() call; the style script keeps its own default
PUB_DIR <- file.path(WURZEL, "figures")
DAT     <- file.path(WURZEL, "figures", "data")
dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)

# pandas writes booleans as "True"/"False", and read.csv reads them as
# character. They are turned into real logicals centrally here; otherwise
# every `!column` and every ifelse() would fail silently.
lies <- function(n) {
  d <- read.csv(file.path(DAT, paste0(n, ".csv")), check.names = FALSE)
  for (k in names(d)) {
    if (is.character(d[[k]])) {
      u <- unique(d[[k]][!is.na(d[[k]])])
      if (length(u) > 0 && all(u %in% c("True", "False")))
        d[[k]] <- d[[k]] == "True"
    }
  }
  d
}

# ------------------------------------------------------------- colour tokens
# Colour codes the status of a calibration or the level -- never the genotype.
EICH <- c("passed" = "#12946B", "failed" = "#C2472A",
          "not calibratable" = GRAU)
ACHSE <- c("osteogenic" = "#E07B12", "adipogenic" = "#7B5EA7",
           "lineage contrast" = GRAU)
KURVE <- c("programme" = "#0E8C8C", "disease genes" = "#1F6FB2")

# Place the panel letters on the finished sheet
tafel <- function(p, marken, name, breite, hoehe) {
  for (typ in c("pdf", "png")) {
    if (typ == "pdf")
      grDevices::cairo_pdf(file.path(PUB_DIR, paste0(name, ".pdf")),
                           width = mm(breite), height = mm(hoehe), family = FONT)
    else
      ragg::agg_png(file.path(PUB_DIR, paste0(name, ".png")),
                    width = mm(breite), height = mm(hoehe), units = "in",
                    res = DPI, background = "white")
    print(p)
    for (i in seq_len(nrow(marken)))
      grid::grid.text(marken$lab[i],
                      x = grid::unit(marken$x[i], "mm"),
                      y = grid::unit(1, "npc") - grid::unit(marken$y[i], "mm"),
                      just = c("left", "top"),
                      gp = grid::gpar(fontfamily = FONT, fontface = "bold",
                                      fontsize = PTL, col = INK))
    invisible(grDevices::dev.off())
  }
  cat(sprintf("  %-4s %3.0f x %3.0f mm\n", name, breite, hoehe))
}

# theme_pub() sets axis.text.x and axis.text.y SEPARATELY (because of the
# margins). A blanket axis.text = element_blank() is overridden by that, which
# left the axis labels standing in the schematic panels. Hence explicitly here.
theme_leer <- function() {
  theme_pub(links = FALSE, unten = FALSE) +
    theme(axis.text = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(2, 2, 2, 2))
}

m <- function(...) do.call(rbind, lapply(list(...), function(v)
  data.frame(lab = v[[1]], x = as.numeric(v[[2]]), y = as.numeric(v[[3]]))))

# =============================================================================
# FIGURE 1 -- material, the screens, and the calibration almost nothing passes
# =============================================================================
f1a_schema <- function() {
  # The two layers as a schematic -- drawn with ggplot on purpose, so that the
  # type and the line weights are identical to those of the data panels.
  kaesten <- data.frame(
    x = c(1, 1), y = c(2, 1), w = c(1.9, 1.9), h = c(0.62, 0.62),
    lab = c("matrix programme\n173 genes", "disease genes\nPanelApp 309"),
    key = c("programme", "disease genes"))
  ggplot() +
    geom_rect(data = kaesten,
              aes(xmin = x - w / 2, xmax = x + w / 2,
                  ymin = y - h / 2, ymax = y + h / 2, colour = key),
              fill = "white", linewidth = LW * 2) +
    geom_text(data = kaesten, aes(x, y, label = lab, colour = key),
              family = FONT, size = gs(PT), lineheight = 0.95) +
    scale_colour_manual(values = KURVE, guide = "none") +
    annotate("segment", x = 1, xend = 1, y = 1.65, yend = 1.35,
             colour = INK, linewidth = LW,
             arrow = arrow(length = unit(1.4, "mm"), ends = "both",
                           type = "closed")) +
    txt(1.08, 1.5, "orthogonal", size = PTS, hjust = 0) +
    txt(0.05, 2.55, "what the cell does", size = PT, hjust = 0) +
    txt(0.05, 0.45, "the machinery that executes it", size = PT, hjust = 0) +
    coord_cartesian(xlim = c(0, 2.1), ylim = c(0.3, 2.7), clip = "off") +
    theme_leer()
}

f1b <- function() {
  d <- lies("F1B_screen_perturbation")
  d$step <- factor(d$step, levels = rev(d$step))
  d$anzeige <- abs(d$n)
  d$type <- ifelse(d$n < 0, "excluded", "kept")
  ggplot(d, aes(anzeige, step, fill = type)) +
    geom_col(width = 0.62, colour = NA) +
    geom_text(aes(anzeige + 3, step, label = abs(n)), hjust = 0,
              family = FONT, size = gs(PTS), colour = INK) +
    scale_fill_manual(values = c(excluded = GRAU, kept = "#1F6FB2"),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = "series", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 8, 1.5, 2))
}

f1c <- function() {
  d <- lies("F1C_screen_diagnoses")
  d$step <- factor(d$step, levels = rev(d$step))
  ggplot(d, aes(n, step, colour = level)) +
    geom_segment(aes(x = 1, xend = n, yend = step), linewidth = LWD) +
    geom_point(size = 1.7) +
  # The count used to sit at the point and cut through it at small values (the
  # "2" of "complete 2x2"). It now stands in a fixed column to the right of
  # the data area -- the same solution as in F4D and F4E.
    geom_text(aes(x = 4200, label = n), hjust = 0, family = FONT,
              size = gs(PTS), colour = INK) +
    scale_colour_manual(values = c(Entity = "#1F6FB2",
                                   Predisposition = "#0E8C8C"),
                        guide = "none") +
    scale_x_continuous(transform = "log10", breaks = c(1, 10, 100, 1000)) +
    labs(x = "series (log scale)", y = NULL) +
    coord_cartesian(xlim = c(1, 3000), clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 13, 1.5, 2))
}

f1d <- function() {
  d <- lies("F1D_calibration_per_dataset")
  d$calibration <- ifelse(!d$calibratable, "not calibratable",
                      ifelse(d$passed, "passed", "failed"))
  d$z_plot <- ifelse(is.na(d$z), 0, d$z)
  d$dataset <- factor(d$dataset, levels = d$dataset)
  zebra <- data.frame(y = seq(1, nrow(d), by = 2))
  ggplot(d, aes(z_plot, dataset)) +
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = SOFT) +
    geom_segment(aes(x = 0, xend = z_plot, yend = dataset, colour = calibration),
                 linewidth = LWD) +
    geom_point(aes(colour = calibration, shape = calibratable), size = 1.5,
               fill = "white", stroke = LW * 2) +
    scale_colour_manual(values = EICH, guide = "none") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21), guide = "none") +
    # The zero line and the threshold line end below the threshold label; as
    # geom_vline they ran straight through it. They stand AFTER geom_point so
    # that the y scale stays discrete.
    annotate("segment", x = c(0, 2), xend = c(0, 2),
             y = 0.45, yend = nrow(d) + 0.5,
             colour = c(LINE, DIFF), linewidth = LW,
             linetype = c("solid", "22")) +
    # The threshold label used to sit on the bottom data row. It now stands
    # ABOVE the top row in the free margin. It is kept short: since the row
    # names carry their accessions the labels of panel E are longer, and a
    # long note from D would run into them.
    txt(-1.55, nrow(d) + 0.9, "dashed: threshold z = 2",
        colour = DIFF, size = PTS, hjust = 0, vjust = 0.5) +
    labs(x = "lineage-marker calibration (z)", y = NULL) +
    coord_cartesian(xlim = c(-1.6, 3.0), ylim = c(0.4, nrow(d) + 1.3),
                    clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))
}

f1e <- function() {
  d <- lies("F1E_calibration_per_cell")
  d$calibration <- ifelse(d$passed, "passed", "failed")
  d$cell <- factor(d$cell, levels = d$cell)
  # The presentation layer writes "engineered" or "not a patient defect";
  # an earlier pattern also matched the negation and drew the SERPINA3 cells
  # as patient lesions by mistake.
  d$echt <- ifelse(grepl("engineered", d$lesion),
                   "engineering", "patient lesion")
  zebra <- data.frame(y = seq(1, nrow(d), by = 2))
  ggplot(d, aes(z, cell)) +
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = SOFT) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = 2, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = z, yend = cell, colour = calibration),
                 linewidth = LWD) +
    geom_point(aes(colour = calibration, shape = echt), size = 1.5,
               fill = "white", stroke = LW * 2) +
    scale_colour_manual(values = EICH, guide = "none") +
    scale_shape_manual(values = c("patient lesion" = 23,
                                  "engineering" = 16), guide = "none") +
    labs(x = "lineage-marker calibration (z)", y = NULL) +
    coord_cartesian(xlim = c(-1.6, 6.2), clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))
}

bau_f1 <- function() {
  # The row names of D and E carry their accessions and are noticeably longer;
  # the sheet gets 8 mm more height for them.
  p <- (f1a_schema() | f1b()) / (f1c() | plot_spacer()) / (f1d() | f1e()) +
    plot_layout(heights = c(1, 0.85, 1.9))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 47), list("D", 1.5, 95), list("E", 89, 95)),
        "F1", SP2, 182)
}

# =============================================================================
# FIGURE 2 -- the main finding
# =============================================================================
MATCHED <- "#7B5EA7"   # the matched (confounder-controlled) null

f2a <- function() {
  # Cross-arm concordance against TWO nulls on one axis: random background genes
  # (grey area) and random genes matched on expression, length and constraint
  # (purple line). The observation is the same needle for both.
  d  <- lies("F2A_cross_arm_concordance")             # unmatched null + rho 0.622
  sm <- lies("F2A_matched_nulls_summary")
  smm <- sm[sm$null_type == "matched", ]
  x <- seq(-0.2, 0.8, length.out = 500)
  yu <- dnorm(x, d$null_mean, d$null_sd)
  ym <- dnorm(x, smm$null_mean, smm$null_sd)
  sc <- max(yu, ym)
  du <- data.frame(x = x, y = yu / sc)
  dm <- data.frame(x = x, y = ym / sc)
  ggplot() +
    geom_area(data = du, aes(x, y), fill = SOFT, colour = GRAU, linewidth = LW) +
    geom_line(data = dm, aes(x, y), colour = MATCHED, linewidth = LWD) +
    annotate("segment", x = d$rho_module_genes, xend = d$rho_module_genes,
             y = 0, yend = 0.80, colour = "#1F6FB2", linewidth = LWD) +
    geom_point(x = d$rho_module_genes, y = 0.80, colour = "#1F6FB2", size = 1.6) +
    txt(d$rho_module_genes, 0.87, sprintf("rho = %.3f", d$rho_module_genes),
        size = PTS, hjust = 0.5, vjust = 0) +
    txt(0.30, 1.10, sprintf("vs random genes: z %+.2f", d$z),
        size = PTS, hjust = 0, colour = INK) +
    txt(0.30, 0.98, sprintf("vs matched set: z %+.2f", smm$z_sd_units),
        size = PTS, hjust = 0, colour = MATCHED) +
    labs(x = "cross-arm concordance of module genes (Spearman rho)",
         y = "null density") +
    coord_cartesian(xlim = c(0, 0.72), ylim = c(0, 1.20),
                    expand = FALSE, clip = "off") +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(6, 4, 1.5, 2))
}

f2b <- function() {
  d <- lies("F2B_module_per_dataset")
  d <- d[order(d$concordance), ]
  d$dataset <- factor(d$dataset, levels = d$dataset)
  zebra <- data.frame(y = seq(1, nrow(d), by = 2))
  ggplot(d, aes(concordance, dataset)) +
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = SOFT) +
    geom_segment(aes(x = concordance_detection_limit, xend = concordance,
                     yend = dataset), colour = GRAU, linewidth = LW) +
    geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.0,
               colour = INK) +
    geom_point(aes(colour = calibration), size = 1.7) +
    scale_colour_manual(values = EICH, guide = "none") +
    # No key inside the field: this sheet has no free area anywhere, and a key
    # sitting on the data rows is worse than none. The colour assignment is in
    # the legend, and panel C carries it as a direct label as well.
    labs(x = "module concordance", y = NULL) +
    coord_cartesian(xlim = c(0.48, 1.02), clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 3, 1.5, 2))
}

f2c <- function() {
  d <- lies("F2C_pooled_by_calibration")
  d$label <- factor(c("calibration\nfailed", "calibration\npassed")[
    match(d$label, c("failed", "passed"))],
    levels = c("calibration\npassed", "calibration\nfailed"))
  d$calibration <- ifelse(grepl("passed", d$label), "passed", "failed")
  ggplot(d, aes(concordance, label)) +
    geom_segment(aes(x = concordance_null, xend = concordance, yend = label),
                 colour = GRAU, linewidth = LW) +
    geom_errorbarh(aes(xmin = concordance_null - concordance_null_sd,
                       xmax = concordance_null + concordance_null_sd),
                   height = 0.12, colour = GRAU, linewidth = LW) +
    geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.4,
               colour = INK) +
    geom_point(aes(colour = calibration), size = 2.0) +
    geom_text(aes(x = concordance + 0.03,
                  label = sprintf("z %+.2f (n %d)", concordance_z,
                                  n_datasets)),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK) +
    scale_colour_manual(values = EICH, guide = "none") +
    labs(x = "pooled module concordance", y = NULL) +
    # The statistic is bounded on [0,1] -- an axis running to 1.26 would
    # suggest a range that does not exist. The z label runs past the right
    # edge (clip = "off" plus a wide margin).
    scale_x_continuous(breaks = c(0.5, 0.75, 1.0)) +
    # The z label used to overhang by 34 mm, down to under panel D. The view
    # now ends BEHIND the label and the margin stays narrow.
    coord_cartesian(xlim = c(0.44, 1.32), clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 12, 1.5, 2))
}

f2d <- function() {
  # Two gene sets per category side by side: the broad, independently curated
  # one (matrisome or Reactome, filled) and the narrow GO set of the first
  # pass (open). Nothing is replaced without both staying visible -- the
  # decision rule is in 06_orthogonal_layers/61_gene_set_enrichment.R. The sensitivity rows (the
  # second Reactome source) are in the panel file and in Table S9, not here.
  d <- lies("F2D_gene_sets_v2")
  d <- d[d$variant %in% c("narrow", "broad"), ]
  reihe <- d[d$variant == "broad", ]
  reihe <- reihe[order(reihe$odds_ratio), ]
  d$category <- factor(d$category, levels = reihe$category)
  d$gene_set <- d$variant
  d$direction <- ifelse(d$odds_ratio > 1, "enriched", "depleted")
  d$y <- as.numeric(d$category) + ifelse(d$gene_set == "broad", 0.16, -0.16)
  bt <- d[d$variant == "broad", ]
  ggplot(d, aes(odds_ratio, y)) +
    geom_vline(xintercept = 1, colour = LINE, linewidth = LW) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, colour = direction),
                   height = 0, linewidth = LW) +
    geom_point(aes(colour = direction, shape = gene_set, alpha = gene_set), size = 1.6) +
    geom_text(aes(x = 40, label = sprintf("k = %d", k)), hjust = 0,
              family = FONT, size = gs(PTS), colour = INK) +
    scale_colour_manual(values = c(enriched = "#12946B",
                                   depleted = "#C2472A"), guide = "none") +
    scale_shape_manual(values = c(broad = 16, narrow = 1), guide = "none") +
    scale_alpha_manual(values = c(broad = 1, narrow = 0.55), guide = "none") +
    scale_y_continuous(breaks = seq_len(nlevels(d$category)),
                       labels = levels(d$category)) +
    scale_x_continuous(transform = "log10",
                       breaks = c(0.1, 0.3, 1, 3, 10, 30)) +
    labs(x = "odds ratio against external gene sets [95% CI]", y = NULL) +
    # The lower end of the TGFb/BMP interval ran into the y-axis labels.
    # clip = "on" now cuts the bars at the edge of the field; the k value
    # still stands in the margin, drawn by its own geom.
    coord_cartesian(xlim = c(0.06, 60),
                    ylim = c(0.5, nlevels(d$category) + 0.5), clip = "on") +
    theme_pub() + theme(plot.margin = margin(2, 12, 1.5, 2))
}

f2e <- function() {
  d <- lies("F2E_atac_per_axis")
  # The axis names are already English in the panel file (presentation
  # layer); only their order is set here.
  d <- d[d$null_model == "H1 baseline-stratified", ]
  d$axis <- factor(d$axis, levels = c("adipogenic", "osteogenic",
                                        "lineage contrast"))
  d$window <- factor(d$window, levels = rev(c("P", "T10", "T50", "GB")))
  d$calibrated <- ifelse(d$calibration_passed %in% TRUE, "calibrated",
                      "NOT calibrated")
  ggplot(d, aes(concordance, window)) +
    geom_segment(aes(x = concordance_detection_limit, xend = concordance, yend = window),
                 colour = GRAU, linewidth = LW) +
    geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.2,
               colour = INK) +
    geom_point(aes(colour = axis, alpha = calibrated), size = 1.8) +
    facet_grid(axis ~ ., switch = "y") +
    scale_colour_manual(values = ACHSE, guide = "none") +
    scale_alpha_manual(values = c(calibrated = 1, `NOT calibrated` = 0.35),
                       guide = "none") +
    # The calibration column stands to the left, BEFORE the data. The view
    # therefore starts at 0.295 rather than 0.40; otherwise the number would
    # sit on the points of the "lineage contrast" row, which reach down to
    # 0.435.
    geom_text(aes(x = 0.300, label = sprintf("cal. z %+.2f", calibration_z)),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK,
              na.rm = TRUE) +
    geom_text(data = data.frame(
                axis = factor("adipogenic",
                               levels = c("adipogenic", "osteogenic",
                                          "lineage contrast")),
                x = 0.300, y = 5.05,
                lab = "faded = axis fails its own calibration"),
              inherit.aes = FALSE, aes(x = x, y = y, label = lab),
              hjust = 0, family = FONT, size = gs(PTS), colour = GRAU) +
    scale_x_continuous(breaks = c(0.4, 0.5, 0.6, 0.7)) +
    labs(x = mit_gse("module concordance, ATAC", "GSE332758"),
         y = "window") +
    coord_cartesian(xlim = c(0.295, 0.73), clip = "off") +
    theme_pub() +
    # The strip titles and the axis title used to sit on top of one another.
    # The strips are now set on one line and moved 3 mm clear.
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, hjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.4, "mm"),
          axis.title.y = element_text(margin = margin(r = 4)),
          plot.margin = margin(2, 4, 1.5, 2))
}

f2f_zerlegung <- function() {
  # The three-way decomposition, across ALL 18 perturbation data sets -- all
  # of them have an undifferentiated arm, so the decomposition is computable
  # in every one. Preregistered in preregistrations/PRAEREG_F2F.md, computed
  # by 04_programme_decomposition/10_decomposition_18_datasets.py. The single case GSE151315 (H3K27ac)
  # is in the supplement as S9C and S9D.
  #
  # Colour codes the CALIBRATION STATUS, not the quantity -- the quantity is
  # on the y axis. That is exactly the statement: the data sets whose
  # calibration fails do not reach their lineage, they leave the
  # undifferentiated state nonetheless, and the module runs in all eighteen.
  #
  # Two fields instead of one, because the module statistic runs to z +13.1
  # while the marker contrasts lie between -4 and +4. A shared view would
  # either compress the markers or clip half of the module -- both would be a
  # distortion. The two x axes are labelled separately.
  d <- lies("F2F_decomposition_eighteen")
  d$calibration <- factor(d$calibration,
                      levels = c("passed", "failed", "not calibratable"))

  mk <- d[d$quantity != "the module", ]
  mk$lab <- factor(mk$quantity,
                   levels = rev(c("undifferentiated markers",
                                  "own-lineage markers",
                                  "adipogenic markers")))
  # Counts per DATA SET (three marker contrasts share the calibration of one
  # data set): 2 passed, 15 failed, 1 not calibratable.
  ds <- mk[!duplicated(mk$dataset), ]
  zebra <- data.frame(y = c(2))
  links <- ggplot(mk, aes(z, lab, colour = calibration)) +
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = SOFT) +
    geom_point(size = 1.4, alpha = 0.85,
               position = position_jitter(width = 0, height = 0.16,
                                          seed = 24)) +
    # The zero line and the threshold lines end above the key; as geom_vline
    # they ran through the two text lines at the foot of the field. The
    # segments stand AFTER geom_point so that the y scale stays discrete.
    annotate("segment", x = 0, xend = 0, y = 0.62, yend = 3.55,
             colour = LINE, linewidth = LW) +
    annotate("segment", x = c(-2, 2), xend = c(-2, 2), y = 0.62, yend = 3.55,
             colour = DIFF, linewidth = LW, linetype = "22") +
    scale_colour_manual(values = EICH, guide = "none") +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
    # The key used to sit on the points of the bottom row. The three data rows
    # lie at y 1 to 3; the key now stands BELOW y = 0.5, in a strip that ylim
    # keeps free for it.
    txt(-4.6, 0.40, "one point per dataset; colour = its own calibration",
        size = PTS, hjust = 0) +
    txt(-4.6, 0.20, sprintf("green passed (%d), red failed (%d), grey not calibratable (%d)",
                            sum(ds$calibration == "passed"),
                            sum(ds$calibration == "failed"),
                            sum(ds$calibration == "not calibratable")),
        size = PTS, hjust = 0) +
    labs(x = "marker-set contrast (z)", y = NULL) +
    coord_cartesian(xlim = c(-4.6, 4.6), ylim = c(0.02, 3.6), clip = "off") +
    theme_pub() + theme(plot.margin = margin(4, 3, 1.5, 2))

  mo <- d[d$quantity == "the module", ]
  mo$lab <- factor("")
  rechts <- ggplot(mo, aes(z, lab, colour = calibration)) +
    geom_point(size = 1.4, alpha = 0.85,
               position = position_jitter(width = 0, height = 0.30,
                                          seed = 25)) +
    scale_colour_manual(values = EICH, guide = "none") +
    scale_x_continuous(breaks = c(6, 9, 12)) +
    # The three lines used to stand far apart to the right of the points.
    # They are now set on two lines above the point cloud, centred.
    txt(9.2, 2.12, "the module: 18 of 18", size = PTS, hjust = 0.5) +
    txt(9.2, 1.88, "above its own limit", size = PTS, hjust = 0.5) +
    labs(x = "module (z)", y = NULL) +
    coord_cartesian(xlim = c(4.8, 13.6), ylim = c(0.30, 2.4), clip = "off") +
    theme_pub(links = FALSE) +
    theme(plot.margin = margin(4, 4, 1.5, 4),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())

  leerstreifen(links) + leerstreifen(rechts) + plot_layout(widths = c(1.9, 1))
}

# --- validation and robustness of the programme (new; 05_programme_validation/10_heldout_and_robustness.py)
HELDOUT <- c("above" = "#12946B", "below" = "#C2472A")

f2g <- function() {
  # Leave-one-study-out: the programme is re-derived from the remaining studies
  # and scored on the held-out dataset(s). One point per held-out dataset;
  # colour = whether it clears its own detection limit (concordance z = 2.8).
  d <- lies("F2G_leave_one_study_out")
  d$status <- ifelse(d$above_mde80 %in% c("True", "TRUE", TRUE), "above", "below")
  n_above <- sum(d$status == "above"); n_tot <- nrow(d)
  ggplot(d, aes(z, 1, colour = status)) +
    annotate("segment", x = 2.8, xend = 2.8, y = 0.55, yend = 1.5,
             colour = DIFF, linewidth = LW, linetype = "22") +
    annotate("segment", x = 0, xend = 0, y = 0.55, yend = 1.5,
             colour = LINE, linewidth = LW) +
    geom_point(size = 1.5, alpha = 0.85,
               position = position_jitter(width = 0, height = 0.28, seed = 35)) +
    scale_colour_manual(values = HELDOUT, guide = "none") +
    txt(2.8, 1.62, "own limit", size = PTS, hjust = 0.5, colour = DIFF) +
    txt(-2.6, 0.44, sprintf("%d of %d held-out datasets above own limit",
                            n_above, n_tot), size = PTS, hjust = 0) +
    labs(x = "held-out concordance (z)", y = NULL) +
    scale_x_continuous(breaks = c(0, 4, 8)) +
    coord_cartesian(xlim = c(-2.8, 10), ylim = c(0.4, 1.75), clip = "off") +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(6, 6, 1.5, 2))
}

f2h <- function() {
  # Robustness: cross-arm rho after dropping the strongest / most-expressed /
  # random gene fractions; the leave-one-gene-out range printed as text.
  d <- lies("F2H_dropout")
  d$pct <- d$removed_frac * 100
  full <- d$rho[d$scheme == "full"]
  jk <- lies("F2H_jackknife")
  nm <- c(drop_top_abs_dwt = "strongest |dWT|",
          drop_top_expression = "most expressed", drop_random = "random")
  d$label <- nm[d$scheme]
  dd <- d[d$scheme != "full", ]
  rnd <- dd[dd$scheme == "drop_random", ]
  COL <- c("strongest |dWT|" = "#C2472A", "most expressed" = "#E07B12",
           "random" = GRAU)
  ggplot(dd, aes(pct, rho, colour = label, group = label)) +
    annotate("segment", x = 0, xend = 20, y = full, yend = full,
             colour = LINE, linewidth = LW, linetype = "22") +
    geom_ribbon(data = rnd, inherit.aes = FALSE, fill = SOFT,
                mapping = aes(x = pct, ymin = rho - rho_sd, ymax = rho + rho_sd)) +
    geom_line(linewidth = LWD) +
    geom_point(size = 1.2) +
    scale_colour_manual(values = COL, name = NULL) +
    txt(0, full + 0.055, sprintf("full = %.3f", full), size = PTS, hjust = 0) +
    txt(0, 0.435, sprintf("leave-one-gene-out:\n%.3f–%.3f",
                          min(jk$rho_without_gene), max(jk$rho_without_gene)),
        size = PTS, hjust = 0, vjust = 1) +
    labs(x = "% of module genes removed", y = "cross-arm rho") +
    coord_cartesian(xlim = c(0, 20), ylim = c(0.25, 0.72), clip = "off") +
    scale_x_continuous(breaks = c(0, 5, 10, 20)) +
    theme_pub() +
    theme(legend.position = c(0.72, 0.20),
          legend.key.height = unit(3.0, "mm"),
          legend.text = element_text(size = PTS),
          legend.background = element_blank(),
          plot.margin = margin(6, 6, 1.5, 2))
}

f2i <- function() {
  # External validation: the locked programme scored on independent GEO
  # differentiation datasets that took no part in its derivation. One row per
  # dataset; colour = whether it clears its own detection limit (z = 2.8).
  d <- lies("F2I_external_validation")
  lab <- c("osteogenic (MSC → osteoblast, GSE37558)" = "osteogenic",
           "adipogenic (hMSC, GSE283759)" = "adipogenic",
           "vascular calcification (VSMC → CVC, GSE37558)" = "vascular calcif.",
           "chondrogenic (iPSC-derived MSC, GSE214987)" = "chondrogenic")
  d$short <- lab[d$dataset]
  d$status <- ifelse(d$above_mde80 %in% c("True", "TRUE", TRUE), "above", "below")
  d <- d[order(d$z), ]
  d$short <- factor(d$short, levels = d$short)
  n_above <- sum(d$status == "above")
  ggplot(d, aes(z, short, colour = status)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = 2.8, colour = DIFF, linewidth = LW, linetype = "22") +
    geom_segment(aes(x = 0, xend = z, yend = short), colour = GRAU, linewidth = LW) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = HELDOUT, guide = "none") +
    txt(2.8, 4.82, "own limit", size = PTS, hjust = 0.5, colour = DIFF) +
    txt(3.4, 1.72, sprintf("%d of %d independent\ndatasets above limit",
                           n_above, nrow(d)), size = PTS, hjust = 0, vjust = 1) +
    labs(x = "concordance of locked programme (z)", y = NULL) +
    scale_x_continuous(breaks = c(0, 4, 8)) +
    coord_cartesian(xlim = c(-0.5, 11.5), ylim = c(0.5, 5.05), clip = "off") +
    theme_pub() + theme(plot.margin = margin(6, 6, 1.5, 2))
}

bau_f2 <- function() {
  # leerstreifen(): F2E is facetted and brings an outer strip column with it
  # on the left. Without compensation the axis of the other panels slides into
  # that column and stands beside its field rather than at it.
  # Rows A-F are the result; G-I are its validation: held-out re-derivation,
  # gene robustness, and external validation on independent data (the
  # matched-null control lives inside panel A).
  p <- (f2a() | f2b()) /
    (leerstreifen(f2c()) | leerstreifen(f2d())) /
    (f2e() | f2f_zerlegung()) /
    (leerstreifen(f2g()) | leerstreifen(f2h()) | leerstreifen(f2i())) +
    plot_layout(heights = c(1.25, 0.95, 1.45, 1.0))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 62), list("D", 89, 62),
             list("E", 1.5, 108), list("F", 89, 108),
             list("G", 1.5, 184), list("H", 60, 184), list("I", 119, 184)),
        "F2", SP2, 230)
}

# =============================================================================
# FIGURE 3 -- in vivo
# =============================================================================
ZONEN <- c("MesCond", "ChondroProg", "RestingChon", "ProlifChon",
           "PrehyperChon", "HyperChon")
ZKURZ <- c("MesCond" = "MesCond", "ChondroProg" = "ChondroProg",
           "RestingChon" = "Resting", "ProlifChon" = "Prolif",
           "PrehyperChon" = "Prehyper", "HyperChon" = "Hyper")

f3a <- function() {
  d <- lies("F3A_atlas_zones")
  # Cell count and sample count in ONE aggregation -- the two earlier
  # aggregate() calls returned different row orders, which assigned the sample
  # counts to the wrong zone.
  a <- do.call(rbind, lapply(split(d, d$zone), function(g)
    data.frame(zone = g$zone[1], zone_rank = g$zone_rank[1],
               n_cells = sum(g$n_cells), np = nrow(g))))
  a <- a[order(a$zone_rank), ]
  a$zone <- factor(ZKURZ[a$zone], levels = ZKURZ[ZONEN])
  # A lollipop rather than bars: on a log axis the area of a bar encodes
  # nothing any more.
  ggplot(a, aes(zone, n_cells)) +
    geom_segment(aes(xend = zone, y = 5, yend = n_cells),
                 colour = "#0E8C8C", linewidth = LWD) +
    geom_point(colour = "#0E8C8C", size = 1.8) +
    geom_text(aes(label = n_cells), vjust = -0.9,
              family = FONT, size = gs(PTS), colour = INK) +
    scale_y_continuous(transform = "log10",
                       breaks = c(10, 100, 1000, 5000)) +
    labs(x = NULL, y = "cells in the atlas (log scale)") +
    coord_cartesian(ylim = c(5, 12000), clip = "off") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 4, 1.5, 2))
}

f3b <- function() {
  d <- lies("F3B_positive_control_per_sample")
  d$zone <- factor(ZKURZ[d$zone], levels = ZKURZ[ZONEN])
  ggplot(d, aes(zone, contrast)) +
    geom_hline(yintercept = 0, colour = LINE, linewidth = LW) +
    geom_point(aes(y = detection_limit), shape = 95, size = 3.2, colour = GRAU) +
    geom_point(colour = "#0E8C8C", size = 1.3,
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    # The note used to sit at the height of the point cloud. It now stands at
    # the top edge of the field, above every point.
    txt(0.6, Inf, "dash = own detection limit", size = PTS, hjust = 0,
        vjust = 1.6, colour = GRAU) +
    labs(x = NULL, y = "chondrogenic minus undifferentiated") +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 4, 1.5, 2))
}

f3c <- function() {
  d <- lies("F3C_module_per_sample")
  t <- lies("F3C_trend_test_donor")
  tm <- t[t$quantity == "module (173 genes)", ]
  d$zone <- factor(ZKURZ[d$zone], levels = ZKURZ[ZONEN])
  ggplot(d, aes(zone, contrast, group = specimen)) +
    geom_hline(yintercept = 0, colour = LINE, linewidth = LW) +
    geom_line(colour = GRAU, linewidth = LW) +
    geom_point(colour = "#0E8C8C", size = 1.2) +
    txt(0.6, 0.235,
        sprintf("rho = %.3f, z = %+.2f, limit rho %.3f\n%d samples, %d specimens",
                tm$rho, tm$z, tm$detection_limit_rho, tm$n_samples, tm$n_specimens),
        size = PTS, hjust = 0, vjust = 1) +
    labs(x = NULL, y = "module contrast (up minus down)") +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 4, 1.5, 8))
}

f3d <- function() {
  d <- lies("F3D_cells_per_zone")
  d$zone <- factor(ZKURZ[d$zone], levels = rev(ZKURZ[ZONEN]))
  d$knapp <- d$n_samples <= 1
  ggplot(d, aes(n_samples, zone)) +
    geom_col(aes(fill = knapp), width = 0.6, colour = NA) +
    geom_text(aes(label = ifelse(n_samples == 1, "1 sample",
                                 sprintf("%d samples", n_samples))),
              hjust = -0.15, family = FONT, size = gs(PTS), colour = INK) +
    scale_fill_manual(values = c("FALSE" = "#0E8C8C", "TRUE" = DIFF),
                      guide = "none") +
    scale_x_continuous(breaks = c(0, 5, 10, 15),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "evaluable samples (>= 5 cells)", y = NULL) +
    # The label on the longest bar (18) needs room to the right -- the view
    # now ends AT the label, not at the bar.
    coord_cartesian(xlim = c(0, 21.5), clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 4, 1.5, 2))
}

bau_f3 <- function() {
  p <- (f3a() | f3b()) / (f3c() | f3d()) + plot_layout(heights = c(1, 1))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 62), list("D", 89, 62)),
        "F3", SP2, 122)
}

# =============================================================================
# FIGURE 4 -- where the disease genes are
# =============================================================================
f4a <- function() {
  lin <- lies("F4A_positive_control_lineage_markers")
  ank <- lies("F4A_positive_control_anchor")
  # The panel files already carry the full panel names (the presentation layer
  # of 10_panel_data_main.py renames them); use them directly.
  d <- rbind(
    data.frame(panel = lin$panel, odds_ratio = lin$odds_ratio_matched,
               detection_limit = lin$odds_ratio_detection_limit, part = "lineage\nmarkers"),
    data.frame(panel = ank$panel, odds_ratio = ank$odds_ratio_raw, detection_limit = NA,
               part = "secretion\nanchor"))
  d$panel <- factor(d$panel, levels = c("Nosology (core)", "Nosology (broad)",
                                        "PanelApp 309"))
  ggplot(d, aes(odds_ratio, panel)) +
    geom_vline(xintercept = 1, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 1, xend = odds_ratio, yend = panel), colour = "#12946B",
                 linewidth = LWD) +
    geom_point(aes(x = detection_limit), shape = 124, size = 2.2, colour = INK,
               na.rm = TRUE) +
    geom_point(colour = "#12946B", size = 1.7) +
    geom_text(aes(x = odds_ratio * 1.25, label = sprintf("%.1f", odds_ratio)), hjust = 0,
              family = FONT, size = gs(PTS), colour = INK) +
    facet_grid(part ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(transform = "log10", breaks = c(1, 3, 10, 30)) +
    labs(x = "odds ratio (positive controls)", y = NULL) +
    coord_cartesian(xlim = c(0.9, 130), clip = "off") +
    theme_pub() +
    # The strip titles used to sit flush left at the top (hjust = 0 from
    # theme_pub) instead of centred beside their rows.
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, hjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(2, 8, 1.5, 2))
}

f4b <- function() {
  d <- lies("F4B_complementarity")
  d$gene_set <- factor(d$gene_set, levels = c("programme (173 genes)",
                                      "Nosology (core)", "Nosology (broad)",
                                      "PanelApp 309"))
  d$compartment <- factor(c("distal secretion" = "distal",
                             "biosynthetic secretion" = "biosynthetic")[
                               d$compartment],
                           levels = c("distal", "biosynthetic"))
  d$direction <- ifelse(d$odds_ratio_matched > 1, "enriched", "depleted")
  ggplot(d, aes(odds_ratio_matched, gene_set)) +
    geom_vline(xintercept = 1, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 1, xend = odds_ratio_matched, yend = gene_set,
                     colour = direction), linewidth = LWD) +
    geom_point(aes(x = odds_ratio_detection_limit), shape = 124, size = 2.2, colour = INK) +
    geom_point(aes(colour = direction), size = 1.7) +
    geom_text(aes(x = 3.9, label = sprintf("z %+.2f", z)), hjust = 1,
              family = FONT, size = gs(PTS), colour = INK) +
    facet_grid(compartment ~ ., switch = "y") +
    scale_colour_manual(values = c(enriched = "#12946B",
                                   depleted = "#C2472A"), guide = "none") +
    labs(x = "odds ratio, matched", y = NULL) +
    coord_cartesian(xlim = c(0.3, 4.0), clip = "off") +
    theme_pub() +
    # "distal" and "biosynthetic" centred on their rows.
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, hjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(2, 4, 1.5, 2))
}

f4c_plot <- function() {
  d <- lies("F4C_mode_of_inheritance")
  # The `expr_rank_med` row is NOT shown: the expression difference between
  # monoallelic and biallelic genes is confounded by study intensity and was
  # explicitly discarded.
  d <- d[d$variable %in% c("loeuf", "dWT_abs"), ]
  vlab <- c(loeuf = "LOEUF\n(gene dosage)",
            dWT_abs = "|dWT|\n(differentiation dynamics)")
  g <- rbind(
    data.frame(variable = d$variable, value = d$median_monoallelic, p_value = d$p_value,
               n = d$n_monoallelic, modus = "monoallelic"),
    data.frame(variable = d$variable, value = d$median_biallelic, p_value = d$p_value,
               n = d$n_biallelic, modus = "biallelic"))
  g$feld <- factor(vlab[g$variable], levels = vlab[c("loeuf", "dWT_abs")])
  g$y <- 1
  lab <- data.frame(feld = factor(vlab[d$variable],
                                  levels = vlab[c("loeuf", "dWT_abs")]),
                    x = 0.5, y = 1.34,
                    lab = ifelse(d$p_value < 0.001,
                                 sprintf("P = %.0e", d$p_value),
                                 sprintf("P = %.2f", d$p_value)))
  # A free x axis per field: LOEUF and |dWT| are different quantities and do
  # not belong on a shared scale.
  ggplot(g, aes(value, y, colour = modus)) +
    geom_line(aes(group = feld), colour = GRAU, linewidth = LW) +
    geom_point(size = 2.0) +
    geom_text(aes(label = sprintf("%s\nn = %d", modus, n),
                  vjust = ifelse(modus == "monoallelic", 1.8, -0.9)),
              family = FONT, size = gs(PTS), lineheight = 0.95,
              show.legend = FALSE) +
    geom_text(data = lab, inherit.aes = FALSE,
              aes(x = x, y = y, label = lab), family = FONT,
              size = gs(PTS), colour = INK, hjust = 0.5) +
    facet_wrap(~ feld, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = c(monoallelic = "#1F6FB2",
                                   biallelic = "#E07B12"), guide = "none") +
    labs(x = "median (PanelApp 309, split by mode of inheritance)",
         y = NULL) +
    coord_cartesian(ylim = c(0.55, 1.45), clip = "off") +
    theme_pub(links = FALSE) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(2, 6, 1.5, 2))
}

achsenpanel <- function(name, xlab, xlim, unit) {
  d <- lies(name)
  lab <- c("PanelApp 309" = "PanelApp 309", "Nosology (core)" = "Nosology (core)",
           "Nosology (broad)" = "Nosology (broad)", "short stature" = "short stature",
           "height GWAS" = "height GWAS",
           "cell cycle (neg. control)" = "cell cycle (neg. control)",
           "programme (173 genes)" = "programme (pos. control)")
  d$satzlab <- lab[d$gene_set]
  d <- d[order(d$z), ]
  d$satzlab <- factor(d$satzlab, levels = d$satzlab)
  d$role <- ifelse(d$gene_set == "programme (173 genes)", "positive control",
                    ifelse(d$gene_set == "cell cycle (neg. control)",
                           "negative control", "disease / trait panel"))
  zebra <- data.frame(y = seq(1, nrow(d), by = 2))
  # The limit column stands in a fixed column on the RIGHT, not at the data
  # point -- at negative z it would otherwise collide with the row name.
  xtext <- xlim[2] + diff(xlim) * 0.03
  ggplot(d, aes(z, satzlab)) +
    # xmin = -Inf points at the minimum of the SCALE, not at the
    # coord_cartesian view -- hence the gap between the axis line and the
    # zebra stripes that looked like a detached y axis.
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = xlim[1], xmax = xlim[2], ymin = y - 0.5,
                  ymax = y + 0.5), fill = SOFT) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 0, xend = z, yend = satzlab, colour = role),
                 linewidth = LWD) +
    geom_point(aes(colour = role), size = 1.7) +
    geom_text(aes(x = xtext, label = sprintf("%.3f", detection_limit_delta)),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK) +
    annotate("text", x = xtext, y = nrow(d) + 0.85,
             label = paste0("limit\n(", unit, ")"), hjust = 0, vjust = 0.5,
             family = FONT, size = gs(PTS), colour = INK, lineheight = 0.95) +
    scale_colour_manual(values = c("positive control" = "#12946B",
                                   "negative control" = "#7B5EA7",
                                   "disease / trait panel" = "#1F6FB2"),
                        guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0)) +
    labs(x = xlab, y = NULL) +
    coord_cartesian(xlim = xlim, ylim = c(0.4, nrow(d) + 1.2),
                    clip = "off") +
    theme_pub() + theme(plot.margin = margin(6, 22, 1.5, 2))
}

f4e <- function() {
  d <- lies("F4F_mechanism_classes")
  d <- d[d$computable, ]
  klab <- c("cilium" = "cilium", "ECM structure" = "ECM structure",
            "glycosylation / linker" = "glycosylation / linker",
            "lysosome" = "lysosome",
            "signalling FGF/WNT/BMP/IHH" = "signalling FGF/WNT/BMP/IHH",
            "transcription factors" = "transcription factors",
            "vesicle / secretion" = "vesicle / secretion")
  d$class <- factor(klab[d$class], levels = rev(sort(unname(klab))))
  ggplot(d, aes(odds_ratio_matched, class)) +
    geom_vline(xintercept = 1, colour = LINE, linewidth = LW) +
    geom_point(aes(x = odds_ratio_detection_limit), shape = 124, size = 2.0, colour = GRAU) +
    geom_point(aes(colour = above_detection_limit), size = 1.3,
               position = position_jitter(width = 0, height = 0.14, seed = 2)) +
    scale_colour_manual(values = c("FALSE" = "#1F6FB2", "TRUE" = DIFF),
                        guide = "none") +
    scale_x_continuous(transform = "log10", breaks = c(0.3, 1, 3, 10)) +
    # The key sits ABOVE the top row (7 classes, so row 7); it used to stand
    # at 6.95 to 8.45 and collided with the cilium points. The block now sits
    # entirely above y = 7.5.
    txt(0.17, 10.55, sprintf("%d of 49 tests computable", nrow(d)),
        size = PTS, hjust = 0) +
    txt(0.17, 10.05, "none above Bonferroni", size = PTS, hjust = 0) +
    txt(0.17, 9.55, "grey bar = own detection limit", size = PTS, hjust = 0,
        colour = GRAU) +
    txt(0.17, 9.05, "red = above own limit (4, all glyco / linker)",
        size = PTS, hjust = 0, colour = DIFF) +
    labs(x = "odds ratio, matched (7 panels x 7 classes)", y = NULL) +
    coord_cartesian(xlim = c(0.15, 14), ylim = c(0.4, 10.9), clip = "off") +
    theme_pub() + theme(plot.margin = margin(11, 4, 1.5, 2))
}

f4g_equivalence <- function() {
  # Effect-size exclusion for the disease-gene negative. The dynamics-axis
  # contrast (|dWT|, publication-matched) in its native units: observed effect
  # with 95% CI against the smallest effect the test would detect at 80% power
  # (MDE80). Where the whole 95% CI lies below the MDE80, an effect of biological
  # size is excluded, not merely non-significant. Reuses F4E_dynamics_axis.csv.
  d <- lies("F4E_dynamics_axis")
  keep <- c("PanelApp 309", "Nosology (core)", "short stature", "height GWAS",
            "cell cycle (neg. control)")
  d <- d[d$gene_set %in% keep, ]
  nm <- c("PanelApp 309" = "dysplasia (PanelApp 309)",
          "Nosology (core)" = "dysplasia (Nosology)",
          "short stature" = "short stature",
          "height GWAS" = "height GWAS",
          "cell cycle (neg. control)" = "cell cycle (control)")
  d$label <- factor(nm[d$gene_set], levels = rev(unname(nm)))
  d$ci_lo <- d$delta_observed - 1.96 * d$null_sd
  d$ci_hi <- d$delta_observed + 1.96 * d$null_sd
  d$excluded <- ifelse(d$above_detection_limit %in% c("True", "TRUE", TRUE),
                       "detected", "effect excluded")
  ggplot(d, aes(delta_observed, label)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_point(aes(x = detection_limit_delta), shape = 124, size = 3.0,
               colour = GRAU) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = excluded),
                   height = 0.16, linewidth = LWD) +
    geom_point(aes(colour = excluded), size = 1.7) +
    scale_colour_manual(values = c("effect excluded" = "#1F6FB2",
                                   "detected" = DIFF), guide = "none") +
    txt(0.175, 5.35, "grey bar = MDE80", size = PTS, hjust = 1, colour = GRAU) +
    txt(0.175, 0.60, "programme, same run: +0.73 (z +18.1)",
        size = PTS, hjust = 1, colour = KURVE[["programme"]]) +
    labs(x = "dynamics-axis effect, publication-matched (|dWT| units)", y = NULL) +
    coord_cartesian(xlim = c(-0.095, 0.175), ylim = c(0.5, 5.6), clip = "off") +
    scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15)) +
    theme_pub() + theme(plot.margin = margin(4, 4, 1.5, 2))
}

bau_f4 <- function() {
  # leerstreifen(): F4A and F4B are facetted and carry an outer strip column
  # on the left. Without compensation patchwork pushes the y axis of the other
  # panels into that strip column -- that was the detached y axis in F4D.
  # Row 4 (G) is the equivalence panel: the disease-gene negative shown as an
  # explicit effect-size exclusion, not merely a non-significant result.
  p <- (f4a() | f4b()) /
    (f4c_plot() | leerstreifen(achsenpanel(
        "F4D_constraint_publication_matched",
        "LOEUF, publication-matched (z)", c(-6.2, 3.2), "LOEUF"))) /
    (leerstreifen(achsenpanel("F4E_dynamics_axis",
                              "|dWT|, publication-matched (z)", c(-2.5, 20),
                              "|dWT|")) | leerstreifen(f4e())) /
    ((leerstreifen(f4g_equivalence()) | patchwork::plot_spacer()) +
       plot_layout(widths = c(2.5, 1))) +
    plot_layout(heights = c(1.05, 1.05, 1.15, 0.62))
  # Row boundaries for heights c(1.05,1.05,1.15,0.62) on a 214 mm sheet:
  # row1 0-58, row2 58-116, row3 116-180, row4 180-214. Letters sit just below
  # each row's top so they never land on the previous row's axis title.
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 59.5), list("D", 89, 59.5),
             list("E", 1.5, 117.5), list("F", 89, 117.5),
             list("G", 1.5, 181)),
        "F4", SP2, 214)
}

# =============================================================================
# FIGURE 5 -- both layers meet at the prehypertrophic transition
# =============================================================================
f5a <- function() {
  d <- lies("F5A_panel_vs_module_per_zone")
  r <- lies("F5A_panel_vs_module_per_sample")
  d$zone <- factor(ZKURZ[d$zone], levels = ZKURZ[ZONEN])
  r$zone <- factor(ZKURZ[r$zone], levels = ZKURZ[ZONEN])
  ggplot(d, aes(zone, contrast_median)) +
    geom_point(data = r, aes(zone, contrast), colour = GRAU, size = 0.9,
               position = position_jitter(width = 0.13, height = 0, seed = 3)) +
    geom_linerange(aes(ymin = contrast_ci_low, ymax = contrast_ci_high),
                   colour = "#1F6FB2", linewidth = LWD) +
    geom_point(colour = "#1F6FB2", size = 1.8) +
    geom_text(aes(y = contrast_ci_high + 0.012,
                  label = sprintf("z %+.2f", z_median)),
              family = FONT, size = gs(PTS), colour = INK, vjust = 0) +
    labs(x = NULL, y = "disease genes minus programme") +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 4, 1.5, 2))
}

f5b <- function() {
  d <- lies("F5B_both_curves")
  d$curve <- ifelse(grepl("^programme", d$curve), "programme",
                    "disease genes")
  d$zone <- factor(ZKURZ[d$zone], levels = ZKURZ[ZONEN])
  ggplot(d, aes(zone, median, colour = curve, group = curve)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = curve), alpha = 0.16,
                colour = NA) +
    geom_line(linewidth = LWD) +
    geom_point(size = 1.5) +
    annotate("segment", x = 5, xend = 5, y = 0.126, yend = 0.168,
             colour = INK, linewidth = LW,
             arrow = arrow(length = unit(1.1, "mm"), ends = "both",
                           type = "closed")) +
    # The arrow marks the two peaks at the prehypertrophic transition; the
    # label stands above the arrow head, in the free field above both bands.
    # It read "widest gap" until the zone medians showed that the vertical
    # distance between the curves is widest at Prolif (0.073 against 0.042
    # here), so that claim was dropped; the defensible statement is the
    # shared peak, as in panel C.
    txt(5.05, 0.187, "both peak here", size = PTS, hjust = 0) +
    scale_colour_manual(values = KURVE, guide = "none") +
    scale_fill_manual(values = KURVE, guide = "none") +
    txt(0.6, 0.20, "disease genes (PanelApp 309)", colour = KURVE[["disease genes"]],
        size = PTS, hjust = 0) +
    txt(0.6, 0.185, "programme (173 genes)", colour = KURVE[["programme"]],
        size = PTS, hjust = 0) +
    labs(x = NULL, y = "contrast against background (median, IQR)") +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 4, 1.5, 2))
}

f5c <- function() {
  k <- data.frame(
    x = c(1, 3, 2), y = c(2.1, 2.1, 0.85),
    w = c(1.7, 1.7, 2.6), h = c(0.66, 0.66, 0.66),
    lab = c("matrix programme\nwhat the cell does",
            "disease genes\nthe machinery that executes it",
            "prehypertrophic transition\nboth peak here"),
    key = c("programme", "disease genes", "meeting"))
  ggplot() +
    geom_rect(data = k, aes(xmin = x - w / 2, xmax = x + w / 2,
                            ymin = y - h / 2, ymax = y + h / 2, colour = key),
              fill = "white", linewidth = LW * 2) +
    geom_text(data = k, aes(x, y, label = lab, colour = key), family = FONT,
              size = gs(PTS), lineheight = 0.95) +
    annotate("segment", x = c(1, 3), xend = c(1.6, 2.4),
             y = c(1.77, 1.77), yend = c(1.18, 1.18), colour = INK,
             linewidth = LW, arrow = arrow(length = unit(1.2, "mm"),
                                           type = "closed")) +
    scale_colour_manual(values = c(KURVE, meeting = "#12946B"),
                        guide = "none") +
    txt(2, 0.28,
        "dysplasia is a failure of the machinery, not of the programme",
        size = PTS, hjust = 0.5) +
    txt(2, 0.06,
        "which is why there is no transcriptional downstream convergence to find",
        size = PTS, hjust = 0.5) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 2.6), clip = "off") +
    theme_leer()
}

bau_f5 <- function() {
  p <- (f5a() | f5b()) / f5c() + plot_layout(heights = c(1.5, 1))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5), list("C", 1.5, 74)),
        "F5", SP2, 120)
}

# =============================================================================
# FIGURE 6 -- levels and detection limits as a forest and range plot
# =============================================================================
# Per level the figure draws the estimate range as a point or a span, and the
# level's OWN detection limit as a cross bar (shape 124) in the same unit --
# the detection-limit idea is the visual signature of the paper. Rows whose
# limit is expressed in a different unit (|dWT|, LOEUF), or that have no
# measurable limit, carry no bar; the limit then stands in the legend and in
# Table S14. Three facets by kind of quantity, so that nothing is forced onto
# a shared axis that does not exist: unit-free statistics (rho, C, K, S1,
# share; on [0,1]), odds ratio (log), and z contrasts. The text version is
# Table S14.
URTEIL <- c("carries" = "#12946B", "good negative" = "#1F6FB2",
            "not measurable" = DIFF, "carries nothing" = GRAU,
            "observation only" = "#E07B12", "fails" = GRAU)

# One field of the levels book: rows as points or spans, the level's OWN
# detection limit as a grey cross bar in the same unit. The detection-limit
# idea is the visual signature of the paper and takes the same form in all
# three fields.
f6_feld <- function(d, xlab, xlim, log = FALSE, null_model = NA) {
  d$level <- factor(d$level, levels = rev(d$level))
  n <- nrow(d)
  zebra <- data.frame(y = seq(1, n, by = 2))
  p <- ggplot(d, aes(y = level)) +
    # Zebra stripes across the WHOLE view (not -Inf: that points at the
    # minimum of the scale and used to leave a gap at the axis)
    geom_rect(data = zebra, inherit.aes = FALSE,
              aes(xmin = xlim[1], xmax = xlim[2], ymin = y - 0.5,
                  ymax = y + 0.5), fill = SOFT)
  if (!is.na(null_model))
    p <- p + geom_vline(xintercept = null_model, colour = LINE, linewidth = LW)
  p +
    # the level's own detection limit: a cross bar where it is a point, a
    # span where it varies across windows or tests
    geom_errorbarh(data = subset(d, !is.na(detection_limit_min)),
                   aes(xmin = detection_limit_min, xmax = detection_limit_max),
                   height = 0.30, colour = GRAU, linewidth = LW * 1.8) +
    # null mean +/- SD where a null distribution exists (donor level)
    geom_errorbarh(data = subset(d, !is.na(null_min)),
                   aes(xmin = null_min, xmax = null_max),
                   height = 0.14, colour = GRAU, linewidth = LW) +
    # the estimate: a point for a single value, a span across windows or tests
    geom_segment(data = subset(d, estimate_max > estimate_min),
                 aes(x = estimate_min, xend = estimate_max, yend = level,
                     colour = verdict_class), linewidth = LWD) +
    geom_point(aes(x = estimate_min, colour = verdict_class), size = 1.8) +
    geom_point(data = subset(d, estimate_max > estimate_min),
               aes(x = estimate_max, colour = verdict_class), size = 1.4) +
    scale_colour_manual(values = URTEIL, guide = "none") +
    (if (log) scale_x_continuous(transform = "log10",
                                 breaks = c(0.6, 1, 1.6, 2.5))
     else scale_x_continuous()) +
    labs(x = xlab, y = NULL) +
    coord_cartesian(xlim = xlim, ylim = c(0.4, n + 0.6), clip = "off") +
    theme_pub() + theme(plot.margin = margin(2, 4, 1.5, 2))
}

bau_f6 <- function() {
  d <- lies("F6_levels_forest")

  # A  unit-free statistics (rho, C, K, S1, share; all on [0,1])
  pa <- f6_feld(d[d$unit == "unit-free statistic", ],
                "estimate and own detection limit (unit-free statistic)",
                c(0.20, 0.78))
  # B  odds ratio on a log axis. A single row -- it stands at the y axis like
  #    all the others, not as free text in the field.
  pb <- f6_feld(d[d$unit == "odds ratio", ],
                "odds ratio (log scale)", c(0.55, 2.6), log = TRUE, null_model = 1)
  # C  the z contrasts, with their own view (the module runs to z +13.1)
  pc <- f6_feld(d[d$unit == "z (contrasts)", ],
                "contrast against own null (z)", c(-3.0, 14.5), null_model = 0)

  nu <- sum(d$unit == "unit-free statistic")
  nz <- sum(d$unit == "z (contrasts)")
  # Stacked rather than side by side: one column of labels, one reading
  # direction. patchwork aligns the three y axes flush.
  p <- pa / pb / pc +
    plot_layout(heights = c(nu + 1.4, 1 + 1.4, nz + 1.4))

  # Panel letters at the height of the top edge of each field. The total
  # height is 150 mm; the proportions follow the heights set above.
  h <- 150
  a1 <- (nu + 1.4); a2 <- (1 + 1.4); a3 <- (nz + 1.4)
  su <- a1 + a2 + a3
  tafel(p, m(list("A", 1.5, 1.5),
             list("B", 1.5, 1.5 + h * a1 / su),
             list("C", 1.5, 1.5 + h * (a1 + a2) / su)),
        "F6", SP2, h)
}

# =============================================================================
main <- function() {
  cat("20_figures_main.R -- main figures to", PUB_DIR, "\n")
  bau_f1(); bau_f2(); bau_f3(); bau_f4(); bau_f5(); bau_f6()
  si <- file.path(WURZEL, "results", "figures_main_sessionInfo.txt")
  capture.output(sessionInfo(), file = si)
  cat("sessionInfo ->", si, "\n")
}

if (!isTRUE(getOption("abbildungen.nur_funktionen"))) main()
