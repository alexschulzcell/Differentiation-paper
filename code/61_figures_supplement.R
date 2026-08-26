# =============================================================================
# 61_figures_supplement.R -- the nine supplementary figures
# =============================================================================
# Purpose  Draws S1 to S9 from the panel CSV files in figures/data/. The style
#          comes unchanged from figure_style/publication_style.R; PUB_DIR is
#          overwritten after the source() call.
#
# Inputs   figures/data/*.csv  (from 50_panel_data.py and 51_supplement_data.py)
# Outputs  figures/S1.pdf to S9.pdf and the same as PNG at 600 dpi
#          results/abbildungen_supplement_sessionInfo.txt
# Runtime  about a minute
#
# S1  the scale critique, in full
# S2  every calibration and every limit
# S3  external triangulation
# S4  patient against control
# S5  robustness and leave-one-out
# S6  the orthogonal levels -- ATAC in full, H3K27ac, and the methylomes as
#     null results WITH their limits
# S7  the screens in detail
# S8  day zero falls, and the publication matching as a methods exhibit
# S9  the lineage contrast across three chromatin cohorts
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
PUB_DIR <- file.path(WURZEL, "figures")
DAT     <- file.path(WURZEL, "figures", "data")

lies <- function(n) {
  d <- read.csv(file.path(DAT, paste0(n, ".csv")), check.names = FALSE)
  for (k in names(d)) if (is.character(d[[k]])) {
    u <- unique(d[[k]][!is.na(d[[k]])])
    if (length(u) > 0 && all(u %in% c("True", "False"))) d[[k]] <- d[[k]] == "True"
  }
  d
}

EICH <- c("passed" = "#12946B", "failed" = "#C2472A",
          "not calibratable" = GRAU)

theme_leer <- function() {
  theme_pub(links = FALSE, unten = FALSE) +
    theme(axis.text = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), plot.margin = margin(2, 2, 2, 2))
}

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
      grid::grid.text(marken$lab[i], x = grid::unit(marken$x[i], "mm"),
                      y = grid::unit(1, "npc") - grid::unit(marken$y[i], "mm"),
                      just = c("left", "top"),
                      gp = grid::gpar(fontfamily = FONT, fontface = "bold",
                                      fontsize = PTL, col = INK))
    invisible(grDevices::dev.off())
  }
  cat(sprintf("  %-4s %3.0f x %3.0f mm\n", name, breite, hoehe))
}

m <- function(...) do.call(rbind, lapply(list(...), function(v)
  data.frame(lab = v[[1]], x = as.numeric(v[[2]]), y = as.numeric(v[[3]]))))

# The global rule: never place text over a plot element. The reliable place
# for a note is the top margin ABOVE the field -- no data element lies there.
# notiz() puts it there (clip = "off" plus a margin of at least 6 mm).
# zeile = 1 is the top line, 2 the one below it.
notiz <- function(text, colour = INK, zeile = 1)
  annotate("text", x = -Inf, y = Inf, label = text, colour = colour,
           family = FONT, size = gs(PTS), hjust = -0.03,
           vjust = -0.55 + (zeile - 1) * 1.25)

# The same for facetted fields: annotate() draws into EVERY strip, which
# repeats the note three times. notiz_oben() binds it to exactly one strip
# value -- the topmost one.
notiz_oben <- function(text, spalte, value, colour = INK, zeile = 1) {
  d <- data.frame(x = -Inf, y = Inf, lab = text)
  d[[spalte]] <- value
  geom_text(data = d, inherit.aes = FALSE,
            aes(x = x, y = y, label = lab), colour = colour, family = FONT,
            size = gs(PTS), hjust = -0.03,
            vjust = -0.55 + (zeile - 1) * 1.25)
}

# =============================================================================
# S1 -- scale critique
# =============================================================================
bau_s1 <- function() {
  # A  why nothing is adjusted on the z scale: the baseline correlates with
  #    the target quantity.
  a <- lies("S1A_baseline_expression")
  pa <- ggplot(a, aes(median_baseline_shift, median_z, colour = artefact)) +
    geom_hline(yintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_colour_manual(values = c("FALSE" = "#1F6FB2", "TRUE" = DIFF),
                        guide = "none") +
    # The text used to stand at the top left, in the middle of the point
    # cloud. It now stands at the top right, where the field is empty.
    txt(Inf, Inf, sprintf("Pearson r = %.3f, %d GO sets",
                          cor(a$median_baseline_shift, a$median_z), nrow(a)),
        size = PTS, hjust = 1.03, vjust = 1.6) +
    labs(x = "median shift of the baseline", y = "median z of the set") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(6, 4, 1.5, 2))

  # B  what additive offset the statistic actually finds
  b <- lies("S1B_sensitivity")
  bs <- aggregate(z_corrected ~ delta + set_size_per_side, b, median)
  bs$set_size_per_side <- factor(bs$set_size_per_side,
                             levels = sort(unique(bs$set_size_per_side)))
  pb <- ggplot(bs, aes(delta, z_corrected, colour = set_size_per_side,
                       group = set_size_per_side)) +
    geom_hline(yintercept = 2, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_line(linewidth = LWD) + geom_point(size = 1.1) +
    # Six set sizes -- the colour ramp is generated from two end colours
    # (rampe() from publication_style.R), so that the scale grows with them.
    scale_colour_manual(
      values = setNames(rampe(seq(0, 1, length.out = nlevels(bs$set_size_per_side)),
                              "#C9CFE9", "#12406B"),
                        levels(bs$set_size_per_side)), guide = "none") +
    # "z = 2" used to sit on the dashed line and on the ends of the curves and
    # was unreadable. It is now a note in the top margin.
    notiz("dashed line: z = 2", colour = DIFF) +
    labs(x = "additive offset applied (z units)",
         y = "recovered statistic (z)") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(7, 4, 1.5, 2))

  # C  self-test: the null rate of the statistic
  cc <- lies("S1C_neutral_contrast")
  rate <- mean(abs(cc$z_corrected) > 2)
  pc <- ggplot(cc, aes(z_corrected)) +
    geom_histogram(bins = 40, fill = SOFT, colour = GRAU, linewidth = LW) +
    geom_vline(xintercept = c(-2, 2), colour = DIFF, linewidth = LW,
               linetype = "22") +
    notiz(sprintf("|z| > 2 in %.1f %% of %d neutral contrasts (nominal 5 %%)",
                  100 * rate, nrow(cc))) +
    labs(x = "statistic under a neutral contrast (z)", y = "count") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(7, 4, 1.5, 2))

  # D  the step function: convergence count against the noise floor
  d <- lies("S1D_convergence_curve")
  d$quantity <- factor(d$quantity)
  pd <- ggplot(d, aes(threshold, observed, colour = quantity,
                      group = quantity)) +
    geom_ribbon(aes(ymin = noise_q05, ymax = noise_q95, fill = quantity),
                alpha = 0.16, colour = NA) +
    geom_line(aes(y = noise_mean), linetype = "22", linewidth = LW) +
    geom_line(linewidth = LWD) + geom_point(size = 1.2) +
    scale_colour_manual(values = c(dWT = "#0E8C8C", iv = "#E07B12"),
                        guide = "none") +
    scale_fill_manual(values = c(dWT = "#0E8C8C", iv = "#E07B12"),
                      guide = "none") +
    scale_y_continuous(transform = "log10") +
    # The key sits at the top right; on the left it covered the first point of
    # the teal curve at threshold 0.80.
    notiz("solid = observed, dashed = noise floor") +
    labs(x = "concordance threshold", y = "convergent genes (log scale)") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(7, 4, 1.5, 2))

  p <- (pa | pb) / (pc | pd)
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 62), list("D", 89, 62)), "S1", SP2, 122)
}

# =============================================================================
# S2 -- every calibration and every limit
# =============================================================================
bau_s2 <- function() {
  d <- lies("S2A_all_calibrations")
  d$calibration <- ifelse(d$status != "ok", "not calibratable",
                      ifelse(d$passed, "passed", "failed"))
  d <- d[order(d$level, d$z), ]
  # The row names used to be drawn as geom_text at x = -1.9 INSIDE the field
  # and therefore lay on the data lines. They now stand on the y axis, where
  # no plot element can cut them. The names already carry their accession from
  # the presentation layer (the accession rule).
  d$reihe <- factor(seq_len(nrow(d)), levels = rev(seq_len(nrow(d))),
                    labels = rev(d$unit))
  ggp <- ggplot(d, aes(z, reihe)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = 2, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = ifelse(is.na(z), 0, z), yend = reihe,
                     colour = calibration), linewidth = LWD) +
    geom_point(aes(colour = calibration), size = 1.4) +
    facet_grid(level ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    scale_colour_manual(values = EICH, guide = "none") +
    # The threshold note now stands below the x axis, not in the field.
    labs(x = "calibration statistic (z); dashed line: threshold z = 2",
         y = NULL) +
    coord_cartesian(xlim = c(-1.6, 6.2), clip = "off") +
    theme_pub() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.4, "mm"),
          plot.margin = margin(4, 4, 1.5, 2))
  tafel(ggp, m(list("", 1.5, 1.5)), "S2", SP2, 168)
}

# =============================================================================
# S3 -- external triangulation
# =============================================================================
bau_s3 <- function() {
  a <- lies("S3A_external_triangulation")
  a$version <- factor(a$version, levels = rev(a$version))
  pa <- ggplot(a, aes(z, version)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 0, xend = z, yend = version, colour = has_detection_limit),
                 linewidth = LWD) +
    geom_point(aes(colour = has_detection_limit), size = 1.8) +
    geom_text(aes(x = z + 0.15,
                  label = sprintf("z %+.2f, P %.4f, n %d", z, p_value, n)),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK) +
    scale_colour_manual(values = c("TRUE" = "#12946B", "FALSE" = GRAU),
                        guide = "none") +
    # The note stands in the top margin (it used to cross the grey line), and
    # the view now ends behind the longest text, so that "n 8" is no longer
    # cut off at the right edge.
    notiz("grey = no detection limit; not reported as the headline",
          colour = GRAU) +
    labs(x = "synthesis statistic (z)", y = NULL) +
    coord_cartesian(xlim = c(0, 5.6), clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))

  b <- lies("S3B_triangulation_per_study")
  b <- b[order(b$z), ]
  b$dataset <- factor(b$dataset, levels = b$dataset)
  pb <- ggplot(b, aes(directional_share, dataset)) +
    geom_vline(xintercept = 0.5, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 0.5, xend = directional_share, yend = dataset),
                 colour = "#1F6FB2", linewidth = LWD) +
    geom_point(aes(x = detection_limit), shape = 124, size = 2.2, colour = INK) +
    geom_point(colour = "#1F6FB2", size = 1.6) +
    geom_text(aes(x = 1.30, label = sprintf("z %+.2f", z)), hjust = 1,
              family = FONT, size = gs(PTS), colour = INK) +
    # The note used to lie on the bottom data row.
    notiz("vertical bar = own detection limit; no single study reaches it") +
    labs(x = "directional share per study", y = NULL) +
    coord_cartesian(xlim = c(0.28, 1.32), clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))

  p <- pa / pb + plot_layout(heights = c(0.65, 1))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 1.5, 44)), "S3", SP2, 108)
}

# =============================================================================
# S4 -- patient vs control, downgraded
# =============================================================================
bau_s4 <- function() {
  a <- lies("S4A_patient_calibration")
  # Three cohorts carry no applicable marker set (whole blood, monocytes,
  # cell-free mRNA). Without a note they stood as empty rows in the field.
  a$calibration <- ifelse(is.na(a$z), "not calibratable",
                      ifelse(a$passed, "passed", "failed"))
  a$lab <- mit_gse(a$entity, a$accession)
  a <- a[order(a$z, na.last = FALSE), ]; a$lab <- factor(a$lab, levels = a$lab)
  pa <- ggplot(a, aes(z, lab)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = 2, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = z, yend = lab, colour = calibration),
                 linewidth = LWD) +
    geom_point(aes(colour = calibration), size = 1.6) +
    geom_text(aes(x = z + 0.2, label = paste0(tolower(marker_set), " markers")),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK,
              na.rm = TRUE) +
    geom_text(data = a[is.na(a$z), ], inherit.aes = FALSE,
              aes(x = 2.25, y = lab), hjust = 0, family = FONT,
              size = gs(PTS), colour = GRAU,
              label = "no marker set applicable to this tissue") +
    scale_colour_manual(values = EICH, guide = "none") +
    scale_x_continuous(breaks = c(0, 2, 4, 6)) +
    labs(x = "calibration of the cohort (z)", y = NULL) +
    # The view reaches to 8 so that the marker set has room behind the point
    # inside the field; it used to run past the right edge of the sheet.
    coord_cartesian(xlim = c(-1, 8), clip = "off") +
    theme_pub() + theme(plot.margin = margin(6, 4, 1.5, 2))

  b <- lies("S4B_patient_concordance")
  b$lab <- mit_gse(b$entity, b$accession)
  b <- b[order(b$concordance), ]; b$lab <- factor(b$lab, levels = b$lab)
  pb <- ggplot(b, aes(concordance, lab)) +
    geom_segment(aes(x = concordance_detection_limit, xend = concordance, yend = lab),
                 colour = GRAU, linewidth = LW) +
    geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.2,
               colour = INK) +
    geom_point(colour = "#1F6FB2", size = 1.6) +
    # The key used to lie on the bottom data row.
    notiz("no cohort reaches its own detection limit; this level carries a limit, not a result") +
    labs(x = "programme concordance, patient vs control", y = NULL) +
    coord_cartesian(xlim = c(0.28, 0.78), clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))

  p <- pa / pb
  tafel(p, m(list("A", 1.5, 1.5), list("B", 1.5, 57)), "S4", SP2, 118)
}

# =============================================================================
# S5 -- robustness and leave-one-out
# =============================================================================
bau_s5 <- function() {
  a <- lies("S5A_leave_one_out")
  a$left_out <- factor(a$left_out, levels = rev(a$left_out))
  pa <- ggplot(a, aes(S1_z, left_out)) +
    geom_vline(xintercept = 2, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = S1_z, yend = left_out), colour = "#12946B",
                 linewidth = LWD) +
    geom_point(colour = "#12946B", size = 1.6) +
    geom_text(aes(x = S1_z + 0.1, label = sprintf("P %.4f", S1_p)), hjust = 0,
              family = FONT, size = gs(PTS), colour = INK) +
    notiz("dashed line: z = 2", colour = DIFF) +
    labs(x = "programme between donors, leave-one-out (z)",
         y = "cell left out") +
    coord_cartesian(xlim = c(0, 5.6), clip = "off") +
    theme_pub() + theme(plot.margin = margin(7, 4, 1.5, 2))

  b <- lies("S5B_self_test")
  # Four statistics, two effect sizes. Without grouping, geom_line() joined
  # all eight points into a zigzag that meant nothing.
  KENN <- c(S1 = "#12946B", S2 = "#C2472A", S3a = "#1F6FB2", S3b = "#7B5EA7")
  b$statistic_id <- factor(b$statistic_id, levels = names(KENN))
  pb <- ggplot(b, aes(effect, fraction_abs_z_above_2, colour = statistic_id,
                      group = statistic_id)) +
    geom_hline(yintercept = 0.05, colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_line(linewidth = LWD) +
    geom_point(size = 1.4) +
    scale_colour_manual(values = KENN, guide = "none") +
    # A direct-labelled block in the free field at the top left: S1 and S3a
    # both end at 1.0, so at the end of the curves the names lay on top of
    # one another.
    txt(0.015, 0.99, "S1", colour = KENN[["S1"]], size = PTS, hjust = 0) +
    txt(0.015, 0.92, "S2", colour = KENN[["S2"]], size = PTS, hjust = 0) +
    txt(0.015, 0.85, "S3a", colour = KENN[["S3a"]], size = PTS, hjust = 0) +
    txt(0.015, 0.78, "S3b", colour = KENN[["S3b"]], size = PTS, hjust = 0) +
    notiz("dashed line: nominal 5 %", colour = DIFF) +
    labs(x = "injected effect", y = "fraction with |z| > 2") +
    coord_cartesian(xlim = c(-0.01, 0.37), clip = "off") + theme_pub() +
    theme(plot.margin = margin(7, 4, 1.5, 2))

  cc <- lies("S5C_calibration_sensitivity")
  nm <- names(cc)
  zv <- grep("^z_", nm, value = TRUE)
  pc <- if (length(zv) >= 2) {
    ggplot(cc, aes(.data[[zv[1]]], .data[[zv[2]]])) +
      geom_abline(slope = 1, intercept = 0, colour = GRAU, linewidth = LW) +
      geom_vline(xintercept = 2, colour = DIFF, linewidth = LW,
                 linetype = "22") +
      geom_hline(yintercept = 2, colour = DIFF, linewidth = LW,
                 linetype = "22") +
      geom_point(colour = "#1F6FB2", size = 1.4) +
      labs(x = "preregistered marker sets (z)",
           y = "in-vitro-reachable markers (z)") +
      coord_cartesian(clip = "off") + theme_pub() +
      theme(plot.margin = margin(6, 4, 1.5, 2))
  } else plot_spacer()

  # This figure has three panels. An earlier draft carried a fourth that was a
  # block of text rather than a data panel; it was removed, and what it said
  # now stands in the exclusion criteria of the Methods.
  p <- (pa | pb) / (pc | plot_spacer())
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 62)), "S5", SP2, 122)
}

# =============================================================================
# S6 -- orthogonal layers
# =============================================================================
bau_s6 <- function() {
  a <- lies("S6A_atac_complete")
  a$axis <- factor(a$axis, levels = c("adipogenic", "osteogenic",
                                        "lineage contrast"))
  a$window <- factor(a$window, levels = rev(c("P", "T10", "T50", "GB")))
  pa <- ggplot(a, aes(concordance, window, colour = null_model)) +
    geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.2,
               colour = INK) +
    geom_point(size = 1.6, position = position_dodge(width = 0.45)) +
    facet_grid(axis ~ ., switch = "y") +
    scale_colour_manual(values = c(background = "#9A9A9A",
                                   `H1 baseline-stratified` = "#1F6FB2"),
                        guide = "none") +
    notiz_oben("grey = background null,", "axis",
               factor("adipogenic", levels = c("adipogenic", "osteogenic",
                                               "lineage contrast"))) +
    notiz_oben("blue = H1 baseline-stratified null", "axis",
               factor("adipogenic", levels = c("adipogenic", "osteogenic",
                                               "lineage contrast")),
               zeile = 2) +
    labs(x = mit_gse("module concordance, chromatin ATAC", "GSE332758"),
         y = "window") +
    coord_cartesian(xlim = c(0.40, 0.78), clip = "off") +
    theme_pub() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(8, 4, 1.5, 2))

  level <- function(name, xlab, title, zweite = NULL) {
    d <- lies(name)
    if (!"concordance" %in% names(d)) return(plot_spacer())
    # The labels have been English in the panel file since the presentation
    # layer; here they are only assembled.
    lab <- if ("axis_name" %in% names(d)) d$axis_name else
      if ("axis" %in% names(d)) d$axis else rep("axis", nrow(d))
    if ("window" %in% names(d)) lab <- paste0(lab, ", ", d$window)
    if ("null_model" %in% names(d))
      lab <- paste0(lab, ifelse(d$null_model == "background",
                                ", background", ", stratified"))
    d$lab <- factor(lab, levels = rev(unique(lab)))
    d$ueber <- d$concordance > d$concordance_detection_limit
    ggplot(d, aes(concordance, lab)) +
      geom_segment(aes(x = concordance_detection_limit, xend = concordance, yend = lab),
                   colour = GRAU, linewidth = LW) +
      geom_point(aes(x = concordance_detection_limit), shape = 124, size = 2.2,
                 colour = INK) +
      geom_point(aes(colour = ueber), size = 1.6) +
      scale_colour_manual(values = c("TRUE" = "#12946B", "FALSE" = DIFF),
                          guide = "none") +
      # FigureFix S6B: the header line lay on the topmost data row.
      # notiz() places it outside the field, in the top margin.
      notiz(title) +
      (if (is.null(zweite)) NULL else notiz(zweite, zeile = 2)) +
      labs(x = xlab, y = NULL) +
      coord_cartesian(clip = "off") + theme_pub() +
      # The heading used to lie on the top data row. notiz() puts it outside
      # the field, in the top margin.
      theme(plot.margin = margin(19, 4, 1.5, 2))
  }

  # Accession rule: the information first, then the accession in parentheses.
  # The headings move 2 mm clear, so that the panel letters C and D get a line
  # of their own above them rather than standing inside them.
  pb <- level("S6B_h3k27ac", "module concordance",
              mit_gse("H3K27ac", "GSE129031"))
  pc <- level("S6C_methylome_27k", "module concordance",
              mit_gse("promoter methylome 27K", "GSE33896"),
              zweite = "null, with its own detection limit")
  pd <- level("S6D_methylome_450k", "module concordance",
              mit_gse("promoter methylome 450K", "GSE129266"),
              zweite = "below its own detection limit")

  # leerstreifen(): S6A is facetted and brings a strip column with it on the
  # left; without compensation the axes of B-D stand beside their fields.
  p <- (pa | leerstreifen(pb)) / (leerstreifen(pc) | leerstreifen(pd)) +
    plot_layout(heights = c(1.35, 1))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 80.5), list("D", 89, 80.5)), "S6", SP2, 142)
}

# =============================================================================
# S7 -- screen detail
# =============================================================================
bau_s7 <- function() {
  a <- lies("S7A_exclusion_codes")
  a$code <- factor(a$code, levels = rev(a$code))
  pa <- ggplot(a, aes(n, code)) +
    geom_col(width = 0.6, fill = "#1F6FB2", colour = NA) +
    # FigureFix S7A: the reason text needed 62 mm of right margin, squeezed
    # the field down to a third of the width ("far too far right") and was
    # still cut off. It stands in full in Table TS2 and in the caption; the
    # panel shows code and count.
    geom_text(aes(label = n), hjust = -0.4, family = FONT,
              size = gs(PTS), colour = INK) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "series excluded (reasons in Table S2)", y = "code") +
    coord_cartesian(xlim = c(0, 86), clip = "off") +
    theme_pub() + theme(plot.margin = margin(4, 4, 1.5, 2))

  b <- lies("S7B_screen_diagnoses_full")
  bs <- as.data.frame(table(b$axis), stringsAsFactors = FALSE)
  names(bs) <- c("axis", "n")
  bs <- bs[order(bs$n), ]
  bs$axis <- factor(bs$axis, levels = bs$axis)
  pb <- ggplot(bs, aes(n, axis)) +
    geom_col(width = 0.6, fill = "#0E8C8C", colour = NA) +
    geom_text(aes(label = n), hjust = -0.3, family = FONT, size = gs(PTS),
              colour = INK) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "hand-checked series", y = "disease axis searched") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(6, 8, 1.5, 2))

  cc <- lies("S7C_screen_by_design")
  cs <- as.data.frame(table(cc$code), stringsAsFactors = FALSE)
  names(cs) <- c("code", "n")
  cs <- cs[order(cs$n), ]
  cs$code <- factor(cs$code, levels = cs$code)
  pc <- ggplot(cs, aes(n, code)) +
    geom_col(width = 0.6, fill = "#7B5EA7", colour = NA) +
    geom_text(aes(label = n), hjust = -0.3, family = FONT, size = gs(PTS),
              colour = INK) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "series screened by design", y = "code") +
    coord_cartesian(clip = "off") + theme_pub() +
    theme(plot.margin = margin(6, 8, 1.5, 2))

  p <- pa / (pb | pc) + plot_layout(heights = c(1, 1.15))
  tafel(p, m(list("A", 1.5, 1.5), list("B", 1.5, 56), list("C", 89, 56)),
        "S7", SP2, 120)
}

# =============================================================================
# S8 -- day-zero falls, and the publication matching
# =============================================================================
bau_s8 <- function() {
  a <- lies("S8A_day_zero_falls")
  a$variant <- factor(a$variant, levels = rev(a$variant))
  pa <- ggplot(a, aes(amplitude_z, variant)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = c(-2, 2), colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = amplitude_z, yend = variant),
                 colour = "#1F6FB2", linewidth = LWD) +
    geom_point(colour = "#1F6FB2", size = 1.6) +
    geom_text(aes(x = amplitude_z - 0.15,
                  label = sprintf("P %.3f (n %d)", amplitude_p, n)),
              hjust = 1, family = FONT, size = gs(PTS), colour = INK) +
    labs(x = "day-0 module value predicting later amplitude (z)", y = NULL) +
    coord_cartesian(xlim = c(-3.4, 2.6), clip = "off") +
    theme_pub() + theme(plot.margin = margin(6, 4, 1.5, 2))

  b <- lies("S8B_cells_day_zero")
  b$calibration <- ifelse(b$passed, "passed", "failed")
  pb <- ggplot(b, aes(calibration, undifferentiated_module_z, colour = calibration)) +
    geom_point(size = 1.5, position = position_jitter(width = 0.14, seed = 4)) +
    scale_colour_manual(values = EICH, guide = "none") +
    # FigureFix S8B: the key lay on the topmost points.
    notiz("cells that later pass and fail their calibration") +
    notiz("do not differ at day 0 (U 23.5, P 0.95)", zeile = 2) +
    labs(x = "later lineage-marker calibration",
         y = "module value at day 0 (z)") +
    coord_cartesian(ylim = c(min(b$undifferentiated_module_z) * 1.12, -0.06),
                    clip = "off") + theme_pub() +
    theme(plot.margin = margin(10, 4, 1.5, 2))

  cc <- lies("S8C_study_intensity")
  cc <- cc[order(cc$factor), ]
  cc$role <- ifelse(cc$gene_set == "programme (173 genes)", "programme",
                     ifelse(cc$gene_set == "cell cycle (neg. control)", "negative control",
                            "panel"))
  # Rule 0: the machine keys of the gene sets are never shown.
  SLAB <- c("PanelApp 309" = "PanelApp 309",
            "Nosology (core)" = "Nosology (core)",
            "Nosology (broad)" = "Nosology (broad)",
            "short stature" = "short stature",
            "height GWAS" = "height GWAS",
            "cell cycle (neg. control)" = "cell cycle (neg. control)",
            "programme (173 genes)" = "programme (173 genes)")
  cc$gene_set <- factor(SLAB[cc$gene_set],
                        levels = SLAB[cc$gene_set][order(cc$factor)])
  pc <- ggplot(cc, aes(factor, gene_set, colour = role)) +
    geom_vline(xintercept = 1, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 1, xend = factor, yend = gene_set), linewidth = LWD) +
    geom_point(size = 1.6) +
    geom_text(aes(x = factor + 0.03,
                  label = sprintf("%.2f x (median %d)", factor,
                                  as.integer(median_publications))),
              hjust = 0, family = FONT, size = gs(PTS), colour = INK) +
    scale_colour_manual(values = c(programme = "#12946B",
                                   `negative control` = "#7B5EA7",
                                   panel = "#1F6FB2"), guide = "none") +
    labs(x = "publications per gene, relative to background", y = NULL) +
    coord_cartesian(xlim = c(0.95, 2.75), clip = "off") +
    theme_pub() + theme(plot.margin = margin(4, 4, 1.5, 2))

  d <- lies("S8D_publication_matching")
  d <- d[d$run %in% c("P1_expression_ohne_pub", "P1_expression_mit_pub"), ]
  d$matching <- ifelse(grepl("log_pub", d$matching), "with publications",
                       "without publications")
  lab <- c("PanelApp 309" = "PanelApp 309", "Nosology (core)" = "Nosology (core)",
           "Nosology (broad)" = "Nosology (broad)", "short stature" = "short stature",
           "height GWAS" = "height GWAS",
           "cell cycle (neg. control)" = "cell cycle (neg. control)",
           "programme (173 genes)" = "programme")
  d$satzlab <- factor(lab[d$gene_set], levels = rev(unname(lab)))
  pd <- ggplot(d, aes(z, satzlab, colour = matching)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_line(aes(group = satzlab), colour = GRAU, linewidth = LW) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = c(`without publications` = "#C2472A",
                                   `with publications` = "#1F6FB2"),
                        guide = "none") +
    # FigureFix S8D: the two-line colour key lay on the topmost two data
    # rows.
    notiz("red = without publication matching", colour = "#C2472A") +
    notiz("blue = with publication matching", colour = "#1F6FB2",
          zeile = 2) +
    labs(x = "absolute expression, matched (z)", y = NULL) +
    coord_cartesian(xlim = c(-3.4, 5.4), ylim = c(0.5, 7.5), clip = "off") +
    theme_pub() + theme(plot.margin = margin(10, 4, 1.5, 2))

  p <- (pa | pb) / (pc | pd)
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 64), list("D", 89, 64)), "S8", SP2, 128)
}

# --- S9: the lineage contrast in three cohorts -------------------------------
bau_s9 <- function() {
  # Accession rule: the assay first, then the accession (the panel files have
  # carried the names in this form since the presentation layer).
  KOH <- c("ATAC (GSE332758)" = "#E07B12", "ATAC (GSE151311)" = "#1F6FB2",
           "H3K27ac (GSE151315)" = "#7B5EA7")
  FW <- c("P", "T10", "T50", "GB")

  # A  calibration L: contrast vs own limit, 3 cohorts x 4 windows
  a <- lies("S9A_calibration_L_three_cohorts")
  a$cohort_label <- factor(a$cohort_label, levels = names(KOH))
  a$window <- factor(a$window, levels = rev(FW))
  pa <- ggplot(a, aes(contrast, window, colour = cohort_label)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_segment(aes(x = 0, xend = contrast, yend = window),
                 linewidth = LWD) +
    geom_point(aes(x = detection_limit), shape = 124, size = 2.4, colour = INK) +
    geom_point(size = 1.7) +
    facet_grid(cohort_label ~ ., switch = "y") +
    scale_colour_manual(values = KOH, guide = "none") +
    notiz_oben("bar = own detection limit; 0 of 12 reach it",
               "cohort_label", factor(names(KOH)[1], levels = names(KOH))) +
    labs(x = "calibration L: osteogenic minus adipogenic markers",
         y = "window") +
    coord_cartesian(xlim = c(-0.42, 0.72), clip = "off") +
    theme_pub() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(8, 4, 1.5, 2))

  # B  where each marker set sits on the difference axis
  b <- lies("S9B_marker_sets_lineage_axis")
  b$cohort_label <- factor(b$cohort_label, levels = names(KOH))
  b$gene_set <- factor(b$gene_set, levels = c("Osteogenic", "Adipogenic",
                                              "Undifferentiated"))
  pb <- ggplot(b, aes(z, gene_set, colour = cohort_label)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = c(-2, 2), colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.16,
                                                      seed = 11)) +
    facet_grid(cohort_label ~ ., switch = "y") +
    scale_colour_manual(values = KOH, guide = "none") +
    scale_y_discrete(limits = rev) +
    notiz_oben("osteogenic markers never above z +2",
               "cohort_label", factor(names(KOH)[1], levels = names(KOH))) +
    labs(x = "marker set on the lineage axis (z)", y = NULL) +
    coord_cartesian(xlim = c(-4.1, 3.4), clip = "off") +
    theme_pub() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(8, 4, 1.5, 2))

  # C  the decoupling, independently replicated
  cc <- lies("S9C_decoupling_second_cohort")
  cc$window <- factor(cc$window, levels = FW)
  g <- rbind(
    data.frame(window = cc$window, axis = cc$axis, quantity = "calibration",
               value = cc$calibration_z, detection_limit = 2),
    data.frame(window = cc$window, axis = cc$axis, quantity = "module",
               value = cc$concordance_z, detection_limit = NA))
  g$quantity <- factor(g$quantity, levels = c("calibration", "module"))
  pc <- ggplot(g, aes(window, value, colour = quantity, shape = axis)) +
    geom_hline(yintercept = 0, colour = LINE, linewidth = LW) +
    geom_hline(yintercept = 2, colour = DIFF, linewidth = LW, linetype = "22") +
    geom_point(size = 1.9, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(calibration = "#C2472A",
                                   module = "#12946B"), guide = "none") +
    scale_shape_manual(values = c(osteogenic = 16, adipogenic = 17),
                       guide = "none") +
    # FigureFix S9C: the three key lines lay on the points.
    notiz("module (green): 8/8 above own limit") +
    notiz("calibration (red): 0/8 pass", zeile = 2) +
    notiz("circle = osteogenic, triangle = adipogenic", zeile = 3) +
    labs(x = mit_gse("window, H3K27ac", "GSE151315"), y = "statistic (z)") +
    coord_cartesian(ylim = c(-0.5, 8), clip = "off") +
    theme_pub() + theme(plot.margin = margin(13, 4, 1.5, 2))

  # D  why the calibration fails: undifferentiated state left, lineage not
  # reached
  d <- lies("S9D_marker_sets_single_axes")
  d$axis <- factor(paste(d$axis, "axis"),
                    levels = c("osteogenic axis", "adipogenic axis"))
  d$gene_set <- factor(d$gene_set, levels = c("Osteogenic", "Adipogenic",
                                              "Undifferentiated"))
  pd <- ggplot(d, aes(z, gene_set, colour = gene_set)) +
    geom_vline(xintercept = 0, colour = LINE, linewidth = LW) +
    geom_vline(xintercept = c(-2, 2), colour = DIFF, linewidth = LW,
               linetype = "22") +
    geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.14,
                                                      seed = 12)) +
    facet_grid(axis ~ ., switch = "y") +
    scale_colour_manual(values = c(Osteogenic = "#E07B12",
                                   Adipogenic = "#7B5EA7",
                                   Undifferentiated = "#1F6FB2"),
                        guide = "none") +
    scale_y_discrete(limits = rev) +
    notiz_oben("undifferentiated state left, lineage not reached",
               "axis", factor("osteogenic axis",
                               levels = c("osteogenic axis", "adipogenic axis"))) +
    labs(x = "marker set on the differentiation axis (z)", y = NULL) +
    coord_cartesian(xlim = c(-3.4, 3.4), clip = "off") +
    theme_pub() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                           margin = margin(r = 3)),
          panel.spacing.y = unit(2.2, "mm"),
          plot.margin = margin(8, 4, 1.5, 2))

  # leerstreifen(): A, B and D are facetted; C evens out the column count.
  p <- (pa | pb) / (leerstreifen(pc) | pd)
  tafel(p, m(list("A", 1.5, 1.5), list("B", 89, 1.5),
             list("C", 1.5, 70), list("D", 89, 70)), "S9", SP2, 138)
}

main <- function() {
  cat("61_figures_supplement.R -- supplementary figures to", PUB_DIR, "\n")
  bau_s1(); bau_s2(); bau_s3(); bau_s4()
  bau_s5(); bau_s6(); bau_s7(); bau_s8(); bau_s9()
  si <- file.path(WURZEL, "results", "figures_supplement_sessionInfo.txt")
  capture.output(sessionInfo(), file = si)
  cat("sessionInfo ->", si, "\n")
}

if (!isTRUE(getOption("abbildungen.nur_funktionen"))) main()
