# =============================================================================
# publication_style.R -- the publication style of the data figures
# =============================================================================
# The page set-up of a journal (Cell Press figure guidelines), because these
# images go into the manuscript as panels:
#
#   * Arial throughout, 7 pt for everything that labels, 6 pt for statistics
#     and secondary annotation. Nothing bold except the panel letter.
#   * Type and axes BLACK. Colour carries data, never text.
#   * Axis line and ticks present, ticks outwards, no grid lines, no frame and
#     no title in the image (title and legend are set by the document).
#   * Widths in mm by column measure: 85 (one column), 114 (1.5), 174 (two).
#   * Lines at least 0.5 pt, points at least 1.2 mm -- so that they stay
#     legible when the figure is reduced.
#   * Every figure as a PDF (vector, Arial embedded) and as a PNG at 600 dpi.
#     No source-script note in the image; that belongs in the legend.
#
# THE COLOUR SEMANTICS are binding: colour codes the programme, never the
# genotype. The genotype runs through the fill (control filled, knockout open).
#
# The rules the figures follow are written out in FIGURE_RULES.md.
# =============================================================================

suppressMessages({
  library(ggplot2); library(grid); library(ragg); library(systemfonts)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
# Output directory. The figure scripts overwrite PUB_DIR after sourcing this
# file, so that this file itself never has to know where it is being used.
.stilordner <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
}
PUB_DIR <- file.path(.stilordner(), "figures")

# ---------------------------------------------------------- colour tokens
APP  <- "#1F6FB2"   # secretory machinery
SEK  <- "#12946B"   # secretome
GRAU <- "#9A9A9A"   # confounders, unclassified categories
DIFF <- "#C2472A"   # exclusively the control - knockout distance
INK  <- "#000000"   # type and axes
MUT  <- "#000000"   # the publication version has no grey type
LINE <- "#000000"
SOFT <- "#F0F0EE"

BED <- c("undifferentiated" = "#7B5EA7", "chondrogenic" = "#0E8C8C",
         "osteogenic" = "#E07B12")

hell <- function(farbe, anteil = 0.45) {
  r <- grDevices::col2rgb(farbe) / 255
  grDevices::rgb(t(r + (1 - r) * anteil))
}
dunkel <- function(farbe, anteil = 0.40) {
  r <- grDevices::col2rgb(farbe) / 255
  grDevices::rgb(t(r * (1 - anteil)))
}

# The six reference states of the bone-marrow atlas
ZUSTAND <- c("Adipo-MSC"  = "#E4CFA6", "THY1+ MSC"  = "#BFD9B8",
             "Fibro-MSC"  = "#E0C6D8", "APOD+ MSCs" = "#C9CFE9",
             "Osteo-MSC"  = "#8FBEDE", "Osteoblast" = "#1F6FB2")

# ------------------------------------------------------- dimensions, type
mm <- function(x) x / 25.4              # mm -> inch
SP1 <- 85; SP15 <- 114; SP2 <- 174      # column widths in mm
DPI <- 600

PT   <- 7.0     # axis text, categories, direct labels
PTS  <- 6.0     # statistics, n counts, secondary axes
PTL  <- 8.0     # panel letter (bold)
LW   <- 0.25    # axes and thin lines (~0.53 pt)
LWD  <- 0.45    # data lines (~0.96 pt)

FONT <- {
  v <- systemfonts::system_fonts()$family
  c("Arial", "Helvetica", "Liberation Sans", "sans")[
    c("Arial", "Helvetica", "Liberation Sans", "sans") %in% c(v, "sans")][1]
}

# text size in ggplot units
gs <- function(pt) pt / .pt

# ------------------------------------------------------------------ the theme
theme_pub <- function(base = PT, links = TRUE, unten = TRUE) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(
      text              = element_text(colour = INK),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      panel.grid        = element_blank(),
      panel.background  = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      axis.line.x       = if (unten) element_line(colour = LINE, linewidth = LW) else element_blank(),
      axis.line.y       = if (links) element_line(colour = LINE, linewidth = LW) else element_blank(),
      axis.ticks        = element_line(colour = LINE, linewidth = LW),
      axis.ticks.length = unit(1.4, "pt"),
      axis.text         = element_text(colour = INK, size = base),
      axis.text.x       = element_text(margin = margin(t = 1.2)),
      axis.text.y       = element_text(margin = margin(r = 1.2)),
      axis.title        = element_text(colour = INK, size = base),
      axis.title.x      = element_text(margin = margin(t = 2.5)),
      axis.title.y      = element_text(margin = margin(r = 2.5)),
      legend.position   = "none",
      strip.text        = element_text(colour = INK, size = base, hjust = 0,
                                       margin = margin(b = 1.5)),
      plot.margin       = margin(2, 2, 1.5, 2)
    )
}

# A label inside the panel: black, plain, 7 pt. Colour only where it names a
# data series directly and thereby replaces the legend.
txt <- function(x, y, label, colour = INK, size = PT, hjust = 0, vjust = 0.5,
                fett = FALSE, ...)
  annotate("text", x = x, y = y, label = label, colour = colour,
           size = gs(size), family = FONT, hjust = hjust, vjust = vjust,
           fontface = if (fett) "bold" else "plain", ...)

# A label with a white halo. For text that unavoidably sits on a line or on a
# field of points: eight white copies in a circle, then the black type on top.
# rx and ry are steps in data units, because the axes here are not square.
saum <- function(d, rx, ry, size = PTS, colour = INK, hjust = 0.5, vjust = 0.5,
                 parse = FALSE, ...) {
  o <- expand.grid(dx = c(-1, 0, 1), dy = c(-1, 0, 1))
  o <- o[!(o$dx == 0 & o$dy == 0), ]
  c(lapply(seq_len(nrow(o)), function(i)
      geom_text(data = d, aes(x + o$dx[i] * rx, y + o$dy[i] * ry, label = lab),
                inherit.aes = FALSE, family = FONT, size = gs(size),
                colour = "white", hjust = hjust, vjust = vjust,
                parse = parse, ...)),
    list(geom_text(data = d, aes(x, y, label = lab), inherit.aes = FALSE,
                   family = FONT, size = gs(size), colour = colour,
                   hjust = hjust, vjust = vjust, parse = parse, ...)))
}

# A colour ramp: v in [0,1] mapped onto two end colours. Needed where the
# colour depth carries a second quantity while the hue still codes the class --
# two ggplot scales of the same aesthetic are not possible.
rampe <- function(v, lo, hi) {
  v <- pmin(pmax(v, 0), 1)
  grDevices::rgb(grDevices::colorRamp(c(lo, hi))(v), maxColorValue = 255)
}

# The panel letter -- only where explicitly wanted; the panels are composed by
# hand and get their letters there.
buchstabe <- function(p, b) p + labs(tag = b) +
  theme(plot.tag = element_text(family = FONT, face = "bold", size = PTL,
                                colour = INK),
        plot.tag.position = c(0, 1))

# The accession rule: the information first, THEN the accession in parentheses
# -- "ATAC (GSE332758)", never "GSE332758 (ATAC)" and never "ATAC GSE332758".
# The form lives only here; every axis, legend and direct label calls
# mit_gse(), so that there is no second version of it. Without an accession
# the label reads "(accession pending)", which is meant to be conspicuous.
mit_gse <- function(info, akzession = NULL) {
  leer <- is.null(akzession) | is.na(akzession) |
    !nzchar(as.character(akzession)) |
    startsWith(as.character(akzession), "GSE_")
  ifelse(leer, paste0(info, " (accession pending)"),
         paste0(info, " (", akzession, ")"))
}

# An empty outer strip on the left -- alignment in mixed sheets.
# patchwork aligns the columns of a sheet by their position, not by their name.
# A facetted panel with strip.placement = "outside" brings an extra strip
# column with it on the left; a non-facetted panel in the same sheet then
# slides its axis into exactly that column and stands visibly BESIDE its data
# field rather than at it. An empty strip matches the column count; it is not
# visible.
leerstreifen <- function(p) {
  if (!inherits(p$facet, "FacetNull")) return(p)
  if (is.null(p$data) || inherits(p$data, "waiver") ||
      !is.data.frame(p$data) || nrow(p$data) == 0) return(p)
  p$data$.streifen <- ""
  p + facet_grid(rows = ggplot2::vars(.streifen), switch = "y") +
    theme(strip.placement = "outside", strip.text.y.left = element_blank(),
          strip.background = element_blank())
}

# significance bracket
klammer <- function(x1, x2, y, label, tick = 0.02, size = PTS) {
  list(annotate("segment", x = x1, xend = x2, y = y, yend = y,
                colour = INK, linewidth = LW),
       annotate("segment", x = c(x1, x2), xend = c(x1, x2), y = y,
                yend = y - tick, colour = INK, linewidth = LW),
       annotate("text", x = (x1 + x2) / 2, y = y + tick * 0.6, label = label,
                colour = INK, family = FONT, size = gs(size), vjust = 0))
}

# p values as a column, in journal setting: 8.5 x 10^-9 or 0.86.
# The return value is plotmath -- always draw it with parse = TRUE.
pexp <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("''")
    if (x >= 0.01) return(sprintf("'%.2f'", x))
    e <- floor(log10(x))
    sprintf("'%.1f'%%*%%10^%d", x / 10^e, e) }, "")
}

# A p value in journal setting: P = 2 x 10^-25
pfmt <- function(p, ziffern = 1) {
  if (is.na(p)) return("")
  if (p >= 0.01) return(sprintf("italic('P')~'= %.2f'", p))
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("italic('P')~'= %.*f'%%*%%10^%d", ziffern, round(m, ziffern), e)
}

# ------------------------------------------------------------------ saving
sichern <- function(p, name, breite = SP1, hoehe = 60) {
  dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(file.path(PUB_DIR, paste0(name, ".pdf")),
                       width = mm(breite), height = mm(hoehe), family = FONT)
  print(p); invisible(grDevices::dev.off())
  ragg::agg_png(file.path(PUB_DIR, paste0(name, ".png")),
                width = mm(breite), height = mm(hoehe), units = "in",
                res = DPI, background = "white")
  print(p); invisible(grDevices::dev.off())
  cat(sprintf("  %-26s %3.0f x %3.0f mm\n", name, breite, hoehe))
}
