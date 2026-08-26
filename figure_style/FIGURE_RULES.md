# Figure rules

The grammar every figure in this repository follows. The style itself lives in
`publication_style.R`; the figures are drawn by `code/60_figures_main.R` and
`code/61_figures_supplement.R`, the graphical abstract by
`code/62_graphical_abstract.py`.

Column widths follow Cell Press: `SP1` 85 mm, `SP15` 114 mm, `SP2` 174 mm.
Every sheet is written as a vector PDF **and** as a PNG at 600 dpi.

## The rules that are not negotiable

1. **Colour codes the programme, never the genotype.** The genotype runs
   through the fill: control filled, knockout open. The condition colours are
   in `BED`, the atlas states in `ZUSTAND`.
2. **Type and axes are black**, Arial, 7 pt (`PT`) for labelling and 6 pt
   (`PTS`) for statistics and secondary annotation. Nothing is bold except the
   panel letter, and that is set on the sheet, not in the plot (`PTL` 8 pt).
3. **No title, no prose and no source-script note inside the image.**
   Everything explanatory belongs in the figure legend
   (`manuscript/CAPTIONS_MAIN.md`, `manuscript/CAPTIONS_SUPPLEMENT.md`), with
   every number and every limitation.
4. **The key is a direct label at data coordinates**, not a ggplot legend: a
   symbol plus black text in a free area of the field.
5. **No `scale_*_continuous(limits =)`.** It discards labels silently. The
   view is always set with `coord_cartesian(..., clip = "off")`.
6. **No colour vector passed to a geom.** Always an aesthetic plus
   `scale_*_manual()`.
7. **Text never crosses a plot element.** The reliable place for a note is the
   margin above the field; `notiz()` and `notiz_oben()` put it there.
8. **The accession always follows the information, in parentheses:**
   `ATAC (GSE332758)`, never `GSE332758 (ATAC)` and never `ATAC GSE332758`.
   `mit_gse()` produces that form, and `code/72_check_language.py` enforces it.

## Layout grammar: what makes a figure calm

* **One row per category, one column per condition**
  (`facet_grid(class ~ column, scales = "free_y", space = "free_y",
  switch = "y")`). Category names read horizontally, and group headings are a
  facet strip rather than a bracket over the field.
* **The comparison of interest sits above one another, not side by side.**
* **One shared axis for all fields.** Free axes per field are not comparable
  and force labelling in every field.
* **Thresholds as a single continuous line**, labelled once.
* **Zebra stripes per row** (`geom_rect` in `SOFT`) lead the eye across the
  columns without introducing grid lines.
* **A different kind of data gets a different symbol shape.** Reference data as
  a diamond beside the circles of the measurements; the fill stays reserved for
  the genotype.
* **n = 3: point = median, line = range.** No box and no error bars, and the
  legend says explicitly that the line is the minimum to maximum.

## Composing panels into one sheet

A journal figure is a sheet with panels A, B, C, flush edges and one legend --
not a staircase of differently wide images stacked on top of one another.

**One target width for every sheet of a document**, here `SP2` = 174 mm. If
every sheet is drawn at that width and placed at exactly that width, Arial 7 pt
is the same size across all sheets, with nothing scaled anywhere.

### The four things that broke here

1. **Do not set `plot.margin` globally.** Several panel functions carry a wide
   margin of their own because their labelling lies outside the field. A
   uniform `plot.margin` cuts exactly that labelling off. Set every margin
   individually.
2. **Panel letters go on the sheet with `grid.text`, not through ggplot.**
   `plot.tag` is placed freely over the drawing area and lands on the rotated
   axis title; `plot.title` hangs on the left edge of the panel block, which
   patchwork stretches, so the letter moves with it. After `print(p)`:

   ```r
   grid::grid.text(m$lab, x = grid::unit(m$x, "mm"),
                   y = grid::unit(1, "npc") - grid::unit(m$y, "mm"),
                   just = c("left", "top"),
                   gp = grid::gpar(fontfamily = FONT, fontface = "bold",
                                   fontsize = PTL, col = INK))
   ```

3. **The letter belongs beside the labelling, not at the edge of the sheet.**
   Measure the leftmost dark pixel per row and allow about 5 mm of air.
4. **patchwork aligns columns, not panels.** To the left of the field it
   creates one shared column each for the margin, the axis title and the axis
   text, each as wide as its maximum across all panels -- and it matches those
   columns by position. A facetted panel with `strip.placement = "outside"`
   brings an extra strip column with it, so a non-facetted panel in the same
   sheet slides its axis into that strip column and ends up standing *beside*
   its field rather than at it. That was the detached y axis in F4D.
   The remedy in this repository is `leerstreifen()` in
   `publication_style.R`: it gives a non-facetted panel an empty outer strip,
   so the column counts match. It is applied in `bau_f2`, `bau_f4`, `bau_s6`
   and `bau_s9`.

### What keeps the block narrow to begin with

* **Two-line category names cost twice the block width.** The parenthetical
  line belongs in the legend.
* **Long rotated axis titles** reach from the middle of the field into the row
  of the panel letter. Set them on two lines.
* **Panel letters are assigned by hand, not with `tag_levels = "A"`.** With
  nested `wrap_plots()` patchwork counts the nesting level, and a
  `plot_spacer()` consumes a letter.

## Pitfalls that have already cost time here

* **A colour vector plus facets gives shifted colours.**
  `geom_point(colour = G$colour)` assigns by position; with facets ggplot
  re-splits the data and a point ends up in the wrong colour in the wrong
  field. Always `aes(colour = key)` plus `scale_colour_manual()`.
* **`element_text()` knows `face`, not `fontface`.** A call with `fontface`
  only warns and is ignored, so the italics go missing silently.
* **The left margin and the row labels belong together.** Row labels drawn with
  `geom_text(x = -0.5, hjust = 1)` lie outside the field and are cut off
  without a matching `plot.margin`; 7 pt Arial needs about 1.4 mm per
  character.
* **Numbers at a data point always go to the right of the widest element of
  the row** (`pmax(value, null + sd) + gap`). Switching sides by direction
  collides with the zero line in one row and with the row label in the next.
* **A value belongs only in a figure that shows the same quantity.** When
  adding a row, check which quantity the figure carries and take the
  comparison value from the *same* run.
* **Look at the result before reporting it.** Read the PNG; do not treat a run
  without an error as a finished figure. For sheets, measure rather than judge
  by eye -- the leftmost dark pixel per row says in one line whether the block
  is too wide:

  ```python
  a = np.array(Image.open(png).convert("L"))
  cols = np.where((a[y0:y1] < 128).any(axis=0))[0]
  ```

## What belongs to a change of a figure

1. the panel itself (PDF and PNG),
2. the call in `code/60_figures_main.R` or `code/61_figures_supplement.R`
   with updated dimensions,
3. the legend in `manuscript/CAPTIONS_MAIN.md` or
   `manuscript/CAPTIONS_SUPPLEMENT.md`, including the changed axis description,
4. the panel table in `README.md`,
5. every sheet in which the panel appears, rebuilt and looked at -- a panel
   that fits on its own can collide in the narrower field of a sheet,
6. if labelling has moved from the image into the legend, the legend must then
   really carry it,
7. `python code/70_check_numbers.py`, which fails if a number in the text or in
   a legend no longer matches its panel file.
