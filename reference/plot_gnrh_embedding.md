# Plot GnRH embedding

Plot GnRH embedding

## Usage

``` r
plot_gnrh_embedding(
  object,
  group_by = "gnrh_status",
  reduction = NULL,
  dims = c(1, 2),
  shuffle = FALSE,
  raster = NULL,
  raster.dpi = c(2048, 2048),
  alpha = 0.9,
  background_alpha = 0.18,
  n.cells = TRUE,
  percentage = FALSE,
  label = FALSE,
  repel = TRUE,
  label.size = 4,
  label.face = "plain",
  cols = NULL,
  axes = TRUE,
  plot.ttl = NULL,
  legend = TRUE,
  leg.ttl = NULL,
  leg.ttl.size = NULL,
  item.size = 3.5,
  leg.pos = "right",
  leg.dir = "vertical",
  leg.size = NULL,
  leg.ncol = NULL,
  item.border = TRUE,
  txtsize = 12,
  pt.size = NULL,
  dark = FALSE,
  total.cells = FALSE,
  style = "classic",
  ...
)
```

## Arguments

- object:

  A Seurat object.

- group_by:

  Metadata column used to colour cells.

- reduction:

  Dimensional reduction; an available reduction is selected when `NULL`.

- dims:

  Two reduction dimensions to display.

- shuffle:

  Randomize plotting order.

- raster:

  Use rasterized points; selected automatically when `NULL`.

- raster.dpi:

  Raster resolution passed to the raster geom.

- alpha, background_alpha:

  Opacity for highlighted and background cells.

- n.cells, percentage:

  Add cell counts or percentages to legend labels.

- label, repel, label.size, label.face:

  Cluster-label controls.

- cols:

  Optional named colour vector.

- axes:

  Show embedding axes.

- plot.ttl:

  Optional plot title.

- legend:

  Show the legend.

- leg.ttl, leg.ttl.size:

  Legend title and title size.

- item.size, item.border:

  Legend-key controls.

- leg.pos, leg.dir, leg.size, leg.ncol:

  Legend layout controls.

- txtsize:

  Base text size.

- pt.size:

  Point size; selected automatically when `NULL`.

- dark:

  Use dark display mode.

- total.cells:

  Include total-cell counts in the legend.

- style:

  Theme style.

- ...:

  Additional graphical arguments.

## Value

A ggplot object.
