# Visualize cells on a Seurat embedding

Creates a customizable two-dimensional or three-dimensional embedding
plot for Seurat objects. This function wraps
[`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)
and adds improved color handling, automatic rasterization, custom
legends, optional cluster labels, figure-panel styling, dark mode, and
support for plotting multiple metadata variables.

## Usage

``` r
cellmap(
  object,
  group.by = NULL,
  reduction = "umap",
  dims = c(1, 2),
  shuffle = FALSE,
  raster = FALSE,
  stroke.size = NULL,
  raster.dpi = c(2048, 2048),
  alpha = 1,
  repel = FALSE,
  n.cells = TRUE,
  label = FALSE,
  label.size = 4,
  label.face = "plain",
  cols = NULL,
  figplot = FALSE,
  axes = TRUE,
  plot.ttl = NULL,
  legend = TRUE,
  leg.ttl = NULL,
  leg.ttl.size = txtsize,
  item.size = 4,
  leg.pos = "right",
  leg.just = "center",
  leg.dir = "vertical",
  leg.size = 10,
  leg.ncol = NULL,
  item.border = TRUE,
  txtsize = 12,
  pt.size = NULL,
  dark = FALSE,
  total.cells = FALSE,
  threeD = FALSE,
  style = "test",
  facet.bg = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- group.by:

  Metadata column used to color cells. If `NULL`, active identities are
  used. A character vector of multiple metadata columns can be supplied
  to generate a patchwork of plots.

- reduction:

  Dimensional reduction to plot. Default is `"umap"`.

- dims:

  Numeric vector specifying dimensions to plot. Use two values for 2D
  plots or three values for 3D plots. Default is `c(1, 2)`.

- shuffle:

  Logical; randomly shuffle plotting order of cells. Default is `FALSE`.

- raster:

  Logical; rasterize points for large datasets. If `NULL`, rasterization
  is enabled automatically for objects with more than 100,000 cells.
  Default is `NULL`.

- stroke.size:

  Optional point stroke size passed to
  [`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html).

- raster.dpi:

  Numeric vector of length 2 specifying rasterization resolution.
  Default is `c(2048, 2048)`.

- alpha:

  Numeric point transparency. Default is `1`.

- repel:

  Logical; repel cluster labels when labels are drawn. Default is
  `FALSE`.

- n.cells:

  Logical; append the number of cells per group to legend labels.
  Default is `TRUE`.

- label:

  Logical; add cluster labels to the embedding. Default is `FALSE`.

- label.size:

  Numeric label text size. Default is `4`.

- label.face:

  Font face for cluster labels. Default is `"plain"`.

- cols:

  Optional named color vector. Missing groups are assigned `"gray70"`.
  If unnamed, colors are matched to group levels in order.

- figplot:

  Logical; generate a minimal figure-style plot with compact arrow axes.
  Default is `FALSE`.

- axes:

  Logical; show embedding axes. Default is `TRUE`.

- plot.ttl:

  Optional plot title.

- legend:

  Logical; show legend. Default is `TRUE`.

- leg.ttl:

  Optional legend title. If `NULL`, `group.by` is used.

- leg.ttl.size:

  Numeric legend title size. Default is `txtsize`.

- item.size:

  Numeric legend item size. Default is `4`.

- leg.pos:

  Legend position passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `"right"`.

- leg.just:

  Legend justification. Default is `"center"`.

- leg.dir:

  Legend direction. Default is `"vertical"`.

- leg.size:

  Numeric legend text size. Default is `10`.

- leg.ncol:

  Number of legend columns. If `NULL`, this is selected automatically
  based on the number of groups.

- item.border:

  Logical; draw borders around legend keys. Default is `TRUE`.

- txtsize:

  Base text size. Default is `12`.

- pt.size:

  Numeric point size. If `NULL`, a point size is chosen automatically
  based on cell number and rasterization mode.

- dark:

  Logical; apply a dark theme. Default is `FALSE`.

- total.cells:

  Logical; append total cell number to the plot title. Default is
  `FALSE`.

- threeD:

  Logical; generate an interactive 3D plot using plotly. Default is
  `FALSE`.

- style:

  Theme style passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `"test"`.

- facet.bg:

  Logical; show facet background in the applied theme. Default is
  `FALSE`.

- ...:

  Additional arguments passed to
  [`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)
  and
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).

## Value

A `ggplot2` object for 2D plots, a patchwork object when multiple
`group.by` variables are supplied, or a plotly object for 3D plots.

## Details

Cell identities are prepared internally to preserve factor-level order
and handle missing values as `"Unknown"`. Colors are strictly matched to
the groups present in the plot.

If `figplot = TRUE`, a minimal publication-panel style is used and a
compact axis-arrow inset is added. If `dark = TRUE`, a dark theme is
applied after the main theme.

## Examples

``` r
if (FALSE) { # \dontrun{
cellmap(obj, group.by = "seurat_clusters")

cellmap(
  obj,
  group.by = c("gnrh_status", "gnrh_stage"),
  reduction = "umap",
  label = TRUE
)

cellmap(
  obj,
  group.by = "gnrh_stage",
  dark = TRUE,
  legend = TRUE
)
} # }
```
