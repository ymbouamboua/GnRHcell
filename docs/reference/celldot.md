# Create an enhanced Seurat dot plot

Generates a customizable dot plot for gene expression patterns across
clusters or metadata-defined groups. This function wraps
[`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html)
and adds improved color control, optional dot outlines, flexible axis
formatting, legend customization, and GnRHcell theme support.

## Usage

``` r
celldot(
  object,
  features,
  group.by = "seurat_clusters",
  th.cols = "RdYlBu",
  rev.th.cols = TRUE,
  dot.scale = 4,
  x.ang = 90,
  vjust.x = NULL,
  hjust.x = NULL,
  flip = FALSE,
  txtsize = 12,
  title = NULL,
  leg.size = 10,
  leg.ttl.size = 10,
  leg.pos = "right",
  leg.just = "bottom",
  leg.hjust = FALSE,
  x.axis.pos = "bottom",
  style = "classic",
  x.face = FALSE,
  y.face = FALSE,
  x.ttl = FALSE,
  y.ttl = FALSE,
  dot.outline = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- features:

  Character vector or named list of features to plot.

- group.by:

  Metadata column used to group cells. Default is `"seurat_clusters"`.

- th.cols:

  RColorBrewer palette name used for the expression gradient. Default is
  `"RdYlBu"`.

- rev.th.cols:

  Logical; reverse the color gradient. Default is `TRUE`.

- dot.scale:

  Numeric dot-size scaling factor passed to
  [`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html).
  Default is `4`.

- x.ang:

  Angle of x-axis labels in degrees. Default is `90`.

- vjust.x, hjust.x:

  Vertical and horizontal justification for x-axis labels.

- flip:

  Logical; flip x and y axes using
  [`ggplot2::coord_flip()`](https://ggplot2.tidyverse.org/reference/coord_flip.html).
  Default is `FALSE`.

- txtsize:

  Base text size passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `12`.

- title:

  Optional plot title.

- leg.size:

  Legend text size. Default is `10`.

- leg.ttl.size:

  Legend title size. Default is `10`.

- leg.pos:

  Legend position. One of `"right"`, `"left"`, `"top"`, `"bottom"`, or
  `NULL`. Default is `"right"`.

- leg.just:

  Legend justification. Default is `"bottom"`.

- leg.hjust:

  Logical; reserved for legend layout customization. Default is `FALSE`.

- x.axis.pos:

  Position of the x-axis. Default is `"bottom"`.

- style:

  Theme style passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `"classic"`.

- x.face, y.face:

  Logical; italicize x- or y-axis labels. Default is `FALSE`.

- x.ttl, y.ttl:

  Logical; show x- or y-axis titles. Default is `FALSE`.

- dot.outline:

  Logical; draw outlines around dots. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html)
  and
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).

## Value

A `ggplot2` object.

## Details

Dot color represents average expression and dot size represents the
percentage of cells expressing each feature, following the standard
[`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html)
convention.

When `dot.outline = TRUE`, dots are drawn with a light outline using
shape 21. This can improve readability in publication figures.

## Examples

``` r
if (FALSE) { # \dontrun{
celldot(pbmc, features = c("MS4A1", "CD3D"))

celldot(
  pbmc,
  features = c("MS4A1", "CD14"),
  th.cols = "Blues",
  dot.outline = TRUE,
  flip = TRUE
)
} # }
```
