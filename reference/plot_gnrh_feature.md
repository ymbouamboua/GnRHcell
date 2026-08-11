# Cell Feature Plot for Seurat Objects

A flexible wrapper around Seurat's `FeaturePlot` to visualize gene
expression or metadata features in a Seurat object. Supports custom
color palettes, viridis, and hotspot/rainbow palettes.

## Usage

``` r
plot_gnrh_feature(
  object,
  features = NULL,
  preset = NULL,
  cols = NULL,
  theme.cols = "gnrh",
  rev.cols = FALSE,
  na.col = "lightgray",
  order = FALSE,
  pt.size = NULL,
  txtsize = 10,
  reduction = NULL,
  na.cutoff = 1e-09,
  raster = NULL,
  raster.dpi = c(512, 512),
  split.by = NULL,
  ncol = NULL,
  layer = "data",
  label = FALSE,
  axes = TRUE,
  combine = TRUE,
  blend = FALSE,
  merge.leg = FALSE,
  style = "classic",
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- features:

  Character vector of features (genes or metadata columns) to plot.

- preset:

  Optional GnRHcell feature preset: `"core"`, `"modules"`, `"staging"`,
  or `"all"`.

- cols:

  Optional character vector of colors for plotting.

- theme.cols:

  Character. Predefined theme color palette (default: `"Reds"`). Options
  include `"Reds"`, `"Blues"`, etc., or custom list palettes
  `"hotspot"`,`"rainbow"`, etc.

- rev.cols:

  Logical. Reverse the color palette (default: `FALSE`).

- na.col:

  Color for NA or below-cutoff expression values (default:
  `"lightgray"`).

- order:

  Logical. If TRUE, plot high-expression cells on top (default:
  `FALSE`).

- pt.size:

  Numeric. Point size. If NULL, automatically calculated based on number
  of cells.

- txtsize:

  Numeric. Base font size for plot titles and axis labels (default: 10).

- reduction:

  Character. Dimensional reduction to use (default: first available in
  Seurat object).

- na.cutoff:

  Numeric. Minimum expression value for coloring; below this will be NA
  if palette requires (default: 1e-9).

- raster:

  Logical. If TRUE, rasterize points for faster plotting of large
  datasets.

- raster.dpi:

  Numeric vector of length 2. DPI for rasterization (default:
  c(512,512)).

- split.by:

  Character. Metadata column to split the plot.

- ncol:

  Numeric. Number of columns when combining multiple plots.

- layer:

  Character. Seurat assay slot to fetch data from (default: `"data"`).

- label:

  Logical. Whether to label clusters (default: FALSE).

- axes:

  Logical. Whether to show axes (default: TRUE).

- combine:

  Logical. Whether to return a single combined plot (default: TRUE).

- blend:

  Logical. Whether to blend exactly two features (default: FALSE).

- merge.leg:

  Logical. Whether to merge multiple legends into one (default: FALSE).

- style:

  Character. ggplot2 theme to apply (default: `"classic"`).

- ...:

  Additional arguments passed to
  [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot.html).

## Value

A `ggplot` object (or `patchwork` object if multiple features).

## Examples

``` r
if (FALSE) { # \dontrun{
path <- system.file("extdata", "hpsc.rds", package = "GnRHcell")
object <- readRDS(path)
plot_gnrh_feature(object, features = "GNRH1")
plot_gnrh_feature(object, preset = "core", ncol = 2)
} # }
```
