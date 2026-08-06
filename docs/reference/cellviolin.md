# Create violin plots for Seurat features

Generates publication-ready violin plots for one or more genes or
metadata features from a Seurat object. This function wraps
[`Seurat::VlnPlot()`](https://satijalab.org/seurat/reference/VlnPlot.html)
and adds shared y-axis scaling, median markers, optional statistical
tests, flexible layout, global axis labels, and GnRHcell theme support.

## Usage

``` r
cellviolin(
  obj,
  features,
  ncol = NULL,
  shared.y = FALSE,
  ttl.pos = c("left", "center", "right"),
  group.by = "seurat_clusters",
  split.by = NULL,
  stack = FALSE,
  assay = "RNA",
  slot = "data",
  log = FALSE,
  cols = NULL,
  med = FALSE,
  med.size = 1,
  pt.size = 0,
  border.size = 0.1,
  style = "classic",
  leg.pos = "none",
  x.ang = 45,
  title = NULL,
  rm.subttl = FALSE,
  flip = FALSE,
  auto.resize = TRUE,
  ylab.global = NULL,
  xlab.global = NULL,
  add.stats = FALSE,
  show.pval = FALSE,
  pairwise = FALSE,
  pval.label = "p.signif",
  p.display = c("stars", "value", "none"),
  txtsize = 10,
  subttl.size = 10,
  ...
)
```

## Arguments

- obj:

  A Seurat object.

- features:

  Character vector of genes or metadata columns to plot.

- ncol:

  Number of columns in the patchwork layout. If `NULL`, selected
  automatically from the number of features.

- shared.y:

  Logical; use the same y-axis range for all features. Default is
  `FALSE`.

- ttl.pos:

  Subplot title position. One of `"left"`, `"center"`, or `"right"`.

- group.by:

  Metadata column used to group cells. Default is `"seurat_clusters"`.

- split.by:

  Optional metadata column used to split violin plots.

- stack:

  Logical; arrange plots in a single column. Default is `FALSE`.

- assay:

  Assay used to retrieve gene expression. Default is `"RNA"`.

- slot:

  Assay slot used for expression values. One of `"data"`, `"counts"`, or
  `"scale.data"`. Default is `"data"`.

- log:

  Logical; if `TRUE` and `slot = "counts"`, use log-transformed counts.
  Default is `FALSE`.

- cols:

  Optional named color vector for groups.

- med:

  Logical; overlay median points on violin plots. Default is `FALSE`.

- med.size:

  Size of median points. Default is `1`.

- pt.size:

  Size of jittered cells. Use `0` to hide points. Default is `0`.

- border.size:

  Violin outline linewidth. Default is `0.1`.

- style:

  Theme style passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `"classic"`.

- leg.pos:

  Legend position. Default is `"none"`.

- x.ang:

  X-axis text angle. Default is `45`.

- title:

  Optional global title for the combined plot.

- rm.subttl:

  Logical; remove individual feature titles. Default is `FALSE`.

- flip:

  Logical; flip x and y axes. Default is `FALSE`.

- auto.resize:

  Logical; add dynamic width and height attributes to the returned
  object. Default is `TRUE`.

- ylab.global:

  Optional global y-axis label. If `NULL`, inferred from `slot`.

- xlab.global:

  Optional global x-axis label. Default is blank.

- add.stats:

  Logical; compute statistical comparisons. Default is `FALSE`.

- show.pval:

  Logical; display p-values or significance labels. Default is `FALSE`.

- pairwise:

  Logical; perform pairwise Wilcoxon tests between groups. If `FALSE`, a
  global Kruskal-Wallis test is used. Default is `FALSE`.

- pval.label:

  Label type passed to statistical annotation utilities. Default is
  `"p.signif"`.

- p.display:

  How to display statistical results. One of `"stars"`, `"value"`, or
  `"none"`. Default is `"stars"`.

- txtsize:

  Base text size. Default is `10`.

- subttl.size:

  Feature subplot title size. Default is `10`.

- ...:

  Additional arguments passed to
  [`Seurat::VlnPlot()`](https://satijalab.org/seurat/reference/VlnPlot.html)
  and
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).

## Value

A patchwork/cowplot object containing the violin plots.

## Details

Features can be either assay genes or metadata columns. Missing features
are silently ignored; an error is raised if none of the requested
features are found.

When `add.stats = TRUE` and `show.pval = TRUE`, statistical annotations
are added if at least two groups are available. Pairwise tests use
Wilcoxon rank-sum tests, while the global comparison uses a
Kruskal-Wallis test.

## Examples

``` r
if (FALSE) { # \dontrun{
cellviolin(
  obj,
  features = c("GNRH1", "ISL1", "DLX5"),
  group.by = "gnrh_stage",
  pt.size = 0.1,
  ncol = 3
)

cellviolin(
  obj,
  features = c("GNRH1", "ISL1"),
  group.by = "gnrh_status",
  add.stats = TRUE,
  show.pval = TRUE,
  pairwise = TRUE
)
} # }
```
