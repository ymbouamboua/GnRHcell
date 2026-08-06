# Plot GnRH feature maps

Wrapper around
[`cellfeature`](https://ymbouamboua.github.io/GnRHcell/reference/cellfeature.md)
for GnRHcell diagnostic features.

## Usage

``` r
plot_gnrh_feature(
  object,
  features = "all",
  feature_type = NULL,
  reduction = "umap",
  ncol = 4,
  style = "classic",
  merge.leg = FALSE,
  split.by = NULL,
  title = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object processed by GnRHcell.

- features:

  Features to plot, or `"all"`.

- feature_type:

  Optional preset: `"core"`, `"modules"`, or `"all"`.

- reduction:

  Dimensional reduction to use.

- ncol:

  Number of columns.

- style:

  Theme style.

- merge.leg:

  Logical; merge legends across panels.

- split.by:

  Optional metadata column for splitting.

- title:

  Optional title.

- ...:

  Additional arguments passed to
  [`cellfeature`](https://ymbouamboua.github.io/GnRHcell/reference/cellfeature.md).

## Value

ggplot2 or patchwork object.
