# Plot GnRHcell metadata on cell embeddings

Visualizes GnRHcell metadata on a Seurat dimensional reduction.

## Usage

``` r
plot_gnrh_embedding(
  object,
  group.by = "all",
  reduction = "umap",
  cols = NULL,
  style = "classic",
  dark = FALSE,
  ncol = NULL
)
```

## Arguments

- object:

  A Seurat object processed by GnRHcell.

- group.by:

  Metadata variable to plot. Use `"all"` to plot `gnrh_status`,
  `gnrh_confident`, and `gnrh_stage`.

- reduction:

  Dimensional reduction to use.

- cols:

  Optional color palette.

- style:

  Theme style.

- dark:

  Logical; use dark theme.

- ncol:

  Number of columns when plotting multiple panels.

## Value

A ggplot2 or patchwork object.
