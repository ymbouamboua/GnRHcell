# Plot metadata distributions

Publication-ready barplot utility for Seurat metadata.

## Usage

``` r
plot_gnrh_distribution(
  object,
  group.by,
  split.by = NULL,
  cols = NULL,
  sort = FALSE,
  decreasing = TRUE,
  proportion = FALSE,
  position = "stack",
  label = TRUE,
  label.size = 3,
  border = TRUE,
  border.col = "black",
  border.size = 0.2,
  width = 0.7,
  x.lab = NULL,
  y.lab = NULL,
  plot.ttl = NULL,
  txtsize = 10,
  x.ang = 45,
  style = "test",
  flip = FALSE,
  ...
)
```

## Arguments

- object:

  Seurat object.

- group.by:

  Metadata column to plot.

- split.by:

  Optional metadata column for grouped plots (e.g. sample, batch,
  condition).

- cols:

  Named color vector.

- sort:

  Sort bars.

- decreasing:

  Sort decreasing.

- proportion:

  Plot proportions instead of counts.

- position:

  Bar position: `"stack"` or `"dodge"`.

- label:

  Add labels.

- label.size:

  Label text size.

- border:

  Draw borders.

- border.col:

  Border color.

- border.size:

  Border linewidth.

- width:

  Bar width.

- x.lab:

  X axis label.

- y.lab:

  Y axis label.

- plot.ttl:

  Plot title.

- txtsize:

  Base text size.

- x.ang:

  X axis angle.

- style:

  Theme style.

- flip:

  Flip coordinates.

- ...:

  Additional arguments passed to
  [`plot_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).

## Value

ggplot object.

## Details

Supports:

- Counts or proportions

- Distribution by sample/group

- Stacked or dodged bars

- Horizontal plots
