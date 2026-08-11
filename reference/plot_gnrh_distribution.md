# Plot GnRHcell metadata distributions

Plot GnRHcell metadata distributions

## Usage

``` r
plot_gnrh_distribution(
  object,
  group.by,
  split.by = NULL,
  cols = NULL,
  proportion = FALSE,
  position = "stack",
  label = TRUE,
  label.size = 3,
  plot.ttl = NULL,
  txtsize = 10,
  x.ang = 45,
  flip = FALSE
)
```

## Arguments

- object:

  A Seurat object.

- group.by:

  Metadata column defining categories.

- split.by:

  Optional metadata column defining bars.

- cols:

  Optional named colour vector.

- proportion:

  Display within-split proportions instead of counts.

- position:

  Bar position, such as `"stack"` or `"dodge"`.

- label:

  Add value labels.

- label.size:

  Label text size.

- plot.ttl:

  Optional plot title.

- txtsize:

  Base text size.

- x.ang:

  X-axis label angle.

- flip:

  Flip coordinates.

## Value

A ggplot object.
