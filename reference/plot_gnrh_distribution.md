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
  x.ang = NULL,
  style = "classic",
  flip = FALSE,
  adaptive = TRUE,
  bar.width = NULL,
  bar.gap = NULL,
  legend.position = "auto"
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

  X-axis label angle. When `NULL`, the angle is selected from

- style:

  Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`. the number
  and length of sample labels.

- flip:

  Flip coordinates.

- adaptive:

  Automatically use compact spacing, an economical legend layout, and
  recommended export dimensions based on the number of samples.

- bar.width:

  Bar width. When `NULL`, it is selected automatically.

- bar.gap:

  Gap between adjacent bar edges, in x-axis units. It is independent of
  `bar.width`; for example, `bar.width = 0.3` and `bar.gap = 0.2`
  produce centres 0.5 units apart. When `NULL`, an adaptive value is
  used.

- legend.position:

  Legend position. Use `"auto"` to place it according to the number of
  samples, or a standard ggplot2 legend position.

## Value

A ggplot object.
