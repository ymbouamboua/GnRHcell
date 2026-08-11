# Plot GnRH module hit distributions

Plot GnRH module hit distributions

## Usage

``` r
plot_gnrh_hits(
  data,
  x,
  fill,
  palette = NULL,
  type = c("count", "fraction"),
  title = NULL,
  txtsize = 10
)
```

## Arguments

- data:

  Data frame containing plotting columns.

- x:

  Column mapped to the x axis.

- fill:

  Column mapped to fill colour.

- palette:

  Optional named colour vector.

- type:

  Display counts or fractions.

- title:

  Optional plot title.

- txtsize:

  Base text size.

## Value

A ggplot object.
