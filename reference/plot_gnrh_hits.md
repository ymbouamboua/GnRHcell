# Plot GnRH module hit distributions

Bar plots of GnRH module hits vs metadata variables

## Usage

``` r
plot_gnrh_hits(
  data,
  x,
  fill,
  palette = NULL,
  type = c("count", "fraction"),
  title = NULL,
  rotate.x = FALSE,
  txtsize = 12,
  style = "classic",
  ...
)
```

## Arguments

- data:

  data.frame (usually object@meta.data)

- x:

  character. grouping variable (e.g. "total_hits_bin")

- fill:

  character. fill variable (e.g. "gnrh_status")

- palette:

  named vector of colors

- type:

  "count" or "fraction"

- title:

  plot title

- rotate.x:

  logical

- txtsize:

  theme text size

- style:

  theme style

- ...:

  Additional arguments passed to
  [`plot_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
