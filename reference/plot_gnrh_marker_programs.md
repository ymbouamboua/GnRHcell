# Plot GnRH marker program results

Plot GnRH marker program results

## Usage

``` r
plot_gnrh_marker_programs(
  programs,
  table = c("summary", "high_confidence", "candidate_table"),
  type = c("bar", "dot", "tile"),
  min_genes = 1,
  mode = "light",
  txtsize = 12,
  x.ang = 45,
  style = c("bw", "test", "classic", "minimal")
)
```

## Arguments

- programs:

  Result returned by
  [`gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md).

- table:

  Result table to visualize.

- type:

  Plot type: bar, dot, or tile.

- min_genes:

  Minimum genes retained for summary plots.

- mode:

  Light or dark display mode.

- txtsize:

  Base text size.

- x.ang:

  X-axis label angle.

- style:

  Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.

## Value

A ggplot object.
