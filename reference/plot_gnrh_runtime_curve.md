# Plot GnRHcell runtime across datasets

Plots runtime across datasets as a line plot. The x-axis displays each
dataset together with the number of detected GnRH-positive cells over
the total number of cells.

## Usage

``` r
plot_gnrh_runtime_curve(
  files = NULL,
  dir = file.path("results", "tables"),
  pattern = "_gnrh_run_info\\.tsv$",
  metric = "total_sec",
  style = "classic",
  x.ang = 45,
  txtsize = 10,
  show_points = TRUE
)
```

## Arguments

- files:

  Optional character vector of GnRH run-info TSV files.

- dir:

  Directory containing run-info tables.

- pattern:

  Regex pattern used to find run-info files.

- metric:

  Runtime column to plot, e.g. `"total_sec"`.

- style:

  Plot theme style passed to
  [`plot_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).

- x.ang:

  X-axis label angle.

- txtsize:

  Base text size.

- show_points:

  Logical; show points on the curve.

## Value

A `ggplot2` object.
