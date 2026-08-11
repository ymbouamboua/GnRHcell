# Plot GnRHcell runtime across datasets

Plot GnRHcell runtime across datasets

## Usage

``` r
plot_gnrh_runtime_curve(
  files = NULL,
  dir = file.path("results", "tables"),
  pattern = "_gnrh_run_info\\.tsv$",
  metric = "total_sec",
  x.ang = 45,
  txtsize = 10,
  show_points = TRUE
)
```

## Arguments

- files:

  Optional named vector of run-information files.

- dir:

  Directory searched when `files` is `NULL`.

- pattern:

  File-selection regular expression.

- metric:

  Runtime column to summarize.

- x.ang:

  X-axis label angle.

- txtsize:

  Base text size.

- show_points:

  Show dataset points.

## Value

A ggplot object.
