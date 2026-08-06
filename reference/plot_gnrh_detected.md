# Plot detected GnRH-positive cells across datasets

Plot detected GnRH-positive cells across datasets

## Usage

``` r
plot_gnrh_detected(
  files = NULL,
  dir = file.path("results", "tables"),
  pattern = "_gnrh_run_info\\.tsv$",
  style = "classic",
  x.ang = 45,
  txtsize = 10,
  debug = TRUE
)
```

## Arguments

- files:

  Optional character vector of run-info TSV files.

- dir:

  Directory containing runtime summary tables.

- pattern:

  Regex pattern used to identify runtime files.

- style:

  Plot style.

- x.ang:

  X-axis text angle.

- txtsize:

  Base text size.

- debug:

  Print loaded table.

## Value

A ggplot2 object.
