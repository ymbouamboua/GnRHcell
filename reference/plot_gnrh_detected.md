# Plot detected GnRH-positive cells across datasets

Plot detected GnRH-positive cells across datasets

## Usage

``` r
plot_gnrh_detected(
  files = NULL,
  dir = file.path("results", "tables"),
  pattern = "_gnrh_run_info\\.tsv$",
  x.ang = 45,
  txtsize = 10,
  debug = FALSE
)
```

## Arguments

- files:

  Optional named vector of run-information files.

- dir:

  Directory searched when `files` is `NULL`.

- pattern:

  File-selection regular expression.

- x.ang:

  X-axis label angle.

- txtsize:

  Base text size.

- debug:

  Print the imported summary columns.

## Value

A ggplot object.
