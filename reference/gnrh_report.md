# Generate a publication-ready GnRHcell diagnostic report

Creates a six-panel quality-control dashboard summarizing detection,
developmental staging, score separation, threshold behavior, detection
classes or ROC discrimination, and marker-program support.

## Usage

``` r
gnrh_report(
  object,
  style = c("test", "classic", "minimal", "bw"),
  mode = c("light", "dark"),
  truth = NULL,
  roc_mode = c("auto", "external", "internal", "none"),
  positive_truth = c("1", "TRUE", "true", "positive", "Positive", "pos", "Pos"),
  txtsize = 10,
  max_points = 100000L,
  seed = 1234L,
  show_class_panel = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object processed with
  [`run_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
  and
  [`gnrh_diagnostics()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md).

- style:

  Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.

- mode:

  Display mode: `"light"` or `"dark"`.

- truth:

  Optional metadata column containing an independent binary reference
  classification for external ROC analysis.

- roc_mode:

  ROC behavior: `"auto"` uses external truth when supplied and otherwise
  restores internal status/stage discrimination curves; `"external"`,
  `"internal"`, and `"none"` force a specific behavior.

- positive_truth:

  Values in `truth` interpreted as positive.

- txtsize:

  Base text size.

- max_points:

  Maximum cells displayed in each scatter panel. All rare GnRH-positive
  cells are retained before negative cells are sampled.

- seed:

  Random seed used for display-only subsampling.

- show_class_panel:

  Show detection-class composition when ROC analysis is disabled or
  unavailable.

- verbose:

  Print progress messages.

## Value

A patchwork object.

## Examples

``` r
if (FALSE) { # \dontrun{
report <- gnrh_report(wang, style = "bw")
report

# External ROC using an independent manual/reference annotation
report <- gnrh_report(
  wang,
  truth = "manual_gnrh",
  roc_mode = "external",
  positive_truth = c("GnRH", "pos", "TRUE"),
  style = "minimal"
)

# Explicit internal score-discrimination curves
report <- gnrh_report(wang, roc_mode = "internal", style = "bw")

ggplot2::ggsave(
  "gnrh_report.pdf", report,
  width = 12, height = 7.5, units = "in",
  device = grDevices::cairo_pdf
)
} # }
```
