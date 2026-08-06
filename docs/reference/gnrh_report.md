# Generate a comprehensive GnRH diagnostic report

Creates a multi-panel quality-control dashboard summarizing GnRH
detection performance, classification behavior, threshold optimization,
ROC performance, module enrichment, developmental staging, and run
parameters from a `GnRHcell` analysis.

## Usage

``` r
gnrh_report(object, style = "test", verbose = TRUE)
```

## Arguments

- object:

  A `Seurat` object processed with
  [`detect_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md)
  and
  [`gnrh_diagnostics()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md).

- style:

  Character scalar specifying the plotting theme style passed to
  [`plot_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `"test"`.

- verbose:

  Logical; if `TRUE`, progress messages are printed.

## Value

A `patchwork` plot object containing the diagnostic dashboard.

## Details

The report automatically extracts diagnostics stored in
`object@misc$gnrh` after running the detection pipeline.

Included panels may include:

- Signal landscape: GNRH1 expression versus composite GnRH score

- Score distribution by GnRH classification

- Threshold optimization curve (sensitivity, specificity, F1)

- ROC curve with AUC estimate

- Module hit burden versus composite score

- Classification enrichment across total module hits

- Developmental stage signal landscape (if staging available)

- Run parameter summary table

Status and stage legends include both cell counts and percentages for
easier interpretation.

ROC curves are generated automatically when valid binary labels are
available in diagnostic metadata.

Diagnostic information is expected in:

- `object@misc$gnrh$diagnostics`

- `object@misc$gnrh$threshold_curve`

- optional staging metadata in `object$gnrh_stage`

- optional run parameters in `object@misc$gnrh_params`

If staging results are unavailable, the stage panel is replaced by a
placeholder panel.

If ROC computation is not possible (for example, only one class
present), the ROC panel is replaced by an informative placeholder.

## See also

[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md),
[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md),
[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md),
[`plot_gnrh_distribution`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_distribution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(hpsc_gnrh)

hpsc_gnrh <- run_gnrh(hpsc_gnrh)

p <- gnrh_report(hpsc_gnrh)
p
} # }
```
