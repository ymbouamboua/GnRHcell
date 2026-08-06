# Plot GnRH coexpression markers

Visualizes genes coexpressed with `GNRH1` and ranks them by coexpression
strength.

## Usage

``` r
plot_gnrh_coexpr(df, coexp_cutoff = 0.25, txtsize = 10)
```

## Arguments

- df:

  A marker table returned by
  [`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md).
  Must contain the columns `gene`, `coexpr`, `score`, and `p_val_adj`.

- coexp_cutoff:

  Minimum coexpression value required for plotting. Default is `0.25`.

- txtsize:

  Base text size passed to
  [`plot_theme`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md).
  Default is `12`.

## Value

A `ggplot2` object showing coexpression strength, marker score, and
adjusted p-value significance.

## Details

`GNRH1` itself is removed before plotting. Genes are ordered by
decreasing coexpression strength.

Point aesthetics:

- x-axis: correlation with `GNRH1`

- y-axis: coexpressed genes

- point size: composite marker score

- point color: \\-log10(adjusted~p-value)\\

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_gnrh_coexpr(
  markers,
  coexp_cutoff = 0.35
)

print(p)
} # }
```
