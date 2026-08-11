# Plot GnRH coexpression markers

Plot GnRH coexpression markers

## Usage

``` r
plot_gnrh_coexpr(df, coexp_cutoff = 0.25, txtsize = 10, style = "bw")
```

## Arguments

- df:

  Marker table containing `gene`, `coexpr`, `score`, and `p_val_adj`.

- coexp_cutoff:

  Minimum GNRH1 co-expression value.

- txtsize:

  Base text size.

- style:

  Theme style.

## Value

A ggplot object.
