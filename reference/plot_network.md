# GnRH gene similarity network

GnRH gene similarity network

## Usage

``` r
plot_network(df, top_n = 25, threshold = 0.4)
```

## Arguments

- df:

  Marker table containing `gene`, `coexpr`, and `score`.

- top_n:

  Maximum number of ranked genes included.

- threshold:

  Minimum co-expression edge weight.

## Value

A ggraph object.
