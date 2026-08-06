# GnRH gene similarity network

Constructs a similarity network based on coexpression structure among
top-ranked genes.

## Usage

``` r
plot_network(df, top_n = 25, threshold = 0.4)
```

## Arguments

- df:

  Data frame containing gene-level scores and coexpression.

- top_n:

  Number of genes to include.

- threshold:

  Minimum edge weight to retain.

## Value

A ggraph object.
