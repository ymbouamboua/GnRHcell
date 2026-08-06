# Detect positive GnRH lineage markers

Detect positive GnRH lineage markers

## Usage

``` r
gnrh_markers(
  object,
  group.by = "gnrh_status",
  ident.1 = "pos",
  ident.2 = NULL,
  assay = NULL,
  layer = "data",
  methods = "wilcox",
  coexpr.method = "spearman",
  min_pct = 0.01,
  min_fc = 0.25,
  max_padj = 0.05,
  coexpr_min = 0.15,
  min_detect = 3,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Seurat object

- group.by:

  Metadata column

- ident.1:

  Positive group

- ident.2:

  Optional comparison group

- assay:

  Assay

- layer:

  Layer

- methods:

  DE methods

- coexpr.method:

  Correlation method

- min_pct:

  Minimum detection fraction

- min_fc:

  Minimum log2FC

- max_padj:

  Maximum adjusted p-value

- coexpr_min:

  Minimum coexpression

- min_detect:

  Minimum detected cells

- verbose:

  Print progress

- ...:

  Additional arguments passed to
  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html).

## Value

Ranked marker table
