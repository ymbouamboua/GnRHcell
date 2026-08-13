# Merge GnRH marker scores across datasets

Builds a gene-by-dataset score matrix from marker tables returned by
[`gnrh_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md).
Marker tables may be supplied as data frames, as one list of data
frames, or as a named vector of TSV files.

## Usage

``` r
merge_gnrh_marker_scores(
  ...,
  dataset_names = NULL,
  dir = ".",
  gene_col = "gene",
  score_col = NULL,
  score_candidates = c("score", "avg_log2FC", "coexpr", "avg_logFC"),
  gene_case = c("upper", "asis"),
  fill = NA_real_
)
```

## Arguments

- ...:

  Marker data frames, one list of marker data frames, or a character
  vector of marker-table paths.

- dataset_names:

  Optional dataset labels.

- dir:

  Directory prepended to relative file paths.

- gene_col:

  Gene-symbol column.

- score_col:

  Score column. `NULL` selects the first available column from
  `score_candidates` independently for each table.

- score_candidates:

  Preferred GnRH marker-score columns.

- gene_case:

  Standardize genes to uppercase or preserve their case.

- fill:

  Value used for genes absent from a dataset. `NA_real_` is safest.

## Value

A numeric gene-by-dataset matrix. Attribute `score_columns` records the
score column selected for each dataset.
