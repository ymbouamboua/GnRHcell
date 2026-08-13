# Plot dataset-specific GnRH markers

Plot dataset-specific GnRH markers

## Usage

``` r
plot_gnrh_dataset_specific_markers(
  files,
  dir = ".",
  gene_col = "gene",
  score_col = NULL,
  dataset_names = names(files),
  gene_case = c("upper", "asis"),
  top_n_per_dataset = 20,
  max_datasets = 2,
  scale_rows = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  fontsize_legend = 10,
  filename = NULL,
  width = 8,
  height = 10
)
```

## Arguments

- files:

  Named character vector of GnRH marker TSV filenames or paths. Names
  are used as dataset labels unless `dataset_names` is supplied.

- dir:

  Directory prepended to relative file paths.

- gene_col:

  Gene-symbol column.

- score_col:

  Score column. `NULL` selects the first available column from
  `score_candidates` independently for each table.

- dataset_names:

  Optional dataset labels.

- gene_case:

  Standardize genes to uppercase or preserve their case.

- top_n_per_dataset:

  Maximum markers selected for each dataset.

- max_datasets:

  Maximum datasets supporting a dataset-specific marker.

- scale_rows:

  Row-standardize scores before plotting.

- cluster_rows, cluster_columns:

  Cluster heatmap rows or columns.

- show_column_names:

  Show dataset labels.

- fontsize_row, fontsize_col, fontsize_legend:

  Text sizes.

- filename:

  Optional PDF filename. Nothing is written when `NULL`.

- width, height:

  PDF dimensions.

## Value

A list with `heatmap`, matrices, specificity scores, and
`specific_table`.
