# Plot conserved GnRH markers across datasets

Plot conserved GnRH markers across datasets

## Usage

``` r
plot_gnrh_conserved_markers(
  files,
  dir = ".",
  gene_col = "gene",
  score_col = NULL,
  dataset_names = names(files),
  gene_case = c("upper", "asis"),
  min_datasets = 2,
  top_n = 50,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  fontsize_row = 8,
  fontsize_col = 9,
  fontsize_legend = 10,
  filename = NULL,
  width = 7,
  height = 9
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

- min_datasets:

  Minimum datasets supporting a gene.

- top_n:

  Maximum genes displayed.

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

A list with `heatmap`, `matrix`, `plot_matrix`, and `conserved`.
