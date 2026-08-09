# Run GnRHcell analysis on a single dataset

Runs the GnRHcell workflow on a Seurat object and generates
dataset-level diagnostic plots, embeddings, feature plots, distribution
plots, marker tables, co-expression plots, and marker networks.

## Usage

``` r
run_gnrh_dataset(
  object,
  dataset_id,
  dataset_label,
  split_by = "orig.ident",
  reduction = "umap",
  output_dir,
  run_markers = TRUE,
  clean_object = TRUE
)
```

## Arguments

- object:

  A Seurat object containing the dataset to analyze.

- dataset_id:

  Character scalar giving a short unique identifier for the dataset.
  This identifier is used in output file names and directories.

- dataset_label:

  Character scalar giving a human-readable dataset name.

- split_by:

  Character scalar giving the metadata column used to split GnRH
  distributions. If unavailable, common alternatives such as
  `"orig.ident"`, `"sample"`, and `"library_id"` are considered.

- reduction:

  Character scalar giving the dimensional reduction used for embedding
  plots. Default is `"umap"`.

- output_dir:

  Character scalar giving the root output directory.

- run_markers:

  Logical. Whether to identify GnRH-associated markers and generate
  marker-based plots. Default is `TRUE`.

- clean_object:

  Logical. Whether to trigger garbage collection after processing the
  dataset. Default is `TRUE`.

## Value

A named list containing:

- object:

  The processed Seurat object.

- run_info:

  Dataset-level GnRHcell run information.

- markers:

  GnRH marker results, or `NULL` when marker analysis was disabled.

- reduction:

  The dimensional reduction used for plotting.

- split_by:

  The metadata column used for distribution plots, or `NULL`.

## Details

The function is primarily used internally by
[`run_gnrh_collection()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh_collection.md)
but can also be called directly for individual datasets.

## See also

[`run_gnrh_collection()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh_collection.md),
[`run_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
