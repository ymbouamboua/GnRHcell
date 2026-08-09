# Run GnRHcell across multiple datasets

Applies the GnRHcell workflow to a collection of Seurat datasets and
optionally performs cross-dataset marker comparisons and conserved
marker program analysis.

## Usage

``` r
run_gnrh_collection(
  datasets,
  output_dir,
  run_markers = TRUE,
  run_comparisons = TRUE,
  run_programs = TRUE,
  clean_objects = TRUE,
  save_objects = FALSE,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A data frame or tibble containing the columns `id`, `label`,
  `species`, `file`, `split_by`, and `reduction`.

- output_dir:

  Character scalar giving the root output directory.

- run_markers:

  Logical. Whether to identify GnRH markers for each dataset. Default is
  `TRUE`.

- run_comparisons:

  Logical. Whether to perform cross-dataset comparisons after all
  datasets have been processed. Default is `TRUE`.

- run_programs:

  Logical. Whether to identify conserved marker programs during
  cross-dataset comparisons. Default is `TRUE`.

- clean_objects:

  Logical. Whether processed Seurat objects should be removed from the
  returned dataset-level results to reduce memory usage. Default is
  `TRUE`.

- save_objects:

  Logical. Whether processed Seurat objects should be saved to disk.
  Default is `FALSE`.

- verbose:

  Logical. Whether to print progress messages. Default is `TRUE`.

## Value

An object of class `"gnrh_collection"` containing:

- datasets:

  The dataset configuration table used for the analysis.

- results:

  Named list of dataset-level GnRHcell results.

- comparisons:

  Cross-dataset comparison results, or `NULL`.

- output_dir:

  Normalized output directory.

## Details

Each row of `datasets` represents one dataset. Dataset files are loaded
sequentially to limit memory usage.

## See also

[`prepare_gnrh_datasets()`](https://ymbouamboua.github.io/GnRHcell/reference/prepare_gnrh_datasets.md),
[`run_gnrh_dataset()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh_dataset.md),
[`compare_gnrh_datasets()`](https://ymbouamboua.github.io/GnRHcell/reference/compare_gnrh_datasets.md)

## Examples

``` r
if (FALSE) { # \dontrun{
datasets <- data.frame(
  id = c("human_hpsc", "mouse_hypomap"),
  label = c("Human hPSC GnRH", "Mouse HypoMap"),
  species = c("Human", "Mouse"),
  file = c("human_hpsc.rds", "mouse_hypomap.rds"),
  split_by = c("orig.ident", "orig.ident"),
  reduction = c("umap", "umap")
)

datasets <- prepare_gnrh_datasets(datasets)

results <- run_gnrh_collection(
  datasets = datasets,
  output_dir = "gnrh_results"
)
} # }
```
