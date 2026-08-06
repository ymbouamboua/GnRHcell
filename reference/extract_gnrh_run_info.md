# Extract GnRHcell run information

Extracts runtime, dataset size, and detection summary statistics from a
Seurat object processed with
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md).

## Usage

``` r
extract_gnrh_run_info(object, dataset_name = NULL)
```

## Arguments

- object:

  A Seurat object containing `object@misc$gnrh$run_info`.

- dataset_name:

  Optional dataset name. If `NULL`, `"dataset"` is used.

## Value

A data frame containing dataset size, GnRH detection counts, runtime
values in seconds, and formatted runtime strings.

## See also

[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md),
[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md)
