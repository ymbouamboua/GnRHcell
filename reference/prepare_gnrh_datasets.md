# Prepare a GnRHcell dataset configuration table

Validates and standardizes a dataset configuration table before running
multi-dataset GnRHcell analyses.

## Usage

``` r
prepare_gnrh_datasets(datasets, check_files = TRUE, remove_missing = FALSE)
```

## Arguments

- datasets:

  A data frame or tibble containing dataset configuration information.

- check_files:

  Logical. Whether to warn when dataset files are missing. Default is
  `TRUE`.

- remove_missing:

  Logical. Whether rows corresponding to missing files should be
  removed. Default is `FALSE`.

## Value

A tibble containing the standardized dataset configuration and an
additional logical column named `exists`.
