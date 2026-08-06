# Validate input Seurat object for GnRH analysis

Validates that the input Seurat object contains the required expression
layers for GnRH analysis.

## Usage

``` r
validate_input(object, assay = NULL, auto_normalize = FALSE, verbose = TRUE)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- assay:

  Assay name to validate. If `NULL`, the default assay of the Seurat
  object is used.

- auto_normalize:

  Logical; automatically normalize the assay if normalized data are
  missing. Default is `FALSE`.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

## Value

A validated Seurat object, optionally normalized if requested.

## Details

This function checks compatibility with both Seurat v5 (`Assay5`) and
earlier assay formats, ensuring that raw counts and normalized
expression data are available.

If normalized data are missing and `auto_normalize = TRUE`,
normalization is performed automatically using
[`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html).

Validation checks:

- presence of raw counts (`counts`)

- presence of normalized expression data (`data`)

- compatibility with Seurat v5 `Assay5` objects

- compatibility with legacy Seurat assay objects

If counts are missing, the function stops with an error.

If normalized data are missing:

- with `auto_normalize = FALSE`: error

- with `auto_normalize = TRUE`: normalization is performed

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`Seurat::NormalizeData`](https://satijalab.org/seurat/reference/NormalizeData.html)
