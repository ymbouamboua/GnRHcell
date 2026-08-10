# Run complete GnRHcell analysis pipeline

Executes the full GnRHcell workflow for identification, developmental
staging, diagnostics, and runtime reporting.

## Usage

``` r
run_gnrh(
  object,
  detect = TRUE,
  stage = TRUE,
  diagnostics = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- detect:

  Logical; run GnRH detection. Default is `TRUE`.

- stage:

  Logical; run developmental staging. Default is `TRUE`.

- diagnostics:

  Logical; run diagnostic analysis. Default is `TRUE`.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md).

## Value

A Seurat object updated with GnRH detection, developmental staging,
diagnostics, and runtime metadata.

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md),
[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md),
[`extract_gnrh_run_info`](https://ymbouamboua.github.io/GnRHcell/reference/extract_gnrh_run_info.md)
