# Stage GnRH lineage cells

Assign developmental states to GnRH-lineage cells using biologically
informed transcriptional programs.

## Usage

``` r
stage_gnrh(object, assay = "RNA", layer = "counts", verbose = TRUE)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- assay:

  Assay used for expression extraction. Default is `"RNA"`.

- layer:

  Expression layer used for staging. Default is `"counts"`.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

## Value

A Seurat object updated with:

- `gnrh_stage`:

  Assigned developmental stage labels.

- `gnrh_identity_score`:

  Identity module score.

- `gnrh_migrating_score`:

  Migration module score.

- `gnrh_mature_score`:

  Mature neuroendocrine score.

- `gnrh_secreting_score`:

  Secretory activity score.

Developmental module definitions are stored in:
`object@misc$gnrh_stage_modules`

## Details

Cells are scored against predefined developmental modules representing
major GnRH neuron states:

- `identity`: lineage specification and early GnRH identity

- `migrating`: migration and axon-guidance programs

- `mature`: neuroendocrine maturation

- `secreting`: secretory and vesicle machinery activation

Each cell is assigned to the stage with the highest module score. If
GnRH classification metadata are present, non-GnRH cells are labeled as
`non-gnrh`.

Stage-specific module scores are stored in object metadata.

Developmental staging is based on predefined marker modules reflecting
known biological programs of GnRH neuron development.

If `gnrh_class` metadata are present from
[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
cells classified as negative are reassigned to `non-gnrh`.

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
