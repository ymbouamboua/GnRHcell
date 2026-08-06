# Detect GnRH neurons from single-cell RNA-seq data

Identifies candidate gonadotropin-releasing hormone (GnRH) neurons in a
Seurat object using a biologically informed multi-signal framework.

## Usage

``` r
detect_gnrh(
  object,
  assay = "RNA",
  layer = "counts",
  reduction = "pca",
  dims = 1:20,
  k = 20,
  min_umi = 2,
  min_counts = 500,
  mad_factor = 2,
  score_q = 0.9,
  scale_factor = 10000,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- assay:

  Assay used for expression extraction. Default is `"RNA"`.

- layer:

  Expression layer used for detection. Default is `"counts"`.

- reduction:

  Dimensional reduction used for kNN neighborhood support. Default is
  `"pca"`.

- dims:

  Dimensions used for neighborhood analysis. Default is `1:20`.

- k:

  Number of nearest neighbors used for neighborhood support. Default is
  20.

- min_umi:

  Minimum raw `GNRH1` UMI count required for detection. Default is 2.

- min_counts:

  Minimum total UMI count required per cell. Default is 500.

- mad_factor:

  Multiplier applied to MAD-based adaptive expression threshold. Default
  is 2.

- score_q:

  Quantile used to define adaptive composite score threshold. Default is
  0.9.

- scale_factor:

  Library normalization scale factor. Default is 10000.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

## Value

A Seurat object updated with GnRH detection metadata, diagnostics, and
stored detection parameters.

## Details

Detection integrates:

- direct `GNRH1` expression

- GnRH-associated marker module scoring

- migration marker support

- neuroendocrine marker support

- ambient RNA correction

- neighborhood enrichment using k-nearest neighbors

- adaptive thresholding of a composite detection score

The method combines transcript abundance, marker co-detection, local
transcriptomic neighborhood structure, and contamination-aware scoring
to improve detection of rare GnRH neurons in sparse single-cell
datasets.

Results are written into object metadata and diagnostics.

Added metadata columns include:

- `gnrh_status`:

  Binary GnRH classification (`neg`, `pos`).

- `gnrh_class`:

  Internal classification labels.

- `gnrh_score`:

  Composite GnRH detection score.

- `gnrh_expr`:

  Normalized `GNRH1` expression.

- `gnrh_raw`:

  Raw `GNRH1` UMI counts.

- `gnrh_core_hits`:

  Number of detected core GnRH markers.

- `gnrh_mig_hits`:

  Number of migration marker hits.

- `gnrh_neuro_hits`:

  Number of neuroendocrine marker hits.

- `gnrh_knn`:

  Neighborhood support score.

- `gnrh_confident`:

  High-confidence gnrh_confident classification.

Detection parameters and classification diagnostics are stored in
`object@misc`.

If `GNRH1` is not found, gene aliases are searched (`GNRH1`, `Gnrh1`,
`gnrh1`).

Composite scoring combines normalized expression, marker module
enrichment, ambient correction, and neighborhood support.

Diagnostic plots and summary outputs are generated via
[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md).

## See also

[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md),
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
