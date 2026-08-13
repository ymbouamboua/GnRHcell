# Detect GnRH neurons from single-cell RNA-seq data

Identifies candidate gonadotropin-releasing hormone (GnRH) neurons in a
Seurat object using direct GNRH1 detection together with independent
transcriptomic evidence.

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
  supported_q = 0.6,
  candidate_q = 0.95,
  scale_factor = 10000,
  max_alternative = 0.75,
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

  Number of nearest neighbors. Default is 20.

- min_umi:

  Minimum raw `GNRH1` UMI count required for direct detection. Default
  is 2.

- min_counts:

  Minimum total UMI count required per cell. Default is 500.

- mad_factor:

  Multiplier applied to the MAD-based adaptive `GNRH1` expression
  threshold. Default is 2.

- supported_q:

  Quantile of the direct-cell transcriptomic support distribution used
  for low-expression supported candidates. Default is 0.60.

- candidate_q:

  Quantile of the direct-cell transcriptomic support distribution used
  only to flag GNRH1-negative transcriptomic candidates. These cells are
  never classified as GnRH-positive. Default is 0.95.

- scale_factor:

  Library normalization scale factor. Default is 10000.

- max_alternative:

  Maximum alternative identity score tolerated for `GNRH1`-dropout
  rescue. Default is 0.75.

- verbose:

  Logical. Whether to print progress messages.

## Value

A Seurat object containing GnRH classifications, scores, diagnostics,
and detection parameters.

## Details

Detection integrates:

- direct `GNRH1` expression;

- GnRH identity and specification programs;

- migration-associated programs;

- neuroendocrine maturation programs;

- hormonal responsiveness;

- alternative neuronal or neuroendocrine identity programs;

- ambient RNA information;

- neighborhood enrichment using k-nearest neighbors; and

- adaptive transcriptomic support thresholds.

Classification uses two GnRH-positive routes: `direct` and `supported`.
Both routes require detectable `GNRH1` expression together with
independent GnRH identity evidence.

Cells without detected `GNRH1` are never classified as GnRH-positive.
However, cells showing strong GnRH-like transcriptomic evidence can be
flagged separately as `gnrh_dropout_candidate` for diagnostic purposes.

Direct candidates require at least `min_umi` raw `GNRH1` counts and
moderate independent GnRH identity evidence. Cells meeting the UMI rule
without identity evidence are retained as `gnrh_direct_isolated` and
`gnrh_direct_signal`, but remain GnRH-negative.

Supported candidates contain detectable but sub-threshold `GNRH1` and
additionally require independent GnRH identity and transcriptomic
support.

Cells with no detected `GNRH1` are not classified as GnRH-positive,
regardless of their transcriptomic support score. Strong GnRH-like
GNRH1-negative cells may instead be flagged as `gnrh_dropout_candidate`
for exploratory or diagnostic analyses.

Migration-associated expression contributes supportive evidence but
cannot independently establish GnRH identity.

The main `gnrh_score` includes direct `GNRH1` information, whereas
`gnrh_support_score` deliberately excludes direct `GNRH1` expression and
ambient-RNA information.

## See also

[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md),
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
