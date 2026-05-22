
<img src="inst/figures/GnRHcell.png" align="right" width="180"/>

# GnRHcell

<!-- badges: start -->

[![R-CMD-check](https://github.com/ymbouamboua/GnRHcell/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ymbouamboua/GnRHcell/actions/workflows/R-CMD-check.yaml)
[![Reproducibility](https://github.com/ymbouamboua/GnRHcell/actions/workflows/reproducibility.yaml/badge.svg)](https://github.com/ymbouamboua/GnRHcell/actions/workflows/reproducibility.yaml)
<!-- badges: end -->

**High-confidence detection, developmental staging, and marker discovery
of GnRH neurons from single-cell RNA-seq data**

`GnRHcell` is an R package for identifying rare **gonadotropin-releasing
hormone (GnRH) neurons** in single-cell transcriptomic datasets.

The package provides an integrated framework for:

- high-confidence GnRH cell detection
- developmental stage assignment
- diagnostic performance evaluation
- marker discovery and coexpression analysis
- network visualization
- publication-ready graphics

Built for **Seurat workflows**, `GnRHcell` is optimized for exploratory
and reproducible analysis of rare GnRH neuron populations.

------------------------------------------------------------------------

# Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ymbouamboua/GnRHcell")
```

# Quick Start

``` r
suppressPackageStartupMessages({
library(GnRHcell)
library(Seurat)
})
```

    ## Warning: package 'Seurat' was built under R version 4.5.2

    ## Warning: package 'SeuratObject' was built under R version 4.5.2

    ## Warning: package 'sp' was built under R version 4.5.2

# Quick Start

## 1. Preprocess Seurat object

``` r
# Example: create a small Seurat object
mat <- matrix(rpois(2000, lambda = 5), nrow = 100)
obj <- CreateSeuratObject(mat)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:20)
```

## 2. Run complete GnRH pipeline

``` r
obj <- run_gnrh(obj)
```

Example console output:

``` text
[GNRH] ==== STARTING GnRHcell PIPELINE ====
[STEP] [1/3] Detecting GnRH cells
[INFO] ==== GNRH DETECTION START ====
[INFO] Using assay: RNA
[INFO] Matrix loaded: 33538 genes by 29708 cells
[INFO] Running diagnostics
[INFO] ==== GNRH DETECTION DONE ====
[DONE] Detection complete. Duration: 6.1s
[STEP] [2/3] Assigning developmental stages
[INFO] ==== GNRH STAGING START ====
[INFO] ==== GNRH STAGING DONE ====
[INFO] [3/3] Running diagnostics
[INFO] Running diagnostics
[DONE] Assigning stages complete. Duration: 6.1s
[INFO] PIPELINE SUMMARY
[INFO] Status:
[INFO]   neg: 26651
[INFO]   pos: 3057
[INFO] Truth:
[INFO]   neg: 27080
[INFO]   pos: 2628
[INFO] Stage:
[INFO]   identity: 1128
[INFO]   migrating: 1857
[INFO]   mature: 69
[INFO]   secreting: 3
[INFO]   non-gnrh: 26651
[DONE] ==== GnRHcell PIPELINE COMPLETE ==== Duration: 7.6s
```

## Visualization

### Embedding visualization

``` r
p <- plot_gnrh_embedding(
  obj,
  group.by = c("gnrh_status", "gnrh_stage")
)

p
```

### Feature expression

``` r
p <- plot_gnrh_feature(
  obj,
  feature_type = "all"
)

p
```

### Distribution plots

GnRH status by sample:

``` r
plot_gnrh_distribution(
  obj,
  group.by = "gnrh_status",
  split.by = "orig.ident",
  proportion = TRUE,
  label = FALSE,
  cols = gnrh_colors("status")
)
```

Developmental stages by sample:

``` r
plot_gnrh_distribution(
  obj,
  group.by = "gnrh_stage",
  split.by = "orig.ident",
  proportion = TRUE,
  label = FALSE,
  cols = gnrh_colors("stage")
)
```

## Diagnostic dashboard

``` r
p <- gnrh_report(obj)
p
```

Includes:

- signal landscape
- score distributions
- threshold performance
- ROC analysis
- module signal summary
- stage composition
- classification diagnostics

## Marker Discovery

Identify GnRH-associated markers:

``` r
markers <- gnrh_markers(obj)
head(markers)
```

Filter coexpressed markers:

``` r
df <- subset(markers, coexpr_flag == TRUE)
```

Coexpression ranking:

``` r
plot_gnrh_coexpr(
  df,
  coexp_cutoff = 0.3
)
```

## Gene network

``` r
plot_network(
  markers,
  top_n = 50,
  threshold = 0.1
)
```

## Core Functions

### Pipeline

- run_gnrh() — complete detection/staging/diagnostics workflow
- detect_gnrh() — GnRH cell detection
- stage_gnrh() — developmental staging
- gnrh_diagnostics() — performance diagnostics

### Visualization

- plot_gnrh_embedding()
- plot_gnrh_feature()
- plot_gnrh_distribution()
- plot_gnrh_coexpr()
- plot_network()
- gnrh_report()

### Marker analysis

- gnrh_markers()

### Utilities

- gnrh_colors()
- cellpal()
- plot_theme()

## Dependencies

Major dependencies:

- Seurat
- ggplot2
- patchwork
- ggrepel
- plotly
- pROC
- igraph
- Matrix

## Development

Run tests:

``` r
devtools::test()
```

Run full package checks:

``` r
devtools::check()
```

Rebuild README:

``` r
devtools::build_readme()
```

## Citation

If you use GnRHcell, please cite:

Yvon Mbouamboua. GnRHcell: High-confidence GnRH neuron detection and
staging from single-cell RNA-seq data.

## License

MIT License

## Status

GnRHcell is under active development. Interfaces may evolve as methods
are refined.
