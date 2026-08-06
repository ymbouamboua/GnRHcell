# Package index

## Main workflow

High-level functions for running the complete GnRH-neuron detection,
staging, diagnostic, and reporting workflow.

- [`run_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
  : Run complete GnRHcell analysis pipeline
- [`gnrh_report()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_report.md)
  : Generate a comprehensive GnRH diagnostic report

## GnRH-neuron detection and staging

Detection, classification, input validation, and developmental-stage
assignment.

- [`detect_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md)
  : Detect GnRH neurons from single-cell RNA-seq data
- [`stage_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md)
  : Stage GnRH lineage cells
- [`validate_input()`](https://ymbouamboua.github.io/GnRHcell/reference/validate_input.md)
  : Validate input Seurat object for GnRH analysis

## Marker discovery and programs

Identification of GnRH-specific genes, marker programs, co-expression
patterns, and cross-dataset overlap.

- [`find_gnrh_genes()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_genes.md)
  : Identify specific and reproducible GnRH-expressed genes
- [`gene_upset()`](https://ymbouamboua.github.io/GnRHcell/reference/gene_upset.md)
  : Gene Set Overlap Analysis and Visualization
- [`gnrh_stage_gene_references()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_gene_references.md)
  : GnRH developmental marker gene references
- [`gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md)
  : Classify GnRH marker candidates into developmental programs
- [`gnrh_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md)
  : Detect positive GnRH lineage markers
- [`plot_gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_marker_programs.md)
  : Plot GnRH Marker Program Results
- [`build_gene_sets()`](https://ymbouamboua.github.io/GnRHcell/reference/build_gene_sets.md)
  : Build gene sets from marker tables
- [`gnrh_stage_modules()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md)
  : Default GnRH developmental marker programs

## Diagnostics and quality control

Diagnostic summaries and quality-control functions.

- [`gnrh_diagnostics()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md)
  : Run GnRH detection diagnostics

## Visualization

Visualization of cellular embeddings, expression, stages, co-expression
patterns, and distributions.

- [`celldot()`](https://ymbouamboua.github.io/GnRHcell/reference/celldot.md)
  : Create an enhanced Seurat dot plot
- [`cellfeature()`](https://ymbouamboua.github.io/GnRHcell/reference/cellfeature.md)
  : Cell Feature Plot for Seurat Objects
- [`cellmap()`](https://ymbouamboua.github.io/GnRHcell/reference/cellmap.md)
  : Visualize cells on a Seurat embedding
- [`cellviolin()`](https://ymbouamboua.github.io/GnRHcell/reference/cellviolin.md)
  : Create violin plots for Seurat features
- [`dark_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/dark_theme.md)
  : Dark ggplot2 theme
- [`plot_gnrh_coexpr()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_coexpr.md)
  : Plot GnRH coexpression markers
- [`plot_gnrh_detected()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_detected.md)
  : Plot detected GnRH-positive cells across datasets
- [`plot_gnrh_distribution()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_distribution.md)
  : Plot metadata distributions
- [`plot_gnrh_embedding()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_embedding.md)
  : Plot GnRHcell metadata on cell embeddings
- [`plot_gnrh_feature()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_feature.md)
  : Plot GnRH feature maps
- [`plot_gnrh_hits()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_hits.md)
  : Plot GnRH module hit distributions
- [`plot_gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_marker_programs.md)
  : Plot GnRH Marker Program Results
- [`plot_gnrh_runtime_curve()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_runtime_curve.md)
  : Plot GnRHcell runtime across datasets
- [`plot_network()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_network.md)
  : GnRH gene similarity network
- [`plot_theme()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_theme.md)
  : GnRHcell ggplot2 theme

## Colors and themes

Color palettes and graphical helpers.

- [`cellpal()`](https://ymbouamboua.github.io/GnRHcell/reference/cellpal.md)
  : Generate a discrete color palette
- [`gnrh_colors()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_colors.md)
  : GnRHcell color palettes

## Example data

Example datasets distributed with GnRHcell.

- [`hpsc_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/hpsc_gnrh.md)
  : Example GnRHcell Seurat object

## Scores and classification

Functions for computing detection scores, classification confidence, and
developmental-stage scores.

## Data extraction and summaries

Extraction and summarization of GnRHcell results.

- [`extract_gnrh_run_info()`](https://ymbouamboua.github.io/GnRHcell/reference/extract_gnrh_run_info.md)
  : Extract GnRHcell run information
