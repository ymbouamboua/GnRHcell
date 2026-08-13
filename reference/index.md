# Package index

## Main workflow

High-level functions for running the complete GnRH-neuron detection,
staging, diagnostic, and reporting workflow.

- [`run_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
  : Run complete GnRHcell analysis pipeline
- [`gnrh_report()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_report.md)
  : Generate a publication-ready GnRHcell diagnostic report

## GnRH-neuron detection and staging

Input validation, detection, classification, and developmental-stage
assignment.

- [`validate_input()`](https://ymbouamboua.github.io/GnRHcell/reference/validate_input.md)
  : Validate input Seurat object for GnRH analysis
- [`detect_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md)
  : Detect GnRH neurons from single-cell RNA-seq data
- [`stage_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md)
  : Stage GnRH lineage cells
- [`gnrh_diagnostics()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md)
  : Run GnRH detection diagnostics

## Marker discovery and biological programs

GnRH-lineage and stage-specific marker discovery, developmental
programs, co-expression, and cross-dataset overlap.

- [`gnrh_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md)
  : Detect positive GnRH lineage markers
- [`find_gnrh_genes()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_genes.md)
  : Identify specific and reproducible GnRH-expressed genes
- [`find_gnrh_stage_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_stage_markers.md)
  : Identify stage-specific markers in GnRH-lineage cells
- [`gnrh_stage_modules()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md)
  : Default GnRH developmental marker programs
- [`gnrh_stage_gene_references()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_gene_references.md)
  : GnRH developmental marker gene references
- [`build_gene_sets()`](https://ymbouamboua.github.io/GnRHcell/reference/build_gene_sets.md)
  : Build gene sets from marker tables
- [`gnrh_gene_upset()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_gene_upset.md)
  : Gene Set Overlap Analysis and Visualization
- [`gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md)
  : Classify GnRH marker candidates into developmental programs
- [`merge_gnrh_marker_scores()`](https://ymbouamboua.github.io/GnRHcell/reference/merge_gnrh_marker_scores.md)
  : Merge GnRH marker scores across datasets
- [`plot_gnrh_conserved_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_conserved_markers.md)
  : Plot conserved GnRH markers across datasets
- [`plot_gnrh_dataset_specific_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_dataset_specific_markers.md)
  : Plot dataset-specific GnRH markers

## Multi-dataset workflows

Preparation, processing, comparison, and validation of single-cell
dataset collections.

- [`prepare_gnrh_datasets()`](https://ymbouamboua.github.io/GnRHcell/reference/prepare_gnrh_datasets.md)
  : Prepare a GnRHcell dataset configuration table
- [`run_gnrh_dataset()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh_dataset.md)
  : Run GnRHcell analysis on a single dataset
- [`run_gnrh_collection()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh_collection.md)
  : Run GnRHcell across multiple datasets
- [`compare_gnrh_datasets()`](https://ymbouamboua.github.io/GnRHcell/reference/compare_gnrh_datasets.md)
  : Compare GnRHcell results across datasets
- [`validate_gnrh_collection()`](https://ymbouamboua.github.io/GnRHcell/reference/validate_gnrh_collection.md)
  : Validate a GnRHcell multi-dataset collection
- [`extract_gnrh_run_info()`](https://ymbouamboua.github.io/GnRHcell/reference/extract_gnrh_run_info.md)
  : Extract GnRHcell run information

## Visualization

Publication-ready cellular embeddings, expression plots, diagnostic
summaries, distributions, networks, and cross-dataset figures.

- [`plot_gnrh_embedding()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_embedding.md)
  : Plot GnRH embedding
- [`plot_gnrh_feature()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_feature.md)
  : Cell Feature Plot for Seurat Objects
- [`plot_gnrh_dot()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_dot.md)
  : Create an enhanced Seurat dot plot
- [`plot_gnrh_distribution()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_distribution.md)
  : Plot GnRHcell metadata distributions
- [`plot_gnrh_hits()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_hits.md)
  : Plot GnRH module hit distributions
- [`plot_gnrh_coexpr()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_coexpr.md)
  : Plot GnRH coexpression markers
- [`plot_gnrh_specificity()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_specificity.md)
  : Plot GnRH specificity landscape
- [`plot_class_counts()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_class_counts.md)
  : Plot GnRH classification counts
- [`plot_network()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_network.md)
  : GnRH gene similarity network
- [`plot_gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_marker_programs.md)
  : Plot GnRH marker program results
- [`plot_gnrh_runtime_curve()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_runtime_curve.md)
  : Plot GnRHcell runtime across datasets
- [`plot_gnrh_detected()`](https://ymbouamboua.github.io/GnRHcell/reference/plot_gnrh_detected.md)
  : Plot detected GnRH-positive cells across datasets

## Colors

Stable color palettes for GnRHcell classifications and stages.

- [`gnrh_colors()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_colors.md)
  : GnRHcell color palettes

## Example data

Sampled Seurat objects and reference data for examples and testing.

- [`hpsc`](https://ymbouamboua.github.io/GnRHcell/reference/hpsc.md) :
  Human hPSC GnRH demo dataset
