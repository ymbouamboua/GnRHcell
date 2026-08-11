# Classify GnRH marker candidates into developmental programs

Combines marker tables, dataset-specific marker overlap results, GNRH1
co-expression evidence, specificity scores, and curated developmental
programs to identify candidate GnRH markers associated with identity,
migration, maturation, or secretion.

## Usage

``` r
gnrh_marker_programs(
  files,
  results,
  outdir,
  modules = gnrh_stage_modules(),
  coexpr_col = "coexpr_flag",
  gene_col = "gene",
  min_medium_score = 0.25,
  min_high_score = 1,
  write_output = TRUE
)
```

## Arguments

- files:

  Named character vector. Names correspond to dataset names and values
  correspond to existing marker-table paths or to file names located in
  `file.path(outdir, "tables")`.

- results:

  A list returned by
  [`gnrh_gene_upset()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_gene_upset.md),
  containing at least `results$unique`, a named list of dataset-specific
  unique genes.

- outdir:

  Character. Output directory containing a `tables/` subdirectory with
  marker tables.

- modules:

  Named list of developmental gene modules. By default, uses
  [`gnrh_stage_modules()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md).

- coexpr_col:

  Character. Name of the column indicating GNRH1 co-expression. Default
  is `"coexpr_flag"`.

- gene_col:

  Character. Name of the gene column in marker tables. Default is
  `"gene"`.

- min_medium_score:

  Numeric. Minimum specificity score required for medium confidence.
  Default is `0.25`.

- min_high_score:

  Numeric. Minimum specificity score required for high confidence.
  Default is `1`.

- write_output:

  Logical. If `TRUE`, writes output tables to
  `file.path(outdir, "tables")`. Default is `TRUE`.

## Value

A named list with four elements:

- `candidate_table`:

  Data frame containing all evaluated marker candidates. Columns are:

  `dataset`

  :   Dataset where the candidate marker was detected.

  `gene`

  :   Candidate marker gene symbol, standardized to uppercase.

  `is_unique`

  :   Logical value indicating whether the gene is unique to one dataset
      according to
      [`gnrh_gene_upset`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_gene_upset.md).

  `known_status`

  :   Biological annotation of the gene. Values include
      `"known_GnRH_or_developmental"`, `"candidate_novel"`,
      `"generic_neuronal"`, and `"unknown_or_context_specific"`.

  `program`

  :   Developmental program assigned using
      [`gnrh_stage_modules`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md).
      Values include `"identity"`, `"migrating"`, `"mature"`,
      `"secreting"`, or `"unassigned"`.

  `coexpr_GNRH1`

  :   Logical value indicating whether the gene is co-expressed with
      `GNRH1` according to the marker table co-expression column.

  `specificity_score`

  :   Numeric marker specificity score. Higher values indicate stronger
      enrichment in GnRH-positive cells.

  `confidence_level`

  :   Final confidence category assigned by the function: `"high"`,
      `"medium"`, or `"low"`.

- `high_confidence`:

  Subset of `candidate_table` classified as high-confidence candidate
  markers.

- `summary`:

  Summary table of candidate genes grouped by dataset, developmental
  program, and confidence level.

- `split_by_dataset_program`:

  List of candidate genes split by dataset and developmental program.

## Details

This function is designed to be used downstream of marker discovery and
overlap analysis. It reads marker tables generated for multiple
datasets, intersects dataset-specific unique markers with
GNRH1-coexpressed genes, assigns candidates to curated GnRH
developmental programs, and reports confidence levels based on
co-expression, specificity, and program membership.

The specificity score is computed as: \$\$ avg\\log2FC \times (pct.1 -
pct.2) \times -log10(p\\val\\adj + \epsilon) \$\$ where \\\epsilon =
1e-300\\. Missing columns are handled safely and replaced with
conservative default values.

Confidence levels are assigned as follows:

- `high`: unique, GNRH1-coexpressed, assigned to a developmental
  program, specificity score greater than or equal to `min_high_score`,
  and not classified as a generic neuronal marker.

- `medium`: unique, GNRH1-coexpressed, specificity score greater than or
  equal to `min_medium_score`, and not classified as generic.

- `low`: all remaining candidates.

Biologically, `candidate_table` separates reference GnRH markers,
generic neuronal markers, potentially novel GnRH-associated markers, and
context-specific candidates. Computationally, genes are prioritized by
dataset specificity, GNRH1 co-expression, developmental program
membership, and marker specificity score.

## See also

[`gnrh_stage_modules`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md),
[`gnrh_gene_upset`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_gene_upset.md),
[`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md)

## Examples

``` r
if (FALSE) { # \dontrun{
files <- c(
  "HuDeCa Nose" = "gnrh_nose_markers.tsv",
  "HuDeCa Hypo" = "gnrh_hudeca_hypo_markers.tsv",
  "HPSC Wang 2022" = "gnrh_wang_markers.tsv",
  "Human HypoMap" = "gnrh_human_hypomap_markers.tsv",
  "Mouse HypoMap" = "gnrh_mouse_hypomap_markers.tsv",
  "Mouse Amato 2024" = "gnrh_mouse_amato_markers.tsv",
  "Mouse POA" = "gnrh_mouse_poa_markers.tsv",
  "Mouse MBH" = "gnrh_mouse_mbh_markers.tsv"
)

marker_programs <- gnrh_marker_programs(
  files = files,
  results = overlap_results,
  outdir = "results/gnrh",
  write_output = TRUE
)

marker_programs$candidate_table
marker_programs$high_confidence
marker_programs$summary
marker_programs$split_by_dataset_program
} # }
```
