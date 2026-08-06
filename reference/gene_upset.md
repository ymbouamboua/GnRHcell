# Gene Set Overlap Analysis and Visualization

Perform overlap analysis between multiple gene sets, export overlap
tables, compute unique/common genes, and generate publication-quality
Venn or UpSet plots depending on the number of gene sets.

## Usage

``` r
gene_upset(
  gene_sets,
  min_size = 1,
  venn_title = "Overlap of Gene Sets",
  outdir = ".",
  save_plot = TRUE,
  plot_width = 10,
  plot_height = 8,
  dpi = 600
)
```

## Arguments

- gene_sets:

  A named list of gene vectors. Each element should contain a character
  vector of gene symbols.

- min_size:

  Integer. Minimum intersection size to display in the UpSet plot.
  Default is `1`.

- venn_title:

  Character string specifying the plot title. Default is
  `"Overlap of Gene Sets"`.

- outdir:

  Output directory where CSV tables and figures will be saved. Default
  is current working directory.

- save_plot:

  Logical indicating whether plots should be exported. Default is
  `TRUE`.

- plot_width:

  Numeric width of exported figures in inches. Default is `10`.

- plot_height:

  Numeric height of exported figures in inches. Default is `8`.

- dpi:

  Numeric resolution for PNG export. Default is `600`.

## Value

A list containing cleaned gene sets, overlap tables, unique genes,
common genes, and the generated plot.

## Details

Prior to overlap analysis, gene symbols are standardized by:

- removing duplicated entries,

- removing missing values,

- removing empty strings,

- converting all gene symbols to uppercase.

This ensures robust overlap comparisons across datasets originating from
different species or annotation conventions.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_sets <- list(
  Dataset_A = c("GNRH1", "KISS1", "TAC3"),
  Dataset_B = c("GNRH1", "TAC3", "PAX6"),
  Dataset_C = c("GNRH1", "DLX1", "DLX2")
)

results <- gene_upset(
  gene_sets = gene_sets,
  venn_title = "GnRH Marker Overlap",
  outdir = "results/gene_overlap"
)

results$venn_plot
} # }
```
