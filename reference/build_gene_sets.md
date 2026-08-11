# Build gene sets from marker tables

Reads marker tables from multiple datasets and extracts unique gene
symbols into a named list suitable for overlap analysis, UpSet plots, or
marker comparison workflows.

## Usage

``` r
build_gene_sets(files, dir, gene_col = "gene")
```

## Arguments

- files:

  Named character vector containing marker table filenames. Names
  correspond to dataset identifiers and values correspond to file names.

- dir:

  Character. Directory containing marker tables.

- gene_col:

  Character. Name of the column containing gene symbols. Default is
  `"gene"`.

## Value

A named list where each element contains a character vector of unique
gene symbols for a dataset.

## Details

Gene names are automatically standardized to uppercase to ensure
consistent comparisons across species and datasets.

For each dataset:

- Marker tables are imported using
  [`read.delim`](https://rdrr.io/r/utils/read.table.html).

- Missing values and empty gene names are removed.

- Duplicate genes are removed.

- Gene symbols are converted to uppercase.

- Gene symbols are sorted alphabetically.

This standardization ensures robust overlap analysis between datasets
originating from different species or annotation conventions.

## See also

[`gnrh_gene_upset`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_gene_upset.md),
[`gnrh_marker_programs`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md)

## Examples

``` r
if (FALSE) { # \dontrun{

files <- c(
  "HuDeCa Nose" = "gnrh_nose_markers.tsv",
  "HPSC Wang 2022" = "gnrh_wang_markers.tsv",
  "Human HypoMap" = "gnrh_human_hypomap_markers.tsv"
)

gene_sets <- build_gene_sets(
  files = files,
  dir = file.path(outdir, "tables")
)

names(gene_sets)
lengths(gene_sets)

} # }
```
