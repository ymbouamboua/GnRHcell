# Identify stage-specific markers in GnRH-lineage cells

Finds genes enriched in each GnRH developmental stage relative to the
other GnRH-positive stages. Unlike
[`gnrh_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md),
this function does not require correlation or co-detection with `GNRH1`,
because those filters favor general GnRH-lineage genes rather than
stage-specific programs.

## Usage

``` r
find_gnrh_stage_markers(
  object,
  stage_col = "gnrh_stage",
  stages = c("identity", "migrating", "mature"),
  class_col = "gnrh_class",
  positive_classes = c("direct", "supported"),
  exclude_stages = c("non-gnrh", "secreting"),
  assay = NULL,
  layer = "data",
  test_use = "wilcox",
  min_cells = 10L,
  min_pct = 0.1,
  min_log2fc = 0.25,
  max_padj = 0.05,
  min_specificity = 0.05,
  donor_col = NULL,
  min_cells_donor = 3L,
  min_donor_support = 0.6,
  min_donors = 2L,
  max_cells_per_ident = Inf,
  exclude_pattern = NULL,
  seed = 1234L,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object processed with
  [`run_gnrh()`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md).

- stage_col:

  Metadata column containing developmental stages.

- stages:

  Stages to test. By default, `identity`, `migrating`, and `mature` are
  tested when present.

- class_col:

  Metadata column defining GnRH detection classes. Set to `NULL` to
  select cells using `stage_col` alone.

- positive_classes:

  Detection classes considered GnRH-lineage positive.

- exclude_stages:

  Stage values excluded from all comparisons.

- assay:

  Assay used for differential expression.

- layer:

  Normalized-expression layer used for donor concordance.

- test_use:

  Differential-expression test passed to
  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html).

- min_cells:

  Minimum number of cells required in both the target stage and its
  reference group.

- min_pct:

  Minimum expression fraction passed to
  [`FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html).

- min_log2fc:

  Minimum positive average log2 fold change.

- max_padj:

  Maximum adjusted p-value.

- min_specificity:

  Minimum `pct.1 - pct.2` detection difference.

- donor_col:

  Optional biological replicate column. When supplied, the direction of
  the stage effect is evaluated independently within donors.

- min_cells_donor:

  Minimum cells required in both comparison groups for a donor to be
  evaluable for a stage.

- min_donor_support:

  Minimum fraction of evaluable donors in which the mean normalized
  expression effect is positive.

- min_donors:

  Minimum number of evaluable donors required before donor support is
  used as a candidate filter.

- max_cells_per_ident:

  Optional maximum cells sampled from each identity by
  [`FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html).
  Useful for very large atlases.

- exclude_pattern:

  Optional regular expression for genes to exclude from the candidate
  table, for example mitochondrial or ribosomal genes.

- seed:

  Random seed used by differential-expression subsampling.

- verbose:

  Display progress messages.

- ...:

  Additional arguments passed to
  [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html).

## Value

A named list containing:

- candidates:

  Filtered and ranked stage-specific markers.

- markers:

  Complete positive marker results for all tested stages.

- donor_effects:

  Long table of within-donor mean-expression effects.

- stage_counts:

  Numbers of selected cells per stage.

- cells:

  Cell barcodes used in the analysis.

- parameters:

  Principal analysis parameters.

## Details

Each stage is tested against the union of the other selected GnRH stages
(one-versus-rest). The reported `specificity` is `pct.1 - pct.2`. The
ranking score combines positive fold change, detection specificity,
statistical significance, and—when available—donor concordance.

Cell-level differential expression is useful for marker discovery but
does not replace a replicate-aware pseudobulk analysis for formal
inference. `donor_support` should therefore be used to prioritize
reproducible markers, while final publication claims should be confirmed
with pseudobulk counts.

## See also

[`gnrh_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md),
[`find_gnrh_genes()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_genes.md)

Other marker discovery:
[`find_gnrh_genes()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_genes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
stage_markers <- find_gnrh_stage_markers(
  object = wang,
  stage_col = "gnrh_stage",
  donor_col = "orig.ident",
  stages = c("identity", "migrating", "mature")
)

head(stage_markers$candidates)
subset(stage_markers$candidates, stage == "migrating")

# Faster discovery in a very large atlas
stage_markers <- find_gnrh_stage_markers(
  object = human_hypomap,
  donor_col = "orig.ident",
  max_cells_per_ident = 5000
)
} # }
```
