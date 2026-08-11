# Identify specific and reproducible GnRH-expressed genes

Identifies genes enriched in GnRH-positive cells relative to matched
control cells, measures their same-cell co-detection with `GNRH1`, and
optionally evaluates their reproducibility across donors.

## Usage

``` r
find_gnrh_genes(
  object,
  status_col = "gnrh_status",
  positive = "pos",
  negative = "neg",
  confidence_col = NULL,
  confidence_value = "pos",
  annotation_col = "ann1",
  control_ident = "Neuronal",
  donor_col = NULL,
  assay = "RNA",
  layer = "data",
  coexpr_gene = "GNRH1",
  methods = "wilcox",
  min_fc = 0.5,
  min_pct = 0.2,
  max_control_pct = 0.1,
  min_coexpr_pct = 20,
  max_padj = 0.05,
  min_donor_pct = 0.1,
  min_cells_donor = 3,
  generic_genes = c("TUBB3", "RBFOX3", "GAPDH", "ACTB", "MALAT1", "GAD1", "GAD2",
    "SLC17A7", "SNAP25", "SYT1"),
  seed = 1234,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing GnRH classification metadata and normalized
  expression data.

- status_col:

  Character. Metadata column containing GnRH-positive and GnRH-negative
  classifications. Default is `"gnrh_status"`.

- positive:

  Character. Value in `status_col` identifying GnRH-positive cells.
  Default is `"pos"`.

- negative:

  Character. Value in `status_col` identifying GnRH-negative cells.
  Default is `"neg"`.

- confidence_col:

  Optional character. Metadata column used to restrict the positive
  population to confident GnRH calls. If `NULL`, no confidence filter is
  applied.

- confidence_value:

  Value or values in `confidence_col` defining confident cells. Default
  is `"pos"`.

- annotation_col:

  Character. Metadata column containing broad cell-type annotations used
  to select matched controls. Default is `"ann1"`.

- control_ident:

  Character vector. Values in `annotation_col` defining the control
  population. Default is `"Neuronal"`.

- donor_col:

  Optional character. Metadata column identifying biological donors or
  samples. If supplied, candidate detection is evaluated independently
  in each donor.

- assay:

  Character. Seurat assay used for marker discovery and co-expression
  analysis. Default is `"RNA"`.

- layer:

  Character. Assay layer containing normalized expression values.
  Default is `"data"`.

- coexpr_gene:

  Character. Reference gene used for same-cell co-expression analysis.
  Default is `"GNRH1"`.

- methods:

  Character vector of differential-expression methods passed to
  [`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md).
  Default is `"wilcox"`.

- min_fc:

  Numeric. Minimum average log2 fold change required for final
  candidates. Default is `0.5`.

- min_pct:

  Numeric. Minimum detection fraction among GnRH-positive cells. Default
  is `0.20`.

- max_control_pct:

  Numeric. Maximum detection fraction allowed among control cells.
  Default is `0.10`.

- min_coexpr_pct:

  Numeric. Minimum percentage of `coexpr_gene`-positive GnRH cells in
  which a candidate must be detected. Expressed on a 0–100 scale.
  Default is `20`.

- max_padj:

  Numeric. Maximum adjusted p-value allowed for candidate markers.
  Default is `0.05`.

- min_donor_pct:

  Numeric. Minimum within-donor detection fraction used to count a
  candidate as supported by a donor. Default is `0.10`.

- min_cells_donor:

  Integer. Minimum number of selected GnRH cells required for a donor to
  be evaluated. Default is `3`.

- generic_genes:

  Character vector of generic neuronal, housekeeping, or broadly
  expressed genes excluded from the final candidate table.

- seed:

  Integer. Random seed used for differential-expression analysis.
  Default is `1234`.

- verbose:

  Logical. Display progress messages. Default is `TRUE`.

## Value

A named list containing:

- `candidates`:

  Filtered and ranked candidate GnRH genes. This table includes
  differential-expression statistics, GnRH-cell and control-cell
  detection fractions, correlation-based co-expression, same-cell
  co-detection counts and percentages, specificity score, and donor
  support when available.

- `markers`:

  The complete marker table returned by
  [`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md),
  augmented with same-cell co-detection metrics.

- `donor_detection`:

  A gene-by-donor matrix containing within-donor detection fractions.
  Returns `NULL` when `donor_col` is not supplied or no donor is
  evaluable.

- `donor_cells`:

  Number of selected GnRH cells per donor, or `NULL`.

- `gnrh_cells`:

  Cell barcodes assigned to the GnRH-positive group.

- `control_cells`:

  Cell barcodes assigned to the matched control group.

- `parameters`:

  Main selection thresholds used for the analysis.

## Details

This function provides a high-level interface to
[`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md).
The latter performs differential-expression testing and
expression-correlation analysis, while `find_gnrh_genes` additionally
selects matched controls, calculates the percentage of `GNRH1`-positive
cells expressing each marker, excludes generic genes, and ranks
candidates by specificity and donor support.

Candidate genes are required to be positively enriched in GnRH cells,
detected in at least `min_pct` of GnRH cells, detected in no more than
`max_control_pct` of matched controls, co-detected with `coexpr_gene` in
at least `min_coexpr_pct` percent of reference gene-positive cells, and
significant at `max_padj`.

The candidate specificity score is: \$\$ logFC \times (pct.1 - pct.2)
\times log(1 + coexpression\\ percentage) \$\$

Correlation-based co-expression and same-cell co-detection measure
distinct properties. Correlation evaluates coordinated expression
variation, whereas co-detection records whether both transcripts are
detected in the same cells. Neither metric alone demonstrates a direct
regulatory or molecular interaction.

Donor support should be interpreted cautiously when few GnRH cells are
available. Donors represented by fewer than `min_cells_donor` selected
cells are excluded from donor-level evaluation.

For Seurat v5 objects containing multiple layers for the selected assay,
layers should be joined with `Seurat::JoinLayers()` before calling this
function.

## See also

[`gnrh_markers`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_markers.md),
[`gnrh_marker_programs`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md),
[`gnrh_stage_modules`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_modules.md)

Other marker discovery:
[`find_gnrh_stage_markers()`](https://ymbouamboua.github.io/GnRHcell/reference/find_gnrh_stage_markers.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- Seurat::JoinLayers(obj, assay = "RNA")

result <- find_gnrh_genes(
  object = obj,
  status_col = "gnrh_status",
  positive = "pos",
  negative = "neg",
  confidence_col = "gnrh_confident",
  confidence_value = "pos",
  annotation_col = "ann1",
  control_ident = "Neuronal",
  donor_col = "status"
)

head(result$candidates, 20)
result$donor_cells
result$donor_detection

# Run without confidence or donor filtering
result <- find_gnrh_genes(
  object = obj,
  annotation_col = "ann1",
  control_ident = "Neuronal"
)
} # }
```
