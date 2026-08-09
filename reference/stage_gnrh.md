# Stage GnRH lineage cells

Assign developmental states to GnRH-lineage cells using biologically
informed transcriptional programs.

## Usage

``` r
stage_gnrh(
  object,
  assay = "RNA",
  layer = "data",
  min_migration_hits = 1L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- assay:

  Assay used for expression extraction. Default is `"RNA"`.

- layer:

  Expression layer used for developmental module scoring. Default is
  `"data"`.

- min_migration_hits:

  Minimum number of expressed migration-core markers required to retain
  a raw `migrating` assignment. Default is `1L`.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

## Value

A Seurat object updated with:

- `gnrh_stage_raw`:

  Developmental stage assigned directly from the maximum module score.

- `gnrh_stage`:

  Final developmental stage after migration validation and masking of
  non-GnRH cells.

- `gnrh_stage_reassigned`:

  Logical indicator specifying whether the raw developmental stage was
  reassigned during migration validation.

- `gnrh_stage_reason`:

  Reason for the final stage assignment.

- `gnrh_migration_core_hits`:

  Number of expressed migration-core markers detected per cell.

- `gnrh_identity_score`:

  Identity module score.

- `gnrh_migrating_score`:

  Migration module score.

- `gnrh_mature_score`:

  Mature neuroendocrine module score.

- `gnrh_secreting_score`:

  Secretory activity module score.

Developmental module definitions, migration-core markers, and staging
parameters are stored in:

- `object@misc$gnrh_stage_modules`

- `object@misc$gnrh_migration_core`

- `object@misc$gnrh_stage_parameters`

## Details

Cells are scored against predefined developmental modules representing
major GnRH neuron states:

- `identity`: lineage specification and early GnRH identity

- `migrating`: migration and axon-guidance programs

- `mature`: neuroendocrine maturation

- `secreting`: secretory and vesicle machinery activation

A raw developmental stage is first assigned from the highest module
score. Migration assignments are then validated using migration-specific
marker evidence. Cells initially assigned as `migrating` but lacking the
required number of migration-core marker hits are reassigned to the
highest-scoring alternative stage.

If GnRH classification metadata are present, non-GnRH cells are labeled
as `non-gnrh` in the final stage assignment.

Developmental staging is based on predefined transcriptional programs
reflecting known biological states of GnRH neuron development.

Raw stage assignments are obtained from the highest developmental module
score. Because general neuronal and axon-guidance genes can produce
elevated migration scores in mature neurons, raw `migrating` assignments
require additional migration-core evidence.

Cells failing this migration criterion are reassigned to the
highest-scoring stage among `identity`, `mature`, and `secreting`.

If `gnrh_status` metadata are present from
[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
cells classified as negative are labeled `non-gnrh` in the final stage
assignment.

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
