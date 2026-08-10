# Stage GnRH lineage cells

Assign developmental states to GnRH-lineage cells using biologically
informed transcriptional programs.

## Usage

``` r
stage_gnrh(
  object,
  assay = "RNA",
  layer = "data",
  min_migration_hits = 2L,
  min_secretory_core_hits = 1L,
  min_secretory_supportive_hits = 2L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- assay:

  Assay used for developmental and secretory module scoring. Default is
  `"RNA"`.

- layer:

  Expression layer used for module scoring. Default is `"data"`.

- min_migration_hits:

  Minimum number of expressed migration-core markers required to retain
  a raw `migrating` assignment. Default is `2L`.

- min_secretory_core_hits:

  Minimum number of core secretory markers required before secretory
  support can be assigned. Default is `1L`.

- min_secretory_supportive_hits:

  Minimum number of supportive secretory markers required when only one
  core secretory marker is detected. Default is `2L`.

- verbose:

  Logical. Retained for API consistency. Progress reporting is normally
  handled by
  [`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md).

## Value

A Seurat object containing developmental stage assignments,
developmental scores, migration-core evidence, and independent secretory
transcriptional support.

Added metadata include:

- `gnrh_stage_raw`:

  Developmental stage assigned directly from the maximum developmental
  module score.

- `gnrh_stage`:

  Final developmental stage after migration-core validation and masking
  of GnRH-negative cells.

- `gnrh_stage_reassigned`:

  Logical indicator specifying whether the raw developmental stage was
  reassigned.

- `gnrh_stage_reason`:

  Reason for the final developmental-stage assignment.

- `gnrh_migration_core_hits`:

  Number of expressed migration-core markers detected per cell.

- `gnrh_secretory_core_hits`:

  Number of expressed core secretory markers detected per cell.

- `gnrh_secretory_supportive_hits`:

  Number of expressed supportive secretory markers detected per cell.

- `gnrh_secretory_hits`:

  Total number of core and supportive secretory markers detected.

- `gnrh_secretory`:

  Independent transcriptional support for secretory machinery,
  classified as `limited`, `supported`, or `non-gnrh`.

## Details

Developmental staging is based on three mutually exclusive states:

- `identity`: lineage specification and early GnRH identity;

- `migrating`: migration and axon-guidance programs; and

- `mature`: neuroendocrine maturation.

Secretory machinery is evaluated independently from developmental stage.
This allows cells to retain a developmental annotation such as
`migrating` or `mature` while independently receiving evidence for a
neuroendocrine secretory transcriptional program.

A raw developmental stage is first assigned from the highest
developmental module score. Raw `migrating` assignments are then
validated using a curated migration-core marker set. Cells lacking
sufficient migration-core evidence are reassigned to the highest-scoring
alternative developmental stage.

Secretory support requires expression of at least one core secretory
marker together with either an additional core marker or sufficient
supportive secretory evidence. The resulting annotation is classified as
`limited` or `supported`.

Cells classified as GnRH-negative by
[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md)
are labeled `non-gnrh` in both developmental-stage and secretory
annotations.

Developmental stage and secretory support are intentionally modeled as
separate dimensions. Secretory-program expression therefore does not
replace or override the developmental-stage assignment.

A raw `migrating` assignment is retained only when at least
`min_migration_hits` migration-core markers are detected. Otherwise, the
cell is reassigned to the highest-scoring alternative developmental
state among `identity` and `mature`.

Secretory support is assigned when the cell expresses at least
`min_secretory_core_hits` core secretory markers and either:

- at least two core secretory markers; or

- at least `min_secretory_supportive_hits` supportive secretory markers.

The secretory annotation reflects transcriptional support for
neuroendocrine secretory machinery and should not be interpreted as a
direct measurement of GnRH peptide release.

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`run_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/run_gnrh.md)
