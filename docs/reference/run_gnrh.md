# Run complete GnRHcell analysis pipeline

Executes the full GnRHcell workflow for identification, developmental
staging, diagnostics, and runtime reporting.

## Usage

``` r
run_gnrh(
  object,
  detect = TRUE,
  stage = TRUE,
  diagnostics = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- detect:

  Logical; run GnRH detection. Default is `TRUE`.

- stage:

  Logical; run developmental staging. Default is `TRUE`.

- diagnostics:

  Logical; run diagnostic analysis. Default is `TRUE`.

- verbose:

  Logical; print progress messages. Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md).

## Value

A Seurat object updated with GnRH detection, developmental staging,
diagnostics, and runtime metadata.

## Details

The pipeline can perform:

- GnRH neuron detection via
  [`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md)

- developmental stage assignment via
  [`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md)

- diagnostic analysis via
  [`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md)

Runtime metrics, dataset summaries, analysis parameters, and session
information are stored in `object@misc$gnrh$run_info`.

Workflow steps:

1.  GnRH detection using transcriptomic and marker-based scoring

2.  developmental state assignment

3.  diagnostic performance evaluation

4.  runtime and summary reporting

Metadata added may include:

- `gnrh_status`:

  Binary GnRH classification.

- `gnrh_class`:

  Internal classification labels.

- `gnrh_truth`:

  High-confidence truth labels.

- `gnrh_stage`:

  Developmental stage assignments.

Stored runtime metadata:

- `detect_sec`:

  Detection runtime in seconds.

- `stage_sec`:

  Staging runtime in seconds.

- `diagnostics_sec`:

  Diagnostics runtime in seconds.

- `total_sec`:

  Total pipeline runtime.

If specific steps are disabled, only selected workflow components are
executed.

## See also

[`detect_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/detect_gnrh.md),
[`stage_gnrh`](https://ymbouamboua.github.io/GnRHcell/reference/stage_gnrh.md),
[`gnrh_diagnostics`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_diagnostics.md),
[`extract_gnrh_run_info`](https://ymbouamboua.github.io/GnRHcell/reference/extract_gnrh_run_info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- run_gnrh(seurat_obj)

obj <- run_gnrh(
  seurat_obj,
  detect = TRUE,
  stage = TRUE,
  diagnostics = TRUE
)
} # }
```
