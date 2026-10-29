# GnRHcell development

- Fixed report legends to retain status/stage names alongside counts and
  percentages, and prevent rare negative proportions from rounding to 100.0%.
- QC zooms now focus on positive cells by default, show their frame on the
  global view, use square distance-preserving limits, and report cells outside
  the central zoom. Added explicit-cell, confidence and full-extent controls.

- Added `gnrh_qc()` with evidence, technical-QC, sample-denominator, stage
  confidence, embedding and candidate-zoom panels, plus exportable data tables.
- Added `gnrh_stability()` to rerun detection across parameter grids, recompute
  adaptive thresholds, and summarize counts, positive-set overlap and per-cell
  frequencies without altering the input object. Both functions read compact
  metadata, and provide printable result objects.

- `run_gnrh()` now returns compact metadata by default, retaining eight essential
  per-cell results. Detailed owned fields are archived by cell name under
  `misc$gnrh$cell_metadata`; `metadata = "full"` restores the full layout.
- Added `gnrh_compact()`, `gnrh_metadata()` and `gnrh_restore()` for migration,
  aligned access after subsetting, and explicit restoration. GnRHcell readers
  and plots support archived fields. User-owned metadata is preserved.
- Stopped creating or overwriting the unprefixed `total_hits` and
  `total_hits_bin` aliases; use the corresponding `gnrh_` fields instead.
