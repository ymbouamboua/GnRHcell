# GnRHcell development

- `run_gnrh()` now returns compact metadata by default, retaining eight essential
  per-cell results. Detailed owned fields are archived by cell name under
  `misc$gnrh$cell_metadata`; `metadata = "full"` restores the full layout.
- Added `gnrh_compact()`, `gnrh_metadata()` and `gnrh_restore()` for migration,
  aligned access after subsetting, and explicit restoration. GnRHcell readers
  and plots support archived fields. User-owned metadata is preserved.
- Stopped creating or overwriting the unprefixed `total_hits` and
  `total_hits_bin` aliases; use the corresponding `gnrh_` fields instead.
