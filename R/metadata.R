#' Metadata fields produced by GnRHcell
#' @keywords internal
#' @noRd
.gnrh_owned_metadata <- function() c(
  "gnrh_alternative_hits",
  "gnrh_alternative_score",
  "gnrh_alternative_strong",
  "gnrh_ambient_ratio",
  "gnrh_class",
  "gnrh_confident",
  "gnrh_core_hits",
  "gnrh_developmental_index",
  "gnrh_direct_isolated",
  "gnrh_direct_signal",
  "gnrh_direct_supported",
  "gnrh_dropout_candidate",
  "gnrh_early_specification_score",
  "gnrh_expr",
  "gnrh_guidance_score",
  "gnrh_hormone_score",
  "gnrh_identity_moderate",
  "gnrh_identity_neuro_hits",
  "gnrh_identity_primary_hits",
  "gnrh_identity_score",
  "gnrh_identity_strong",
  "gnrh_identity_supportive_hits",
  "gnrh_independent_support",
  "gnrh_knn",
  "gnrh_lineage_identity_score",
  "gnrh_mig_hits",
  "gnrh_migration_core_hits",
  "gnrh_migration_primary_hits",
  "gnrh_migration_score",
  "gnrh_migration_support",
  "gnrh_migration_supportive_hits",
  "gnrh_neuro_hits",
  "gnrh_neuro_primary_hits",
  "gnrh_neuro_score",
  "gnrh_neuro_support",
  "gnrh_neuro_supportive_hits",
  "gnrh_raw",
  "gnrh_reference_positive",
  "gnrh_review_candidate",
  "gnrh_score",
  "gnrh_score_raw",
  "gnrh_secretory",
  "gnrh_secretory_core_hits",
  "gnrh_secretory_hits",
  "gnrh_secretory_supported",
  "gnrh_secretory_supportive_hits",
  "gnrh_signal_class",
  "gnrh_signal_status",
  "gnrh_stage",
  "gnrh_stage_confident",
  "gnrh_stage_early_expression",
  "gnrh_stage_early_fraction",
  "gnrh_stage_early_score",
  "gnrh_stage_margin",
  "gnrh_stage_mature_expression",
  "gnrh_stage_mature_fraction",
  "gnrh_stage_mature_score",
  "gnrh_stage_migrating_expression",
  "gnrh_stage_migrating_fraction",
  "gnrh_stage_migrating_score",
  "gnrh_stage_raw",
  "gnrh_stage_reason",
  "gnrh_stage_reassigned",
  "gnrh_stage_resolution",
  "gnrh_stage_score",
  "gnrh_stage_second_score",
  "gnrh_status",
  "gnrh_support_score",
  "gnrh_support_score_raw",
  "gnrh_total_hits",
  "gnrh_total_hits_bin",
  "gnrh_transcriptomic_candidate"
)

#' Essential GnRHcell per-cell results
#' @keywords internal
#' @noRd
.gnrh_essential_metadata <- function() c(
  "gnrh_status", "gnrh_class", "gnrh_score", "gnrh_confident",
  "gnrh_stage", "gnrh_stage_confident", "gnrh_stage_resolution", "gnrh_secretory"
)

#' Validate requested retained metadata
#' @keywords internal
#' @noRd
.gnrh_check_keep <- function(keep) {
  if (!is.null(keep) && (!is.character(keep) || anyNA(keep) ||
      any(!nzchar(keep)) || anyDuplicated(keep) ||
      any(!keep %in% .gnrh_owned_metadata())))
    stop("`keep` / `metadata_keep` must contain unique package-owned GnRH metadata names.", call.=FALSE)
  invisible(keep)
}

#' Read GnRHcell results from metadata and miscellaneous storage
#'
#' @param object A Seurat object.
#' @param columns Optional character vector of columns to return. `NULL` returns
#'   ordinary metadata together with all archived GnRH fields.
#' @return A data frame indexed and ordered by current cell names. Current
#'   `meta.data` values take precedence over archived values. Cells not found in
#'   the archive have missing archived values. Factor levels are preserved.
#' @details Seurat subsetting does not necessarily subset `misc`. This accessor
#'   aligns archived rows by cell name, never by position. If cell names are
#'   changed, restore metadata before renaming, then compact again.
#' @export
gnrh_metadata <- function(object, columns=NULL) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call.=FALSE)
  md <- object[[]]
  archived <- object@misc$gnrh$cell_metadata
  if (!is.null(archived)) {
    if (!is.data.frame(archived) || is.null(rownames(archived)) || anyDuplicated(rownames(archived)))
      stop("Invalid archived GnRH cell metadata: expected a data frame with unique cell names.", call.=FALSE)
    add <- setdiff(colnames(archived), colnames(md))
    if (!is.null(columns)) add <- intersect(add, columns)
    if (length(add)) md[add] <- archived[match(rownames(md), rownames(archived)), add, drop=FALSE]
  }
  if (!is.null(columns)) {
    if (!is.character(columns) || anyNA(columns)) stop("`columns` must be a character vector.", call.=FALSE)
    missing <- setdiff(columns, colnames(md))
    if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse=", "), call.=FALSE)
    md <- md[, columns, drop=FALSE]
  }
  md
}

#' Internal metadata view
#' @keywords internal
#' @noRd
.gnrh_metadata <- function(object) gnrh_metadata(object)

#' Keep only essential GnRHcell columns in metadata
#'
#' @param object A Seurat object, including an existing fully annotated object.
#' @param keep Optional character vector of package-owned GnRH columns to retain.
#'   `NULL` keeps `gnrh_status`, `gnrh_class`, `gnrh_score`, `gnrh_confident`,
#'   `gnrh_stage`, `gnrh_stage_confident`, `gnrh_stage_resolution` and
#'   `gnrh_secretory`, when available. `character()` archives all owned fields.
#' @return Updated Seurat object. Detailed per-cell columns are moved to
#'   `object@misc$gnrh$cell_metadata`. User metadata, including custom `gnrh_`
#'   columns, is left unchanged. Existing diagnostics and run information remain
#'   in `misc`. No detection or staging is rerun.
#' @details This simplifies the metadata layout; moving data to `misc` does not
#'   by itself substantially reduce the object's total memory footprint.
#' @examples
#' \dontrun{
#' obj <- gnrh_compact(obj)
#' details <- gnrh_metadata(obj)
#' obj_full <- gnrh_restore(obj)
#' }
#' @export
gnrh_compact <- function(object, keep=NULL) {
  .gnrh_check_keep(keep)
  md <- gnrh_metadata(object)
  object <- .init_gnrh_misc(object)
  if (is.null(keep)) keep <- .gnrh_essential_metadata()
  owned <- intersect(.gnrh_owned_metadata(), colnames(md))
  extra <- setdiff(owned, keep)
  object@misc$gnrh$cell_metadata <- md[, extra, drop=FALSE]
  for (nm in intersect(extra, colnames(object[[]]))) object[[nm]] <- NULL
  # A custom keep selection can promote a previously archived field.
  for (nm in intersect(keep, colnames(md))) object[[nm]] <- md[[nm]]
  object@misc$gnrh$metadata_layout <- list(
    mode="compact", kept=intersect(keep, owned), archived=extra,
    schema_version=1L
  )
  object
}

#' Restore detailed GnRHcell columns to metadata
#' @param object A Seurat object with archived GnRHcell metadata.
#' @return Updated Seurat object with archived columns restored in current cell
#'   order. The per-cell archive is cleared to avoid duplicate storage. Other
#'   `misc` entries remain unchanged. Existing metadata values take precedence.
#' @seealso [gnrh_compact()], [gnrh_metadata()]
#' @export
gnrh_restore <- function(object) {
  md <- gnrh_metadata(object)
  add <- setdiff(colnames(md), colnames(object[[]]))
  for (nm in add) object[[nm]] <- md[[nm]]
  if (is.list(object@misc$gnrh)) {
    object@misc$gnrh$cell_metadata <- NULL
    object@misc$gnrh$metadata_layout <- list(mode="full", schema_version=1L)
  }
  object
}

#' Restore only columns needed by a Seurat plotting operation
#' @keywords internal
#' @noRd
.gnrh_restore_columns <- function(object, columns) {
  md <- gnrh_metadata(object)
  add <- intersect(setdiff(columns, colnames(object[[]])), colnames(md))
  for (nm in add) object[[nm]] <- md[[nm]]
  object
}
