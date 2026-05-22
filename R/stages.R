#' Assign dominant developmental stage
#'
#' Internal helper that assigns each cell to the developmental
#' stage with the highest module score.
#'
#' Ties are resolved by selecting the first maximum.
#'
#' @param scores Numeric matrix or data frame of per-cell stage
#' module scores, with stages in columns.
#'
#' @return Character vector of assigned stage labels.
#'
#' @keywords internal
#' @noRd
assign_stage <- function(scores) {

  colnames(scores)[
    max.col(
      scores,
      ties.method = "first"
    )
  ]
}


#' Stage GnRH lineage cells
#'
#' Assign developmental states to GnRH-lineage cells using
#' biologically informed transcriptional programs.
#'
#' Cells are scored against predefined developmental modules
#' representing major GnRH neuron states:
#' \itemize{
#'   \item \code{identity}: lineage specification and early GnRH identity
#'   \item \code{migrating}: migration and axon-guidance programs
#'   \item \code{mature}: neuroendocrine maturation
#'   \item \code{secreting}: secretory and vesicle machinery activation
#' }
#'
#' Each cell is assigned to the stage with the highest module score.
#' If GnRH classification metadata are present, non-GnRH cells are
#' labeled as \code{non-gnrh}.
#'
#' Stage-specific module scores are stored in object metadata.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay used for expression extraction.
#' Default is \code{"RNA"}.
#' @param layer Expression layer used for staging.
#' Default is \code{"counts"}.
#' @param verbose Logical; print progress messages.
#' Default is \code{TRUE}.
#'
#' @return A Seurat object updated with:
#' \describe{
#'   \item{\code{gnrh_stage}}{Assigned developmental stage labels.}
#'   \item{\code{gnrh_identity_score}}{Identity module score.}
#'   \item{\code{gnrh_migrating_score}}{Migration module score.}
#'   \item{\code{gnrh_mature_score}}{Mature neuroendocrine score.}
#'   \item{\code{gnrh_secreting_score}}{Secretory activity score.}
#' }
#'
#' Developmental module definitions are stored in:
#' \code{object@misc$gnrh_stage_modules}
#'
#' @details
#' Developmental staging is based on predefined marker modules
#' reflecting known biological programs of GnRH neuron development.
#'
#' If \code{gnrh_class} metadata are present from
#' \code{\link{detect_gnrh}}, cells classified as negative are
#' reassigned to \code{non-gnrh}.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{run_gnrh}}
#'
#' @export
stage_gnrh <- function(
    object,
    assay = "RNA",
    layer = "counts",
    verbose = TRUE
) {

  stopifnot(
    inherits(object, "Seurat")
  )

  log <- .msg(verbose)

  log("==== GNRH STAGING START ====")

  expr <- .get_expr(
    object,
    assay = assay,
    layer = layer
  )

  modules <- build_stage_modules(
    rownames(expr)
  )

  mod <- .score_modules(
    expr = expr,
    modules = modules
  )

  scores <- mod$score

  stage <- assign_stage(scores)

  if ("gnrh_class" %in% colnames(object[[]])) {

    neg <- object$gnrh_class %in% c("neg")
    stage[neg] <- "non-gnrh"
  }

  object$gnrh_stage <- factor(
    stage,
    levels = c(
      "identity",
      "migrating",
      "mature",
      "secreting",
      "non-gnrh"
    )
  )

  # store scores
  for (nm in colnames(scores)) {

    object[[paste0(
      "gnrh_",
      nm,
      "_score"
    )]] <- scores[[nm]]
  }

  object@misc$gnrh_stage_modules <- modules

  log(
    paste(
      names(table(object$gnrh_stage)),
      table(object$gnrh_stage),
      collapse = " | "
    )
  )

  log("==== GNRH STAGING DONE ====")

  object
}
