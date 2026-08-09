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
assign_stage <- function(
    scores,
    migration_core_hits = NULL,
    min_migration_hits = 1L
) {

  scores <- as.matrix(scores)

  stage <- colnames(scores)[
    max.col(
      scores,
      ties.method = "first"
    )
  ]

  if (is.null(migration_core_hits)) {
    return(stage)
  }

  if (length(migration_core_hits) != nrow(scores)) {
    stop(
      "`migration_core_hits` must have one value per cell.",
      call. = FALSE
    )
  }

  weak_migration <-
    stage == "migrating" &
    migration_core_hits < min_migration_hits

  if (any(weak_migration)) {

    alternatives <- scores[
      weak_migration,
      c(
        "identity",
        "mature",
        "secreting"
      ),
      drop = FALSE
    ]

    stage[weak_migration] <-
      colnames(alternatives)[
        max.col(
          alternatives,
          ties.method = "first"
        )
      ]
  }

  stage
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
#' A raw developmental stage is first assigned from the highest
#' module score. Migration assignments are then validated using
#' migration-specific marker evidence. Cells initially assigned as
#' \code{migrating} but lacking the required number of migration-core
#' marker hits are reassigned to the highest-scoring alternative stage.
#'
#' If GnRH classification metadata are present, non-GnRH cells are
#' labeled as \code{non-gnrh} in the final stage assignment.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay used for expression extraction.
#' Default is \code{"RNA"}.
#' @param layer Expression layer used for developmental module scoring.
#' Default is \code{"data"}.
#' @param min_migration_hits Minimum number of expressed migration-core
#' markers required to retain a raw \code{migrating} assignment.
#' Default is \code{1L}.
#' @param verbose Logical; print progress messages.
#' Default is \code{TRUE}.
#'
#' @return A Seurat object updated with:
#' \describe{
#'   \item{\code{gnrh_stage_raw}}{
#'   Developmental stage assigned directly from the maximum module score.
#'   }
#'   \item{\code{gnrh_stage}}{
#'   Final developmental stage after migration validation and masking
#'   of non-GnRH cells.
#'   }
#'   \item{\code{gnrh_stage_reassigned}}{
#'   Logical indicator specifying whether the raw developmental stage
#'   was reassigned during migration validation.
#'   }
#'   \item{\code{gnrh_stage_reason}}{
#'   Reason for the final stage assignment.
#'   }
#'   \item{\code{gnrh_migration_core_hits}}{
#'   Number of expressed migration-core markers detected per cell.
#'   }
#'   \item{\code{gnrh_identity_score}}{
#'   Identity module score.
#'   }
#'   \item{\code{gnrh_migrating_score}}{
#'   Migration module score.
#'   }
#'   \item{\code{gnrh_mature_score}}{
#'   Mature neuroendocrine module score.
#'   }
#'   \item{\code{gnrh_secreting_score}}{
#'   Secretory activity module score.
#'   }
#' }
#'
#' Developmental module definitions, migration-core markers, and
#' staging parameters are stored in:
#' \itemize{
#'   \item \code{object@misc$gnrh_stage_modules}
#'   \item \code{object@misc$gnrh_migration_core}
#'   \item \code{object@misc$gnrh_stage_parameters}
#' }
#'
#' @details
#' Developmental staging is based on predefined transcriptional
#' programs reflecting known biological states of GnRH neuron
#' development.
#'
#' Raw stage assignments are obtained from the highest developmental
#' module score. Because general neuronal and axon-guidance genes can
#' produce elevated migration scores in mature neurons, raw
#' \code{migrating} assignments require additional migration-core
#' evidence.
#'
#' Cells failing this migration criterion are reassigned to the
#' highest-scoring stage among \code{identity}, \code{mature}, and
#' \code{secreting}.
#'
#' If \code{gnrh_status} metadata are present from
#' \code{\link{detect_gnrh}}, cells classified as negative are labeled
#' \code{non-gnrh} in the final stage assignment.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{run_gnrh}}
#'
#' @export
stage_gnrh <- function(
    object,
    assay = "RNA",
    layer = "data",
    min_migration_hits = 1L,
    verbose = TRUE
) {

  stopifnot(
    inherits(object, "Seurat")
  )

  if (
    length(min_migration_hits) != 1L ||
    is.na(min_migration_hits) ||
    min_migration_hits < 0
  ) {
    stop(
      "`min_migration_hits` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  min_migration_hits <- as.integer(min_migration_hits)

  log <- .msg(verbose)

  log("==== GNRH STAGING START ====")

  # ---------------------------------------------------------------------------
  # Expression for module scoring
  # ---------------------------------------------------------------------------

  expr <- .get_expr(
    object,
    assay = assay,
    layer = layer
  )

  modules <- .build_stage_modules(
    rownames(expr)
  )

  mod <- .score_modules(
    expr = expr,
    modules = modules
  )

  scores <- mod$score

  if (!all(
    c(
      "identity",
      "migrating",
      "mature",
      "secreting"
    ) %in% colnames(scores)
  )) {
    stop(
      "Developmental module scores are incomplete.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Raw stage assignment
  # ---------------------------------------------------------------------------

  raw_stage <- colnames(scores)[
    max.col(
      as.matrix(scores),
      ties.method = "first"
    )
  ]


  # ---------------------------------------------------------------------------
  # Migration-specific evidence
  # ---------------------------------------------------------------------------

  migration_core <- .build_migration_core(
    rownames(expr)
  )

  if (length(migration_core) > 0L) {

    migration_core_hits <- Matrix::colSums(
      expr[
        migration_core,
        ,
        drop = FALSE
      ] > 0
    )

  } else {

    migration_core_hits <- rep(
      0L,
      ncol(expr)
    )
  }

  migration_core_hits <- as.integer(
    migration_core_hits
  )


  # ---------------------------------------------------------------------------
  # Final developmental stage
  # ---------------------------------------------------------------------------

  stage <- assign_stage(
    scores = scores,
    migration_core_hits = migration_core_hits,
    min_migration_hits = min_migration_hits
  )


  # ---------------------------------------------------------------------------
  # Track reassignment
  # ---------------------------------------------------------------------------

  stage_reassigned <- raw_stage != stage

  stage_reason <- rep(
    "max_score",
    length(stage)
  )

  stage_reason[
    raw_stage == "migrating" &
      migration_core_hits < min_migration_hits
  ] <- "migration_core_filter"


  # ---------------------------------------------------------------------------
  # Mask non-GnRH cells
  # ---------------------------------------------------------------------------

  if ("gnrh_status" %in% colnames(object[[]])) {

    neg <- as.character(object$gnrh_status) == "neg"

    stage[neg] <- "non-gnrh"

    # Masking is not considered developmental-stage reassignment.
    stage_reassigned[neg] <- FALSE

    stage_reason[neg] <- "non_gnrh"
  }


  # ---------------------------------------------------------------------------
  # Store stage assignments
  # ---------------------------------------------------------------------------

  object$gnrh_stage_raw <- factor(
    raw_stage,
    levels = c(
      "identity",
      "migrating",
      "mature",
      "secreting"
    )
  )

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

  object$gnrh_stage_reassigned <-
    as.logical(stage_reassigned)

  object$gnrh_stage_reason <- factor(
    stage_reason,
    levels = c(
      "max_score",
      "migration_core_filter",
      "non_gnrh"
    )
  )


  # ---------------------------------------------------------------------------
  # Store module scores
  # ---------------------------------------------------------------------------

  for (nm in colnames(scores)) {

    object[[
      paste0(
        "gnrh_",
        nm,
        "_score"
      )
    ]] <- scores[[nm]]
  }


  # ---------------------------------------------------------------------------
  # Store migration evidence
  # ---------------------------------------------------------------------------

  object$gnrh_migration_core_hits <-
    migration_core_hits


  # ---------------------------------------------------------------------------
  # Store staging metadata
  # ---------------------------------------------------------------------------

  object@misc$gnrh_stage_modules <-
    modules

  object@misc$gnrh_migration_core <-
    migration_core

  object@misc$gnrh_stage_parameters <- list(
    assay = assay,
    layer = layer,
    min_migration_hits = min_migration_hits
  )


  # ---------------------------------------------------------------------------
  # Report
  # ---------------------------------------------------------------------------

  if (verbose) {

    positive <- if (
      "gnrh_status" %in% colnames(object[[]])
    ) {
      as.character(object$gnrh_status) != "neg"
    } else {
      rep(TRUE, ncol(object))
    }

    raw_tab <- table(
      object$gnrh_stage_raw[positive],
      useNA = "no"
    )

    final_tab <- table(
      object$gnrh_stage[positive],
      useNA = "no"
    )

    n_reassigned <- sum(
      object$gnrh_stage_reassigned[positive],
      na.rm = TRUE
    )

    log("Raw stage assignment:")

    for (nm in names(raw_tab)) {
      log(
        sprintf(
          "  %-10s %d",
          nm,
          raw_tab[[nm]]
        )
      )
    }

    log("Final stage assignment:")

    for (nm in names(final_tab)) {
      log(
        sprintf(
          "  %-10s %d",
          nm,
          final_tab[[nm]]
        )
      )
    }

    log(
      sprintf(
        "Migration-filter reassigned: %d",
        n_reassigned
      )
    )
  }

  log("==== GNRH STAGING DONE ====")

  object
}
