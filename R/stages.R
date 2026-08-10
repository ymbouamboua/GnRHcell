#' Assign dominant developmental stage
#'
#' Internal helper that assigns each cell to the developmental stage with the
#' highest module score.
#'
#' Developmental stages are restricted to \code{identity},
#' \code{migrating}, and \code{mature}. Migration assignments are additionally
#' validated using migration-core marker evidence.
#'
#' @param scores Numeric matrix or data frame of per-cell developmental stage
#'   scores, with stages in columns.
#' @param migration_core_hits Optional integer vector containing the number of
#'   expressed migration-core markers per cell.
#' @param min_migration_hits Minimum number of migration-core hits required to
#'   retain a raw \code{migrating} assignment. Default is \code{2L}.
#'
#' @return Character vector containing developmental stage labels.
#'
#' @keywords internal
#' @noRd
assign_stage <- function(
    scores,
    migration_core_hits = NULL,
    min_migration_hits = 2L
) {

  scores <- as.matrix(scores)

  required_stages <- c(
    "identity",
    "migrating",
    "mature"
  )

  missing_stages <- setdiff(
    required_stages,
    colnames(scores)
  )

  if (length(missing_stages) > 0L) {
    stop(
      "Missing developmental stage score(s): ",
      paste(
        missing_stages,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  scores <- scores[
    ,
    required_stages,
    drop = FALSE
  ]

  # --------------------------------------------------------------------------- #
  # Initial assignment
  # --------------------------------------------------------------------------- #

  stage <- colnames(scores)[
    max.col(
      scores,
      ties.method = "first"
    )
  ]

  if (is.null(migration_core_hits)) {
    return(stage)
  }

  if (
    length(migration_core_hits) !=
    nrow(scores)
  ) {
    stop(
      "`migration_core_hits` must have one value per cell.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------- #
  # Validate migrating assignments
  # --------------------------------------------------------------------------- #

  weak_migration <-
    stage == "migrating" &
    migration_core_hits < min_migration_hits

  if (any(weak_migration)) {

    alternatives <- scores[
      weak_migration,
      c(
        "identity",
        "mature"
      ),
      drop = FALSE
    ]

    stage[
      weak_migration
    ] <- colnames(alternatives)[
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
#' Assign developmental states to GnRH-lineage cells using biologically
#' informed transcriptional programs.
#'
#' Developmental staging is based on three mutually exclusive states:
#' \itemize{
#'   \item \code{identity}: lineage specification and early GnRH identity;
#'   \item \code{migrating}: migration and axon-guidance programs; and
#'   \item \code{mature}: neuroendocrine maturation.
#' }
#'
#' Secretory machinery is evaluated independently from developmental stage.
#' This allows cells to retain a developmental annotation such as
#' \code{migrating} or \code{mature} while independently receiving evidence
#' for a neuroendocrine secretory transcriptional program.
#'
#' A raw developmental stage is first assigned from the highest developmental
#' module score. Raw \code{migrating} assignments are then validated using a
#' curated migration-core marker set. Cells lacking sufficient migration-core
#' evidence are reassigned to the highest-scoring alternative developmental
#' stage.
#'
#' Secretory support requires expression of at least one core secretory marker
#' together with either an additional core marker or sufficient supportive
#' secretory evidence. The resulting annotation is classified as
#' \code{limited} or \code{supported}.
#'
#' Cells classified as GnRH-negative by \code{\link{detect_gnrh}} are labeled
#' \code{non-gnrh} in both developmental-stage and secretory annotations.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay used for developmental and secretory module scoring.
#'   Default is \code{"RNA"}.
#' @param layer Expression layer used for module scoring.
#'   Default is \code{"data"}.
#' @param min_migration_hits Minimum number of expressed migration-core
#'   markers required to retain a raw \code{migrating} assignment.
#'   Default is \code{2L}.
#' @param min_secretory_core_hits Minimum number of core secretory markers
#'   required before secretory support can be assigned.
#'   Default is \code{1L}.
#' @param min_secretory_supportive_hits Minimum number of supportive secretory
#'   markers required when only one core secretory marker is detected.
#'   Default is \code{2L}.
#' @param verbose Logical. Retained for API consistency. Progress reporting is
#'   normally handled by \code{\link{run_gnrh}}.
#'
#' @return A Seurat object containing developmental stage assignments,
#' developmental scores, migration-core evidence, and independent secretory
#' transcriptional support.
#'
#' Added metadata include:
#' \describe{
#'   \item{\code{gnrh_stage_raw}}{
#'   Developmental stage assigned directly from the maximum developmental
#'   module score.}
#'   \item{\code{gnrh_stage}}{
#'   Final developmental stage after migration-core validation and masking
#'   of GnRH-negative cells.}
#'   \item{\code{gnrh_stage_reassigned}}{
#'   Logical indicator specifying whether the raw developmental stage was
#'   reassigned.}
#'   \item{\code{gnrh_stage_reason}}{
#'   Reason for the final developmental-stage assignment.}
#'   \item{\code{gnrh_migration_core_hits}}{
#'   Number of expressed migration-core markers detected per cell.}
#'   \item{\code{gnrh_secretory_core_hits}}{
#'   Number of expressed core secretory markers detected per cell.}
#'   \item{\code{gnrh_secretory_supportive_hits}}{
#'   Number of expressed supportive secretory markers detected per cell.}
#'   \item{\code{gnrh_secretory_hits}}{
#'   Total number of core and supportive secretory markers detected.}
#'   \item{\code{gnrh_secretory}}{
#'   Independent transcriptional support for secretory machinery, classified
#'   as \code{limited}, \code{supported}, or \code{non-gnrh}.}
#' }
#'
#' @details
#' Developmental stage and secretory support are intentionally modeled as
#' separate dimensions. Secretory-program expression therefore does not
#' replace or override the developmental-stage assignment.
#'
#' A raw \code{migrating} assignment is retained only when at least
#' \code{min_migration_hits} migration-core markers are detected. Otherwise,
#' the cell is reassigned to the highest-scoring alternative developmental
#' state among \code{identity} and \code{mature}.
#'
#' Secretory support is assigned when the cell expresses at least
#' \code{min_secretory_core_hits} core secretory markers and either:
#' \itemize{
#'   \item at least two core secretory markers; or
#'   \item at least \code{min_secretory_supportive_hits} supportive secretory
#'   markers.
#' }
#'
#' The secretory annotation reflects transcriptional support for
#' neuroendocrine secretory machinery and should not be interpreted as a
#' direct measurement of GnRH peptide release.
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
    min_migration_hits = 2L,
    min_secretory_core_hits = 1L,
    min_secretory_supportive_hits = 2L,
    verbose = TRUE
) {

  stopifnot(
    inherits(
      object,
      "Seurat"
    )
  )

  # --------------------------------------------------------------------------- #
  # Validate arguments
  # --------------------------------------------------------------------------- #

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

  if (
    length(min_secretory_core_hits) != 1L ||
    is.na(min_secretory_core_hits) ||
    min_secretory_core_hits < 0
  ) {
    stop(
      "`min_secretory_core_hits` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  if (
    length(min_secretory_supportive_hits) != 1L ||
    is.na(min_secretory_supportive_hits) ||
    min_secretory_supportive_hits < 0
  ) {
    stop(
      "`min_secretory_supportive_hits` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  min_migration_hits <- as.integer(
    min_migration_hits
  )

  min_secretory_core_hits <- as.integer(
    min_secretory_core_hits
  )

  min_secretory_supportive_hits <- as.integer(
    min_secretory_supportive_hits
  )

  # --------------------------------------------------------------------------- #
  # Expression matrix
  # --------------------------------------------------------------------------- #

  expr <- .get_expr(
    object,
    assay = assay,
    layer = layer
  )

  genes <- rownames(
    expr
  )

  # --------------------------------------------------------------------------- #
  # Developmental modules
  # --------------------------------------------------------------------------- #

  modules <- .build_stage_modules(
    genes
  )

  required_stages <- c(
    "identity",
    "migrating",
    "mature"
  )

  missing_modules <- setdiff(
    required_stages,
    names(modules)
  )

  if (length(missing_modules) > 0L) {
    stop(
      "Developmental stage modules are incomplete. Missing: ",
      paste(
        missing_modules,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  developmental_modules <- modules[
    required_stages
  ]

  mod <- .score_modules(
    expr = expr,
    modules = developmental_modules
  )

  scores <- as.data.frame(
    mod$score
  )

  missing_scores <- setdiff(
    required_stages,
    colnames(scores)
  )

  if (length(missing_scores) > 0L) {
    stop(
      "Developmental module scores are incomplete. Missing: ",
      paste(
        missing_scores,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------- #
  # Raw developmental-stage assignment
  # --------------------------------------------------------------------------- #

  raw_stage <- colnames(scores)[
    max.col(
      as.matrix(
        scores[
          ,
          required_stages,
          drop = FALSE
        ]
      ),
      ties.method = "first"
    )
  ]

  # --------------------------------------------------------------------------- #
  # Migration-core evidence
  # --------------------------------------------------------------------------- #

  migration_core <- .build_migration_core(
    genes
  )

  migration_core_hits <- if (
    length(migration_core) > 0L
  ) {

    as.integer(
      Matrix::colSums(
        expr[
          migration_core,
          ,
          drop = FALSE
        ] > 0
      )
    )

  } else {

    integer(
      ncol(expr)
    )
  }

  # --------------------------------------------------------------------------- #
  # Refined developmental-stage assignment
  # --------------------------------------------------------------------------- #

  stage <- assign_stage(
    scores = scores,
    migration_core_hits = migration_core_hits,
    min_migration_hits = min_migration_hits
  )

  stage_reassigned <-
    raw_stage != stage

  stage_reason <- rep(
    "max_score",
    length(stage)
  )

  migration_filtered <-
    raw_stage == "migrating" &
    migration_core_hits < min_migration_hits

  stage_reason[
    migration_filtered
  ] <- "migration_core_filter"

  # --------------------------------------------------------------------------- #
  # GnRH-positive mask
  # --------------------------------------------------------------------------- #

  positive <- if (
    "gnrh_status" %in%
    colnames(
      object[[]]
    )
  ) {

    as.character(
      object$gnrh_status
    ) == "pos"

  } else {

    rep(
      TRUE,
      ncol(object)
    )
  }

  negative <- !positive

  # --------------------------------------------------------------------------- #
  # Mask GnRH-negative cells
  # --------------------------------------------------------------------------- #

  stage[
    negative
  ] <- "non-gnrh"

  stage_reassigned[
    negative
  ] <- FALSE

  stage_reason[
    negative
  ] <- "non_gnrh"

  # --------------------------------------------------------------------------- #
  # Secretory module
  #
  # Secretory activity is evaluated independently from developmental stage.
  # --------------------------------------------------------------------------- #

  secretory_module <- .build_secretory_module(
    genes
  )

  if (
    !all(
      c(
        "core",
        "supportive"
      ) %in%
      names(secretory_module)
    )
  ) {
    stop(
      "Secretory module must contain `core` and `supportive` gene sets.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------- #
  # Secretory core hits
  # --------------------------------------------------------------------------- #

  secretory_core_hits <- if (
    length(secretory_module$core) > 0L
  ) {

    as.integer(
      Matrix::colSums(
        expr[
          secretory_module$core,
          ,
          drop = FALSE
        ] > 0
      )
    )

  } else {

    integer(
      ncol(expr)
    )
  }

  # --------------------------------------------------------------------------- #
  # Secretory supportive hits
  # --------------------------------------------------------------------------- #

  secretory_supportive_hits <- if (
    length(secretory_module$supportive) > 0L
  ) {

    as.integer(
      Matrix::colSums(
        expr[
          secretory_module$supportive,
          ,
          drop = FALSE
        ] > 0
      )
    )

  } else {

    integer(
      ncol(expr)
    )
  }

  # --------------------------------------------------------------------------- #
  # Secretory total hits
  # --------------------------------------------------------------------------- #

  secretory_hits <-
    secretory_core_hits +
    secretory_supportive_hits

  # --------------------------------------------------------------------------- #
  # Secretory activity
  #
  # High secretory activity requires:
  #
  # - at least min_secretory_core_hits core markers; and
  # - either a second core marker or sufficient supportive evidence.
  #
  # With the defaults, this means:
  #
  # core >= 1 AND (core >= 2 OR supportive >= 2)
  # --------------------------------------------------------------------------- #

  secretory_supported <-
    positive &
    secretory_core_hits >= min_secretory_core_hits &
    (
      secretory_core_hits >=
        max(
          2L,
          min_secretory_core_hits
        ) |
        secretory_supportive_hits >=
        min_secretory_supportive_hits
    )

  # --------------------------------------------------------------------------- #
  # Store developmental-stage assignments
  # --------------------------------------------------------------------------- #

  object$gnrh_stage_raw <- factor(
    raw_stage,
    levels = required_stages
  )

  object$gnrh_stage <- factor(
    stage,
    levels = c(
      required_stages,
      "non-gnrh"
    )
  )

  object$gnrh_stage_reassigned <-
    as.logical(
      stage_reassigned
    )

  object$gnrh_stage_reason <- factor(
    stage_reason,
    levels = c(
      "max_score",
      "migration_core_filter",
      "non_gnrh"
    )
  )

  # --------------------------------------------------------------------------- #
  # Store developmental scores
  # --------------------------------------------------------------------------- #

  for (nm in required_stages) {

    object[[
      paste0(
        "gnrh_",
        nm,
        "_score"
      )
    ]] <- scores[[nm]]
  }

  # --------------------------------------------------------------------------- #
  # Store migration evidence
  # --------------------------------------------------------------------------- #

  object$gnrh_migration_core_hits <-
    migration_core_hits

  # --------------------------------------------------------------------------- #
  # Store secretory evidence
  # --------------------------------------------------------------------------- #

  object$gnrh_secretory_core_hits <-
    secretory_core_hits

  object$gnrh_secretory_supportive_hits <-
    secretory_supportive_hits

  object$gnrh_secretory_hits <-
    secretory_hits

  object$gnrh_secretory <- factor(
    ifelse(
      negative,
      "non-gnrh",
      ifelse(
        secretory_supported,
        "supported",
        "limited"
      )
    ),
    levels = c(
      "limited",
      "supported",
      "non-gnrh"
    )
  )

  # --------------------------------------------------------------------------- #
  # Store staging metadata
  # --------------------------------------------------------------------------- #

  object@misc$gnrh_stage_modules <-
    developmental_modules

  object@misc$gnrh_migration_core <-
    migration_core

  object@misc$gnrh_secretory_module <-
    secretory_module

  object@misc$gnrh_stage_parameters <- list(
    assay = assay,
    layer = layer,
    min_migration_hits = min_migration_hits,
    min_secretory_core_hits =
      min_secretory_core_hits,
    min_secretory_supportive_hits =
      min_secretory_supportive_hits
  )

  object
}
