#' Assign dominant developmental stage
#'
#' Internal helper assigning the dominant developmental state from curated
#' GnRH developmental-program scores.
#'
#' Raw assignments are based on the highest score among \code{identity},
#' \code{migrating}, and \code{mature}. Migrating assignments can optionally
#' be rejected when insufficient migration-core evidence is present.
#'
#' @param scores Numeric matrix or data frame containing developmental-stage
#'   scores in columns.
#' @param migration_core_hits Optional integer vector containing migration-core
#'   marker hits.
#' @param min_migration_hits Minimum migration-core hits required to retain a
#'   migrating assignment.
#'
#' @return A list containing final stage, raw stage, dominant score,
#'   second-best score, score margin, and migration-filter status.
#'
#' @keywords internal
#' @noRd
assign_stage <- function(
    scores,
    migration_core_hits = NULL,
    min_migration_hits = 2L
) {
  scores <- as.matrix(scores)

  required <- c(
    "identity",
    "migrating",
    "mature"
  )

  missing <- setdiff(
    required,
    colnames(scores)
  )

  if (length(missing)) {
    stop(
      "Missing developmental-stage score(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  scores <- scores[
    ,
    required,
    drop = FALSE
  ]

  if (!nrow(scores)) {
    return(
      list(
        stage = character(0),
        raw_stage = character(0),
        top_score = numeric(0),
        second_score = numeric(0),
        margin = numeric(0),
        migration_filtered = logical(0)
      )
    )
  }

  finite <- is.finite(scores)

  scores_work <- scores
  scores_work[!finite] <- -Inf

  no_information <- rowSums(finite) == 0L

  top_index <- max.col(
    scores_work,
    ties.method = "first"
  )

  raw_stage <- colnames(scores_work)[
    top_index
  ]

  raw_stage[
    no_information
  ] <- NA_character_

  top_score <- scores_work[
    cbind(
      seq_len(nrow(scores_work)),
      top_index
    )
  ]

  second_score <- vapply(
    seq_len(nrow(scores_work)),
    function(i) {
      x <- sort(
        scores_work[i, ],
        decreasing = TRUE
      )

      if (
        length(x) < 2L ||
        !is.finite(x[[2L]])
      ) {
        return(NA_real_)
      }

      x[[2L]]
    },
    numeric(1)
  )

  top_score[
    no_information
  ] <- NA_real_

  margin <- top_score - second_score

  stage <- raw_stage

  migration_filtered <- rep(
    FALSE,
    nrow(scores_work)
  )

  if (!is.null(migration_core_hits)) {
    if (
      length(migration_core_hits) !=
      nrow(scores_work)
    ) {
      stop(
        "`migration_core_hits` must contain one value per cell.",
        call. = FALSE
      )
    }

    migration_core_hits <- suppressWarnings(
      as.numeric(
        migration_core_hits
      )
    )

    migration_core_hits[
      !is.finite(migration_core_hits)
    ] <- 0

    migration_filtered <-
      !is.na(stage) &
      stage == "migrating" &
      migration_core_hits <
      min_migration_hits

    if (any(migration_filtered)) {
      alternative <- scores_work[
        migration_filtered,
        c(
          "identity",
          "mature"
        ),
        drop = FALSE
      ]

      alternative_index <- max.col(
        alternative,
        ties.method = "first"
      )

      stage[
        migration_filtered
      ] <- colnames(alternative)[
        alternative_index
      ]
    }
  }

  list(
    stage = stage,
    raw_stage = raw_stage,
    top_score = top_score,
    second_score = second_score,
    margin = margin,
    migration_filtered = migration_filtered
  )
}


#' Stage GnRH lineage cells
#'
#' Assigns developmental states to GnRH-positive cells using curated
#' transcriptional programs representing lineage identity, migration, and
#' neuroendocrine maturation.
#'
#' Developmental stage and secretory machinery are modeled independently.
#'
#' @param object A Seurat object previously processed with
#'   \code{\link{detect_gnrh}}.
#' @param assay Assay used for developmental and secretory scoring.
#' @param layer Expression layer used for module scoring.
#' @param min_migration_hits Minimum migration-core hits required to retain a
#'   migrating assignment.
#' @param min_secretory_core_hits Minimum core secretory marker hits.
#' @param min_secretory_supportive_hits Minimum supportive secretory hits.
#' @param expression_weight Weight assigned to normalized module expression.
#' @param detection_weight Weight assigned to module detection fraction.
#' @param stage_margin Minimum difference between best and second-best
#'   developmental score required for \code{gnrh_stage_confident = TRUE}.
#' @param verbose Print progress messages.
#'
#' @return A Seurat object containing developmental-stage scores,
#' assignments, confidence metrics, migration evidence, and secretory support.
#'
#' @details
#' Developmental stage is derived from integrated module activity combining
#' normalized expression and the fraction of detected module genes.
#'
#' The continuous developmental-program scores are retained independently from
#' the categorical dominant-stage assignment.
#'
#' \code{gnrh_stage_confident} indicates whether the dominant program exceeds
#' the second-best program by at least \code{stage_margin}.
#'
#' \code{gnrh_stage_resolution} stores the corresponding categorical
#' interpretation: \code{"resolved"}, \code{"transitional"}, or
#' \code{"non-gnrh"}.
#'
#' Developmental scores are stored using the \code{gnrh_stage_*} namespace and
#' therefore do not overwrite detection scores such as
#' \code{gnrh_identity_score}.
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
    expression_weight = 0.5,
    detection_weight = 0.5,
    stage_margin = 0.10,
    verbose = TRUE
) {
  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  log <- .msg(verbose)

  validate_integer <- function(x, name) {
    if (
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      x != floor(x)
    ) {
      stop(
        "`",
        name,
        "` must be a single non-negative integer.",
        call. = FALSE
      )
    }

    as.integer(x)
  }

  min_migration_hits <- validate_integer(
    min_migration_hits,
    "min_migration_hits"
  )

  min_secretory_core_hits <- validate_integer(
    min_secretory_core_hits,
    "min_secretory_core_hits"
  )

  min_secretory_supportive_hits <- validate_integer(
    min_secretory_supportive_hits,
    "min_secretory_supportive_hits"
  )

  if (
    length(stage_margin) != 1L ||
    !is.finite(stage_margin) ||
    stage_margin < 0
  ) {
    stop(
      "`stage_margin` must be a single non-negative number.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------------- #
  # Input
  # ------------------------------------------------------------------------- #

  object <- validate_input(
    object,
    assay = assay,
    required_layers = layer,
    auto_normalize = FALSE,
    verbose = FALSE
  )

  md <- object[[]]

  if (!"gnrh_status" %in% colnames(md)) {
    stop(
      "`gnrh_status` is missing. Run `detect_gnrh()` first.",
      call. = FALSE
    )
  }

  positive <-
    !is.na(md$gnrh_status) &
    as.character(md$gnrh_status) == "pos"

  negative <- !positive

  log(
    sprintf(
      "Staging %s GnRH-positive cells",
      format(
        sum(positive),
        big.mark = ","
      )
    )
  )

  # ------------------------------------------------------------------------- #
  # Expression
  # ------------------------------------------------------------------------- #

  expr <- .get_expr(
    object,
    assay = assay,
    layer = layer
  )

  genes <- rownames(expr)

  # ------------------------------------------------------------------------- #
  # Developmental programs
  # ------------------------------------------------------------------------- #

  modules <- .build_stage_modules(
    genes
  )

  required_stages <- c(
    "identity",
    "migrating",
    "mature"
  )

  missing <- setdiff(
    required_stages,
    names(modules)
  )

  if (length(missing)) {
    stop(
      "Developmental-stage modules are incomplete. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  developmental_modules <- modules[
    required_stages
  ]

  mod <- .score_modules(
    expr = expr,
    modules = developmental_modules,
    expression_weight = expression_weight,
    detection_weight = detection_weight
  )

  expression_scores <- as.data.frame(
    mod$score,
    check.names = FALSE
  )

  detection_fractions <- as.data.frame(
    mod$fraction,
    check.names = FALSE
  )

  stage_scores <- as.data.frame(
    mod$integrated,
    check.names = FALSE
  )

  # ------------------------------------------------------------------------- #
  # Migration-core evidence
  # ------------------------------------------------------------------------- #

  migration_core <- .build_migration_core(
    genes
  )

  migration_core_hits <- if (
    length(migration_core)
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

  # ------------------------------------------------------------------------- #
  # Stage assignment
  # ------------------------------------------------------------------------- #

  assignment <- assign_stage(
    scores = stage_scores,
    migration_core_hits = migration_core_hits,
    min_migration_hits = min_migration_hits
  )

  raw_stage <- assignment$raw_stage
  stage <- assignment$stage

  stage_reassigned <-
    !is.na(raw_stage) &
    !is.na(stage) &
    raw_stage != stage

  stage_reason <- rep(
    "max_score",
    length(stage)
  )

  stage_reason[
    is.na(raw_stage)
  ] <- "insufficient_stage_signal"

  stage_reason[
    assignment$migration_filtered
  ] <- "migration_core_filter"

  # ------------------------------------------------------------------------- #
  # Stage confidence
  # ------------------------------------------------------------------------- #

  stage_confident <-
    positive &
    !is.na(stage) &
    is.finite(assignment$margin) &
    assignment$margin >= stage_margin

  stage_resolution <- ifelse(
    negative,
    "non-gnrh",
    ifelse(
      stage_confident,
      "resolved",
      "transitional"
    )
  )

  # ------------------------------------------------------------------------- #
  # Mask negative cells
  # ------------------------------------------------------------------------- #

  stage[
    negative
  ] <- "non-gnrh"

  stage_reassigned[
    negative
  ] <- FALSE

  stage_confident[
    negative
  ] <- FALSE

  stage_reason[
    negative
  ] <- "non_gnrh"

  # ------------------------------------------------------------------------- #
  # Secretory program
  # ------------------------------------------------------------------------- #

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

  secretory_core_hits <- if (
    length(secretory_module$core)
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

  secretory_supportive_hits <- if (
    length(secretory_module$supportive)
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

  secretory_hits <-
    secretory_core_hits +
    secretory_supportive_hits

  secretory_supported <-
    positive &
    secretory_core_hits >=
    min_secretory_core_hits &
    (
      secretory_core_hits >=
        max(
          2L,
          min_secretory_core_hits
        ) |
        secretory_supportive_hits >=
        min_secretory_supportive_hits
    )

  secretory_state <- ifelse(
    negative,
    "non-gnrh",
    ifelse(
      secretory_supported,
      "supported",
      "limited"
    )
  )

  # ------------------------------------------------------------------------- #
  # Metadata: developmental stage
  # ------------------------------------------------------------------------- #

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
      "insufficient_stage_signal",
      "non_gnrh"
    )
  )

  object$gnrh_stage_score <-
    assignment$top_score

  object$gnrh_stage_second_score <-
    assignment$second_score

  object$gnrh_stage_margin <-
    assignment$margin

  object$gnrh_stage_confident <-
    stage_confident

  object$gnrh_stage_resolution <- factor(
    stage_resolution,
    levels = c(
      "resolved",
      "transitional",
      "non-gnrh"
    )
  )

  # ------------------------------------------------------------------------- #
  # Metadata: continuous developmental programs
  # ------------------------------------------------------------------------- #

  for (nm in required_stages) {
    object[[
      paste0(
        "gnrh_stage_",
        nm,
        "_score"
      )
    ]] <- stage_scores[[nm]]

    object[[
      paste0(
        "gnrh_stage_",
        nm,
        "_expression"
      )
    ]] <- expression_scores[[nm]]

    object[[
      paste0(
        "gnrh_stage_",
        nm,
        "_fraction"
      )
    ]] <- detection_fractions[[nm]]
  }

  object$gnrh_migration_core_hits <-
    migration_core_hits

  # ------------------------------------------------------------------------- #
  # Metadata: secretory program
  # ------------------------------------------------------------------------- #

  object$gnrh_secretory_core_hits <-
    secretory_core_hits

  object$gnrh_secretory_supportive_hits <-
    secretory_supportive_hits

  object$gnrh_secretory_hits <-
    secretory_hits

  object$gnrh_secretory_supported <-
    secretory_supported

  object$gnrh_secretory <- factor(
    secretory_state,
    levels = c(
      "limited",
      "supported",
      "non-gnrh"
    )
  )

  # ------------------------------------------------------------------------- #
  # Stored information
  # ------------------------------------------------------------------------- #

  object@misc$gnrh_stage_modules <-
    developmental_modules

  object@misc$gnrh_migration_core <-
    migration_core

  object@misc$gnrh_secretory_module <-
    secretory_module

  object@misc$gnrh_stage_parameters <- list(
    assay = assay,
    layer = layer,

    expression_weight =
      expression_weight,

    detection_weight =
      detection_weight,

    stage_margin =
      stage_margin,

    min_migration_hits =
      min_migration_hits,

    min_secretory_core_hits =
      min_secretory_core_hits,

    min_secretory_supportive_hits =
      min_secretory_supportive_hits
  )

  if (is.null(object@misc$gnrh)) {
    object@misc$gnrh <- list()
  }

  object@misc$gnrh$stage <- list(
    modules =
      developmental_modules,

    migration_core =
      migration_core,

    secretory_module =
      secretory_module,

    parameters =
      object@misc$gnrh_stage_parameters
  )

  object
}
