#' Run GnRH detection diagnostics
#'
#' Generates descriptive and, when an independent truth annotation is
#' available, validation diagnostics for GnRH detection results.
#'
#' Descriptive diagnostics are always generated. Performance metrics such as
#' ROC, precision-recall curves, AUC, AUPRC, sensitivity, specificity,
#' precision, F1 score, and false-positive rate are computed only when an
#' independent binary truth annotation is supplied.
#'
#' @param object A Seurat object previously processed with
#'   \code{\link{detect_gnrh}}.
#' @param truth Optional metadata column name or vector containing independent
#'   binary truth labels.
#' @param positive Optional value identifying the positive truth class.
#'   If \code{NULL}, logical truth uses \code{TRUE}; otherwise the second
#'   observed class is used.
#' @param predictor Metadata column used for threshold validation.
#'   Default is \code{"gnrh_support_score_raw"}.
#' @param n_thresholds Number of thresholds evaluated. Default is 200.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A Seurat object with diagnostic results stored in
#'   \code{object@misc$gnrh}.
#'
#' @details
#' Without an independent truth annotation, this function stores per-cell
#' diagnostics and classification summaries but does not compute ROC or
#' threshold-performance statistics.
#'
#' When \code{truth} is supplied, threshold metrics include sensitivity,
#' specificity, precision, recall, F1 score, balanced accuracy, false-positive
#' rate, and false positives per 100,000 truth-negative cells.
#'
#' Precision-recall metrics are emphasized because GnRH neurons are typically
#' rare and ROC AUC alone may overstate performance in highly imbalanced data.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{run_gnrh}}
#'
#' @export
gnrh_diagnostics <- function(
    object,
    truth = NULL,
    positive = NULL,
    predictor = "gnrh_support_score_raw",
    n_thresholds = 200L,
    verbose = TRUE
) {
  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  log <- .msg(verbose)
  log("Running diagnostics")

  md <- object[[]]

  # ------------------------------------------------------------------------- #
  # Validate detection outputs
  # ------------------------------------------------------------------------- #

  if (is.null(object@misc$gnrh_params)) {
    stop(
      "Missing `gnrh_params`; run `detect_gnrh()` first.",
      call. = FALSE
    )
  }

  required_cols <- c(
    "gnrh_raw",
    "gnrh_expr",
    "gnrh_score",
    "gnrh_status",
    "gnrh_class",
    "gnrh_core_hits",
    "gnrh_mig_hits",
    "gnrh_neuro_hits",
    "gnrh_confident"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(md)
  )

  if (length(missing_cols)) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  object <- .init_gnrh_misc(object)

  # ------------------------------------------------------------------------- #
  # Add total marker evidence
  # ------------------------------------------------------------------------- #

  object <- .make_total_hits(object)
  md <- object[[]]

  # ------------------------------------------------------------------------- #
  # Helper for optional metadata
  # ------------------------------------------------------------------------- #

  optional <- function(x, default = NA) {
    if (x %in% colnames(md)) {
      md[[x]]
    } else {
      rep(default, nrow(md))
    }
  }

  # ------------------------------------------------------------------------- #
  # Unified per-cell diagnostics
  # ------------------------------------------------------------------------- #

  diagnostics_table <- data.frame(
    raw_expr = md$gnrh_raw,
    expr = md$gnrh_expr,

    score = md$gnrh_score,
    score_raw = optional(
      "gnrh_score_raw",
      NA_real_
    ),

    support_score = optional(
      "gnrh_support_score",
      NA_real_
    ),

    support_score_raw = optional(
      "gnrh_support_score_raw",
      NA_real_
    ),

    status = as.character(
      md$gnrh_status
    ),

    class = as.character(
      md$gnrh_class
    ),

    confident = as.logical(
      md$gnrh_confident
    ),

    direct_signal = as.logical(
      optional(
        "gnrh_direct_signal",
        FALSE
      )
    ),

    direct_supported = as.logical(
      optional(
        "gnrh_direct_supported",
        FALSE
      )
    ),

    direct_isolated = as.logical(
      optional(
        "gnrh_direct_isolated",
        FALSE
      )
    ),

    reference_positive = as.logical(
      optional(
        "gnrh_reference_positive",
        FALSE
      )
    ),

    transcriptomic_candidate = as.logical(
      optional(
        "gnrh_transcriptomic_candidate",
        optional(
          "gnrh_dropout_candidate",
          FALSE
        )
      )
    ),

    identity_moderate = as.logical(
      optional(
        "gnrh_identity_moderate",
        FALSE
      )
    ),

    identity_strong = as.logical(
      optional(
        "gnrh_identity_strong",
        FALSE
      )
    ),

    independent_support = as.logical(
      optional(
        "gnrh_independent_support",
        FALSE
      )
    ),

    neuro_support = as.logical(
      optional(
        "gnrh_neuro_support",
        FALSE
      )
    ),

    migration_support = as.logical(
      optional(
        "gnrh_migration_support",
        FALSE
      )
    ),

    alternative_strong = as.logical(
      optional(
        "gnrh_alternative_strong",
        FALSE
      )
    ),

    knn = optional(
      "gnrh_knn",
      NA_real_
    ),

    core_hits = md$gnrh_core_hits,
    mig_hits = md$gnrh_mig_hits,
    neuro_hits = md$gnrh_neuro_hits,

    total_hits = optional(
      "gnrh_total_hits",
      optional(
        "total_hits",
        NA_real_
      )
    ),

    stringsAsFactors = FALSE
  )

  rownames(
    diagnostics_table
  ) <- rownames(md)

  # ------------------------------------------------------------------------- #
  # Add optional staging / secretory metadata
  # ------------------------------------------------------------------------- #

  optional_cols <- c(
    "gnrh_stage",
    "gnrh_stage_raw",
    "gnrh_stage_resolution",
    "gnrh_stage_confident",
    "gnrh_stage_score",
    "gnrh_stage_second_score",
    "gnrh_stage_margin",
    "gnrh_stage_reason",
    "gnrh_stage_reassigned",

    "gnrh_stage_identity_score",
    "gnrh_stage_migrating_score",
    "gnrh_stage_mature_score",

    "gnrh_stage_identity_expression",
    "gnrh_stage_migrating_expression",
    "gnrh_stage_mature_expression",

    "gnrh_stage_identity_fraction",
    "gnrh_stage_migrating_fraction",
    "gnrh_stage_mature_fraction",

    "gnrh_migration_core_hits",

    "gnrh_secretory",
    "gnrh_secretory_supported",
    "gnrh_secretory_core_hits",
    "gnrh_secretory_supportive_hits",
    "gnrh_secretory_hits"
  )

  for (
    column in
    intersect(
      optional_cols,
      colnames(md)
    )
  ) {
    diagnostics_table[[column]] <-
      md[[column]]
  }

  object@misc$gnrh$diagnostics <-
    diagnostics_table

  # ------------------------------------------------------------------------- #
  # Descriptive classification summary
  # ------------------------------------------------------------------------- #

  count_value <- function(x, value) {
    sum(
      as.character(x) == value,
      na.rm = TRUE
    )
  }

  count_true <- function(x) {
    sum(
      x %in% TRUE,
      na.rm = TRUE
    )
  }

  object@misc$gnrh$diagnostic_summary <- list(
    n_cells = nrow(md),

    neg = count_value(
      md$gnrh_status,
      "neg"
    ),

    pos = count_value(
      md$gnrh_status,
      "pos"
    ),

    direct = count_value(
      md$gnrh_class,
      "direct"
    ),

    supported = count_value(
      md$gnrh_class,
      "supported"
    ),

    confident = count_true(
      md$gnrh_confident
    ),

    direct_signal = count_true(
      optional(
        "gnrh_direct_signal",
        FALSE
      )
    ),

    direct_supported = count_true(
      optional(
        "gnrh_direct_supported",
        FALSE
      )
    ),

    direct_isolated = count_true(
      optional(
        "gnrh_direct_isolated",
        FALSE
      )
    ),

    reference_positive = count_true(
      optional(
        "gnrh_reference_positive",
        FALSE
      )
    ),

    transcriptomic_candidate = count_true(
      optional(
        "gnrh_transcriptomic_candidate",
        optional(
          "gnrh_dropout_candidate",
          FALSE
        )
      )
    ),

    stage_resolved = count_true(
      optional(
        "gnrh_stage_confident",
        FALSE
      )
    )
  )

  # ------------------------------------------------------------------------- #
  # No independent truth: descriptive diagnostics only
  # ------------------------------------------------------------------------- #

  if (is.null(truth)) {
    object@misc$gnrh$threshold_curve <- NULL
    object@misc$gnrh$best_threshold <- NA_real_

    object@misc$gnrh$roc <- NULL
    object@misc$gnrh$auc <- NA_real_

    object@misc$gnrh$pr_curve <- NULL
    object@misc$gnrh$auprc <- NA_real_

    object@misc$gnrh$validation <- NULL

    log(
      "No independent truth supplied; validation metrics skipped.",
      type = "info"
    )

    return(object)
  }

  # ------------------------------------------------------------------------- #
  # Resolve truth
  # ------------------------------------------------------------------------- #

  truth_values <- if (
    is.character(truth) &&
    length(truth) == 1L
  ) {
    if (!truth %in% colnames(md)) {
      stop(
        "Truth column `",
        truth,
        "` not found in object metadata.",
        call. = FALSE
      )
    }

    md[[truth]]
  } else {
    truth
  }

  if (
    length(truth_values) !=
    nrow(md)
  ) {
    stop(
      "`truth` must contain one value per cell.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------------- #
  # Resolve predictor
  # ------------------------------------------------------------------------- #

  predictor_name <- as.character(
    predictor
  )

  if (
    length(predictor_name) != 1L ||
    is.na(predictor_name) ||
    !nzchar(predictor_name)
  ) {
    stop(
      "`predictor` must be one metadata column name.",
      call. = FALSE
    )
  }

  if (
    !predictor_name %in%
    colnames(md)
  ) {
    stop(
      "Predictor column `",
      predictor_name,
      "` not found.",
      call. = FALSE
    )
  }

  predictor_values <- suppressWarnings(
    as.numeric(
      md[[predictor_name]]
    )
  )

  # ------------------------------------------------------------------------- #
  # Convert truth to logical
  # ------------------------------------------------------------------------- #

  if (is.logical(truth_values)) {
    truth_positive <-
      truth_values

    positive_label <-
      TRUE

  } else {
    truth_character <- as.character(
      truth_values
    )

    classes <- unique(
      truth_character[
        !is.na(
          truth_character
        )
      ]
    )

    if (length(classes) != 2L) {
      stop(
        "`truth` must contain exactly two non-missing classes.",
        call. = FALSE
      )
    }

    positive_label <- positive

    if (is.null(positive_label)) {
      positive_label <-
        classes[[2L]]
    }

    if (
      !as.character(
        positive_label
      ) %in%
      classes
    ) {
      stop(
        "`positive` value not present in truth labels.",
        call. = FALSE
      )
    }

    truth_positive <-
      truth_character ==
      as.character(
        positive_label
      )
  }

  # ------------------------------------------------------------------------- #
  # Remove missing observations
  # ------------------------------------------------------------------------- #

  ok <-
    !is.na(
      truth_positive
    ) &
    is.finite(
      predictor_values
    )

  predictor_values <-
    predictor_values[
      ok
    ]

  truth_positive <-
    truth_positive[
      ok
    ]

  if (!length(predictor_values)) {
    stop(
      "No usable observations remain after removing missing values.",
      call. = FALSE
    )
  }

  if (
    length(
      unique(
        truth_positive
      )
    ) != 2L
  ) {
    stop(
      "Validation requires both positive and negative truth classes.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------------- #
  # Threshold grid
  # ------------------------------------------------------------------------- #

  if (
    length(n_thresholds) != 1L ||
    is.na(n_thresholds) ||
    !is.finite(n_thresholds) ||
    n_thresholds < 2 ||
    n_thresholds != floor(n_thresholds)
  ) {
    stop(
      "`n_thresholds` must be an integer >= 2.",
      call. = FALSE
    )
  }

  n_thresholds <-
    as.integer(
      n_thresholds
    )

  score_range <- range(
    predictor_values,
    finite = TRUE
  )

  if (
    !all(
      is.finite(
        score_range
      )
    ) ||
    diff(
      score_range
    ) <= 0
  ) {
    stop(
      "Predictor has no usable range.",
      call. = FALSE
    )
  }

  thresholds <- sort(
    unique(
      c(
        Inf,

        seq(
          score_range[[1L]],
          score_range[[2L]],
          length.out =
            n_thresholds
        ),

        -Inf
      )
    ),
    decreasing = TRUE
  )

  # ------------------------------------------------------------------------- #
  # Metric helper
  # ------------------------------------------------------------------------- #

  safe_div <- function(
    numerator,
    denominator
  ) {
    if (
      length(denominator) != 1L ||
      denominator <= 0
    ) {
      return(
        NA_real_
      )
    }

    numerator /
      denominator
  }

  # ------------------------------------------------------------------------- #
  # Threshold metrics
  # ------------------------------------------------------------------------- #

  curve <- lapply(
    thresholds,
    function(threshold) {
      predicted <-
        predictor_values >=
        threshold

      TP <- sum(
        predicted &
          truth_positive
      )

      FP <- sum(
        predicted &
          !truth_positive
      )

      FN <- sum(
        !predicted &
          truth_positive
      )

      TN <- sum(
        !predicted &
          !truth_positive
      )

      sensitivity <- safe_div(
        TP,
        TP + FN
      )

      specificity <- safe_div(
        TN,
        TN + FP
      )

      precision <- safe_div(
        TP,
        TP + FP
      )

      F1 <- if (
        is.finite(
          precision
        ) &&
        is.finite(
          sensitivity
        ) &&
        precision +
        sensitivity > 0
      ) {
        2 *
          precision *
          sensitivity /
          (
            precision +
              sensitivity
          )
      } else {
        NA_real_
      }

      balanced_accuracy <- if (
        is.finite(sensitivity) &&
        is.finite(specificity)
      ) {
        mean(
          c(
            sensitivity,
            specificity
          )
        )
      } else {
        NA_real_
      }

      fpr <- safe_div(
        FP,
        FP + TN
      )

      fp_per_100k <- safe_div(
        FP * 1e5,
        FP + TN
      )

      data.frame(
        threshold = threshold,

        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,

        sensitivity = sensitivity,
        recall = sensitivity,
        specificity = specificity,
        precision = precision,
        F1 = F1,

        balanced_accuracy =
          balanced_accuracy,

        fpr = fpr,
        tpr = sensitivity,

        fp_per_100k =
          fp_per_100k,

        stringsAsFactors = FALSE
      )
    }
  )

  curve <- do.call(
    rbind,
    curve
  )

  rownames(curve) <- NULL

  # ------------------------------------------------------------------------- #
  # Best threshold
  # ------------------------------------------------------------------------- #

  valid_f1 <- which(
    is.finite(
      curve$F1
    )
  )

  best_index <- if (
    length(valid_f1)
  ) {
    valid_f1[
      which.max(
        curve$F1[
          valid_f1
        ]
      )
    ]
  } else {
    NA_integer_
  }

  best_threshold <- if (
    !is.na(best_index)
  ) {
    curve$threshold[
      best_index
    ]
  } else {
    NA_real_
  }

  # ------------------------------------------------------------------------- #
  # ROC curve
  # ------------------------------------------------------------------------- #

  roc <- curve[
    is.finite(
      curve$fpr
    ) &
      is.finite(
        curve$tpr
      ),
    c(
      "threshold",
      "fpr",
      "tpr",
      "sensitivity",
      "specificity"
    ),
    drop = FALSE
  ]

  roc <- roc[
    order(
      roc$fpr,
      roc$tpr
    ),
    ,
    drop = FALSE
  ]

  roc <- roc[
    !duplicated(
      roc[
        ,
        c(
          "fpr",
          "tpr"
        ),
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  auc <- if (
    nrow(roc) >= 2L
  ) {
    sum(
      diff(
        roc$fpr
      ) *
        (
          utils::head(
            roc$tpr,
            -1L
          ) +
            utils::tail(
              roc$tpr,
              -1L
            )
        ) /
        2
    )
  } else {
    NA_real_
  }

  # ------------------------------------------------------------------------- #
  # Precision-recall curve
  # ------------------------------------------------------------------------- #

  pr <- curve[
    is.finite(
      curve$recall
    ) &
      is.finite(
        curve$precision
      ),
    c(
      "threshold",
      "recall",
      "precision"
    ),
    drop = FALSE
  ]

  if (nrow(pr)) {
    pr <- stats::aggregate(
      precision ~ recall,
      data = pr,
      FUN = max
    )

    pr <- pr[
      order(
        pr$recall
      ),
      ,
      drop = FALSE
    ]
  }

  auprc <- if (
    nrow(pr) >= 2L
  ) {
    sum(
      diff(
        pr$recall
      ) *
        (
          utils::head(
            pr$precision,
            -1L
          ) +
            utils::tail(
              pr$precision,
              -1L
            )
        ) /
        2
    )
  } else {
    NA_real_
  }

  # ------------------------------------------------------------------------- #
  # Performance of actual GnRHcell classification
  # ------------------------------------------------------------------------- #

  actual_call <-
    as.character(
      md$gnrh_status[
        ok
      ]
    ) ==
    "pos"

  TP <- sum(
    actual_call &
      truth_positive
  )

  FP <- sum(
    actual_call &
      !truth_positive
  )

  FN <- sum(
    !actual_call &
      truth_positive
  )

  TN <- sum(
    !actual_call &
      !truth_positive
  )

  sensitivity <- safe_div(
    TP,
    TP + FN
  )

  specificity <- safe_div(
    TN,
    TN + FP
  )

  precision <- safe_div(
    TP,
    TP + FP
  )

  F1 <- if (
    is.finite(
      precision
    ) &&
    is.finite(
      sensitivity
    ) &&
    precision +
    sensitivity > 0
  ) {
    2 *
      precision *
      sensitivity /
      (
        precision +
          sensitivity
      )
  } else {
    NA_real_
  }

  balanced_accuracy <- if (
    is.finite(sensitivity) &&
    is.finite(specificity)
  ) {
    mean(
      c(
        sensitivity,
        specificity
      )
    )
  } else {
    NA_real_
  }

  # ------------------------------------------------------------------------- #
  # Store validation summary
  # ------------------------------------------------------------------------- #

  validation <- list(
    predictor =
      predictor_name,

    positive =
      positive_label,

    n =
      length(
        truth_positive
      ),

    n_positive =
      sum(
        truth_positive
      ),

    n_negative =
      sum(
        !truth_positive
      ),

    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,

    sensitivity =
      sensitivity,

    recall =
      sensitivity,

    specificity =
      specificity,

    precision =
      precision,

    F1 =
      F1,

    balanced_accuracy =
      balanced_accuracy,

    false_positive_rate =
      safe_div(
        FP,
        FP + TN
      ),

    false_positives_per_100k =
      safe_div(
        FP * 1e5,
        FP + TN
      ),

    auc =
      auc,

    auprc =
      auprc,

    best_threshold =
      best_threshold,

    prevalence =
      mean(
        truth_positive
      )
  )

  # ------------------------------------------------------------------------- #
  # Store validation diagnostics
  # ------------------------------------------------------------------------- #

  object@misc$gnrh$threshold_curve <-
    curve

  object@misc$gnrh$best_threshold <-
    best_threshold

  object@misc$gnrh$roc <-
    roc

  object@misc$gnrh$auc <-
    auc

  object@misc$gnrh$pr_curve <-
    pr

  object@misc$gnrh$auprc <-
    auprc

  object@misc$gnrh$validation <-
    validation

  log(
    sprintf(
      "Validation: AUC %.3f | AUPRC %.3f | F1 %.3f",
      auc,
      auprc,
      F1
    ),
    type = "info"
  )

  object
}
