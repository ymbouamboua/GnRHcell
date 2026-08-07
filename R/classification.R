#' Classify candidate GnRH cells
#'
#' Internal helper for assigning GnRH-positive or GnRH-negative status using
#' direct \code{GNRH1} expression and an independent transcriptomic GnRH
#' support score.
#'
#' Classification uses three complementary routes:
#' \itemize{
#'   \item \strong{Direct}: raw \code{GNRH1} expression greater than or equal
#'   to \code{min_umi};
#'   \item \strong{Supported}: detectable but sub-threshold \code{GNRH1}
#'   expression together with independent GnRH transcriptomic support; and
#'   \item \strong{Dropout rescue}: undetected \code{GNRH1} together with
#'   strong GnRH identity, developmental or neuroendocrine support, strong
#'   local neighborhood support, and no dominant alternative neuronal
#'   identity.
#' }
#'
#' The transcriptomic support score is intentionally independent of direct
#' \code{GNRH1} expression. Thresholds for supported and dropout-rescue
#' classification are estimated from cells classified through the direct
#' route whenever enough direct cells are available.
#'
#' Alternative neuronal programs are used only for the dropout-rescue route.
#' They do not reject cells showing direct \code{GNRH1} evidence.
#'
#' @param raw Numeric vector containing raw \code{GNRH1} UMI counts.
#' @param norm Numeric vector containing normalized \code{GNRH1} expression.
#' @param score Numeric vector containing the composite GnRH detection score.
#' @param support_score Numeric vector containing the transcriptomic GnRH
#'   support score calculated independently of direct \code{GNRH1}
#'   expression.
#' @param hits List containing marker hit vectors. Required elements are
#'   \code{core}, \code{mig}, and \code{neuro}. Additional elements such as
#'   \code{identity_primary} are used when available.
#' @param lib Numeric vector containing total library sizes.
#' @param min_umi Minimum raw \code{GNRH1} UMI count required for direct
#'   detection. Default is 2.
#' @param min_counts Minimum total UMI count required for classification.
#'   Default is 500.
#' @param supported_q Quantile of the direct-cell transcriptomic support
#'   distribution used for supported low-expression candidates.
#'   Default is 0.25.
#' @param dropout_q Quantile of the direct-cell transcriptomic support
#'   distribution used for dropout rescue. Default is 0.90.
#' @param expr_thr Adaptive normalized \code{GNRH1} expression threshold.
#'   Retained for diagnostics.
#' @param knn Optional numeric vector containing k-nearest-neighbor support.
#' @param alternative_score Optional numeric vector containing the strongest
#'   alternative neuronal or neuroendocrine identity score for each cell.
#' @param max_alternative Maximum alternative identity score tolerated for
#'   dropout-rescue classification. Default is 0.75.
#'
#' @return A list containing classification vectors, adaptive support
#'   thresholds, individual rules, and rule summaries.
#'
#' @keywords internal
#' @noRd
.classify <- function(
    raw,
    norm,
    score,
    support_score,
    hits,
    lib,
    min_umi = 2,
    min_counts = 500,
    supported_q = 0.60,
    dropout_q = 0.95,
    expr_thr,
    knn = NULL,
    alternative_score = NULL,
    max_alternative = 0.75
) {

  # --------------------------------------------------------------------------- #
  # Input checks
  # --------------------------------------------------------------------------- #

  n <- length(score)

  required_vectors <- list(
    raw = raw,
    norm = norm,
    score = score,
    support_score = support_score,
    lib = lib
  )

  bad_length <- vapply(
    required_vectors,
    length,
    integer(1)
  ) != n

  if (any(bad_length)) {
    stop(
      paste0(
        "raw, norm, score, support_score, and lib ",
        "must have identical lengths."
      ),
      call. = FALSE
    )
  }

  if (
    length(supported_q) != 1L ||
    !is.finite(supported_q) ||
    supported_q < 0 ||
    supported_q > 1
  ) {
    stop(
      "supported_q must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    length(dropout_q) != 1L ||
    !is.finite(dropout_q) ||
    dropout_q < 0 ||
    dropout_q > 1
  ) {
    stop(
      "dropout_q must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (dropout_q < supported_q) {
    stop(
      "dropout_q must be greater than or equal to supported_q.",
      call. = FALSE
    )
  }

  required_hits <- c(
    "core",
    "mig",
    "neuro"
  )

  missing_hits <- setdiff(
    required_hits,
    names(hits)
  )

  if (length(missing_hits)) {
    stop(
      "Missing marker hit vectors: ",
      paste(
        missing_hits,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------- #
  # Basic evidence
  # --------------------------------------------------------------------------- #

  lib_ok <- is.finite(lib) &
    lib >= min_counts

  umi_any <- is.finite(raw) &
    raw > 0

  umi_ok <- is.finite(raw) &
    raw >= min_umi

  expr_ok <- is.finite(norm) &
    norm >= expr_thr

  score_ok <- is.finite(score)

  support_score_ok <- is.finite(
    support_score
  )

  # --------------------------------------------------------------------------- #
  # Neighborhood support
  # --------------------------------------------------------------------------- #

  if (is.null(knn)) {

    knn <- rep(
      0,
      n
    )

  } else if (length(knn) != n) {

    stop(
      "knn must have the same length as score.",
      call. = FALSE
    )
  }

  knn[
    !is.finite(knn)
  ] <- 0

  knn_ok <- knn > 0.05
  knn_strong <- knn > 0.15

  # --------------------------------------------------------------------------- #
  # Marker support
  # --------------------------------------------------------------------------- #

  core_hits <- hits$core
  mig_hits <- hits$mig
  neuro_hits <- hits$neuro

  identity_primary_hits <- if (
    "identity_primary" %in% names(hits)
  ) {

    hits$identity_primary

  } else {

    rep(
      0L,
      n
    )
  }

  core_supported <- core_hits >= 2
  core_strong <- core_hits >= 3

  migration_supported <- mig_hits >= 1
  migration_strong <- mig_hits >= 2

  neuro_supported <- neuro_hits >= 1
  neuro_strong <- neuro_hits >= 2

  primary_identity <- identity_primary_hits >= 1
  strong_primary_identity <- identity_primary_hits >= 2

  # --------------------------------------------------------------------------- #
  # Alternative neuronal identities
  # --------------------------------------------------------------------------- #

  if (is.null(alternative_score)) {

    alternative_score <- rep(
      0,
      n
    )

  } else if (length(alternative_score) != n) {

    stop(
      "alternative_score must have the same length as score.",
      call. = FALSE
    )
  }

  alternative_score[
    !is.finite(alternative_score)
  ] <- 0

  alternative_low <- alternative_score <= max_alternative

  # --------------------------------------------------------------------------- #
  # Biological support
  # --------------------------------------------------------------------------- #

  identity_support <- (
    primary_identity |
      core_supported
  )

  broad_support <- (
    identity_support &
      (
        migration_strong |
          neuro_supported |
          knn_strong
      )
  )

  strong_program_support <- (
    (
      strong_primary_identity &
        core_strong
    ) |
      (
        core_strong &
          migration_strong &
          neuro_supported
      ) |
      (
        primary_identity &
          migration_strong &
          neuro_strong
      )
  )

  # --------------------------------------------------------------------------- #
  # Route 1: direct GNRH1 detection
  # --------------------------------------------------------------------------- #

  direct <- (
    lib_ok &
      umi_ok
  )

  # --------------------------------------------------------------------------- #
  # Reference transcriptomic support distribution
  #
  # Direct cells define the biological support distribution against which
  # weaker-expression candidates are evaluated.
  # --------------------------------------------------------------------------- #

  reference_support <- support_score[
    direct &
      support_score_ok
  ]

  finite_support <- support_score[
    support_score_ok
  ]

  # --------------------------------------------------------------------------- #
  # Supported threshold
  # --------------------------------------------------------------------------- #

  if (length(reference_support) >= 20L) {

    supported_thr <- as.numeric(
      stats::quantile(
        reference_support,
        probs = supported_q,
        na.rm = TRUE,
        names = FALSE
      )
    )

  } else if (length(finite_support)) {

    supported_thr <- as.numeric(
      stats::quantile(
        finite_support,
        probs = 0.90,
        na.rm = TRUE,
        names = FALSE
      )
    )

  } else {

    supported_thr <- Inf
  }

  # --------------------------------------------------------------------------- #
  # Dropout-rescue threshold
  # --------------------------------------------------------------------------- #

  if (length(reference_support) >= 20L) {

    dropout_thr <- as.numeric(
      stats::quantile(
        reference_support,
        probs = dropout_q,
        na.rm = TRUE,
        names = FALSE
      )
    )

  } else {

    # Dropout rescue is intentionally disabled when there are not enough
    # direct cells to define a reliable reference distribution.
    dropout_thr <- Inf
  }

  support_supported <- (
    support_score_ok &
      support_score >= supported_thr
  )

  support_dropout <- (
    support_score_ok &
      support_score >= dropout_thr
  )

  # --------------------------------------------------------------------------- #
  # Route 2: detectable but sub-threshold GNRH1
  # --------------------------------------------------------------------------- #

  supported_candidate <- (
    lib_ok &
      umi_any &
      !umi_ok &
      primary_identity &
      broad_support
  )

  supported <- (
    supported_candidate &
      support_supported
  )

  # --------------------------------------------------------------------------- #
  # Route 3: GNRH1 dropout rescue
  # --------------------------------------------------------------------------- #

  dropout_candidate <- (
    lib_ok &
      !umi_any &
      primary_identity &
      strong_program_support &
      knn_strong &
      alternative_low
  )

  dropout_rescue <- (
    dropout_candidate &
      support_dropout
  )

  # --------------------------------------------------------------------------- #
  # Final classification
  # --------------------------------------------------------------------------- #

  cls <- rep("neg", n)
  cls[dropout_rescue] <- "dropout_rescue"
  cls[supported] <- "supported"
  cls[direct] <- "direct"
  status <- ifelse(cls == "neg", "neg", "pos")

  # --------------------------------------------------------------------------- #
  # High-confidence candidates
  # --------------------------------------------------------------------------- #

  keep <- (
    direct |
      supported |
      dropout_rescue
    ) &
    lib_ok

  # --------------------------------------------------------------------------- #
  # Diagnostic rules
  # --------------------------------------------------------------------------- #

  rules <- data.frame(
    library_ok = lib_ok,
    gnrh_detected = umi_any,
    gnrh_min_umi = umi_ok,
    gnrh_expr_high = expr_ok,

    primary_identity = primary_identity,
    strong_primary_identity = strong_primary_identity,

    identity_support = identity_support,
    broad_support = broad_support,
    strong_program_support = strong_program_support,

    migration_support = migration_supported,
    migration_strong = migration_strong,

    neuro_support = neuro_supported,
    neuro_strong = neuro_strong,

    knn_support = knn_ok,
    knn_strong = knn_strong,

    alternative_low = alternative_low,

    support_supported = support_supported,
    support_dropout = support_dropout,

    direct = direct,
    supported_candidate = supported_candidate,
    supported = supported,

    dropout_candidate = dropout_candidate,
    dropout_rescue = dropout_rescue,

    stringsAsFactors = FALSE
  )

  rule_summary <- vapply(
    rules,
    function(x) {
      sum(
        x,
        na.rm = TRUE
      )
    },
    integer(1)
  )

  list(
    status = status,
    class = cls,
    keep = keep,

    # `thr` retained for limited backward compatibility.
    thr = supported_thr,

    supported_thr = supported_thr,
    dropout_thr = dropout_thr,

    rules = rules,
    rule_summary = rule_summary
  )
}
