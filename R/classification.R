#' Classify candidate GnRH cells
#'
#' Internal helper for assigning GnRH-positive or GnRH-negative status using
#' direct \code{GNRH1} expression together with independent transcriptomic
#' evidence of GnRH neuronal identity.
#'
#' Classification uses two positive routes:
#' \itemize{
#'   \item \strong{Direct}: raw \code{GNRH1} expression greater than or equal
#'   to \code{min_umi}, together with independent GnRH identity evidence; and
#'   \item \strong{Supported}: detectable but sub-threshold \code{GNRH1}
#'   expression together with independent GnRH identity and transcriptomic
#'   support.
#' }
#'
#' Cells without detected \code{GNRH1} are never classified as GnRH-positive.
#' However, cells with strong GnRH-like transcriptomic evidence can be flagged
#' separately as \code{dropout_candidate} for diagnostic purposes.
#'
#' Cells with direct \code{GNRH1} signal but no independent GnRH identity
#' evidence are retained as \code{direct_isolated} diagnostic signals but are
#' not classified as GnRH-positive.
#'
#' @param raw Numeric vector containing raw \code{GNRH1} UMI counts.
#' @param norm Numeric vector containing normalized \code{GNRH1} expression.
#' @param score Numeric vector containing the composite GnRH detection score.
#' @param support_score Numeric vector containing the independent
#'   transcriptomic GnRH support score.
#' @param hits List containing marker-hit vectors.
#' @param lib Numeric vector containing total library sizes.
#' @param min_umi Minimum raw \code{GNRH1} UMI count required for direct
#'   detection. Default is 2.
#' @param min_counts Minimum total UMI count required for classification.
#'   Default is 500.
#' @param supported_q Quantile of the direct-cell support distribution used
#'   for supported low-expression candidates. Default is 0.60.
#' @param candidate_q Quantile of the direct-cell support distribution used
#'   only to flag transcriptomic dropout candidates. Default is 0.95.
#' @param expr_thr Adaptive normalized \code{GNRH1} expression threshold.
#' @param knn Optional numeric vector containing k-nearest-neighbor support.
#' @param alternative_score Optional numeric vector containing the strongest
#'   alternative neuronal identity score.
#' @param max_alternative Maximum alternative identity score tolerated for
#'   transcriptomic dropout candidates. Default is 0.75.
#' @param identity_strong Logical vector indicating strong GnRH identity.
#' @param identity_moderate Logical vector indicating moderate GnRH identity.
#' @param independent_support Logical vector indicating sufficient independent
#'   GnRH biological support.
#'
#' @return A list containing classification vectors, diagnostic subclasses,
#'   adaptive thresholds, rules, and rule summaries.
#'
#' Classify candidate GnRH cells
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
    candidate_q = 0.95,
    min_reference_cells = 20L,
    expr_thr,
    knn = NULL,
    alternative_score = NULL,
    alternative_strong = NULL,
    max_alternative = 0.75,
    identity_strong,
    identity_moderate,
    independent_support
) {
  n <- length(raw)

  required <- list(
    raw = raw, norm = norm, score = score, support_score = support_score,
    lib = lib, identity_strong = identity_strong,
    identity_moderate = identity_moderate,
    independent_support = independent_support
  )

  if (any(vapply(required, length, integer(1)) != n))
    stop("Classification vectors must have identical lengths.", call. = FALSE)

  for (nm in c("supported_q", "candidate_q")) {
    x <- get(nm)
    if (length(x) != 1L || !is.finite(x) || x < 0 || x > 1)
      stop(sprintf("`%s` must be between 0 and 1.", nm), call. = FALSE)
  }

  if (candidate_q < supported_q)
    stop("`candidate_q` must be >= `supported_q`.", call. = FALSE)

  min_reference_cells <- as.integer(min_reference_cells)

  if (!is.finite(min_reference_cells) || min_reference_cells < 1L)
    stop("`min_reference_cells` must be >= 1.", call. = FALSE)

  required_hits <- c("core", "mig", "neuro")
  missing_hits <- setdiff(required_hits, names(hits))

  if (length(missing_hits))
    stop("Missing marker-hit vectors: ",
         paste(missing_hits, collapse = ", "), call. = FALSE)

  if (any(vapply(hits, length, integer(1)) != n))
    stop("All marker-hit vectors must have length `n`.", call. = FALSE)

  # ------------------------------------------------------------------------- #
  # Basic evidence
  # ------------------------------------------------------------------------- #

  lib_ok <- is.finite(lib) & lib >= min_counts
  umi_any <- is.finite(raw) & raw > 0
  umi_ok <- is.finite(raw) & raw >= min_umi
  umi_low <- is.finite(raw) & raw > 0 & raw < min_umi

  expr_ok <- is.finite(norm) & norm >= expr_thr
  support_score_ok <- is.finite(support_score)

  identity_strong <- !is.na(identity_strong) & identity_strong
  identity_moderate <- !is.na(identity_moderate) & identity_moderate
  independent_support <- !is.na(independent_support) & independent_support

  # ------------------------------------------------------------------------- #
  # kNN identity support
  # ------------------------------------------------------------------------- #

  if (is.null(knn)) knn <- rep(0, n)

  if (length(knn) != n)
    stop("`knn` must have length `n`.", call. = FALSE)

  knn[!is.finite(knn)] <- 0
  knn <- pmin(pmax(knn, 0), 1)

  knn_support_thr <- 0.05
  knn_strong_thr <- 0.15

  knn_ok <- knn > knn_support_thr
  knn_strong <- knn > knn_strong_thr

  # ------------------------------------------------------------------------- #
  # Marker evidence
  # ------------------------------------------------------------------------- #

  core_hits <- hits$core
  mig_hits <- hits$mig

  identity_primary_hits <- hits$identity_primary %||% integer(n)
  neuro_primary_hits <- hits$neuro_primary %||% integer(n)
  neuro_supportive_hits <- hits$neuro_supportive %||% integer(n)

  core_supported <- core_hits >= 2L
  core_strong <- core_hits >= 3L

  migration_supported <- mig_hits >= 1L
  migration_strong <- mig_hits >= 2L

  neuro_supported <-
    neuro_primary_hits >= 1L |
    neuro_supportive_hits >= 2L

  neuro_strong <-
    neuro_primary_hits >= 2L |
    (neuro_primary_hits >= 1L & neuro_supportive_hits >= 2L)

  primary_identity <- identity_primary_hits >= 1L
  strong_primary_identity <- identity_primary_hits >= 2L

  # ------------------------------------------------------------------------- #
  # Alternative identities
  # ------------------------------------------------------------------------- #

  if (is.null(alternative_score))
    alternative_score <- rep(0, n)

  if (length(alternative_score) != n)
    stop("`alternative_score` must have length `n`.", call. = FALSE)

  alternative_score[!is.finite(alternative_score)] <- 0

  if (is.null(alternative_strong))
    alternative_strong <- alternative_score > max_alternative

  if (length(alternative_strong) != n)
    stop("`alternative_strong` must have length `n`.", call. = FALSE)

  alternative_strong <-
    !is.na(alternative_strong) &
    alternative_strong

  alternative_low <-
    !alternative_strong &
    alternative_score <= max_alternative

  # ------------------------------------------------------------------------- #
  # Route 1: direct detection
  # ------------------------------------------------------------------------- #

  direct_signal <- lib_ok & umi_ok

  direct_supported <-
    direct_signal &
    identity_moderate

  direct_isolated <-
    direct_signal &
    !identity_moderate

  direct <- direct_supported

  # ------------------------------------------------------------------------- #
  # High-specificity calibration population
  # ------------------------------------------------------------------------- #

  reference_positive <-
    direct_signal &
    identity_strong &
    independent_support

  reference_support <-
    support_score[
      reference_positive &
        support_score_ok
    ]

  reference_n <- length(reference_support)

  # ------------------------------------------------------------------------- #
  # Supported threshold
  #
  # No global-cell fallback: without enough trusted reference cells,
  # supported calls are disabled.
  # ------------------------------------------------------------------------- #

  supported_thr <- if (reference_n >= min_reference_cells) {
    as.numeric(stats::quantile(
      reference_support,
      probs = supported_q,
      na.rm = TRUE,
      names = FALSE
    ))
  } else {
    Inf
  }

  support_supported <-
    support_score_ok &
    support_score >= supported_thr

  # ------------------------------------------------------------------------- #
  # Route 2: low-UMI GNRH1
  #
  # Normalized GNRH1 is retained diagnostically but not used as a required
  # gate because one UMI is diluted in high-depth libraries.
  # ------------------------------------------------------------------------- #

  supported_candidate <-
    lib_ok &
    umi_low &
    identity_moderate &
    independent_support

  supported <-
    supported_candidate &
    support_supported

  # ------------------------------------------------------------------------- #
  # GNRH1-negative transcriptomic candidate
  #
  # Diagnostic only.
  # ------------------------------------------------------------------------- #

  candidate_thr <- if (reference_n >= min_reference_cells) {
    as.numeric(stats::quantile(
      reference_support,
      probs = candidate_q,
      na.rm = TRUE,
      names = FALSE
    ))
  } else {
    Inf
  }

  candidate_support <-
    support_score_ok &
    support_score >= candidate_thr

  transcriptomic_candidate <-
    lib_ok &
    !umi_any &
    identity_strong &
    neuro_supported &
    independent_support &
    knn_strong &
    alternative_low &
    candidate_support

  # ------------------------------------------------------------------------- #
  # Final classes
  # ------------------------------------------------------------------------- #

  cls <- rep("neg", n)
  cls[supported] <- "supported"
  cls[direct] <- "direct"

  status <- ifelse(cls == "neg", "neg", "pos")
  keep <- cls != "neg"

  # ------------------------------------------------------------------------- #
  # Diagnostics
  # ------------------------------------------------------------------------- #

  rules <- data.frame(
    library_ok = lib_ok,
    gnrh_detected = umi_any,
    gnrh_min_umi = umi_ok,
    gnrh_low_umi = umi_low,
    gnrh_expr_high = expr_ok,

    primary_identity = primary_identity,
    strong_primary_identity = strong_primary_identity,
    identity_moderate = identity_moderate,
    identity_strong = identity_strong,
    independent_support = independent_support,

    core_supported = core_supported,
    core_strong = core_strong,

    migration_support = migration_supported,
    migration_strong = migration_strong,

    neuro_support = neuro_supported,
    neuro_strong = neuro_strong,

    knn_support = knn_ok,
    knn_strong = knn_strong,

    alternative_low = alternative_low,
    alternative_strong = alternative_strong,

    reference_positive = reference_positive,
    support_supported = support_supported,

    direct_signal = direct_signal,
    direct_supported = direct_supported,
    direct_isolated = direct_isolated,
    direct = direct,

    supported_candidate = supported_candidate,
    supported = supported,

    transcriptomic_candidate = transcriptomic_candidate,
    stringsAsFactors = FALSE
  )

  rule_summary <- vapply(
    rules,
    function(x) sum(x, na.rm = TRUE),
    integer(1)
  )

  list(
    status = status,
    class = cls,
    keep = keep,

    direct_signal = direct_signal,
    direct_supported = direct_supported,
    direct_isolated = direct_isolated,

    transcriptomic_candidate = transcriptomic_candidate,

    # Backward compatibility
    dropout_candidate = transcriptomic_candidate,

    reference_positive = reference_positive,
    reference_n = reference_n,

    thr = supported_thr,
    supported_thr = supported_thr,
    candidate_thr = candidate_thr,

    knn_support_thr = knn_support_thr,
    knn_strong_thr = knn_strong_thr,

    rules = rules,
    rule_summary = rule_summary
  )
}
