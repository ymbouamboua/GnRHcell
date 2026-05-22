#' Classify candidate GnRH cells
#'
#' Internal helper for assigning GnRH-positive/negative status
#' based on expression thresholds, marker support, library size,
#' and composite scoring.
#'
#' A cell is considered marker-supported if it satisfies one of:
#' \itemize{
#'   \item at least 2 core markers with sufficient raw UMI
#'   \item at least 2 core markers and at least 1 migration marker
#'   \item at least 1 core marker, at least 2 neuronal markers, and sufficient raw UMI
#'   \item at least 2 core markers with positive kNN support
#' }
#'
#' Cells passing initial filtering are used to estimate an adaptive
#' score threshold based on the specified score quantile.
#'
#' @param raw Numeric vector of raw GnRH UMI counts.
#' @param norm Numeric vector of normalized GnRH expression values.
#' @param score Numeric vector of composite detection scores.
#' @param hits List containing marker hit counts with elements:
#' \code{core}, \code{mig}, and \code{neuro}.
#' @param lib Numeric vector of total library sizes (UMI counts).
#' @param min_umi Minimum raw GnRH UMI required for expression support.
#' Default is 2.
#' @param min_counts Minimum total library size required.
#' Default is 500.
#' @param score_q Quantile used to define adaptive score threshold.
#' Default is 0.9.
#' @param expr_thr Minimum normalized expression threshold.
#' @param knn Optional numeric vector of k-nearest-neighbor support scores.
#' If \code{NULL}, neighborhood support is ignored.
#'
#' @return A list with:
#' \describe{
#'   \item{class}{Character vector of classification labels
#'   (\code{"neg"} or \code{"pos"}).}
#'   \item{keep}{Logical vector indicating high-confidence retained cells.}
#'   \item{thr}{Numeric adaptive score threshold used for filtering.}
#' }
#'
#' @keywords internal
#' @noRd
.classify <- function(
    raw,
    norm,
    score,
    hits,
    lib,
    min_umi = 2,
    min_counts = 500,
    score_q = 0.9,
    expr_thr,
    knn = NULL
) {

  umi_ok <- raw >= min_umi
  lib_ok <- lib >= min_counts
  expr_ok <- norm >= expr_thr

  knn_ok <- if (!is.null(knn)) knn > 0.05 else FALSE

  marker_ok <- (
    (hits$core >= 2 & umi_ok) |
      (hits$core >= 2 & hits$mig >= 1) |
      (hits$core >= 1 & hits$neuro >= 2 & umi_ok) |
      (hits$core >= 2 & knn_ok)
  )

  keep0 <- umi_ok & lib_ok & expr_ok & marker_ok

  thr <- if (any(keep0)) {
    stats::quantile(score[keep0], score_q, na.rm = TRUE)
  } else {
    stats::median(score, na.rm = TRUE)
  }

  keep <- keep0 & score >= thr & score > 1

  cls <- rep("neg", length(score))
  cls[umi_ok & lib_ok & marker_ok] <- "pos"
  cls[keep] <- "pos"

  list(
    status = cls,
    keep = keep,
    thr = thr
  )
}
