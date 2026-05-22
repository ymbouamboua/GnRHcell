#' Compute module scores and detection metrics
#'
#' Internal helper for scoring predefined gene modules across cells.
#'
#' For each module, this function computes:
#' \itemize{
#'   \item normalized average expression score
#'   \item number of detected genes per cell
#'   \item fraction of module genes detected per cell
#' }
#'
#' Module scores are library-size normalized and log-transformed.
#'
#' @param expr Sparse gene expression matrix (genes x cells).
#' @param modules Named list of gene vectors defining biological modules.
#' @param scale_factor Library normalization scale factor.
#' Default is \code{10000}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{score}}{Data frame of normalized module activity scores.}
#'   \item{\code{hits}}{Data frame of detected gene counts per module.}
#'   \item{\code{fraction}}{Data frame of per-cell detection fractions.}
#'   \item{\code{genes}}{Matched module gene definitions.}
#' }
#'
#' @keywords internal
#' @noRd
.score_modules <- function(
    expr,
    modules,
    scale_factor = 10000
) {

  lib <- Matrix::colSums(expr)

  # Normalize module score
  score_fun <- function(g) {

    g <- intersect(
      g,
      rownames(expr)
    )

    if (length(g) == 0) {
      return(rep(0, ncol(expr)))
    }

    mat <- expr[g, , drop = FALSE]

    avg_expr <- Matrix::colMeans(mat)

    log1p(
      (avg_expr / lib) * scale_factor
    )
  }

  # Detection counts
  hit_fun <- function(g) {

    g <- intersect(
      g,
      rownames(expr)
    )

    if (length(g) == 0) {
      return(rep(0, ncol(expr)))
    }

    Matrix::colSums(
      expr[g, , drop = FALSE] > 0
    )
  }

  # Detection fraction
  frac_fun <- function(g) {

    g <- intersect(
      g,
      rownames(expr)
    )

    if (length(g) == 0) {
      return(rep(0, ncol(expr)))
    }

    Matrix::colSums(
      expr[g, , drop = FALSE] > 0
    ) / length(g)
  }

  # Compute
  score <- lapply(
    modules,
    score_fun
  )

  hits <- lapply(
    modules,
    hit_fun
  )

  fraction <- lapply(
    modules,
    frac_fun
  )

  # Convert to data.frames
  to_df <- function(x) {

    x <- as.data.frame(
      x,
      check.names = FALSE
    )

    rownames(x) <- colnames(expr)

    x
  }

  list(
    score = to_df(score),
    hits = to_df(hits),
    fraction = to_df(fraction),
    genes = modules
  )
}



#' Estimate ambient RNA background signal
#'
#' Internal helper for estimating ambient RNA contamination
#' from low-complexity cells.
#'
#' Ambient signal is estimated as the mean expression of the
#' target gene across cells with very low library size
#' (\code{lib < 100}).
#'
#' @param expr Sparse gene expression matrix (genes x cells).
#' @param gene Character vector or gene identifier used for ambient estimation.
#' @param lib Numeric vector of per-cell library sizes.
#'
#' @return Numeric ambient RNA estimate for the target gene.
#'
#' @keywords internal
#' @noRd
.ambient <- function(expr, gene, lib) {

  low <- lib < 100

  if (!any(low)) {
    return(0)
  }

  x <- Matrix::rowMeans(
    expr[, low, drop = FALSE]
  )[gene]

  ifelse(is.na(x), 0, x)
}



#' Compute neighborhood support signal
#'
#' Internal helper for estimating local transcriptomic support
#' using k-nearest neighbors in a reduced-dimensional embedding.
#'
#' For each cell, the mean signal across its nearest neighbors
#' is computed as a neighborhood support metric.
#'
#' @param object A Seurat object containing dimensional reductions.
#' @param signal Numeric per-cell signal vector.
#' @param reduction Dimensional reduction used for neighborhood search.
#' Default is \code{"pca"}.
#' @param dims Dimensions used from the selected reduction.
#' Default is \code{1:20}.
#' @param k Number of nearest neighbors. Default is 20.
#'
#' @return Numeric vector of neighborhood support scores.
#' Returns zeros if the requested reduction is unavailable.
#'
#' @keywords internal
#' @noRd
.knn_signal <- function(
    object,
    signal,
    reduction = "pca",
    dims = 1:20,
    k = 20
) {

  if (!reduction %in%
      names(object@reductions)) {

    return(rep(0, length(signal)))
  }

  emb <- Seurat::Embeddings(
    object,
    reduction = reduction
  )[, dims, drop = FALSE]

  nn <- FNN::get.knn(
    emb,
    k = k
  )

  vapply(
    seq_along(signal),
    function(i) {
      mean(signal[
        nn$nn.index[i, ]
      ])
    },
    numeric(1)
  )
}



