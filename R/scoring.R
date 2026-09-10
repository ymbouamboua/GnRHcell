#' Compute module scores and detection metrics
#'
#' Computes expression, detection, and integrated activity metrics for
#' predefined gene modules.
#'
#' @param expr Sparse gene-by-cell expression matrix.
#' @param modules Named list of gene vectors defining biological modules.
#' @param scale_factor Library normalization scale factor.
#' @param expression_weight Weight assigned to normalized expression.
#' @param detection_weight Weight assigned to the fraction of module genes
#'   detected.
#'
#' @return A list containing expression scores, marker hits, detection
#' fractions, integrated program scores, and matched genes.
#'
#' @keywords internal
#' @noRd
.score_modules <- function(
    expr,
    modules,
    scale_factor = 10000,
    expression_weight = 0.5,
    detection_weight = 0.5
) {
  if (!length(modules))
    stop(
      "`modules` must contain at least one gene set.",
      call. = FALSE
    )

  if (
    length(expression_weight) != 1L ||
    !is.finite(expression_weight) ||
    expression_weight < 0
  ) {
    stop(
      "`expression_weight` must be a non-negative number.",
      call. = FALSE
    )
  }

  if (
    length(detection_weight) != 1L ||
    !is.finite(detection_weight) ||
    detection_weight < 0
  ) {
    stop(
      "`detection_weight` must be a non-negative number.",
      call. = FALSE
    )
  }

  if (
    expression_weight +
    detection_weight <= 0
  ) {
    stop(
      "At least one module-score weight must be positive.",
      call. = FALSE
    )
  }

  weight_sum <-
    expression_weight +
    detection_weight

  expression_weight <-
    expression_weight /
    weight_sum

  detection_weight <-
    detection_weight /
    weight_sum

  lib <- pmax(
    Matrix::colSums(expr),
    1
  )

  matched_modules <- lapply(
    modules,
    intersect,
    y = rownames(expr)
  )

  score_fun <- function(g) {
    if (!length(g))
      return(
        rep(
          0,
          ncol(expr)
        )
      )

    avg_expr <- Matrix::colMeans(
      expr[
        g,
        ,
        drop = FALSE
      ]
    )

    log1p(
      avg_expr /
        lib *
        scale_factor
    )
  }

  hit_fun <- function(g) {
    if (!length(g))
      return(
        integer(
          ncol(expr)
        )
      )

    as.integer(
      Matrix::colSums(
        expr[
          g,
          ,
          drop = FALSE
        ] > 0
      )
    )
  }

  fraction_fun <- function(g) {
    if (!length(g))
      return(
        rep(
          0,
          ncol(expr)
        )
      )

    Matrix::colSums(
      expr[
        g,
        ,
        drop = FALSE
      ] > 0
    ) / length(g)
  }

  to_df <- function(x) {
    x <- as.data.frame(
      x,
      check.names = FALSE
    )

    rownames(x) <-
      colnames(expr)

    x
  }

  score <- to_df(
    lapply(
      matched_modules,
      score_fun
    )
  )

  hits <- to_df(
    lapply(
      matched_modules,
      hit_fun
    )
  )

  fraction <- to_df(
    lapply(
      matched_modules,
      fraction_fun
    )
  )

  integrated <- as.data.frame(
    expression_weight *
      as.matrix(score) +
      detection_weight *
      as.matrix(fraction),
    check.names = FALSE
  )

  rownames(integrated) <-
    colnames(expr)

  list(
    score = score,
    hits = hits,
    fraction = fraction,
    integrated = integrated,
    genes = matched_modules,
    weights = c(
      expression = expression_weight,
      detection = detection_weight
    )
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
#'
#' @keywords internal
#' @noRd
.ambient <- function(
    expr,
    gene,
    lib,
    max_library = 100
) {
  low <-
    is.finite(lib) &
    lib > 0 &
    lib < max_library

  if (!any(low))
    return(0)

  x <- Matrix::rowMeans(
    expr[, low, drop = FALSE]
  )[gene]

  x <- as.numeric(x)

  if (!length(x) || !is.finite(x))
    return(0)

  x
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
  n <- length(signal)

  if (!reduction %in% names(object@reductions))
    return(rep(0, n))

  emb <- Seurat::Embeddings(
    object,
    reduction = reduction
  )

  dims <- dims[dims >= 1L & dims <= ncol(emb)]

  if (!length(dims))
    return(rep(0, n))

  emb <- emb[, dims, drop = FALSE]

  if (nrow(emb) != n)
    stop("Signal length does not match reduction.", call. = FALSE)

  if (n <= 1L)
    return(rep(0, n))

  k <- min(
    as.integer(k),
    n - 1L
  )

  if (k < 1L)
    return(rep(0, n))

  signal[!is.finite(signal)] <- 0

  nn <- FNN::get.knn(
    emb,
    k = k
  )

  rowMeans(
    matrix(
      signal[nn$nn.index],
      nrow = n,
      ncol = k
    )
  )
}


