# =========================================================
# utils-core.R
# Clean + factorized utilities
# =========================================================

# ---------------------------------------------------------
# Logger
# ---------------------------------------------------------

#' Lightweight logger
#'
#' @keywords internal
#' @noRd
log_msg <- function(verbose = TRUE) {

  start_time <- Sys.time()

  function(...) {

    if (!isTRUE(verbose)) {
      return(invisible(NULL))
    }

    msg <- paste(..., collapse = " ")

    if (grepl("DONE", msg, fixed = TRUE)) {

      elapsed <- difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )

      msg <- sprintf(
        "%s (%.2fs)",
        msg,
        as.numeric(elapsed)
      )
    }

    cat("INFO ", msg, "\n", sep = "")

    invisible(msg)
  }
}


# ---------------------------------------------------------
# Gene utilities
# ---------------------------------------------------------

#' Match gene symbols ignoring case
#'
#' @param target Character vector of target genes.
#' @param pool Character vector of available genes.
#'
#' @return Matched gene symbols.
#'
#' @export
match_gene_symbols <- function(target, pool) {

  pool[
    match(
      toupper(target),
      toupper(pool)
    )
  ]
}


#' Detect first matching gene
#'
#' @keywords internal
#' @noRd
detect_gene <- function(
    pool,
    candidates
) {

  hit <- match_gene_symbols(
    target = candidates,
    pool = pool
  )

  hit <- hit[!is.na(hit)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


# ---------------------------------------------------------
# Seurat helpers
# ---------------------------------------------------------

#' Extract expression matrix
#'
#' @param object Seurat object.
#' @param assay Assay name.
#' @param layer Data layer.
#'
#' @return Sparse expression matrix.
#'
#' @export
extract_expression_matrix <- function(
    object,
    assay = "RNA",
    layer = "counts"
) {

  stopifnot(inherits(object, "Seurat"))

  assays <- SeuratObject::Assays(object)

  if (!assay %in% assays) {

    stop(
      sprintf(
        "Assay '%s' not found",
        assay
      ),
      call. = FALSE
    )
  }

  SeuratObject::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )
}


#' Convert matrix to sparse dgCMatrix
#'
#' @keywords internal
#' @noRd
as_sparse_matrix <- function(x) {

  if (inherits(x, "dgCMatrix")) {
    return(x)
  }

  methods::as(x, "dgCMatrix")
}


# ---------------------------------------------------------
# Normalization
# ---------------------------------------------------------

#' Library size normalization
#'
#' @keywords internal
#' @noRd
library_normalize <- function(
    x,
    lib,
    scale_factor = GNRH_SCALE_FACTOR
) {

  norm <- scale_factor / pmax(lib, 1)

  log1p(x * norm)
}


# ---------------------------------------------------------
# Safe statistics
# ---------------------------------------------------------

#' Safe standardization
#'
#' @keywords internal
#' @noRd
safe_scale <- function(x) {

  s <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}


#' Safe median
#'
#' @keywords internal
#' @noRd
safe_median <- function(x) {

  stats::median(
    x,
    na.rm = TRUE
  )
}


#' Safe quantile
#'
#' @keywords internal
#' @noRd
safe_quantile <- function(x, q = 0.5) {

  stats::quantile(
    x,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )
}


#' Safe MAD
#'
#' @keywords internal
#' @noRd
safe_mad <- function(x) {

  stats::mad(
    x,
    na.rm = TRUE
  )
}


#' Safe SD
#'
#' @keywords internal
#' @noRd
safe_sd <- function(x) {

  stats::sd(
    x,
    na.rm = TRUE
  )
}


# ---------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------

#' Replace NULL with default
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {

  if (is.null(x)) y else x
}


#' Clamp numeric vector
#'
#' @keywords internal
#' @noRd
clamp <- function(
    x,
    lower = -Inf,
    upper = Inf
) {

  pmin(
    pmax(x, lower),
    upper
  )
}


#' Remove invalid rows
#'
#' @keywords internal
#' @noRd
drop_invalid_rows <- function(
    df,
    cols
) {

  keep <- Reduce(
    `&`,
    lapply(cols, function(cl) {

      is.finite(df[[cl]]) &
        !is.na(df[[cl]])
    })
  )

  df[keep, , drop = FALSE]
}
