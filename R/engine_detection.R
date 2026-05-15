# =========================================================
# detection_modules.R
# =========================================================

# ---------------------------------------------------------
# Build detection gene modules
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
build_detection_modules <- function(genes) {

  map_genes <- function(features) {
    na.omit(match_gene_symbols(features, genes))
  }

  list(

    core = map_genes(c(
      "GNRH1", "FEZF1", "ISL1",
      "OTX2", "SIX3", "SIX6",
      "DLX1", "DLX2", "DLX5", "DLX6"
    )),

    migration = map_genes(c(
      "ANOS1", "PROKR2", "PROK2",
      "NRP1", "NRP2", "SEMA3A",
      "SEMA3C", "SEMA3F", "ROBO1",
      "ROBO2", "L1CAM", "DCX"
    )),

    neuroendo = map_genes(c(
      "KISS1R", "TAC3", "TACR3",
      "GNRHR", "PCSK1", "PCSK2",
      "SCG2", "CHGA", "CHGB",
      "CPE", "VGF", "SYP", "RAB3A"
    )),

    hormone = map_genes(c(
      "ESR1", "PGR", "AR"
    ))
  )
}


# =========================================================
# detection_scoring.R
# =========================================================

# ---------------------------------------------------------
# Module score computation
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
compute_module_scores <- function(expr, modules) {

  mean_score <- function(gset) {

    if (length(gset) == 0) {
      return(rep(0, ncol(expr)))
    }

    mat <- expr[gset, , drop = FALSE]

    detected <- Matrix::colSums(mat > 0)
    detected[detected == 0] <- 1

    Matrix::colSums(mat) / detected
  }

  hit_score <- function(gset) {

    if (length(gset) == 0) {
      return(rep(0, ncol(expr)))
    }

    Matrix::colSums(
      expr[gset, , drop = FALSE] > 0
    )
  }

  scores <- lapply(modules, mean_score)
  hits   <- lapply(modules, hit_score)

  list(
    scores = scores,
    hits = hits
  )
}


# ---------------------------------------------------------
# Ambient RNA estimation
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
estimate_ambient_expression <- function(
    expr,
    gene,
    lib
) {

  low_depth <- lib < 100

  if (!any(low_depth)) {
    return(0)
  }

  ambient <- Matrix::rowMeans(
    expr[, low_depth, drop = FALSE]
  )[gene]

  ifelse(is.na(ambient), 0, ambient)
}


# ---------------------------------------------------------
# Detection score
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
compute_detection_score <- function(
    gnrh,
    scores,
    ambient_ratio
) {

  (
    2.5 * safe_scale(log1p(gnrh)) +
      2.0 * safe_scale(scores$core) +
      1.0 * safe_scale(scores$migration) +
      1.5 * safe_scale(scores$neuroendo) +
      1.0 * safe_scale(log1p(ambient_ratio))
  )
}


# =========================================================
# detection_classification.R
# =========================================================

# ---------------------------------------------------------
# Cell classification
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
classify_gnrh_cells <- function(
    gnrh_raw,
    gnrh,
    score,
    hits,
    min_umi,
    min_counts,
    gnrh_threshold,
    score_quantile,
    lib
) {

  umi_ok   <- gnrh_raw >= min_umi
  depth_ok <- lib >= min_counts
  expr_ok  <- gnrh >= gnrh_threshold

  marker_ok <-
    hits$core >= 2 &
    (
      hits$migration >= 1 |
        hits$neuroendo >= 1
    )

  base_cells <-
    umi_ok &
    depth_ok &
    expr_ok &
    marker_ok

  score_threshold <- if (any(base_cells)) {

    stats::quantile(
      score[base_cells],
      probs = score_quantile,
      na.rm = TRUE
    )

  } else {

    stats::median(score, na.rm = TRUE)

  }

  keep <-
    base_cells &
    score >= score_threshold &
    score > 1

  label <- rep("gnrh_neg", length(score))

  label[
    umi_ok &
      marker_ok &
      depth_ok
  ] <- "gnrh_low"

  label[keep] <- "gnrh_high"

  list(
    label = label,
    keep = keep,
    threshold = score_threshold
  )
}



# =========================================================
# run_detection_diagnostics.R
# =========================================================

# ---------------------------------------------------------
# Diagnostic metrics
# ---------------------------------------------------------

#' Build diagnostic metrics for GnRH detection
#'
#' Computes diagnostic statistics, threshold sweeps, and
#' truth-proxy evaluation for GnRH neuron detection results.
#'
#' This function is typically called internally by
#' \code{\link{detect_gnrh_cells}}, but may also be run manually
#' after detection to regenerate diagnostics.
#'
#' Diagnostic outputs are stored in the Seurat object's
#' \code{@misc} slot:
#'
#' \itemize{
#'   \item \code{object@misc$gnrh_threshold_curve}
#'   \item \code{object@misc$gnrh_diag}
#'   \item \code{object@misc$gnrh_params}
#' }
#'
#' @param object A Seurat object containing GnRH metadata.
#' @param verbose Logical; print progress messages.
#'
#' @return
#' A Seurat object with diagnostic metrics and threshold
#' optimization results stored in \code{@misc}.
#'
#' @export
run_detection_diagnostics <- function(object, verbose = TRUE) {

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  log <- log_msg(verbose)
  log("[10] GnRH diagnostic layer...")

  meta <- object[[]]

  required <- c(
    "gnrh_expr",
    "gnrh_raw",
    "gnrh_score",
    "gnrh_class",
    "gnrh_core_hits",
    "gnrh_migration_hits",
    "gnrh_neuroendo_hits"
  )

  missing <- setdiff(required, colnames(meta))

  if (length(missing) > 0) {
    stop(sprintf(
      "Missing metadata: %s (run detect_gnrh_cells first)",
      paste(missing, collapse = ", ")
    ))
  }

  params <- object@misc$gnrh_params

  if (is.null(params)) {
    stop("Missing gnrh_params (run detect_gnrh_cells first)")
  }

  df <- data.frame(
    cell = rownames(meta),
    expr = meta$gnrh_expr,
    raw_expr = meta$gnrh_raw,
    score = meta$gnrh_score,
    class = meta$gnrh_class,
    core_hits = meta$gnrh_core_hits,
    migration_hits = meta$gnrh_migration_hits,
    neuroendo_hits = meta$gnrh_neuroendo_hits
  )

  gnrh_thr <- params$expr_threshold_strict
  score_thr <- params$score_threshold_strict
  min_umi <- params$min_umi

  df$above_threshold <-
    df$expr >= gnrh_thr &
    df$score >= score_thr

  df$keep_any <- df$class %in% c("gnrh_low", "gnrh_high")

  truth <-
    (df$core_hits >= 2) &
    ((df$migration_hits >= 1) | (df$neuroendo_hits >= 1)) &
    (df$raw_expr >= min_umi)

  object[["gnrh_truth_proxy"]] <- truth

  if (length(unique(df$expr)) > 1 &&
      sum(truth, na.rm = TRUE) > 0) {

    thresholds <- seq(
      min(df$expr, na.rm = TRUE),
      max(df$expr, na.rm = TRUE),
      length.out = 50
    )

    res <- data.frame(
      threshold = thresholds,
      sensitivity = NA_real_,
      specificity = NA_real_,
      precision = NA_real_,
      F1 = NA_real_
    )

    for (i in seq_along(thresholds)) {

      pred <- df$expr >= thresholds[i]

      TP <- sum(pred & truth, na.rm = TRUE)
      FP <- sum(pred & !truth, na.rm = TRUE)
      FN <- sum(!pred & truth, na.rm = TRUE)
      TN <- sum(!pred & !truth, na.rm = TRUE)

      sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
      prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

      res$F1[i] <- if (is.finite(prec) && is.finite(sens) && (prec + sens) > 0) {
        2 * prec * sens / (prec + sens)
      } else NA_real_

      res$sensitivity[i] <- sens
      res$specificity[i] <- spec
      res$precision[i] <- prec
    }

    best <- which.max(res$F1)

    params$best_expr_threshold <- res$threshold[best]

    object@misc$gnrh_threshold_curve <- res
    object@misc$gnrh_diag <- df
    object@misc$gnrh_params <- params
  }

  object
}



# =========================================================
# detection_pipeline.R
# =========================================================

# ---------------------------------------------------------
# Detect GnRH-expressing cells
# ---------------------------------------------------------

#' Detect GnRH-expressing cells
#'
#' Identifies GnRH-positive cells using:
#' \itemize{
#'   \item normalized GNRH1 expression
#'   \item developmental marker modules
#'   \item co-detection signatures
#'   \item ambient RNA correction
#' }
#'
#' @param object A Seurat object.
#' @param assay Assay name.
#' @param layer Data layer.
#' @param min_umi Minimum raw GNRH UMI count.
#' @param mad_factor MAD multiplier for adaptive thresholding.
#' @param score_quantile Quantile used for score cutoff.
#' @param scale_factor Library normalization scale factor.
#' @param min_counts Minimum total counts per cell.
#' @param verbose Print progress messages.
#'
#' @return Updated Seurat object.
#'
#' @export
detect_gnrh_cells <- function(
    object,
    assay = "RNA",
    layer = "counts",
    min_umi = GNRH_DEFAULT_MIN_UMI,
    mad_factor = 2,
    score_quantile = GNRH_DEFAULT_SCORE_Q,
    scale_factor = GNRH_SCALE_FACTOR,
    min_counts = 500,
    verbose = TRUE
) {

  stopifnot(inherits(object, "Seurat"))

  log <- log_msg(verbose)

  log("==== GNRH DETECTION START ====")

  # -------------------------------------------------------
  # Expression matrix
  # -------------------------------------------------------

  expr <- extract_expression_matrix(
    object = object,
    assay = assay,
    layer = layer
  )

  genes <- rownames(expr)

  # -------------------------------------------------------
  # Detection modules
  # -------------------------------------------------------

  modules <- build_detection_modules(genes)

  # -------------------------------------------------------
  # Detect GNRH gene
  # -------------------------------------------------------

  gnrh_gene <- match_gene_symbols(
    "GNRH1",
    genes
  )

  if (is.na(gnrh_gene)) {
    stop("GNRH1 not found")
  }

  # -------------------------------------------------------
  # Library sizes
  # -------------------------------------------------------

  lib <- Matrix::colSums(expr)

  # -------------------------------------------------------
  # Raw expression
  # -------------------------------------------------------

  gnrh_raw <- as.numeric(
    expr[gnrh_gene, ]
  )

  # -------------------------------------------------------
  # Normalized expression
  # -------------------------------------------------------

  gnrh <- library_normalize(
    x = gnrh_raw,
    lib = lib,
    scale_factor = scale_factor
  )

  # -------------------------------------------------------
  # Ambient correction
  # -------------------------------------------------------

  ambient <- estimate_ambient_expression(
    expr = expr,
    gene = gnrh_gene,
    lib = lib
  )

  ambient_ratio <-
    (gnrh_raw + 1) /
    (ambient + 1)

  # -------------------------------------------------------
  # Module scores
  # -------------------------------------------------------

  module_res <- compute_module_scores(
    expr = expr,
    modules = modules
  )

  # -------------------------------------------------------
  # Detection score
  # -------------------------------------------------------

  score <- compute_detection_score(
    gnrh = gnrh,
    scores = module_res$scores,
    ambient_ratio = ambient_ratio
  )

  # -------------------------------------------------------
  # Adaptive threshold
  # -------------------------------------------------------

  nonzero_expr <- gnrh[gnrh > 0]

  gnrh_threshold <- if (
    length(nonzero_expr) > 20 &&
    stats::mad(nonzero_expr) > 0
  ) {

    stats::median(nonzero_expr) +
      mad_factor * stats::mad(nonzero_expr)

  } else if (length(nonzero_expr) > 0) {

    stats::quantile(nonzero_expr, 0.99)

  } else {

    Inf

  }

  # -------------------------------------------------------
  # Cell classification
  # -------------------------------------------------------

  cls <- classify_gnrh_cells(
    gnrh_raw = gnrh_raw,
    gnrh = gnrh,
    score = score,
    hits = module_res$hits,
    min_umi = min_umi,
    min_counts = min_counts,
    gnrh_threshold = gnrh_threshold,
    score_quantile = score_quantile,
    lib = lib
  )

  # -------------------------------------------------------
  # Metadata
  # -------------------------------------------------------

  object$gnrh_class <- factor(
    cls$label,
    levels = c(
      "gnrh_low",
      "gnrh_high",
      "gnrh_neg"
    )
  )

  object$gnrh_status <- factor(
    ifelse(
      cls$label != "gnrh_neg",
      "pos",
      "neg"
    ),
    levels = c("pos", "neg")
  )

  object$gnrh_score <- score
  object$gnrh_expr  <- gnrh
  object$gnrh_raw   <- gnrh_raw

  # -------------------------------------------------------
  # Parameters
  # -------------------------------------------------------

  object@misc$gnrh_params <- list(
    expr_threshold_strict  = gnrh_threshold,
    score_threshold_strict = cls$threshold,
    min_umi = min_umi,
    min_counts = min_counts,
    score_quantile = score_quantile,
    mad_factor = mad_factor,
    scale_factor = scale_factor
  )

  # -------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------

  object <- run_detection_diagnostics(object, verbose)

  log("==== GNRH DETECTION DONE ====")

  object
}
