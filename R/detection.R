#' Detect GnRH neurons from single-cell RNA-seq data
#'
#' Identifies candidate gonadotropin-releasing hormone (GnRH) neurons
#' in a Seurat object using a biologically informed multi-signal framework.
#'
#' Detection integrates:
#' \itemize{
#'   \item direct \code{GNRH1} expression
#'   \item GnRH-associated marker module scoring
#'   \item migration marker support
#'   \item neuroendocrine marker support
#'   \item ambient RNA correction
#'   \item neighborhood enrichment using k-nearest neighbors
#'   \item adaptive thresholding of a composite detection score
#' }
#'
#' The method combines transcript abundance, marker co-detection,
#' local transcriptomic neighborhood structure, and contamination-aware
#' scoring to improve detection of rare GnRH neurons in sparse
#' single-cell datasets.
#'
#' Results are written into object metadata and diagnostics.
#'
#' Added metadata columns include:
#' \describe{
#'   \item{\code{gnrh_status}}{Binary GnRH classification (\code{neg}, \code{pos}).}
#'   \item{\code{gnrh_class}}{Internal classification labels.}
#'   \item{\code{gnrh_score}}{Composite GnRH detection score.}
#'   \item{\code{gnrh_expr}}{Normalized \code{GNRH1} expression.}
#'   \item{\code{gnrh_raw}}{Raw \code{GNRH1} UMI counts.}
#'   \item{\code{gnrh_core_hits}}{Number of detected core GnRH markers.}
#'   \item{\code{gnrh_mig_hits}}{Number of migration marker hits.}
#'   \item{\code{gnrh_neuro_hits}}{Number of neuroendocrine marker hits.}
#'   \item{\code{gnrh_knn}}{Neighborhood support score.}
#'   \item{\code{gnrh_truth}}{High-confidence truth classification.}
#' }
#'
#' Detection parameters and classification diagnostics are stored in
#' \code{object@misc}.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay used for expression extraction. Default is \code{"RNA"}.
#' @param layer Expression layer used for detection. Default is \code{"counts"}.
#' @param reduction Dimensional reduction used for kNN neighborhood support.
#' Default is \code{"pca"}.
#' @param dims Dimensions used for neighborhood analysis. Default is \code{1:20}.
#' @param k Number of nearest neighbors used for neighborhood support.
#' Default is 20.
#' @param min_umi Minimum raw \code{GNRH1} UMI count required for detection.
#' Default is 2.
#' @param min_counts Minimum total UMI count required per cell.
#' Default is 500.
#' @param mad_factor Multiplier applied to MAD-based adaptive expression threshold.
#' Default is 2.
#' @param score_q Quantile used to define adaptive composite score threshold.
#' Default is 0.9.
#' @param scale_factor Library normalization scale factor.
#' Default is 10000.
#' @param verbose Logical; print progress messages. Default is \code{TRUE}.
#'
#' @return A Seurat object updated with GnRH detection metadata,
#' diagnostics, and stored detection parameters.
#'
#' @details
#' If \code{GNRH1} is not found, gene aliases are searched
#' (\code{GNRH1}, \code{Gnrh1}, \code{gnrh1}).
#'
#' Composite scoring combines normalized expression, marker module
#' enrichment, ambient correction, and neighborhood support.
#'
#' Diagnostic plots and summary outputs are generated via
#' \code{\link{gnrh_diagnostics}}.
#'
#' @seealso
#' \code{\link{gnrh_diagnostics}},
#' \code{\link{stage_gnrh}},
#' \code{\link{run_gnrh}}
#'
#' @export
detect_gnrh <- function(
    object,
    assay = "RNA",
    layer = "counts",
    reduction = "pca",
    dims = 1:20,
    k = 20,
    min_umi = 2,
    min_counts = 500,
    mad_factor = 2,
    score_q = 0.9,
    scale_factor = 10000,
    verbose = TRUE
) {

  log <- .msg(verbose)
  log("==== GNRH DETECTION START ====")

  object <- validate_input(object, assay = assay, verbose = verbose)

  expr <- .get_expr(object, assay = assay, layer = layer)
  genes <- rownames(expr)

  log(sprintf("Matrix loaded: %d genes by %d cells", nrow(expr), ncol(expr)))

  gene <- .match_genes(c("GNRH1", "Gnrh1", "gnrh1"), genes)

  if (!length(gene)) {
    stop("GnRH gene not found. Tried: GNRH1, Gnrh1")
  }

  gnrh_gene <- gene[1]

  modules <- .gnrh_modules(genes)

  # core signals
  lib  <- Matrix::colSums(expr)
  raw  <- as.numeric(expr[gnrh_gene, ])
  norm <- log1p((raw / lib) * scale_factor)

  amb <- .ambient(expr, gene, lib)
  amb_ratio <- (raw + 1) / (amb + 1)

  mod  <- .score_modules(expr, modules)
  knn  <- .knn_signal(object, raw, reduction, dims, k)

  # unified score
  score <- .scale0(
    2.5 * log1p(norm) +
      2.0 * mod$score$core +
      1.0 * mod$score$mig +
      1.5 * mod$score$neuro +
      1.0 * log1p(amb_ratio) +
      0.75 * log1p(knn)
  )

  # expression threshold
  nz <- norm[norm > 0]

  expr_thr <- if (length(nz) > 20 && stats::mad(nz) > 0) {
    stats::median(nz) + mad_factor * stats::mad(nz)
  } else if (length(nz)) {
    stats::quantile(nz, 0.99)
  } else {
    Inf
  }

  # classification
  cls <- .classify(
    raw = raw,
    norm = norm,
    score = score,
    hits = mod$hits,
    lib = lib,
    min_umi = min_umi,
    min_counts = min_counts,
    score_q = score_q,
    expr_thr = expr_thr,
    knn = knn
  )

  # metadata assignment (clean + consistent)
  object$gnrh_status <- factor(
    ifelse(cls$status != "neg", "pos", "neg"),
    levels = c("neg", "pos")
  )

  object$gnrh_score <- score
  object$gnrh_expr  <- norm
  object$gnrh_raw   <- raw
  object@misc$gnrh_gene <- gnrh_gene

  object$gnrh_core_hits  <- mod$hits$core
  object$gnrh_mig_hits   <- mod$hits$mig
  object$gnrh_neuro_hits <- mod$hits$neuro
  object$gnrh_knn        <- knn

  # truth (single source of truth)
  object$gnrh_truth <- .compute_gnrh_truth(object[[]], min_umi)

  # params + diagnostics
  object@misc$gnrh_params <- list(
    min_umi = min_umi,
    min_counts = min_counts,
    score_q = score_q,
    expr_thr = expr_thr,
    score_thr = cls$thr,
    marker_rules = c(
      "core >= 2 & GNRH UMI >= min_umi",
      "core >= 2 & migration >= 1",
      "core >= 1 & neuroendocrine >= 2 & GNRH UMI >= min_umi",
      "core >= 2 & KNN support > 0.05"
    )
  )

  object@misc$gnrh$classify_rules <- cls$rules
  object@misc$gnrh$classify_summary <- cls$rule_summary

  object <- gnrh_diagnostics(object, verbose = verbose)

  log("==== GNRH DETECTION DONE ====")
  object
}
