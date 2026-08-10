#' Detect GnRH neurons from single-cell RNA-seq data
#'
#' Identifies candidate gonadotropin-releasing hormone (GnRH) neurons in a
#' Seurat object using direct GNRH1 detection together with independent
#' transcriptomic evidence.
#'
#' Detection integrates:
#' \itemize{
#'   \item direct \code{GNRH1} expression;
#'   \item GnRH identity and specification programs;
#'   \item migration-associated programs;
#'   \item neuroendocrine maturation programs;
#'   \item hormonal responsiveness;
#'   \item alternative neuronal or neuroendocrine identity programs;
#'   \item ambient RNA information;
#'   \item neighborhood enrichment using k-nearest neighbors; and
#'   \item adaptive transcriptomic support thresholds.
#' }
#'
#' Classification uses two GnRH-positive routes: \code{direct} and
#' \code{supported}. Both routes require detectable \code{GNRH1} expression.
#'
#' Cells without detected \code{GNRH1} are never classified as GnRH-positive.
#' However, cells showing strong GnRH-like transcriptomic evidence can be
#' flagged separately as \code{gnrh_dropout_candidate} for diagnostic
#' purposes.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay used for expression extraction. Default is \code{"RNA"}.
#' @param layer Expression layer used for detection. Default is
#'   \code{"counts"}.
#' @param reduction Dimensional reduction used for kNN neighborhood support.
#'   Default is \code{"pca"}.
#' @param dims Dimensions used for neighborhood analysis. Default is
#'   \code{1:20}.
#' @param k Number of nearest neighbors. Default is 20.
#' @param min_umi Minimum raw \code{GNRH1} UMI count required for direct
#'   detection. Default is 2.
#' @param min_counts Minimum total UMI count required per cell. Default is 500.
#' @param mad_factor Multiplier applied to the MAD-based adaptive
#'   \code{GNRH1} expression threshold. Default is 2.
#' @param supported_q Quantile of the direct-cell transcriptomic support
#'   distribution used for low-expression supported candidates.
#'   Default is 0.60.
#' @param candidate_q Quantile of the direct-cell transcriptomic support
#'   distribution used only to flag GNRH1-negative transcriptomic candidates.
#'   These cells are never classified as GnRH-positive. Default is 0.95.
#' @param scale_factor Library normalization scale factor. Default is 10000.
#' @param max_alternative Maximum alternative identity score tolerated for
#'   \code{GNRH1}-dropout rescue. Default is 0.75.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A Seurat object containing GnRH classifications, scores,
#'   diagnostics, and detection parameters.
#'
#' @details
#' Direct candidates require at least \code{min_umi} raw \code{GNRH1}
#' counts.
#'
#' Supported candidates contain detectable but sub-threshold \code{GNRH1}
#' and additionally require independent GnRH identity and transcriptomic
#' support.
#'
#' Cells with no detected \code{GNRH1} are not classified as GnRH-positive,
#' regardless of their transcriptomic support score. Strong GnRH-like
#' GNRH1-negative cells may instead be flagged as
#' \code{gnrh_dropout_candidate} for exploratory or diagnostic analyses.
#'
#' Migration-associated expression contributes supportive evidence but cannot
#' independently establish GnRH identity.
#'
#' The main \code{gnrh_score} includes direct \code{GNRH1} information,
#' whereas \code{gnrh_support_score} deliberately excludes direct
#' \code{GNRH1} expression and ambient-RNA information.
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
    supported_q = 0.60,
    candidate_q = 0.95,
    scale_factor = 10000,
    max_alternative = 0.75,
    verbose = TRUE
) {

  log <- .msg(verbose)

  log(
    "==== GNRH DETECTION START ===="
  )

  # --------------------------------------------------------------------------- #
  # Validate input
  # --------------------------------------------------------------------------- #

  object <- validate_input(
    object,
    assay = assay,
    verbose = verbose
  )

  expr <- .get_expr(
    object,
    assay = assay,
    layer = layer
  )

  genes <- rownames(expr)

  log(
    sprintf(
      "Matrix loaded: %d genes by %d cells",
      nrow(expr),
      ncol(expr)
    )
  )

  # --------------------------------------------------------------------------- #
  # Resolve GNRH1
  # --------------------------------------------------------------------------- #

  gene <- .match_genes(
    c(
      "GNRH1",
      "Gnrh1",
      "gnrh1"
    ),
    genes
  )

  if (!length(gene)) {
    stop(
      "GnRH gene not found. Tried: GNRH1, Gnrh1, gnrh1",
      call. = FALSE
    )
  }

  gnrh_gene <- gene[[1]]

  log(
    "Using GnRH gene: ",
    gnrh_gene
  )

  # --------------------------------------------------------------------------- #
  # Curated modules
  # --------------------------------------------------------------------------- #

  modules <- .gnrh_modules(
    genes
  )

  alternative_modules <- .gnrh_alternative_modules(
    genes
  )

  # --------------------------------------------------------------------------- #
  # Basic expression signals
  # --------------------------------------------------------------------------- #

  lib <- Matrix::colSums(
    expr
  )

  raw <- as.numeric(
    expr[
      gnrh_gene,
      ,
      drop = TRUE
    ]
  )

  safe_lib <- pmax(
    lib,
    1
  )

  norm <- log1p(
    (raw / safe_lib) *
      scale_factor
  )

  # --------------------------------------------------------------------------- #
  # Ambient RNA
  # --------------------------------------------------------------------------- #

  amb <- .ambient(
    expr,
    gene,
    lib
  )

  amb_ratio <-
    (raw + 1) /
    (amb + 1)

  # --------------------------------------------------------------------------- #
  # Module helpers
  # --------------------------------------------------------------------------- #

  score_gene_set <- function(gene_set) {

    if (!length(gene_set)) {
      return(
        rep(
          0,
          ncol(expr)
        )
      )
    }

    Matrix::colMeans(
      expr[
        gene_set,
        ,
        drop = FALSE
      ] > 0
    )
  }

  count_gene_hits <- function(gene_set) {

    if (!length(gene_set)) {
      return(
        integer(
          ncol(expr)
        )
      )
    }

    as.integer(
      Matrix::colSums(
        expr[
          gene_set,
          ,
          drop = FALSE
        ] > 0
      )
    )
  }

  # --------------------------------------------------------------------------- #
  # Module scores
  # --------------------------------------------------------------------------- #

  identity_primary <- score_gene_set(
    modules$identity$primary
  )

  identity_supportive <- score_gene_set(
    modules$identity$supportive
  )

  migration_primary <- score_gene_set(
    modules$migration$primary
  )

  migration_supportive <- score_gene_set(
    modules$migration$supportive
  )

  neuro_primary <- score_gene_set(
    modules$neuroendocrine$primary
  )

  neuro_supportive <- score_gene_set(
    modules$neuroendocrine$supportive
  )

  hormone_supportive <- score_gene_set(
    modules$hormone$supportive
  )

  guidance_environment <- score_gene_set(
    modules$guidance_environment$supportive
  )

  # --------------------------------------------------------------------------- #
  # Weighted biological scores
  # --------------------------------------------------------------------------- #

  identity_score <-
    3.0 * identity_primary +
    1.5 * identity_supportive

  migration_score <-
    1.0 * migration_primary +
    0.5 * migration_supportive

  neuro_score <-
    1.5 * neuro_primary +
    0.75 * neuro_supportive

  hormone_score <-
    0.25 * hormone_supportive

  # --------------------------------------------------------------------------- #
  # Marker hits
  # --------------------------------------------------------------------------- #

  identity_primary_hits <- count_gene_hits(
    modules$identity$primary
  )

  identity_supportive_hits <- count_gene_hits(
    modules$identity$supportive
  )

  migration_primary_hits <- count_gene_hits(
    modules$migration$primary
  )

  migration_supportive_hits <- count_gene_hits(
    modules$migration$supportive
  )

  neuro_primary_hits <- count_gene_hits(
    modules$neuroendocrine$primary
  )

  neuro_supportive_hits <- count_gene_hits(
    modules$neuroendocrine$supportive
  )

  core_hits <-
    identity_primary_hits +
    identity_supportive_hits

  mig_hits <-
    migration_primary_hits +
    migration_supportive_hits

  neuro_hits <-
    neuro_primary_hits +
    neuro_supportive_hits

  # --------------------------------------------------------------------------- #
  # Alternative neuronal / neuroendocrine programs
  # --------------------------------------------------------------------------- #

  alternative <- .score_gnrh_alternatives(
    expr = expr,
    modules = alternative_modules
  )

  alternative_score <-
    alternative$score

  alternative_hits <-
    alternative$hits

  # --------------------------------------------------------------------------- #
  # Neighborhood support
  # --------------------------------------------------------------------------- #

  knn <- .knn_signal(
    object,
    raw,
    reduction,
    dims,
    k
  )

  # --------------------------------------------------------------------------- #
  # Main GnRH score
  #
  # Includes direct GNRH1 information.
  # --------------------------------------------------------------------------- #

  score_raw <-
    3.0 * log1p(norm) +
    identity_score +
    migration_score +
    neuro_score +
    hormone_score +
    1.0 * log1p(amb_ratio) +
    0.75 * log1p(knn)

  score <- .scale0(
    score_raw
  )

  # --------------------------------------------------------------------------- #
  # Independent transcriptomic support score
  #
  # IMPORTANT:
  #
  # - excludes direct GNRH1 expression;
  # - excludes ambient GNRH1 ratio;
  # - heavily favors GnRH identity;
  # - migration and generic neuroendocrine programs are supportive only.
  # --------------------------------------------------------------------------- #

  support_score_raw <-
    4.0 * identity_score +
    0.35 * migration_score +
    0.75 * neuro_score +
    0.10 * hormone_score +
    0.35 * log1p(knn)

  support_score <- .scale0(
    support_score_raw
  )

  # --------------------------------------------------------------------------- #
  # Adaptive GNRH1 expression threshold
  # --------------------------------------------------------------------------- #

  nz <- norm[
    is.finite(norm) &
      norm > 0
  ]

  expr_thr <- if (
    length(nz) > 20L &&
    stats::mad(nz) > 0
  ) {

    stats::median(nz) +
      mad_factor *
      stats::mad(nz)

  } else if (length(nz) > 0L) {

    as.numeric(
      stats::quantile(
        nz,
        probs = 0.99,
        names = FALSE
      )
    )

  } else {

    Inf
  }

  # --------------------------------------------------------------------------- #
  # Marker structure
  # --------------------------------------------------------------------------- #

  hits <- list(
    core = core_hits,
    mig = mig_hits,
    neuro = neuro_hits,

    identity_primary =
      identity_primary_hits,

    identity_supportive =
      identity_supportive_hits,

    migration_primary =
      migration_primary_hits,

    migration_supportive =
      migration_supportive_hits,

    neuro_primary =
      neuro_primary_hits,

    neuro_supportive =
      neuro_supportive_hits
  )

  # --------------------------------------------------------------------------- #
  # GnRH identity gates
  #
  # These gates prevent migration or generic neuroendocrine programs from
  # independently defining a GnRH candidate.
  # --------------------------------------------------------------------------- #

  identity_strong <-
    identity_primary_hits >= 2L |
    (
      identity_primary_hits >= 1L &
        identity_supportive_hits >= 2L
    )

  identity_moderate <-
    identity_primary_hits >= 1L |
    identity_supportive_hits >= 2L

  neuro_support <-
    neuro_primary_hits >= 1L |
    neuro_supportive_hits >= 2L

  migration_support <-
    migration_primary_hits >= 1L

  independent_support <-
    identity_strong |
    (
      identity_moderate &
        neuro_support
    )

  # --------------------------------------------------------------------------- #
  # Classification
  # --------------------------------------------------------------------------- #

  cls <- .classify(
    raw = raw,
    norm = norm,
    score = score,
    support_score = support_score,
    hits = hits,
    lib = lib,

    min_umi = min_umi,
    min_counts = min_counts,

    supported_q = supported_q,
    candidate_q = candidate_q,

    expr_thr = expr_thr,

    knn = knn,

    alternative_score = alternative_score,
    max_alternative = max_alternative,

    identity_strong = identity_strong,
    identity_moderate = identity_moderate,
    independent_support = independent_support
  )

  # --------------------------------------------------------------------------- #
  # Metadata: classification
  # --------------------------------------------------------------------------- #

  object$gnrh_status <- factor(
    cls$status,
    levels = c(
      "neg",
      "pos"
    )
  )

  object$gnrh_class <- factor(
    cls$class,
    levels = c(
      "neg",
      "supported",
      "direct"
    )
  )

  # Direct-cell diagnostic subclass.
  object$gnrh_direct_supported <-
    cls$direct_supported

  object$gnrh_direct_isolated <-
    cls$direct_isolated

  # --------------------------------------------------------------------------- #
  # Metadata: classification
  # --------------------------------------------------------------------------- #

  object$gnrh_status <- factor(
    cls$status,
    levels = c(
      "neg",
      "pos"
    )
  )

  object$gnrh_class <- factor(
    cls$class,
    levels = c(
      "neg",
      "supported",
      "direct"
    )
  )

  # Direct-cell diagnostic subclasses.
  object$gnrh_direct_supported <-
    cls$direct_supported

  object$gnrh_direct_isolated <-
    cls$direct_isolated

  # GNRH1-negative transcriptomic candidate.
  # Diagnostic only: these cells remain gnrh_status == "neg".
  object$gnrh_dropout_candidate <-
    cls$dropout_candidate

  # --------------------------------------------------------------------------- #
  # Metadata: scores
  # --------------------------------------------------------------------------- #

  object$gnrh_score <-
    score

  object$gnrh_support_score <-
    support_score

  # --------------------------------------------------------------------------- #
  # Metadata: scores
  # --------------------------------------------------------------------------- #

  object$gnrh_score <-
    score

  object$gnrh_support_score <-
    support_score

  object$gnrh_expr <-
    norm

  object$gnrh_raw <-
    raw

  object$gnrh_identity_score <- .scale0(
    identity_score
  )

  object$gnrh_migration_score <- .scale0(
    migration_score
  )

  object$gnrh_neuro_score <- .scale0(
    neuro_score
  )

  object$gnrh_hormone_score <- .scale0(
    hormone_score
  )

  object$gnrh_guidance_score <- .scale0(
    guidance_environment
  )

  object$gnrh_alternative_score <-
    alternative_score

  # --------------------------------------------------------------------------- #
  # Metadata: marker hits
  # --------------------------------------------------------------------------- #

  object$gnrh_core_hits <-
    core_hits

  object$gnrh_identity_primary_hits <-
    identity_primary_hits

  object$gnrh_identity_supportive_hits <-
    identity_supportive_hits

  object$gnrh_mig_hits <-
    mig_hits

  object$gnrh_migration_primary_hits <-
    migration_primary_hits

  object$gnrh_migration_supportive_hits <-
    migration_supportive_hits

  object$gnrh_neuro_hits <-
    neuro_hits

  object$gnrh_neuro_primary_hits <-
    neuro_primary_hits

  object$gnrh_neuro_supportive_hits <-
    neuro_supportive_hits

  object$gnrh_alternative_hits <-
    alternative_hits

  object$gnrh_knn <-
    knn

  # --------------------------------------------------------------------------- #
  # Metadata: biological gates
  # --------------------------------------------------------------------------- #

  object$gnrh_identity_moderate <-
    identity_moderate

  object$gnrh_identity_strong <-
    identity_strong

  object$gnrh_neuro_support <-
    neuro_support

  object$gnrh_migration_support <-
    migration_support

  object$gnrh_independent_support <-
    independent_support

  # --------------------------------------------------------------------------- #
  # Store detected GnRH gene
  # --------------------------------------------------------------------------- #

  object@misc$gnrh_gene <-
    gnrh_gene

  # --------------------------------------------------------------------------- #
  # High-confidence classification
  # --------------------------------------------------------------------------- #

  object$gnrh_confident <- .compute_gnrh_confident(
    object[[]],
    min_umi
  )

  # --------------------------------------------------------------------------- #
  # Store parameters
  # --------------------------------------------------------------------------- #

  object@misc$gnrh_params <- list(
    assay = assay,
    layer = layer,

    reduction = reduction,
    dims = dims,
    k = k,

    min_umi = min_umi,
    min_counts = min_counts,

    mad_factor = mad_factor,

    supported_q = supported_q,
    candidate_q = candidate_q,

    scale_factor = scale_factor,

    expr_thr = expr_thr,

    supported_thr =
      cls$supported_thr,

    dropout_thr =
      cls$dropout_thr,

    knn_support_thr =
      cls$knn_support_thr,

    knn_strong_thr =
      cls$knn_strong_thr,

    max_alternative =
      max_alternative,

    module_weights = c(
      identity_primary = 3.0,
      identity_supportive = 1.5,

      migration_primary = 1.0,
      migration_supportive = 0.5,

      neuroendocrine_primary = 1.5,
      neuroendocrine_supportive = 0.75,

      hormone_supportive = 0.25,

      guidance_environment = 0
    ),

    support_weights = c(
      identity = 4.0,
      migration = 0.35,
      neuroendocrine = 0.75,
      hormone = 0.10,
      knn = 0.35
    ),

    identity_rules = list(
      moderate = paste0(
        "identity_primary_hits >= 1 OR ",
        "identity_supportive_hits >= 2"
      ),

      strong = paste0(
        "identity_primary_hits >= 2 OR ",
        "(identity_primary_hits >= 1 AND ",
        "identity_supportive_hits >= 2)"
      ),

      independent_support = paste0(
        "identity_strong OR ",
        "(identity_moderate AND neuro_support)"
      )
    )
  )

  # --------------------------------------------------------------------------- #
  # Store modules
  # --------------------------------------------------------------------------- #

  object@misc$gnrh_modules <-
    modules

  object@misc$gnrh_alternative_modules <-
    alternative_modules

  object@misc$gnrh_alternative_scores <-
    alternative$scores

  # --------------------------------------------------------------------------- #
  # Classification diagnostics
  # --------------------------------------------------------------------------- #

  if (is.null(object@misc$gnrh)) {
    object@misc$gnrh <- list()
  }

  object@misc$gnrh$classify_rules <-
    cls$rules

  object@misc$gnrh$classify_summary <-
    cls$rule_summary

  # --------------------------------------------------------------------------- #
  # Diagnostics
  # --------------------------------------------------------------------------- #

  object <- gnrh_diagnostics(
    object,
    verbose = verbose
  )

  log(
    "==== GNRH DETECTION DONE ===="
  )

  object
}
