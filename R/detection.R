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
#' \code{supported}. Both routes require detectable \code{GNRH1} expression
#' together with independent GnRH identity evidence.
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
#' @param min_reference_cells Minimum number of direct-signal cells required
#'   to construct empirical reference distributions used for detection
#'   thresholds. Default is \code{20L}.
#' @param max_alternative Maximum alternative identity score tolerated for
#'   \code{GNRH1}-dropout rescue. Default is 0.75.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A Seurat object containing GnRH classifications, scores,
#'   diagnostics, and detection parameters.
#'
#' @details
#' Direct candidates require at least \code{min_umi} raw \code{GNRH1} counts
#' and moderate independent GnRH identity evidence. Cells meeting the UMI rule
#' without identity evidence are retained as \code{gnrh_direct_isolated} and
#' \code{gnrh_direct_signal}, but remain GnRH-negative.
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
#' Detect GnRH neurons from single-cell RNA-seq data
#'
#' @export
detect_gnrh <- function(object,
                        assay="RNA",
                        layer="counts",
                        reduction="pca",
                        dims=1:20,
                        k=20,
                        min_umi=2,
                        min_counts=500,
                        mad_factor=2,
                        supported_q=0.60,
                        candidate_q=0.95,
                        min_reference_cells=20L,
                        scale_factor=10000,
                        max_alternative=0.75,
                        verbose=TRUE) {
  log <- .msg(verbose)
  object <- validate_input(object,assay=assay,required_layers=layer,verbose=verbose)
  expr <- .get_expr(object,assay=assay,layer=layer)
  genes <- rownames(expr)
  log(sprintf("Matrix loaded: %d genes by %d cells",nrow(expr),ncol(expr)))
  gene <- .match_genes(c("GNRH1","Gnrh1","gnrh1"),genes)
  if (!length(gene)) stop("GnRH gene not found. Tried: GNRH1, Gnrh1, gnrh1",call.=FALSE)
  gnrh_gene <- gene[[1L]]
  log("Using GnRH gene:",gnrh_gene)
  modules <- .gnrh_modules(genes)
  alternative_modules <- .gnrh_alternative_modules(genes)
  lib <- Matrix::colSums(expr)
  raw <- as.numeric(expr[gnrh_gene,,drop=TRUE])
  norm <- log1p(raw/pmax(lib,1)*scale_factor)
  ambient <- .ambient(expr,gnrh_gene,lib)
  ambient_ratio <- (raw+1)/(ambient+1)
  score_set <- function(g) if (!length(g)) rep(0,ncol(expr)) else Matrix::colMeans(expr[g,,drop=FALSE]>0)
  count_hits <- function(g) if (!length(g)) integer(ncol(expr)) else as.integer(Matrix::colSums(expr[g,,drop=FALSE]>0))
  identity_primary <- score_set(modules$identity$primary)
  identity_supportive <- score_set(modules$identity$supportive)
  migration_primary <- score_set(modules$migration$primary)
  migration_supportive <- score_set(modules$migration$supportive)
  neuro_primary <- score_set(modules$neuroendocrine$primary)
  neuro_supportive <- score_set(modules$neuroendocrine$supportive)
  hormone_supportive <- score_set(modules$hormone$supportive)
  guidance_environment <- score_set(modules$guidance_environment$supportive)
  identity_score <- 3*identity_primary+1.5*identity_supportive
  migration_score <- migration_primary+0.5*migration_supportive
  neuro_score <- 1.5*neuro_primary+0.75*neuro_supportive
  hormone_score <- 0.25*hormone_supportive
  identity_primary_hits <- count_hits(modules$identity$primary)
  identity_supportive_hits <- count_hits(modules$identity$supportive)
  migration_primary_hits <- count_hits(modules$migration$primary)
  migration_supportive_hits <- count_hits(modules$migration$supportive)
  neuro_primary_hits <- count_hits(modules$neuroendocrine$primary)
  neuro_supportive_hits <- count_hits(modules$neuroendocrine$supportive)
  core_hits <- identity_primary_hits+identity_supportive_hits
  mig_hits <- migration_primary_hits+migration_supportive_hits
  neuro_hits <- neuro_primary_hits+neuro_supportive_hits
  identity_moderate <- identity_primary_hits>=1L | identity_supportive_hits>=2L
  identity_strong <- identity_primary_hits>=2L | (identity_primary_hits>=1L & identity_supportive_hits>=2L)
  neuro_support <- neuro_primary_hits>=1L | neuro_supportive_hits>=2L
  migration_support <- migration_primary_hits>=1L
  independent_support <- identity_strong | (identity_moderate & neuro_support)
  alternative <- .score_gnrh_alternatives(expr=expr,modules=alternative_modules)
  alternative_score <- alternative$score
  alternative_hits <- alternative$hits
  alternative_strong <- alternative$strong
  knn_identity <- .knn_signal(object,as.numeric(identity_moderate),reduction=reduction,dims=dims,k=k)
  score_raw <- 3*norm+identity_score+migration_score+neuro_score+hormone_score+0.75*knn_identity
  score <- .scale0(score_raw)
  support_score_raw <- 4*identity_score+0.35*migration_score+0.75*neuro_score+0.10*hormone_score+0.35*knn_identity
  support_score <- .scale0(support_score_raw)
  nz <- norm[is.finite(norm) & norm>0]
  expr_thr <- if (length(nz)>20L && stats::mad(nz)>0) stats::median(nz)+mad_factor*stats::mad(nz) else if (length(nz)) as.numeric(stats::quantile(nz,0.99,names=FALSE)) else Inf
  hits <- list(
    core=core_hits,mig=mig_hits,neuro=neuro_hits,
    identity_primary=identity_primary_hits,identity_supportive=identity_supportive_hits,
    migration_primary=migration_primary_hits,migration_supportive=migration_supportive_hits,
    neuro_primary=neuro_primary_hits,neuro_supportive=neuro_supportive_hits
  )
  cls <- .classify(
    raw=raw,norm=norm,score=score,support_score=support_score_raw,hits=hits,lib=lib,
    min_umi=min_umi,min_counts=min_counts,supported_q=supported_q,candidate_q=candidate_q,
    min_reference_cells=min_reference_cells,expr_thr=expr_thr,knn=knn_identity,
    alternative_score=alternative_score,alternative_strong=alternative_strong,
    max_alternative=max_alternative,identity_strong=identity_strong,
    identity_moderate=identity_moderate,independent_support=independent_support
  )
  object$gnrh_status <- factor(cls$status,levels=c("neg","pos"))
  object$gnrh_class <- factor(cls$class,levels=c("neg","supported","direct"))
  object$gnrh_direct_supported <- cls$direct_supported
  object$gnrh_direct_isolated <- cls$direct_isolated
  object$gnrh_direct_signal <- cls$direct_signal
  object$gnrh_signal_status <- factor(ifelse(cls$direct_signal,"signal","no_signal"),levels=c("no_signal","signal"))
  object$gnrh_transcriptomic_candidate <- cls$transcriptomic_candidate
  object$gnrh_dropout_candidate <- cls$transcriptomic_candidate
  object$gnrh_reference_positive <- cls$reference_positive
  object$gnrh_score <- score
  object$gnrh_score_raw <- score_raw
  object$gnrh_support_score <- support_score
  object$gnrh_support_score_raw <- support_score_raw
  object$gnrh_expr <- norm
  object$gnrh_raw <- raw
  object$gnrh_identity_score <- .scale0(identity_score)
  object$gnrh_migration_score <- .scale0(migration_score)
  object$gnrh_neuro_score <- .scale0(neuro_score)
  object$gnrh_hormone_score <- .scale0(hormone_score)
  object$gnrh_guidance_score <- .scale0(guidance_environment)
  object$gnrh_ambient_ratio <- ambient_ratio
  object$gnrh_core_hits <- core_hits
  object$gnrh_identity_primary_hits <- identity_primary_hits
  object$gnrh_identity_supportive_hits <- identity_supportive_hits
  object$gnrh_mig_hits <- mig_hits
  object$gnrh_migration_primary_hits <- migration_primary_hits
  object$gnrh_migration_supportive_hits <- migration_supportive_hits
  object$gnrh_neuro_hits <- neuro_hits
  object$gnrh_neuro_primary_hits <- neuro_primary_hits
  object$gnrh_neuro_supportive_hits <- neuro_supportive_hits
  object$gnrh_alternative_score <- alternative_score
  object$gnrh_alternative_hits <- alternative_hits
  object$gnrh_alternative_strong <- alternative_strong
  object$gnrh_knn <- knn_identity
  object$gnrh_identity_moderate <- identity_moderate
  object$gnrh_identity_strong <- identity_strong
  object$gnrh_neuro_support <- neuro_support
  object$gnrh_migration_support <- migration_support
  object$gnrh_independent_support <- independent_support
  object$gnrh_confident <- .compute_gnrh_confident(object[[]],min_umi)
  object@misc$gnrh_gene <- gnrh_gene
  object@misc$gnrh_params <- list(
    assay=assay,layer=layer,reduction=reduction,dims=dims,k=k,
    min_umi=min_umi,min_counts=min_counts,mad_factor=mad_factor,
    supported_q=supported_q,candidate_q=candidate_q,min_reference_cells=min_reference_cells,
    scale_factor=scale_factor,expr_thr=expr_thr,
    supported_thr=cls$supported_thr,candidate_thr=cls$candidate_thr,
    dropout_thr=cls$candidate_thr,reference_n=cls$reference_n,
    knn_support_thr=cls$knn_support_thr,knn_strong_thr=cls$knn_strong_thr,
    max_alternative=max_alternative,
    module_weights=.gnrh_module_weights(),
    support_weights=c(identity=4,migration=0.35,neuroendocrine=0.75,hormone=0.10,knn_identity=0.35),
    identity_rules=list(
      moderate="identity_primary_hits >= 1 OR identity_supportive_hits >= 2",
      strong="identity_primary_hits >= 2 OR (identity_primary_hits >= 1 AND identity_supportive_hits >= 2)",
      independent_support="identity_strong OR (identity_moderate AND neuro_support)"
    )
  )
  object@misc$gnrh_modules <- modules
  object@misc$gnrh_alternative_modules <- alternative_modules
  object@misc$gnrh_alternative_scores <- alternative$scores
  object <- .init_gnrh_misc(object)
  object@misc$gnrh$classify_rules <- cls$rules
  object@misc$gnrh$classify_summary <- cls$rule_summary
  object
}
