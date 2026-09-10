#' Run GnRH detection diagnostics
#' @param object A Seurat object previously processed with `detect_gnrh()`.
#' @param truth Optional metadata column name or vector with binary truth labels.
#' @param positive Value identifying the positive truth class.
#' @param predictor Metadata column used for threshold validation.
#' @param n_thresholds Number of thresholds evaluated.
#' @param verbose Print progress messages.
#' @return A Seurat object with diagnostics stored in `object@misc$gnrh`.
#' @export
gnrh_diagnostics <- function(object,truth=NULL,positive=NULL,predictor="gnrh_support_score_raw",n_thresholds=200L,verbose=TRUE) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  log <- .msg(verbose)
  log("Running diagnostics")
  md <- object[[]]
  if (is.null(object@misc$gnrh_params)) stop("Missing `gnrh_params`; run `detect_gnrh()` first.",call.=FALSE)
  req <- c("gnrh_raw","gnrh_expr","gnrh_score","gnrh_status","gnrh_class","gnrh_core_hits","gnrh_mig_hits","gnrh_neuro_hits","gnrh_confident")
  miss <- setdiff(req,colnames(md))
  if (length(miss)) stop("Missing required metadata columns: ",paste(miss,collapse=", "),call.=FALSE)
  object <- .init_gnrh_misc(object)
  object <- .make_total_hits(object)
  md <- object[[]]
  opt <- function(x,default=NA) if (x %in% colnames(md)) md[[x]] else rep(default,nrow(md))
  diagnostics <- data.frame(
    raw_expr=md$gnrh_raw,
    expr=md$gnrh_expr,
    score=md$gnrh_score,
    score_raw=opt("gnrh_score_raw",NA_real_),
    support_score=opt("gnrh_support_score",NA_real_),
    support_score_raw=opt("gnrh_support_score_raw",NA_real_),
    status=as.character(md$gnrh_status),
    class=as.character(md$gnrh_class),
    confident=as.logical(md$gnrh_confident),
    direct_signal=as.logical(opt("gnrh_direct_signal",FALSE)),
    direct_supported=as.logical(opt("gnrh_direct_supported",FALSE)),
    direct_isolated=as.logical(opt("gnrh_direct_isolated",FALSE)),
    reference_positive=as.logical(opt("gnrh_reference_positive",FALSE)),
    transcriptomic_candidate=as.logical(opt("gnrh_transcriptomic_candidate",opt("gnrh_dropout_candidate",FALSE))),
    identity_moderate=as.logical(opt("gnrh_identity_moderate",FALSE)),
    identity_strong=as.logical(opt("gnrh_identity_strong",FALSE)),
    independent_support=as.logical(opt("gnrh_independent_support",FALSE)),
    neuro_support=as.logical(opt("gnrh_neuro_support",FALSE)),
    migration_support=as.logical(opt("gnrh_migration_support",FALSE)),
    alternative_strong=as.logical(opt("gnrh_alternative_strong",FALSE)),
    knn=opt("gnrh_knn",NA_real_),
    core_hits=md$gnrh_core_hits,
    mig_hits=md$gnrh_mig_hits,
    neuro_hits=md$gnrh_neuro_hits,
    total_hits=opt("gnrh_total_hits",opt("total_hits",NA_real_)),
    stringsAsFactors=FALSE
  )
  rownames(diagnostics) <- rownames(md)
  optional_cols <- c(
    "gnrh_stage","gnrh_stage_raw","gnrh_stage_resolution","gnrh_stage_confident",
    "gnrh_stage_score","gnrh_stage_second_score","gnrh_stage_margin","gnrh_stage_reason",
    "gnrh_stage_reassigned","gnrh_stage_early_score","gnrh_stage_migrating_score",
    "gnrh_stage_mature_score","gnrh_stage_early_expression","gnrh_stage_migrating_expression",
    "gnrh_stage_mature_expression","gnrh_stage_early_fraction","gnrh_stage_migrating_fraction",
    "gnrh_stage_mature_fraction","gnrh_migration_core_hits","gnrh_secretory",
    "gnrh_secretory_supported","gnrh_secretory_core_hits","gnrh_secretory_supportive_hits",
    "gnrh_secretory_hits"
  )
  for (nm in intersect(optional_cols,colnames(md))) diagnostics[[nm]] <- md[[nm]]
  object@misc$gnrh$diagnostics <- diagnostics
  count_value <- function(x,v) sum(as.character(x)==v,na.rm=TRUE)
  count_true <- function(x) sum(x %in% TRUE,na.rm=TRUE)
  object@misc$gnrh$diagnostic_summary <- list(
    n_cells=nrow(md),
    neg=count_value(md$gnrh_status,"neg"),
    pos=count_value(md$gnrh_status,"pos"),
    direct=count_value(md$gnrh_class,"direct"),
    supported=count_value(md$gnrh_class,"supported"),
    confident=count_true(md$gnrh_confident),
    direct_signal=count_true(opt("gnrh_direct_signal",FALSE)),
    direct_supported=count_true(opt("gnrh_direct_supported",FALSE)),
    direct_isolated=count_true(opt("gnrh_direct_isolated",FALSE)),
    reference_positive=count_true(opt("gnrh_reference_positive",FALSE)),
    transcriptomic_candidate=count_true(opt("gnrh_transcriptomic_candidate",opt("gnrh_dropout_candidate",FALSE))),
    stage_resolved=count_true(opt("gnrh_stage_confident",FALSE))
  )
  if (is.null(truth)) {
    object@misc$gnrh$threshold_curve <- NULL
    object@misc$gnrh$best_threshold <- NA_real_
    object@misc$gnrh$roc <- NULL
    object@misc$gnrh$auc <- NA_real_
    object@misc$gnrh$pr_curve <- NULL
    object@misc$gnrh$auprc <- NA_real_
    object@misc$gnrh$validation <- NULL
    log("No independent truth supplied; validation metrics skipped.",type="info")
    return(object)
  }
  truth_values <- if (is.character(truth) && length(truth)==1L) {
    if (!truth %in% colnames(md)) stop("Truth column `",truth,"` not found in object metadata.",call.=FALSE)
    md[[truth]]
  } else truth
  if (length(truth_values)!=nrow(md)) stop("`truth` must contain one value per cell.",call.=FALSE)
  predictor_name <- as.character(predictor)
  if (length(predictor_name)!=1L || is.na(predictor_name) || !nzchar(predictor_name)) stop("`predictor` must be one metadata column name.",call.=FALSE)
  if (!predictor_name %in% colnames(md)) stop("Predictor column `",predictor_name,"` not found.",call.=FALSE)
  predictor_values <- suppressWarnings(as.numeric(md[[predictor_name]]))
  if (is.logical(truth_values)) {
    truth_positive <- truth_values
    positive_label <- TRUE
  } else {
    truth_chr <- as.character(truth_values)
    classes <- unique(truth_chr[!is.na(truth_chr)])
    if (length(classes)!=2L) stop("`truth` must contain exactly two non-missing classes.",call.=FALSE)
    positive_label <- if (is.null(positive)) classes[[2L]] else positive
    if (!as.character(positive_label) %in% classes) stop("`positive` value not present in truth labels.",call.=FALSE)
    truth_positive <- truth_chr==as.character(positive_label)
  }
  ok <- !is.na(truth_positive) & is.finite(predictor_values)
  predictor_values <- predictor_values[ok]
  truth_positive <- truth_positive[ok]
  if (!length(predictor_values)) stop("No usable observations remain after removing missing values.",call.=FALSE)
  if (length(unique(truth_positive))!=2L) stop("Validation requires both positive and negative truth classes.",call.=FALSE)
  if (length(n_thresholds)!=1L || is.na(n_thresholds) || !is.finite(n_thresholds) || n_thresholds<2 || n_thresholds!=floor(n_thresholds)) stop("`n_thresholds` must be an integer >= 2.",call.=FALSE)
  n_thresholds <- as.integer(n_thresholds)
  r <- range(predictor_values,finite=TRUE)
  if (!all(is.finite(r)) || diff(r)<=0) stop("Predictor has no usable range.",call.=FALSE)
  thresholds <- sort(unique(c(Inf,seq(r[1],r[2],length.out=n_thresholds),-Inf)),decreasing=TRUE)
  safe_div <- function(a,b) if (length(b)!=1L || b<=0) NA_real_ else a/b
  metric_row <- function(threshold) {
    pred <- predictor_values>=threshold
    TP <- sum(pred & truth_positive)
    FP <- sum(pred & !truth_positive)
    FN <- sum(!pred & truth_positive)
    TN <- sum(!pred & !truth_positive)
    sens <- safe_div(TP,TP+FN)
    spec <- safe_div(TN,TN+FP)
    prec <- safe_div(TP,TP+FP)
    F1 <- if (is.finite(prec) && is.finite(sens) && prec+sens>0) 2*prec*sens/(prec+sens) else NA_real_
    bal <- if (is.finite(sens) && is.finite(spec)) mean(c(sens,spec)) else NA_real_
    data.frame(
      threshold=threshold,TP=TP,FP=FP,FN=FN,TN=TN,
      sensitivity=sens,recall=sens,specificity=spec,precision=prec,F1=F1,
      balanced_accuracy=bal,fpr=safe_div(FP,FP+TN),tpr=sens,
      fp_per_100k=safe_div(FP*1e5,FP+TN)
    )
  }
  curve <- do.call(rbind,lapply(thresholds,metric_row))
  rownames(curve) <- NULL
  valid <- which(is.finite(curve$F1))
  best_i <- if (length(valid)) valid[which.max(curve$F1[valid])] else NA_integer_
  best_threshold <- if (!is.na(best_i)) curve$threshold[best_i] else NA_real_
  roc <- curve[is.finite(curve$fpr) & is.finite(curve$tpr),c("threshold","fpr","tpr","sensitivity","specificity"),drop=FALSE]
  roc <- roc[order(roc$fpr,roc$tpr),,drop=FALSE]
  roc <- roc[!duplicated(roc[,c("fpr","tpr"),drop=FALSE]),,drop=FALSE]
  auc <- if (nrow(roc)>=2L) sum(diff(roc$fpr)*(head(roc$tpr,-1L)+tail(roc$tpr,-1L))/2) else NA_real_
  pr <- curve[is.finite(curve$recall) & is.finite(curve$precision),c("threshold","recall","precision"),drop=FALSE]
  if (nrow(pr)) {
    pr <- stats::aggregate(precision~recall,data=pr,FUN=max)
    pr <- pr[order(pr$recall),,drop=FALSE]
  }
  auprc <- if (nrow(pr)>=2L) sum(diff(pr$recall)*(head(pr$precision,-1L)+tail(pr$precision,-1L))/2) else NA_real_
  actual <- as.character(md$gnrh_status[ok])=="pos"
  TP <- sum(actual & truth_positive)
  FP <- sum(actual & !truth_positive)
  FN <- sum(!actual & truth_positive)
  TN <- sum(!actual & !truth_positive)
  sensitivity <- safe_div(TP,TP+FN)
  specificity <- safe_div(TN,TN+FP)
  precision <- safe_div(TP,TP+FP)
  F1 <- if (is.finite(precision) && is.finite(sensitivity) && precision+sensitivity>0) 2*precision*sensitivity/(precision+sensitivity) else NA_real_
  balanced_accuracy <- if (is.finite(sensitivity) && is.finite(specificity)) mean(c(sensitivity,specificity)) else NA_real_
  validation <- list(
    predictor=predictor_name,
    positive=positive_label,
    n=length(truth_positive),
    n_positive=sum(truth_positive),
    n_negative=sum(!truth_positive),
    TP=TP,FP=FP,FN=FN,TN=TN,
    sensitivity=sensitivity,
    recall=sensitivity,
    specificity=specificity,
    precision=precision,
    F1=F1,
    balanced_accuracy=balanced_accuracy,
    false_positive_rate=safe_div(FP,FP+TN),
    false_positives_per_100k=safe_div(FP*1e5,FP+TN),
    auc=auc,
    auprc=auprc,
    best_threshold=best_threshold,
    prevalence=mean(truth_positive)
  )
  object@misc$gnrh$threshold_curve <- curve
  object@misc$gnrh$best_threshold <- best_threshold
  object@misc$gnrh$roc <- roc
  object@misc$gnrh$auc <- auc
  object@misc$gnrh$pr_curve <- pr
  object@misc$gnrh$auprc <- auprc
  object@misc$gnrh$validation <- validation
  log(sprintf("Validation: AUC %.3f | AUPRC %.3f | F1 %.3f",auc,auprc,F1),type="info")
  object
}
