# ============================================================================= #
# GnRHcell pipeline
# ============================================================================= #
#' Run complete GnRHcell analysis pipeline
#' @param object A Seurat object.
#' @param detect Run GnRH detection.
#' @param stage Run developmental staging.
#' @param diagnostics Run diagnostics.
#' @param detect_args Named arguments passed to `detect_gnrh()`.
#' @param stage_args Named arguments passed to `stage_gnrh()`.
#' @param diagnostic_args Named arguments passed to `gnrh_diagnostics()`.
#' @param summary_level Summary level: concise, detailed, or none.
#' @param verbose Print progress messages.
#' @param ... Additional arguments passed to `detect_gnrh()`.
#' @return Updated Seurat object.
#' @export
run_gnrh <- function(object,detect=TRUE,stage=TRUE,diagnostics=TRUE,detect_args=list(),stage_args=list(),diagnostic_args=list(),summary_level=c("concise","detailed","none"),verbose=TRUE,...) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  summary_level <- match.arg(summary_level)
  args <- list(detect_args=detect_args,stage_args=stage_args,diagnostic_args=diagnostic_args)
  for (nm in names(args)) {
    x <- args[[nm]]
    if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x)))))) stop("`",nm,"` must be a named list.",call.=FALSE)
  }
  object <- .init_gnrh_misc(object)
  log <- .msg(verbose)
  step_times <- list()
  total_start <- Sys.time()
  log("==== STARTING GnRHcell PIPELINE ====",type="header")
  legacy <- list(...)
  dup <- intersect(names(legacy),names(detect_args))
  if (length(dup)) legacy[dup] <- NULL
  detect_args <- c(legacy,detect_args)
  detect_args$verbose <- NULL
  stage_args$verbose <- NULL
  diagnostic_args$verbose <- NULL
  active <- c(detect=isTRUE(detect),stage=isTRUE(stage),diagnostics=isTRUE(diagnostics))
  n_steps <- sum(active)
  i <- 0L
  step_label <- function(x) {
    i <<- i+1L
    sprintf("[%d/%d] %s",i,n_steps,x)
  }
  time_step <- function(fun,args,label,key,verbose_fun=FALSE) {
    log(step_label(label),type="step")
    t0 <- Sys.time()
    object <<- do.call(fun,c(list(object=object,verbose=verbose_fun),args))
    step_times[[key]] <<- as.numeric(difftime(Sys.time(),t0,units="secs"))
    log(paste0(sub("ing.*$","",label)," complete."),type="done",duration=step_times[[key]])
  }
  if (isTRUE(detect)) {
    log(step_label("Detecting GnRH cells"),type="step")
    t0 <- Sys.time()
    object <- do.call(detect_gnrh,c(list(object=object,verbose=verbose),detect_args))
    step_times$detect_sec <- as.numeric(difftime(Sys.time(),t0,units="secs"))
    log("Detection complete.",type="done",duration=step_times$detect_sec)
  }
  if (isTRUE(stage)) {
    if (!"gnrh_status" %in% colnames(object[[]])) stop("Developmental staging requires GnRH detection metadata.",call.=FALSE)
    log(step_label("Assigning developmental stages"),type="step")
    t0 <- Sys.time()
    object <- do.call(stage_gnrh,c(list(object=object,verbose=FALSE),stage_args))
    step_times$stage_sec <- as.numeric(difftime(Sys.time(),t0,units="secs"))
    log("Staging complete.",type="done",duration=step_times$stage_sec)
  }
  if (isTRUE(diagnostics)) {
    if (!"gnrh_status" %in% colnames(object[[]])) stop("Diagnostics require GnRH detection metadata.",call.=FALSE)
    log(step_label("Running diagnostics"),type="step")
    t0 <- Sys.time()
    object <- do.call(gnrh_diagnostics,c(list(object=object,verbose=FALSE),diagnostic_args))
    step_times$diagnostics_sec <- as.numeric(difftime(Sys.time(),t0,units="secs"))
    log("Diagnostics complete.",type="done",duration=step_times$diagnostics_sec)
  }
  step_times$total_sec <- as.numeric(difftime(Sys.time(),total_start,units="secs"))
  params <- list(
    detect=isTRUE(detect),
    stage=isTRUE(stage),
    diagnostics=isTRUE(diagnostics),
    detect_args=detect_args,
    stage_args=stage_args,
    diagnostic_args=diagnostic_args
  )
  object <- .add_run_info(object,step_times=step_times,params=params)
  md <- object[[]]
  count_true <- function(x) if (x %in% colnames(md)) sum(md[[x]] %in% TRUE,na.rm=TRUE) else NA_integer_
  count_value <- function(x,v) if (x %in% colnames(md)) sum(as.character(md[[x]])==v,na.rm=TRUE) else 0L
  print_table <- function(x,title,exclude=NULL,order=NULL) {
    x <- as.character(x)
    x <- x[!is.na(x) & (is.null(exclude) || !x %in% exclude)]
    if (!length(x)) return(invisible(NULL))
    tab <- table(x)
    if (!is.null(order)) tab <- tab[c(intersect(order,names(tab)),setdiff(names(tab),order))]
    log(paste0(title,":"))
    for (nm in names(tab)) log(sprintf("  %s: %d",nm,tab[[nm]]))
    invisible(tab)
  }
  if (summary_level!="none") log("PIPELINE SUMMARY",type="info")
  # ------------------------------------------------------------------------- #
  # Concise summary
  # ------------------------------------------------------------------------- #
  if (summary_level=="concise") {
    total <- nrow(md)
    pos <- count_value("gnrh_status","pos")
    conf <- count_true("gnrh_confident")
    log(sprintf("  Cells: %s | GnRH+: %s (%.1f%%) | high confidence: %s",format(total,big.mark=","),format(pos,big.mark=","),if (total) 100*pos/total else 0,format(conf,big.mark=",")))
    direct <- if ("gnrh_class" %in% colnames(md)) count_value("gnrh_class","direct") else count_true("gnrh_direct_supported")
    supported <- count_value("gnrh_class","supported")
    isolated <- count_true("gnrh_direct_isolated")
    log(sprintf("  Evidence: direct %s | transcriptomic %s | isolated GNRH1 signal %s",format(direct,big.mark=","),format(supported,big.mark=","),format(isolated,big.mark=",")))
    if ("gnrh_stage" %in% colnames(md) && pos>0L) {
      lev <- c("early","migrating","post-migratory","mature","transitional")
      st <- table(factor(as.character(md$gnrh_stage),levels=lev))
      log(sprintf("  Developmental state: early %d | migrating %d | post-migratory %d | mature %d | transitional %d",st["early"],st["migrating"],st["post-migratory"],st["mature"],st["transitional"]))
    }
    if ("gnrh_stage_resolution" %in% colnames(md) && pos>0L) {
      positive <- as.character(md$gnrh_status)=="pos"
      resolution <- as.character(md$gnrh_stage_resolution)
      resolved <- sum(positive & resolution=="resolved",na.rm=TRUE)
      transitional <- sum(positive & resolution=="transitional",na.rm=TRUE)
      log(sprintf("  Stage resolution: resolved %d (%.1f%%) | transitional %d (%.1f%%)",resolved,100*resolved/pos,transitional,100*transitional/pos))
    }
  }
  # ------------------------------------------------------------------------- #
  # Detailed summary
  # ------------------------------------------------------------------------- #
  if (summary_level=="detailed") {
    if ("gnrh_status" %in% colnames(md)) print_table(md$gnrh_status,"Status")
    if ("gnrh_class" %in% colnames(md)) print_table(md$gnrh_class,"Class")
    evidence <- c(
      gnrh_direct_signal="GNRH1 direct signal",
      gnrh_direct_supported="identity-supported direct",
      gnrh_direct_isolated="isolated GNRH1 signal",
      gnrh_reference_positive="high-specificity reference",
      gnrh_transcriptomic_candidate="transcriptomic candidates",
      gnrh_confident="high-confidence GnRH"
    )
    available <- intersect(names(evidence),colnames(md))
    if (length(available)) {
      log("Detection evidence:")
      for (nm in available) log(sprintf("  %s: %d",evidence[[nm]],count_true(nm)))
    }
    if ("gnrh_stage" %in% colnames(md)) print_table(md$gnrh_stage,"Stage",exclude="non-gnrh",order=c("early","migrating","post-migratory","mature","transitional","undetermined"))
    if ("gnrh_stage_resolution" %in% colnames(md)) {
      x <- as.character(md$gnrh_stage_resolution)
      if ("gnrh_status" %in% colnames(md)) x <- x[as.character(md$gnrh_status)=="pos"] else x <- x[x!="non-gnrh"]
      x <- x[!is.na(x) & x!="non-gnrh"]
      if (length(x)) {
        tab <- table(x)
        n <- sum(tab)
        log("Stage resolution:")
        for (nm in names(tab)) log(sprintf("  %s: %d (%.1f%%)",nm,tab[[nm]],100*tab[[nm]]/n))
      }
    }
    if (all(c("gnrh_stage","gnrh_stage_resolution") %in% colnames(md))) {
      st <- as.character(md$gnrh_stage)
      rs <- as.character(md$gnrh_stage_resolution)
      ok <- !is.na(st) & !is.na(rs) & st!="non-gnrh" & rs!="non-gnrh"
      if (any(ok)) {
        tab <- table(st[ok],rs[ok])
        resolved <- if ("resolved" %in% colnames(tab)) tab[,"resolved"] else setNames(rep(0,nrow(tab)),rownames(tab))
        totals <- rowSums(tab)
        order <- c("early","migrating","post-migratory","mature","transitional","undetermined")
        order <- c(intersect(order,rownames(tab)),setdiff(rownames(tab),order))
        log("Stage confidence:")
        for (nm in order) log(sprintf("  %s: %d/%d resolved (%.1f%%)",nm,resolved[[nm]],totals[[nm]],100*resolved[[nm]]/totals[[nm]]))
      }
    }
    if ("gnrh_secretory" %in% colnames(md)) print_table(md$gnrh_secretory,"Secretory",exclude="non-gnrh")
    validation <- object@misc$gnrh$validation
    if (is.list(validation)) {
      metric <- function(x) {
        z <- validation[[x]]
        if (is.null(z) || length(z)!=1L || !is.finite(z)) NA_real_ else as.numeric(z)
      }
      vals <- c(auc=metric("auc"),auprc=metric("auprc"),precision=metric("precision"),recall=metric("recall"),F1=metric("F1"),fp100k=metric("false_positives_per_100k"))
      if (any(is.finite(vals))) {
        log("Validation:")
        if (is.finite(vals["auc"])) log(sprintf("  ROC AUC: %.3f",vals["auc"]))
        if (is.finite(vals["auprc"])) log(sprintf("  PR AUC: %.3f",vals["auprc"]))
        if (is.finite(vals["precision"])) log(sprintf("  Precision: %.3f",vals["precision"]))
        if (is.finite(vals["recall"])) log(sprintf("  Recall: %.3f",vals["recall"]))
        if (is.finite(vals["F1"])) log(sprintf("  F1: %.3f",vals["F1"]))
        if (is.finite(vals["fp100k"])) log(sprintf("  False positives / 100k: %.1f",vals["fp100k"]))
      }
    }
  }
  log("==== GnRHcell PIPELINE COMPLETE ====",type="done",duration=step_times$total_sec)
  object
}
