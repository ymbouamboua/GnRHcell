#' Evaluate sensitivity of GnRH detection to parameter choices
#'
#' @param object A Seurat object previously processed by [detect_gnrh()] or
#'   [run_gnrh()]. Raw detection inputs must still be available.
#' @param grid Data frame with numeric columns chosen from `min_umi`,
#'   `supported_q`, and `k`. Each row is a complete trial; unlisted parameters
#'   use the recorded detection settings. `NULL` varies the recorded `min_umi`
#'   by -1/0/+1 (minimum 1) and `supported_q` by -0.1/0/+0.1, bounded by 0 and
#'   the recorded `candidate_q`. Duplicate rows are removed. Include `k` to
#'   assess neighborhood sensitivity. Detection is rerun for every trial.
#' @param detect_args Named list overriding recorded detection arguments for
#'   the baseline and every trial, except arguments varied in `grid`.
#' @param sample.by Optional donor/sample metadata field for per-trial summaries.
#' @param base_size Plot font size.
#' @param verbose Print progress. Computation is sequential; no Seurat objects
#'   are retained per trial. Runtime scales with the number of grid rows.
#' @return A list with `plot`, named `plots`, `runs` (counts, Jaccard overlap,
#'   class changes and reference-cell counts per trial), `cells` (positive and
#'   class-agreement frequencies across trials), `samples` (when requested),
#'   `baseline` (recomputed classes and comparison to stored labels), and
#'   `parameters`. The input object is not modified.
#' @details The baseline is recomputed on current cells using the recorded
#'   detection settings plus `detect_args`, ensuring correct behavior after
#'   subsetting. Adaptive thresholds are recalibrated separately in every run.
#'   Jaccard is defined as 1 when both positive sets are empty; this does not
#'   demonstrate biological accuracy. Frequencies measure sensitivity to this
#'   particular grid, not probabilities, confidence intervals or independent
#'   validation. Original labels are compared to the recomputed baseline and
#'   disagreements are returned explicitly. Every trial must succeed: invalid
#'   or failed runs stop instead of silently changing the frequency denominator.
#' @examples
#' \dontrun{
#' stability <- gnrh_stability(obj, sample.by="donor")
#' stability$plot
#' stability$runs
#' grid <- expand.grid(min_umi=c(1,2,3), supported_q=c(.5,.6,.8), k=c(15,30))
#' stability <- gnrh_stability(obj, grid=grid)
#' }
#' @export
gnrh_stability <- function(object,grid=NULL,detect_args=list(),sample.by=NULL,
                           base_size=11,verbose=TRUE) {
  md <- gnrh_metadata(object)
  if (!all(c("gnrh_status","gnrh_class") %in% names(md)) || !is.list(object@misc$gnrh_params))
    stop("Run GnRH detection first; stored detection settings are required.",call.=FALSE)
  if (!is.list(detect_args) || (length(detect_args) && (is.null(names(detect_args)) ||
      anyNA(names(detect_args)) || any(!nzchar(names(detect_args))) || anyDuplicated(names(detect_args)))))
    stop("`detect_args` must be a uniquely named list.",call.=FALSE)
  allowed <- setdiff(names(formals(detect_gnrh)),c("object","verbose"))
  if (any(!names(detect_args) %in% allowed)) stop("Unknown `detect_args`.",call.=FALSE)
  params <- object@misc$gnrh_params[intersect(allowed,names(object@misc$gnrh_params))]
  params[names(detect_args)] <- detect_args
  # Resolve defaults for older saved objects.
  for (nm in setdiff(allowed,names(params))) params[[nm]] <- eval(formals(detect_gnrh)[[nm]],envir=environment(detect_gnrh))
  if (!is.null(sample.by) && (!is.character(sample.by) || length(sample.by)!=1L || is.na(sample.by) || !sample.by %in% names(md)))
    stop("`sample.by` must name one metadata column.",call.=FALSE)
  if (!is.logical(verbose) || length(verbose)!=1L || is.na(verbose)) stop("`verbose` must be TRUE or FALSE.",call.=FALSE)
  if (!is.numeric(base_size) || length(base_size)!=1L || !is.finite(base_size) || base_size<=0)
    stop("`base_size` must be positive.",call.=FALSE)
  if (is.null(grid)) grid <- expand.grid(
    min_umi=unique(pmax(1,params$min_umi+c(-1,0,1))),
    supported_q=unique(pmin(params$candidate_q,pmax(0,params$supported_q+c(-.1,0,.1)))))
  if (!is.data.frame(grid) || !nrow(grid) || !ncol(grid) || anyDuplicated(names(grid)) ||
      any(!names(grid) %in% c("min_umi","supported_q","k")))
    stop("`grid` must be a nonempty data frame with columns from min_umi, supported_q, k.",call.=FALSE)
  for (nm in names(grid)) {
    x <- grid[[nm]]
    if (!is.numeric(x) || any(!is.finite(x))) stop("Grid columns must be finite numeric values.",call.=FALSE)
    if (nm %in% c("min_umi","k") && any(x<1 | x!=floor(x) | x>.Machine$integer.max))
      stop("Grid `",nm,"` must contain positive integers.",call.=FALSE)
    if (nm=="supported_q" && any(x<0 | x>params$candidate_q))
      stop("Grid `supported_q` must lie between 0 and candidate_q.",call.=FALSE)
  }
  grid <- unique(grid)
  rownames(grid) <- NULL
  execute <- function(args,label) tryCatch(do.call(detect_gnrh,c(list(object=object,verbose=FALSE),args)),
    error=function(e) stop(label,": ",conditionMessage(e),call.=FALSE))
  if (verbose) message("Recomputing baseline and ",nrow(grid)," sensitivity trials (sequential).")
  fit <- execute(params,"Baseline failed")
  base_md <- gnrh_metadata(fit,c("gnrh_status","gnrh_class"))
  baseline_pos <- as.character(base_md$gnrh_status)=="pos"
  baseline_class <- as.character(base_md$gnrh_class)
  baseline <- data.frame(cell=rownames(md),class=baseline_class,positive=baseline_pos,
    stored_class=as.character(md$gnrh_class),stored_positive=as.character(md$gnrh_status)=="pos",
    differs_from_stored=is.na(md$gnrh_class) | is.na(md$gnrh_status) |
      baseline_class!=as.character(md$gnrh_class) | baseline_pos!=(as.character(md$gnrh_status)=="pos"),stringsAsFactors=FALSE)
  rm(fit)
  positive_frequency <- class_agreement <- numeric(nrow(md))
  runs <- samples <- vector("list",nrow(grid))
  for (i in seq_len(nrow(grid))) {
    trial <- params
    trial[names(grid)] <- lapply(grid,function(x) x[[i]])
    if (verbose) message("  Trial ",i,"/",nrow(grid))
    fit <- execute(trial,paste("Trial",i,"failed"))
    d <- gnrh_metadata(fit,c("gnrh_status","gnrh_class","gnrh_confident"))
    pos <- as.character(d$gnrh_status)=="pos"
    cls <- as.character(d$gnrh_class)
    positive_frequency <- positive_frequency+pos
    class_agreement <- class_agreement+(cls==baseline_class)
    union <- sum(pos | baseline_pos)
    runs[[i]] <- data.frame(run=i,min_umi=trial$min_umi,supported_q=trial$supported_q,k=trial$k,
      n_positive=sum(pos),n_direct=sum(cls=="direct"),n_supported=sum(cls=="supported"),
      n_confident=sum(d$gnrh_confident %in% TRUE),
      gained=sum(pos & !baseline_pos),lost=sum(!pos & baseline_pos),
      class_changes=sum(cls!=baseline_class),jaccard=if (union) sum(pos & baseline_pos)/union else 1,
      reference_n=fit@misc$gnrh_params$reference_n,
      supported_thr=fit@misc$gnrh_params$supported_thr,stringsAsFactors=FALSE)
    if (!is.null(sample.by)) {
      smp <- as.character(md[[sample.by]])
      samples[[i]] <- do.call(rbind,lapply(unique(smp),function(s) {
        ii <- if (is.na(s)) is.na(smp) else !is.na(smp)&smp==s
        data.frame(run=i,sample=s,n_cells=sum(ii),n_positive=sum(pos[ii]),percent_positive=100*mean(pos[ii]))
      }))
    }
    rm(fit)
  }
  runs <- do.call(rbind,runs)
  cells <- data.frame(cell=rownames(md),baseline_positive=baseline_pos,
    positive_frequency=positive_frequency/nrow(grid),
    class_agreement=class_agreement/nrow(grid),stringsAsFactors=FALSE)
  cells$stability <- ifelse(cells$positive_frequency==0,"never positive",
    ifelse(cells$positive_frequency==1,"always positive","parameter-sensitive"))
  labels <- paste0("R",runs$run,"\nUMI ",runs$min_umi," \u00b7 q ",signif(runs$supported_q,3),"\nk ",runs$k)
  runs$setting <- factor(labels,levels=labels)
  theme <- function(p) p+ggplot2::theme_bw(base_size=base_size)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=0,hjust=.5,size=base_size*.8),plot.title=ggplot2::element_text(face="bold"))
  plots <- list()
  plots$counts <- theme(ggplot2::ggplot(runs,ggplot2::aes(x=.data$setting,y=.data$n_positive))+
    ggplot2::geom_col(fill="#0072B2")+ggplot2::geom_hline(yintercept=sum(baseline_pos),linetype=2)+
    ggplot2::labs(title="Detection sensitivity to parameters",subtitle="Dashed line: recomputed baseline",x=NULL,y="Positive cells"))
  plots$agreement <- theme(ggplot2::ggplot(runs,ggplot2::aes(x=.data$setting,y=.data$jaccard))+
    ggplot2::geom_point(size=2,colour="#009E73")+ggplot2::scale_y_continuous(limits=c(0,1))+
    ggplot2::labs(title="Positive-set overlap with baseline",x=NULL,y="Jaccard overlap",
      subtitle="Both sets empty gives 1; agreement is not accuracy"))
  selected <- cells[cells$positive_frequency>0 | cells$baseline_positive,,drop=FALSE]
  if (nrow(selected)) plots$frequency <- theme(ggplot2::ggplot(selected,ggplot2::aes(x=.data$positive_frequency))+
    ggplot2::geom_histogram(breaks=seq(0,1,length.out=min(nrow(grid)+1L,21L)),fill="#E69F00",colour="white",closed="right")+
    ggplot2::scale_x_continuous(limits=c(0,1))+
    ggplot2::labs(title="Per-cell robustness across the grid",x="Fraction of trials called positive",y="Cells",
      subtitle="Union of baseline and trial positives; not a probability"))
  out <- list(plot=patchwork::wrap_plots(plots,ncol=1),plots=plots,runs=runs,cells=cells,
    samples=if (is.null(sample.by)) NULL else do.call(rbind,samples),baseline=baseline,
    parameters=list(grid=grid,baseline=params,sample.by=sample.by,n_runs=nrow(grid),
      baseline_disagreements=sum(baseline$differs_from_stored,na.rm=TRUE)))
  class(out) <- c("gnrh_stability_result","list")
  out
}

#' Print the sensitivity composite without dumping its data tables
#' @param x A GnRHcell sensitivity result.
#' @param ... Passed to the plot print method.
#' @keywords internal
#' @export
print.gnrh_stability_result <- function(x, ...) {
  print(x$plot, ...)
  invisible(x)
}
