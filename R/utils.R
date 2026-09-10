# ============================================================================= #
# GnRHcell internal utilities
# ============================================================================= #
# ============================================================================= #
# Runtime / logging
# ============================================================================= #
#' Format elapsed runtime into human-readable text
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {
  if (is.null(seconds) || !length(seconds) || !is.finite(seconds)) return(NA_character_)
  seconds <- max(as.numeric(seconds),0)
  hrs <- floor(seconds/3600)
  mins <- floor((seconds%%3600)/60)
  secs <- seconds%%60
  if (hrs>0) sprintf("%dh %dm %.1fs",hrs,mins,secs) else if (mins>0) sprintf("%dm %.1fs",mins,secs) else sprintf("%.1fs",secs)
}


#' Format runtime duration
#' @keywords internal
#' @noRd
.format_runtime <- function(seconds) .format_duration(seconds)


#' Styled GnRHcell logger
#' @keywords internal
#' @noRd
.msg <- function(verbose=TRUE,color=interactive()) {
  t0 <- Sys.time()
  c0 <- if (color) "\033[0m" else ""
  bold <- if (color) "\033[1m" else ""
  cols <- if (color) c(info="\033[90m",step="\033[36m",done="\033[32m",warn="\033[33m",error="\033[31m",header="\033[35m") else setNames(rep("",6),c("info","step","done","warn","error","header"))
  function(...,type=c("info","step","done","warn","error","header"),duration=NULL) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    type <- match.arg(type)
    txt <- paste(...,collapse=" ")
    if (type=="done") {
      if (is.null(duration)) duration <- as.numeric(difftime(Sys.time(),t0,units="secs"))
      txt <- sprintf("%s Duration: %s",txt,.format_duration(duration))
    }
    prefix <- c(info="[INFO]",step="[STEP]",done="[DONE]",warn="[WARN]",error="[ERROR]",header="[GNRH]")[[type]]
    col <- cols[[type]]
    if (type=="header") col <- paste0(col,bold)
    cat(col,prefix,c0," ",col,txt,c0,"\n",sep="")
    invisible(NULL)
  }
}


# ============================================================================= #
# Numeric / gene helpers
# ============================================================================= #
#' Standardize numeric vector
#' @keywords internal
#' @noRd
.scale0 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (!any(ok)) return(rep(0,length(x)))
  s <- stats::sd(x[ok])
  if (!is.finite(s) || s==0) return(rep(0,length(x)))
  out <- rep(NA_real_,length(x))
  out[ok] <- (x[ok]-mean(x[ok]))/s
  out
}


#' Match gene symbols case-insensitively
#' @keywords internal
#' @noRd
.match_genes <- function(features,genes) {
  if (!length(features) || !length(genes)) return(character())
  idx <- match(toupper(as.character(features)),toupper(as.character(genes)))
  unique(genes[idx[!is.na(idx)]])
}


#' Extract expression matrix from a Seurat object
#' @keywords internal
#' @noRd
.get_expr <- function(object,assay=NULL,layer="counts") {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  assay <- assay %||% Seurat::DefaultAssay(object)
  if (!assay %in% names(object@assays)) stop("Assay `",assay,"` not found.",call.=FALSE)
  Seurat::GetAssayData(object,assay=assay,layer=layer)
}


# ============================================================================= #
# GnRH confidence
# ============================================================================= #
#' Compute high-confidence GnRH classification
#' @keywords internal
#' @noRd
.compute_gnrh_confident <- function(metadata,min_umi=2) {
  req <- c("gnrh_class","gnrh_raw","gnrh_identity_strong","gnrh_independent_support")
  miss <- setdiff(req,colnames(metadata))
  if (length(miss)) stop("Missing metadata columns required for confidence classification: ",paste(miss,collapse=", "),call.=FALSE)
  cls <- as.character(metadata$gnrh_class)
  raw <- as.numeric(metadata$gnrh_raw)
  strong <- !is.na(metadata$gnrh_identity_strong) & metadata$gnrh_identity_strong
  support <- !is.na(metadata$gnrh_independent_support) & metadata$gnrh_independent_support
  out <- (cls=="direct" & raw>=min_umi & strong) | (cls=="supported" & strong & support)
  out[is.na(out)] <- FALSE
  out
}


# ============================================================================= #
# Marker-hit summaries
# ============================================================================= #
#' Add total GnRH marker hit summaries
#' @keywords internal
#' @noRd
.make_total_hits <- function(object) {
  md <- object[[]]
  req <- c("gnrh_core_hits","gnrh_mig_hits","gnrh_neuro_hits")
  miss <- setdiff(req,colnames(md))
  if (length(miss)) stop("Missing GnRH marker-hit columns: ",paste(miss,collapse=", "),call.=FALSE)
  identity_neuro <- md$gnrh_core_hits+md$gnrh_neuro_hits
  total <- identity_neuro+md$gnrh_mig_hits
  bin_hits <- function(x) factor(cut(x,breaks=c(-Inf,0,1,2,3,4,Inf),labels=c("0","1","2","3","4","5+"),right=TRUE),levels=c("0","1","2","3","4","5+"))
  object$gnrh_identity_neuro_hits <- identity_neuro
  object$gnrh_total_hits <- total
  object$gnrh_total_hits_bin <- bin_hits(total)
  object$total_hits <- total
  object$total_hits_bin <- object$gnrh_total_hits_bin
  object
}


# ============================================================================= #
# ROC validation
# ============================================================================= #
#' Build GnRH ROC curve data
#' @keywords internal
#' @noRd
.gnrh_build_roc <- function(object,truth,predictor="gnrh_support_score_raw",positive=NULL) {
  if (!requireNamespace("pROC",quietly=TRUE)) stop("Package 'pROC' is required.",call.=FALSE)
  md <- object[[]]
  if (missing(truth) || is.null(truth)) stop("`truth` must contain an independent binary reference annotation.",call.=FALSE)
  y <- if (is.character(truth) && length(truth)==1L) {
    if (!truth %in% colnames(md)) stop("Truth column `",truth,"` not found.",call.=FALSE)
    md[[truth]]
  } else truth
  if (!predictor %in% colnames(md)) stop("Predictor column `",predictor,"` not found.",call.=FALSE)
  x <- suppressWarnings(as.numeric(md[[predictor]]))
  if (length(y)!=length(x)) stop("`truth` must contain one value per cell.",call.=FALSE)
  keep <- !is.na(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  classes <- unique(as.character(y))
  if (length(classes)!=2L) stop("ROC requires exactly two truth classes; found: ",paste(classes,collapse=", "),call.=FALSE)
  if (!is.null(positive)) {
    positive <- as.character(positive)
    if (!positive %in% classes) stop("`positive` is not present in `truth`.",call.=FALSE)
    y <- factor(as.character(y),levels=c(setdiff(classes,positive),positive))
  } else {
    y <- factor(y)
    if (nlevels(y)!=2L) stop("Unable to determine binary truth levels.",call.=FALSE)
  }
  roc <- pROC::roc(response=y,predictor=x,levels=levels(y),direction="<",quiet=TRUE)
  data.frame(threshold=roc$thresholds,fpr=1-roc$specificities,tpr=roc$sensitivities,auc=as.numeric(pROC::auc(roc)),predictor=predictor,stringsAsFactors=FALSE)
}


# ============================================================================= #
# Misc storage
# ============================================================================= #
#' Initialize GnRHcell miscellaneous storage
#' @keywords internal
#' @noRd
.init_gnrh_misc <- function(object) {
  if (is.null(object@misc$gnrh) || !is.list(object@misc$gnrh)) object@misc$gnrh <- list()
  object
}


#' Count values
#' @keywords internal
#' @noRd
.count_factor <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.character(x)
  as.list(table(x[!is.na(x)]))
}


#' Add GnRHcell run information
#' @keywords internal
#' @noRd
.add_run_info <- function(object,step_times=list(),params=list()) {
  object <- .init_gnrh_misc(object)
  md <- object[[]]
  count_if <- function(x) if (x %in% colnames(md)) .count_factor(md[[x]]) else NULL
  sum_if <- function(x) if (x %in% colnames(md)) sum(md[[x]] %in% TRUE,na.rm=TRUE) else NULL
  object@misc$gnrh$run_info <- list(
    date=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),
    dataset=list(n_genes=nrow(object),n_cells=ncol(object)),
    summary=list(
      status=count_if("gnrh_status"),
      class=count_if("gnrh_class"),
      confident=count_if("gnrh_confident"),
      stage=count_if("gnrh_stage"),
      secretory=count_if("gnrh_secretory"),
      reference_positive=sum_if("gnrh_reference_positive"),
      transcriptomic_candidate=sum_if("gnrh_transcriptomic_candidate"),
      direct_isolated=sum_if("gnrh_direct_isolated"),
      direct_signal=sum_if("gnrh_direct_signal")
    ),
    timing=list(seconds=step_times,human=lapply(step_times,.format_duration)),
    params=params,
    session=list(r_version=R.version.string,platform=R.version$platform)
  )
  object
}


# ============================================================================= #
# Run statistics
# ============================================================================= #
#' Load GnRHcell run statistics
#' @keywords internal
#' @noRd
.load_gnrh_stats <- function(files=NULL,dir=".",pattern="_gnrh_run_info\\.tsv$") {
  if (is.null(files)) files <- list.files(dir,pattern=pattern,full.names=TRUE)
  if (!length(files)) stop("No GnRH run info files found.",call.=FALSE)
  miss <- files[!file.exists(files)]
  if (length(miss)) stop("Missing files: ",paste(miss,collapse=", "),call.=FALSE)
  stats <- do.call(rbind,lapply(files,function(f) {
    x <- utils::read.delim(f,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
    x$source_file <- basename(f)
    x
  }))
  num <- c("n_genes","n_cells","neg","pos","direct","supported","confident","reference_positive","transcriptomic_candidate","direct_isolated","direct_signal","detect_sec","stage_sec","diagnostics_sec","total_sec")
  for (nm in intersect(num,colnames(stats))) stats[[nm]] <- suppressWarnings(as.numeric(stats[[nm]]))
  zero <- intersect(c("neg","pos","direct","supported","confident","reference_positive","transcriptomic_candidate","direct_isolated","direct_signal"),colnames(stats))
  for (nm in zero) stats[[nm]][is.na(stats[[nm]])] <- 0
  stats
}


#' Extract GnRHcell run information
#' @param object A Seurat object processed with `run_gnrh()`.
#' @param dataset_name Optional dataset name.
#' @return One-row data frame of GnRHcell run statistics.
#' @export
extract_gnrh_run_info <- function(object,dataset_name=NULL) {
  if (is.null(object@misc$gnrh$run_info)) stop("Missing object@misc$gnrh$run_info. Run run_gnrh() first.",call.=FALSE)
  info <- object@misc$gnrh$run_info
  sec <- info$timing$seconds %||% list()
  hum <- info$timing$human %||% list()
  md <- object[[]]
  count_value <- function(col,val) if (col %in% colnames(md)) sum(as.character(md[[col]])==val,na.rm=TRUE) else 0L
  count_true <- function(col) if (col %in% colnames(md)) sum(md[[col]] %in% TRUE,na.rm=TRUE) else 0L
  data.frame(
    dataset=dataset_name %||% "dataset",
    n_genes=as.integer(info$dataset$n_genes),
    n_cells=as.integer(info$dataset$n_cells),
    neg=count_value("gnrh_status","neg"),
    pos=count_value("gnrh_status","pos"),
    direct=count_value("gnrh_class","direct"),
    supported=count_value("gnrh_class","supported"),
    confident=count_true("gnrh_confident"),
    reference_positive=count_true("gnrh_reference_positive"),
    transcriptomic_candidate=count_true("gnrh_transcriptomic_candidate"),
    direct_signal=count_true("gnrh_direct_signal"),
    direct_isolated=count_true("gnrh_direct_isolated"),
    detect_sec=as.numeric(sec$detect_sec %||% NA_real_),
    stage_sec=as.numeric(sec$stage_sec %||% NA_real_),
    diagnostics_sec=as.numeric(sec$diagnostics_sec %||% NA_real_),
    total_sec=as.numeric(sec$total_sec %||% NA_real_),
    detect_time=hum$detect_sec %||% .format_duration(sec$detect_sec %||% NA_real_),
    stage_time=hum$stage_sec %||% .format_duration(sec$stage_sec %||% NA_real_),
    diagnostics_time=hum$diagnostics_sec %||% .format_duration(sec$diagnostics_sec %||% NA_real_),
    total_time=hum$total_sec %||% .format_duration(sec$total_sec %||% NA_real_),
    stringsAsFactors=FALSE
  )
}


# ============================================================================= #
# Build marker gene sets
# ============================================================================= #
#' Build gene sets from marker tables
#' @param files Named character vector of marker files.
#' @param dir Directory containing files.
#' @param gene_col Gene column.
#' @return Named list of gene sets.
#' @export
build_gene_sets <- function(files,dir=".",gene_col="gene") {
  if (!is.character(files) || !length(files)) stop("`files` must be a non-empty character vector.",call.=FALSE)
  if (is.null(names(files)) || anyNA(names(files)) || any(!nzchar(names(files))) || anyDuplicated(names(files))) stop("`files` must have unique non-empty names.",call.=FALSE)
  out <- lapply(files,function(f) {
    path <- if (file.exists(f)) f else file.path(dir,f)
    if (!file.exists(path)) stop("Marker file not found: ",path,call.=FALSE)
    x <- utils::read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
    if (!gene_col %in% colnames(x)) stop("Column `",gene_col,"` not found in ",basename(path),".",call.=FALSE)
    g <- toupper(trimws(as.character(x[[gene_col]])))
    unique(g[!is.na(g) & nzchar(g)])
  })
  names(out) <- names(files)
  out
}


# ============================================================================= #
# Gene-set overlap
# ============================================================================= #
#' Gene set overlap analysis and visualization
#' @export
gnrh_gene_upset <- function(gene_sets,min_size=1,venn_title="Overlap of Gene Sets",outdir=".",save_plot=TRUE,max_intersections=Inf,plot_width=NULL,plot_height=NULL,dpi=600) {
  for (pkg in c("ggplot2","ComplexUpset","patchwork")) if (!requireNamespace(pkg,quietly=TRUE)) stop("Package '",pkg,"' is required.",call.=FALSE)
  if (!is.list(gene_sets) || length(gene_sets)<2L) stop("`gene_sets` must contain at least two sets.",call.=FALSE)
  labels <- names(gene_sets)
  if (is.null(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) stop("`gene_sets` must have unique non-empty names.",call.=FALSE)
  min_size <- as.integer(min_size)
  if (!is.finite(min_size) || min_size<1L) stop("`min_size` must be >= 1.",call.=FALSE)
  if (!is.null(max_intersections) && (length(max_intersections)!=1L || is.na(max_intersections) || max_intersections<1)) stop("`max_intersections` must be NULL, positive, or Inf.",call.=FALSE)
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  clean <- function(x) sort(unique(toupper(trimws(as.character(x)))[!is.na(x) & nzchar(trimws(as.character(x)))]))
  gene_sets <- lapply(gene_sets,clean)
  if (!any(lengths(gene_sets))) stop("All gene sets are empty.",call.=FALSE)
  summary_df <- data.frame(dataset=labels,n_genes=lengths(gene_sets),stringsAsFactors=FALSE)
  utils::write.csv(summary_df,file.path(outdir,"gene_set_sizes.csv"),row.names=FALSE)
  pairwise_results <- list()
  pairwise_summary <- do.call(rbind,lapply(utils::combn(labels,2,simplify=FALSE),function(z) {
    ov <- intersect(gene_sets[[z[1]]],gene_sets[[z[2]]])
    nm <- paste(z,collapse="_vs_")
    pairwise_results[[nm]] <<- ov
    utils::write.csv(data.frame(gene=ov),file.path(outdir,paste0("genes_",gsub("[^A-Za-z0-9_\\-]+","_",nm),".csv")),row.names=FALSE)
    data.frame(dataset1=z[1],dataset2=z[2],n_overlap=length(ov),stringsAsFactors=FALSE)
  }))
  utils::write.csv(pairwise_summary,file.path(outdir,"pairwise_overlap_summary.csv"),row.names=FALSE)
  common_all <- Reduce(intersect,gene_sets)
  utils::write.csv(data.frame(gene=common_all),file.path(outdir,"genes_common_all.csv"),row.names=FALSE)
  unique_results <- setNames(vector("list",length(labels)),labels)
  unique_summary <- do.call(rbind,lapply(labels,function(lbl) {
    u <- setdiff(gene_sets[[lbl]],Reduce(union,gene_sets[names(gene_sets)!=lbl]))
    unique_results[[lbl]] <<- u
    utils::write.csv(data.frame(gene=u),file.path(outdir,paste0("genes_unique_",gsub("[^A-Za-z0-9_\\-]+","_",lbl),".csv")),row.names=FALSE)
    data.frame(dataset=lbl,n_unique=length(u))
  }))
  utils::write.csv(unique_summary,file.path(outdir,"unique_gene_summary.csv"),row.names=FALSE)
  all_genes <- sort(unique(unlist(gene_sets,use.names=FALSE)))
  membership <- data.frame(gene=all_genes,check.names=FALSE)
  for (lbl in labels) membership[[lbl]] <- all_genes %in% gene_sets[[lbl]]
  utils::write.csv(membership,file.path(outdir,"gene_set_membership.csv"),row.names=FALSE)
  sig <- apply(membership[,labels,drop=FALSE],1,function(x) paste(as.integer(x),collapse=""))
  ints <- sort(table(sig),decreasing=TRUE)
  ints <- ints[ints>=min_size]
  n_available <- length(ints)
  if (!n_available) stop("No intersections satisfy `min_size = ",min_size,"`.",call.=FALSE)
  n_sets <- length(gene_sets)
  if (is.null(max_intersections)) max_intersections <- if (n_sets<=4) 25L else if (n_sets<=6) 35L else if (n_sets<=10) 50L else 60L
  n_displayed <- if (is.infinite(max_intersections)) n_available else min(n_available,as.integer(max_intersections))
  if (is.null(plot_width)) plot_width <- max(9,min(18,7+0.16*n_displayed+0.25*n_sets))
  if (is.null(plot_height)) plot_height <- max(6,min(11,4.6+0.45*n_sets))
  width_ratio <- min(0.32,max(0.19,0.17+max(nchar(labels))/250))
  count_size <- if (n_displayed<=20) 3.6 else if (n_displayed<=40) 3 else 2.5
  text_size <- if (n_sets<=6) 9 else if (n_sets<=10) 8 else 7
  plot <- ComplexUpset::upset(
    membership,intersect=labels,min_size=min_size,n_intersections=n_displayed,width_ratio=width_ratio,
    base_annotations=list("Intersection size"=ComplexUpset::intersection_size(counts=TRUE,bar_number_threshold=0.82,text=list(size=count_size,fontface="bold"))),
    set_sizes=ComplexUpset::upset_set_size(geom=ggplot2::geom_bar(width=0.68,fill="#4D4D4D")),
    sort_sets="descending",sort_intersections_by="cardinality"
  ) +
    patchwork::plot_annotation(title=venn_title,theme=ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5,face="bold",size=max(12,min(17,19-0.4*n_sets)),margin=ggplot2::margin(b=8)))) &
    ggplot2::theme(axis.text=ggplot2::element_text(size=text_size))
  if (isTRUE(save_plot)) for (ext in c("pdf","png")) ggplot2::ggsave(file.path(outdir,paste0("gene_overlap_plot.",ext)),plot,width=plot_width,height=plot_height,dpi=dpi,bg="white")
  list(gene_sets=gene_sets,summary=summary_df,pairwise=pairwise_results,pairwise_summary=pairwise_summary,common_all=common_all,unique=unique_results,unique_summary=unique_summary,membership=membership,intersections_available=n_available,intersections_displayed=n_displayed,plot_width=plot_width,plot_height=plot_height,plot=plot,upset_plot=plot,plot_type="upset")
}


# ============================================================================= #
# GnRH module weighting
# ============================================================================= #
#' GnRH module weights
#' @keywords internal
#' @noRd
.gnrh_module_weights <- function() c(
  identity_primary=3,
  identity_supportive=1.5,
  migration_primary=1,
  migration_supportive=0.5,
  neuroendocrine_primary=1.5,
  neuroendocrine_supportive=0.75,
  hormone_supportive=0.25,
  guidance_environment=0
)


#' Compute a weighted GnRH module score
#' @keywords internal
#' @noRd
.weight_gnrh_scores <- function(scores,weights=.gnrh_module_weights(),normalize=FALSE) {
  if ((!is.list(scores) && !is.data.frame(scores)) || is.null(names(scores))) stop("`scores` must be a named list or data frame.",call.=FALSE)
  if (!is.numeric(weights) || is.null(names(weights))) stop("`weights` must be a named numeric vector.",call.=FALSE)
  common <- intersect(names(weights),names(scores))
  if (!length(common)) stop("No GnRH module scores matched the supplied weights.",call.=FALSE)
  w <- weights[common]
  keep <- is.finite(w) & w!=0
  common <- common[keep]
  w <- w[keep]
  if (!length(common)) stop("All matched GnRH module weights are zero or non-finite.",call.=FALSE)
  mat <- do.call(cbind,lapply(scores[common],as.numeric))
  if (is.null(dim(mat))) mat <- matrix(mat,ncol=1L)
  out <- rowSums(sweep(mat,2,w,"*"),na.rm=TRUE)
  if (isTRUE(normalize)) out <- out/sum(abs(w))
  out
}
