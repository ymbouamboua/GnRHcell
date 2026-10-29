#' Diagnostic evidence categories
#' @keywords internal
#' @noRd
.gnrh_qc_categories <- function(md) {
  flag <- function(nm) if (nm %in% names(md)) md[[nm]] %in% TRUE else rep(FALSE,nrow(md))
  cls <- as.character(md$gnrh_class)
  out <- rep("negative",nrow(md))
  out[flag("gnrh_transcriptomic_candidate") | flag("gnrh_dropout_candidate")] <- "transcriptomic candidate"
  out[flag("gnrh_direct_isolated")] <- "isolated signal"
  out[flag("gnrh_review_candidate")] <- "manual review"
  out[!is.na(cls) & cls=="supported"] <- "supported"
  out[!is.na(cls) & cls=="direct"] <- "direct"
  factor(out,levels=c("negative","isolated signal","manual review","transcriptomic candidate","supported","direct"))
}

#' Diagnostic colors
#' @keywords internal
#' @noRd
.gnrh_qc_colors <- function() c(
  negative="#B8B8B8", `isolated signal`="#CC79A7", `manual review`="#E69F00",
  `transcriptomic candidate`="#56B4E9", supported="#009E73", direct="#0072B2"
)

#' Validate a positive integer scalar
#' @keywords internal
#' @noRd
.gnrh_qc_integer <- function(x, name, minimum=1) {
  if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<minimum ||
      x!=floor(x) || x>.Machine$integer.max)
    stop("`",name,"` must be a finite integer >= ",minimum,".",call.=FALSE)
  invisible(x)
}

#' Select background points without changing the global random state
#' @keywords internal
#' @noRd
.gnrh_qc_points <- function(data,max_background,seed) {
  .gnrh_qc_integer(max_background,"max_background")
  .gnrh_qc_integer(seed,"seed",0)
  bg <- which(data$category=="negative")
  if (length(bg)>max_background) {
    had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if (had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else
      if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv))
    set.seed(seed)
    bg <- bg[sample.int(length(bg),max_background)]
  }
  rbind(data[bg,,drop=FALSE],data[data$category!="negative",,drop=FALSE])
}

#' Summarize detection evidence and technical QC
#'
#' @param object A Seurat object processed by [run_gnrh()] or [detect_gnrh()].
#' @param sample.by Metadata column identifying donors or samples. `NULL`
#'   selects the first available of `donor`, `sample`, `orig.ident`. The selected
#'   name is returned in `parameters`. With none available, all cells form one
#'   group. Missing sample IDs form a separate explicitly labeled group.
#' @param qc.columns Numeric technical-QC columns to display. `NULL` selects
#'   available RNA count, feature, mitochondrial and doublet-score fields.
#'   `character()` skips this panel. Values are shown without normalization.
#' @param reduction Optional reduction for context and candidate-zoom panels.
#'   Default `"umap"`; an unavailable reduction is reported in `notes` and skipped.
#'   Use `NULL` to omit embedding panels.
#' @param dims Two reduction dimensions to display.
#' @param max_background Maximum negative cells drawn in scatter plots. All
#'   candidate cells are retained; tables and summaries always use all cells.
#' @param seed Seed for background sampling; the global RNG state is preserved.
#' @param base_size Plot font size.
#' @param ncol Columns in the composite figure.
#' @param zoom.by Cells defining the zoom: `"positive"` (default),
#'   `"confident"`, or `"candidates"` (all non-negative evidence categories).
#'   All categories remain drawn inside the selected view.
#' @param zoom.cells Optional explicit cell names defining a local zoom,
#'   overriding `zoom.by`. Names must identify cells in the current object.
#' @param zoom.padding Non-negative fraction added around the zoom extent.
#'   Default `0.15`. Both embedding panels preserve distances with square
#'   coordinate limits. The zoom rectangle is drawn on the context panel.
#' @param zoom.quantile Central interval fraction used on each axis for the
#'   zoom frame. Default `0.95` limits expansion by dispersed outliers; `1`
#'   includes the full selected extent. Padding is then applied and limits are
#'   made square. Actual selected-cell counts inside/outside are shown in the
#'   subtitle and returned in `data$zoom`; no cells are filtered from the object
#'   or the global view. This is not a biological selection criterion.
#' @return A list with `plot` (patchwork), named `plots`, and `data` containing
#'   cell evidence, long technical-QC values, sample counts and proportions,
#'   stage confidence, and plotted cells. Also returns `parameters` and `notes`.
#'   The input object and annotations are never modified.
#' @details These are descriptive diagnostics, not validation against independent
#'   truth. The evidence panel uses the raw support score (excluding direct
#'   GNRH1 expression). Reference lines show individual gates, not the complete
#'   multi-rule decision boundary. Sample rates use all analyzed cells in that
#'   sample as the denominator. They are not estimates of tissue abundance.
#'   All archived metadata is aligned by cell name via [gnrh_metadata()].
#' @examples
#' \dontrun{
#' qc <- gnrh_qc(obj, sample.by="donor")
#' qc$plot
#' qc$data$samples
#' ggplot2::ggsave("gnrh_qc.png", qc$plot, width=16, height=16, dpi=300)
#' }
#' @export
gnrh_qc <- function(object,sample.by=NULL,qc.columns=NULL,reduction="umap",dims=1:2,
                    max_background=30000L,seed=1234L,base_size=11,ncol=2,
                    zoom.by=c("positive","confident","candidates"),zoom.cells=NULL,zoom.padding=.15,zoom.quantile=.95) {
  zoom.by <- match.arg(zoom.by)
  if (!is.numeric(zoom.padding) || length(zoom.padding)!=1L || !is.finite(zoom.padding) || zoom.padding<0)
    stop("`zoom.padding` must be a finite non-negative number.",call.=FALSE)
  if (!is.numeric(zoom.quantile) || length(zoom.quantile)!=1L || !is.finite(zoom.quantile) || zoom.quantile<=0 || zoom.quantile>1)
    stop("`zoom.quantile` must be in (0, 1].",call.=FALSE)
  md <- gnrh_metadata(object)
  if (!is.null(zoom.cells) && (!is.character(zoom.cells) || !length(zoom.cells) || anyNA(zoom.cells) || any(!zoom.cells %in% rownames(md))))
    stop("`zoom.cells` must contain current cell names.",call.=FALSE)
  if (!is.null(zoom.cells)) zoom.cells <- unique(zoom.cells)
  req <- c("gnrh_raw","gnrh_support_score_raw","gnrh_class","gnrh_status")
  missing <- setdiff(req,names(md))
  if (length(missing)) stop("Missing GnRH results: ",paste(missing,collapse=", "),". Run detection first.",call.=FALSE)
  .gnrh_qc_integer(max_background,"max_background")
  .gnrh_qc_integer(seed,"seed",0)
  .gnrh_qc_integer(ncol,"ncol")
  if (!is.numeric(base_size) || length(base_size)!=1L || !is.finite(base_size) || base_size<=0)
    stop("`base_size` must be positive.",call.=FALSE)
  notes <- character()
  if (anyNA(md$gnrh_status) || anyNA(md$gnrh_class) ||
      any(!as.character(md$gnrh_status) %in% c("neg","pos")) ||
      any(!as.character(md$gnrh_class) %in% c("neg","supported","direct")) ||
      any((as.character(md$gnrh_status)=="pos") != (as.character(md$gnrh_class)!="neg")))
    stop("GnRH status/class contains missing or inconsistent values. Rerun detection on the current cells.",call.=FALSE)
  if (is.null(sample.by)) {
    hits <- intersect(c("donor","sample","orig.ident"),names(md))
    sample.by <- if (length(hits)) hits[1L] else NULL
  }
  if (!is.null(sample.by) && (!is.character(sample.by) || length(sample.by)!=1L || is.na(sample.by) || !sample.by %in% names(md)))
    stop("`sample.by` must name one metadata column.",call.=FALSE)
  sample <- if (is.null(sample.by)) rep("all cells",nrow(md)) else as.character(md[[sample.by]])
  missing_label <- "<missing sample>"
  while (missing_label %in% sample[!is.na(sample)]) missing_label <- paste0(missing_label,"*")
  sample[is.na(sample)] <- missing_label
  cell <- data.frame(cell=rownames(md),sample=sample,category=.gnrh_qc_categories(md),
    class=as.character(md$gnrh_class),positive=as.character(md$gnrh_status)=="pos",
    gnrh_raw=as.numeric(md$gnrh_raw),support_score=as.numeric(md$gnrh_support_score_raw),
    stringsAsFactors=FALSE)
  plotted <- .gnrh_qc_points(cell,max_background,seed)
  if (!is.null(zoom.cells)) {
    extra <- cell[cell$cell %in% zoom.cells & !cell$cell %in% plotted$cell,,drop=FALSE]
    plotted <- rbind(extra,plotted)
  }
  colors <- .gnrh_qc_colors()
  style <- function(p) {
    wrap <- function(x,width) if (is.null(x)) NULL else paste(strwrap(x,width=width),collapse="\n")
    p+ggplot2::labs(title=wrap(p$labels$title,42),subtitle=wrap(p$labels$subtitle,52))+
    ggplot2::theme_bw(base_size=base_size)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),legend.position="bottom",
      axis.text.x=ggplot2::element_text(angle=30,hjust=1))
  }
  color_scale <- function() ggplot2::scale_colour_manual(values=colors,drop=FALSE,guide=ggplot2::guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=2,alpha=1)))
  plots <- list()
  plots$evidence <- style(ggplot2::ggplot(plotted,
    ggplot2::aes(x=log1p(.data$gnrh_raw),y=.data$support_score,colour=.data$category))+
    ggplot2::geom_point(size=.7,alpha=.5,na.rm=TRUE)+color_scale()+
    ggplot2::labs(title="Detection evidence",x="log1p(GNRH1 counts)",
      y="Independent support score (raw)",colour=NULL,
      subtitle="Reference lines are individual gates; classification uses additional rules"))
  pars <- object@misc$gnrh_params
  if (length(pars$min_umi)==1L && is.finite(pars$min_umi))
    plots$evidence <- plots$evidence+ggplot2::geom_vline(xintercept=log1p(pars$min_umi),linetype=2,colour="#333333")
  if (length(pars$supported_thr)==1L && is.finite(pars$supported_thr))
    plots$evidence <- plots$evidence+ggplot2::geom_hline(yintercept=pars$supported_thr,linetype=2,colour="#333333")
  if (is.null(qc.columns)) {
    mito <- intersect(c("percent_mito","percent.mt","percent_mitochondrial","percent_mt"),names(md))
    dbl <- intersect(c("scDblFinder_score","scDblFinder.score","doublet_score","demux_scds_score"),names(md))
    qc.columns <- c(intersect(c("nCount_RNA","nFeature_RNA"),names(md)),utils::head(mito,1),utils::head(dbl,1))
  }
  if (!is.character(qc.columns) || anyNA(qc.columns) || anyDuplicated(qc.columns) || any(!qc.columns %in% names(md)))
    stop("`qc.columns` must contain unique existing numeric metadata fields.",call.=FALSE)
  if (any(!vapply(md[qc.columns],is.numeric,logical(1))))
    stop("All `qc.columns` must be numeric.",call.=FALSE)
  technical <- data.frame(cell=character(),sample=character(),category=factor(levels=levels(cell$category)),metric=character(),value=numeric())
  if (length(qc.columns)) {
    technical <- do.call(rbind,lapply(qc.columns,function(nm)
      data.frame(cell=cell$cell,sample=sample,category=cell$category,metric=nm,value=md[[nm]])))
    technical <- technical[is.finite(technical$value),,drop=FALSE]
    if (nrow(technical)) plots$technical <- style(ggplot2::ggplot(technical,
      ggplot2::aes(x=.data$category,y=.data$value,fill=.data$category))+
      ggplot2::geom_boxplot(outlier.shape=NA,width=.65)+
      ggplot2::facet_wrap(~metric,scales="free_y",ncol=2)+
      ggplot2::scale_fill_manual(values=colors,guide="none")+
      ggplot2::scale_x_discrete(labels=function(x) ifelse(x=="transcriptomic candidate","transcript.\ncandidate",gsub(" ","\n",x,fixed=TRUE)))+
      ggplot2::theme(axis.text.x=ggplot2::element_text(size=base_size*.7,angle=0,hjust=.5))+
      ggplot2::labs(title="Technical QC by evidence category",x=NULL,y="Observed value",
        subtitle="Boxes use all finite values; metrics are not normalized"))
    else notes <- c(notes,"No finite technical-QC values; panel omitted.")
  } else notes <- c(notes,"No technical-QC columns selected or available.")
  if (!is.null(plots$technical)) plots$technical <- plots$technical+ggplot2::theme(axis.text.x=ggplot2::element_text(size=base_size*.7,angle=0,hjust=.5))
  samples <- do.call(rbind,lapply(unique(sample),function(nm) {
    d <- cell[cell$sample==nm,,drop=FALSE]
    data.frame(sample=nm,n_cells=nrow(d),n_positive=sum(d$positive),
      n_direct=sum(d$class=="direct"),n_supported=sum(d$class=="supported"),
      n_review=sum(d$category=="manual review"),n_isolated=sum(d$category=="isolated signal"),
      n_transcriptomic_candidate=sum(d$category=="transcriptomic candidate"),
      percent_positive=100*mean(d$positive),stringsAsFactors=FALSE)
  }))
  rates <- rbind(data.frame(sample=samples$sample,category="direct",percent=100*samples$n_direct/samples$n_cells),
                 data.frame(sample=samples$sample,category="supported",percent=100*samples$n_supported/samples$n_cells))
  rates$sample <- factor(rates$sample,levels=samples$sample)
  plots$samples <- style(ggplot2::ggplot(rates,ggplot2::aes(x=.data$sample,y=.data$percent,fill=.data$category))+
    ggplot2::geom_col(width=.7)+ggplot2::scale_fill_manual(values=colors,breaks=c("direct","supported"))+
    ggplot2::geom_text(data=samples,ggplot2::aes(x=.data$sample,y=.data$percent_positive,
      label=paste0(.data$n_positive,"/",.data$n_cells)),inherit.aes=FALSE,vjust=-.5,size=3)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.18)))+
    ggplot2::labs(title="Detected cells by sample",x=sample.by %||% "Sample",y="Positive cells (%)",fill=NULL,
      subtitle="Labels: positive / all analyzed cells; descriptive proportions"))
  stage <- data.frame()
  zoom_info <- NULL
  if (all(c("gnrh_stage","gnrh_stage_margin","gnrh_stage_resolution") %in% names(md))) {
    stage <- data.frame(cell=cell$cell,stage=md$gnrh_stage,margin=md$gnrh_stage_margin,
      resolution=md$gnrh_stage_resolution,positive=cell$positive)
    stage <- stage[stage$positive & is.finite(stage$margin),,drop=FALSE]
    if (nrow(stage)) plots$stage <- style(ggplot2::ggplot(stage,
      ggplot2::aes(x=.data$stage,y=.data$margin,fill=.data$resolution))+
      ggplot2::geom_boxplot(outlier.shape=NA)+
      ggplot2::scale_fill_manual(values=c(resolved="#0072B2",transitional="#E69F00",`non-gnrh`="#B8B8B8"))+
      ggplot2::labs(title="Developmental-state confidence",x=NULL,y="Top minus second program score",fill="Resolution",
        subtitle="GnRH-positive cells only; expression programs do not establish a trajectory"))
  }
  if (!is.null(reduction)) {
    if (!is.character(reduction) || length(reduction)!=1L || is.na(reduction)) stop("Invalid `reduction`.",call.=FALSE)
    if (!reduction %in% SeuratObject::Reductions(object)) notes <- c(notes,paste0("Reduction '",reduction,"' unavailable; embedding panels omitted."))
    else {
      emb <- SeuratObject::Embeddings(object,reduction=reduction)
      if (!is.numeric(dims) || length(dims)!=2L || any(!is.finite(dims)) ||
          any(dims!=floor(dims)) || any(dims<1 | dims>ncol(emb)) || anyDuplicated(dims))
        stop("`dims` must select two distinct reduction dimensions.",call.=FALSE)
      ii <- match(plotted$cell,rownames(emb))
      embedding <- plotted
      embedding$x <- emb[ii,dims[1]];embedding$y <- emb[ii,dims[2]]
      embedding <- embedding[is.finite(embedding$x)&is.finite(embedding$y),,drop=FALSE]
      n_missing <- sum(!rownames(md) %in% rownames(emb))
      if (n_missing) notes <- c(notes,paste(n_missing,"cells lack embedding coordinates."))
      if (nrow(embedding)) {
        embedding_base <- style(ggplot2::ggplot(embedding,ggplot2::aes(x=.data$x,y=.data$y,colour=.data$category))+
          ggplot2::geom_point(size=.6,alpha=.6)+color_scale()+
          ggplot2::labs(title="Embedding context",x=paste0(reduction,"_",dims[1]),y=paste0(reduction,"_",dims[2]),colour=NULL))
        square_limits <- function(x,y,padding,coverage=1) {
          probs <- c((1-coverage)/2,1-(1-coverage)/2)
          xr <- as.numeric(stats::quantile(x,probs=probs,names=FALSE))
          yr <- as.numeric(stats::quantile(y,probs=probs,names=FALSE))
          span <- max(diff(xr),diff(yr),.1)*(1+2*padding)
          list(x=mean(xr)+c(-.5,.5)*span,y=mean(yr)+c(-.5,.5)*span)
        }
        ci <- match(cell$cell,rownames(emb))
        cx <- emb[ci,dims[1]];cy <- emb[ci,dims[2]]
        cv <- is.finite(cx)&is.finite(cy)
        context_limits <- square_limits(cx[cv],cy[cv],.05)
        plots$embedding <- embedding_base+ggplot2::coord_equal(xlim=context_limits$x,ylim=context_limits$y,expand=FALSE)
        focus_cells <- if (!is.null(zoom.cells)) zoom.cells else switch(zoom.by,
          positive=cell$cell[cell$positive],
          confident={
            if (!"gnrh_confident" %in% names(md)) stop("`zoom.by = 'confident'` requires gnrh_confident.",call.=FALSE)
            cell$cell[md$gnrh_confident %in% TRUE]
          },
          candidates=cell$cell[cell$category!="negative"])
        # Zoom selection uses every requested cell, independently of background sampling.
        fi <- match(focus_cells,rownames(emb))
        fx <- emb[fi,dims[1]];fy <- emb[fi,dims[2]]
        valid_focus <- is.finite(fx)&is.finite(fy)
        if (any(valid_focus)) {
          lim <- square_limits(fx[valid_focus],fy[valid_focus],zoom.padding,zoom.quantile)
          inside <- valid_focus & fx>=lim$x[1] & fx<=lim$x[2] & fy>=lim$y[1] & fy<=lim$y[2]
          inside[is.na(inside)] <- FALSE
          zoom_label <- if (!is.null(zoom.cells)) "Selected-cell zoom" else switch(zoom.by,
            positive="GnRH-positive-cell zoom",confident="High-confidence-cell zoom",candidates="Candidate-region zoom")
          zoom_info <- list(selection=if (is.null(zoom.cells)) zoom.by else "explicit",
            cells=focus_cells[valid_focus],cells_inside=focus_cells[inside],
            cells_outside=focus_cells[valid_focus & !inside],
            n_missing_coordinates=sum(!valid_focus),xlim=lim$x,ylim=lim$y,
            context_xlim=context_limits$x,context_ylim=context_limits$y)
          plots$zoom <- embedding_base+ggplot2::coord_equal(xlim=lim$x,ylim=lim$y,expand=FALSE)+
            ggplot2::labs(title=zoom_label,subtitle=paste0(sum(inside),"/",sum(valid_focus)," selected cells in frame; ",sum(valid_focus & !inside)," outside\nOutside cells remain in the global view"))
          plots$embedding <- plots$embedding+ggplot2::annotate("rect",
            xmin=lim$x[1],xmax=lim$x[2],ymin=lim$y[1],ymax=lim$y[2],
            fill=NA,colour="#333333",linetype=2,linewidth=.4)
          if (diff(lim$x)>=.9*diff(context_limits$x))
            notes <- c(notes,"Selected cells span most of the embedding; use zoom.cells for a specific local population.")
        } else notes <- c(notes,"No selected zoom cells have valid coordinates; zoom omitted.")
      }
    }
  }
  composite_plots <- plots
  composite_plots$samples <- composite_plots$samples+ggplot2::guides(fill="none")
  out <- list(plot=patchwork::wrap_plots(composite_plots,ncol=ncol,widths=rep(1,ncol),
      heights=if (length(plots)==6L && ncol==2L) c(1,1,1.35) else NULL,guides="collect") & ggplot2::theme(legend.position="bottom"),plots=plots,
    data=list(cells=cell,technical=technical,samples=samples,stage=stage,plotted=plotted,zoom=zoom_info),
    parameters=list(sample.by=sample.by,qc.columns=qc.columns,reduction=reduction,dims=dims,
      max_background=max_background,seed=seed,zoom.by=zoom.by,zoom.cells=zoom.cells,zoom.padding=zoom.padding,zoom.quantile=zoom.quantile),notes=notes)
  class(out) <- c("gnrh_qc_result","list")
  out
}

#' Print the QC composite without dumping its data tables
#' @param x A GnRHcell QC result.
#' @param ... Passed to the plot print method.
#' @keywords internal
#' @export
print.gnrh_qc_result <- function(x, ...) {
  print(x$plot, ...)
  invisible(x)
}
