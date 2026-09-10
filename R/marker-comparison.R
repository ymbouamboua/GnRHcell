# ============================================================================= #
# GnRHcell cross-dataset marker visualization
# ============================================================================= #
# ============================================================================= #
# Adaptive helpers
# ============================================================================= #
.gnrh_text_size <- function(
    n,
    base=11,
    min_size=8,
    max_size=12
) {
  n <- max(1L,as.integer(n))
  size <- base*sqrt(40/max(40,n))
  max(min_size,min(max_size,size))
}
.gnrh_plot_dims <- function(
    n_datasets=1L,
    n_genes=1L,
    type=c("heatmap","upset","program")
) {
  type <- match.arg(type)
  n_datasets <- max(1L,as.integer(n_datasets))
  n_genes <- max(1L,as.integer(n_genes))
  if (type=="heatmap") {
    width <- max(7,min(16,4.5+0.85*n_datasets))
    height <- max(5,min(22,3.5+0.14*n_genes))
  } else if (type=="upset") {
    width <- max(9,min(22,7+0.15*n_genes+0.35*n_datasets))
    height <- max(6,min(12,4.8+0.45*n_datasets))
  } else {
    width <- max(8,min(20,6+1.5*n_datasets))
    height <- max(6,min(20,4+0.20*n_genes))
  }
  c(width=width,height=height)
}
# ============================================================================= #
# Merge marker scores
# ============================================================================= #
#' Merge GnRH marker scores across datasets
#'
#' @param ... Marker data frames, marker lists, or marker table paths.
#' @param dataset_names Optional dataset labels.
#' @param dir Directory prepended to relative file paths.
#' @param gene_col Gene-symbol column.
#' @param score_col Optional score column used in every dataset.
#' @param score_candidates Preferred score columns when \code{score_col = NULL}.
#' @param gene_case Standardize genes to uppercase or preserve original case.
#' @param fill Value used for genes absent from a dataset.
#' @param require_same_score Require the same score column across datasets.
#'
#' @return Numeric gene-by-dataset score matrix.
#' @export
merge_gnrh_marker_scores <- function(
    ...,
    dataset_names=NULL,
    dir=".",
    gene_col="gene",
    score_col=NULL,
    score_candidates=c(
      "final_score",
      "marker_score",
      "specificity_score",
      "avg_log2FC",
      "avg_logFC",
      "coexpr_cor",
      "coexpr"
    ),
    gene_case=c("upper","asis"),
    fill=NA_real_,
    require_same_score=TRUE
) {
  gene_case <- match.arg(gene_case)
  inputs <- list(...)
  if (length(inputs)==1L && is.character(inputs[[1L]])) inputs <- as.list(inputs[[1L]])
  if (length(inputs)==1L && is.list(inputs[[1L]]) && !is.data.frame(inputs[[1L]])) inputs <- inputs[[1L]]
  if (!length(inputs)) stop("Supply at least one marker table.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Dataset names
  # ------------------------------------------------------------------------- #
  input_names <- names(inputs)
  if (is.null(dataset_names)) dataset_names <- input_names
  if (
    is.null(dataset_names) ||
    length(dataset_names)!=length(inputs) ||
    anyNA(dataset_names) ||
    any(!nzchar(dataset_names))
  ) {
    dataset_names <- paste0("dataset_",seq_along(inputs))
  }
  if (anyDuplicated(dataset_names)) stop("`dataset_names` must be unique.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Read files
  # ------------------------------------------------------------------------- #
  is_path <- vapply(
    inputs,
    function(x) is.character(x) && length(x)==1L,
    logical(1)
  )
  if (all(is_path)) {
    paths <- vapply(
      inputs,
      function(x) if (file.exists(x)) x else file.path(dir,x),
      character(1)
    )
    missing <- paths[!file.exists(paths)]
    if (length(missing)) stop("Missing marker files: ",paste(missing,collapse=", "),call.=FALSE)
    inputs <- lapply(
      paths,
      utils::read.delim,
      check.names=FALSE,
      stringsAsFactors=FALSE
    )
  }
  # ------------------------------------------------------------------------- #
  # Select score columns
  # ------------------------------------------------------------------------- #
  selected_columns <- vapply(seq_along(inputs),function(i) {
    tab <- inputs[[i]]
    if (!is.data.frame(tab)) stop("Every marker input must be a data frame or TSV file.",call.=FALSE)
    if (!gene_col %in% colnames(tab)) {
      stop("Column `",gene_col,"` is missing from ",dataset_names[[i]],".",call.=FALSE)
    }
    if (!is.null(score_col)) {
      if (!score_col %in% colnames(tab)) {
        stop("Score column `",score_col,"` is missing from ",dataset_names[[i]],".",call.=FALSE)
      }
      return(score_col)
    }
    available <- score_candidates[score_candidates %in% colnames(tab)]
    if (!length(available)) {
      stop(
        "No supported score column found in ",dataset_names[[i]],
        ". Available columns: ",paste(colnames(tab),collapse=", "),
        call.=FALSE
      )
    }
    available[[1L]]
  },character(1))
  if (isTRUE(require_same_score) && length(unique(selected_columns))>1L) {
    stop(
      "Different score columns would be used across datasets: ",
      paste(paste0(dataset_names,"=",selected_columns),collapse=", "),
      ". Supply `score_col` explicitly or use `require_same_score = FALSE`.",
      call.=FALSE
    )
  }
  # ------------------------------------------------------------------------- #
  # Build score vectors
  # ------------------------------------------------------------------------- #
  scores <- lapply(seq_along(inputs),function(i) {
    tab <- inputs[[i]]
    selected <- selected_columns[[i]]
    gene <- trimws(as.character(tab[[gene_col]]))
    if (gene_case=="upper") gene <- toupper(gene)
    value <- suppressWarnings(as.numeric(tab[[selected]]))
    keep <- !is.na(gene) & nzchar(gene) & is.finite(value)
    gene <- gene[keep]
    value <- value[keep]
    if (!length(gene)) return(stats::setNames(numeric(0),character(0)))
    split_values <- split(value,gene)
    vapply(split_values,max,numeric(1),na.rm=TRUE)
  })
  # ------------------------------------------------------------------------- #
  # Build matrix
  # ------------------------------------------------------------------------- #
  genes <- sort(unique(unlist(lapply(scores,names),use.names=FALSE)))
  if (!length(genes)) stop("No usable marker scores were found.",call.=FALSE)
  mat <- vapply(scores,function(x) {
    out <- rep(fill,length(genes))
    names(out) <- genes
    hit <- intersect(genes,names(x))
    out[hit] <- x[hit]
    out
  },numeric(length(genes)))
  if (is.null(dim(mat))) mat <- matrix(mat,ncol=1L)
  rownames(mat) <- genes
  colnames(mat) <- dataset_names
  attr(mat,"score_columns") <- stats::setNames(selected_columns,dataset_names)
  attr(mat,"score_type") <- if (length(unique(selected_columns))==1L) unique(selected_columns) else "mixed"
  mat
}
# ============================================================================= #
# Heatmap helper
# ============================================================================= #
#' Build GnRH marker heatmap
#'
#' @keywords internal
#' @noRd
.gnrh_marker_heatmap <- function(
    matrix,
    row_split=NULL,
    scale_rows=TRUE,
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    show_column_names=TRUE,
    fontsize_row=NULL,
    fontsize_col=NULL,
    fontsize_legend=10,
    title=NULL
) {
  if (!requireNamespace("ComplexHeatmap",quietly=TRUE) || !requireNamespace("circlize",quietly=TRUE)) {
    stop("Packages `ComplexHeatmap` and `circlize` are required.",call.=FALSE)
  }
  if (!is.matrix(matrix) || !nrow(matrix) || !ncol(matrix)) {
    stop("`matrix` must be a non-empty numeric matrix.",call.=FALSE)
  }
  # ------------------------------------------------------------------------- #
  # Adaptive text
  # ------------------------------------------------------------------------- #
  if (is.null(fontsize_row)) {
    fontsize_row <- .gnrh_text_size(
      n=nrow(matrix),
      base=11,
      min_size=8,
      max_size=12
    )
  }
  if (is.null(fontsize_col)) {
    fontsize_col <- .gnrh_text_size(
      n=ncol(matrix),
      base=11,
      min_size=9,
      max_size=12
    )
  }
  plot_matrix <- matrix
  # ------------------------------------------------------------------------- #
  # Row scaling
  # ------------------------------------------------------------------------- #
  if (isTRUE(scale_rows)) {
    scaled <- matrix(
      NA_real_,
      nrow=nrow(matrix),
      ncol=ncol(matrix),
      dimnames=dimnames(matrix)
    )
    for (i in seq_len(nrow(matrix))) {
      values <- matrix[i,]
      observed <- is.finite(values)
      n_observed <- sum(observed)
      if (n_observed>=2L) {
        s <- stats::sd(values[observed])
        if (is.finite(s) && s>0) {
          scaled[i,observed] <- (values[observed]-mean(values[observed]))/s
        } else {
          scaled[i,observed] <- 0
        }
      } else if (n_observed==1L) {
        scaled[i,observed] <- 0
      }
    }
    plot_matrix <- scaled
    colour <- circlize::colorRamp2(
      c(-2,0,2),
      c("#2166AC","white","#B2182B")
    )
    legend_title <- "Row z-score"
  } else {
    observed <- plot_matrix[is.finite(plot_matrix)]
    if (!length(observed)) stop("Heatmap matrix contains no finite values.",call.=FALSE)
    rng <- range(observed)
    if (diff(rng)==0) rng <- c(rng[[1L]]-1,rng[[1L]]+1)
    midpoint <- stats::median(observed)
    colour <- circlize::colorRamp2(
      c(rng[[1L]],midpoint,rng[[2L]]),
      c("#F7FBFF","#6BAED6","#08306B")
    )
    score_type <- attr(matrix,"score_type")
    legend_title <- if (!is.null(score_type) && score_type!="mixed") score_type else "Marker score"
  }
  # ------------------------------------------------------------------------- #
  # Clustering matrix
  # ------------------------------------------------------------------------- #
  clustering_matrix <- plot_matrix
  for (i in seq_len(nrow(clustering_matrix))) {
    missing <- !is.finite(clustering_matrix[i,])
    if (any(missing)) {
      observed <- clustering_matrix[i,!missing]
      replacement <- if (length(observed) && any(is.finite(observed))) mean(observed,na.rm=TRUE) else 0
      clustering_matrix[i,missing] <- replacement
    }
  }
  # ------------------------------------------------------------------------- #
  # Dendrograms
  # ------------------------------------------------------------------------- #
  row_clustering <- FALSE
  if (isTRUE(cluster_rows) && nrow(clustering_matrix)>1L) {
    row_clustering <- stats::hclust(stats::dist(clustering_matrix))
  }
  column_clustering <- FALSE
  if (isTRUE(cluster_columns) && ncol(clustering_matrix)>1L) {
    column_clustering <- stats::hclust(stats::dist(t(clustering_matrix)))
  }
  # ------------------------------------------------------------------------- #
  # Heatmap
  # ------------------------------------------------------------------------- #
  heatmap <- ComplexHeatmap::Heatmap(
    plot_matrix,
    name=legend_title,
    col=colour,
    row_split=row_split,
    row_names_gp=grid::gpar(
      fontsize=fontsize_row,
      fontface="italic"
    ),
    column_names_gp=grid::gpar(
      fontsize=fontsize_col
    ),
    cluster_rows=row_clustering,
    cluster_columns=column_clustering,
    show_column_names=show_column_names,
    show_row_names=TRUE,
    column_title=title,
    na_col="#F2F2F2",
    border=TRUE,
    heatmap_legend_param=list(
      title_gp=grid::gpar(
        fontsize=fontsize_legend,
        fontface="bold"
      ),
      labels_gp=grid::gpar(
        fontsize=max(fontsize_legend-1,1)
      )
    )
  )
  list(
    heatmap=heatmap,
    plot_matrix=plot_matrix,
    clustering_matrix=clustering_matrix
  )
}
# ============================================================================= #
# Plot export helper
# ============================================================================= #
#' Save ComplexHeatmap object
#'
#' @keywords internal
#' @noRd
.save_gnrh_heatmap <- function(
    heatmap,
    filename,
    width,
    height
) {
  if (is.null(filename)) return(invisible(NULL))
  ext <- tolower(tools::file_ext(filename))
  if (!nzchar(ext)) {
    ext <- "pdf"
    filename <- paste0(filename,".pdf")
  }
  if (!(ext %in% c("pdf","png"))) {
    stop("`filename` must end in .pdf or .png.",call.=FALSE)
  }
  dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE)
  if (ext=="pdf") {
    grDevices::pdf(
      filename,
      width=width,
      height=height,
      useDingbats=FALSE
    )
  } else {
    grDevices::png(
      filename,
      width=width,
      height=height,
      units="in",
      res=300
    )
  }
  on.exit(grDevices::dev.off(),add=TRUE)
  ComplexHeatmap::draw(heatmap)
  invisible(filename)
}
# ============================================================================= #
# Conserved markers
# ============================================================================= #
#' Plot conserved GnRH markers across datasets
#'
#' @inheritParams merge_gnrh_marker_scores
#' @param files Named marker TSV filenames or paths.
#' @param min_datasets Minimum number of datasets supporting a gene.
#' @param top_n Maximum number of genes displayed.
#' @param exclude_genes Optional genes excluded from visualization.
#' @param scale_rows Row-standardize marker scores.
#' @param cluster_rows,cluster_columns Cluster heatmap rows or columns.
#' @param show_column_names Show dataset labels.
#' @param fontsize_row,fontsize_col,fontsize_legend Font sizes.
#' @param filename Optional PDF or PNG filename.
#' @param width,height Optional figure dimensions in inches.
#'
#' @return List containing heatmap matrices and conserved-marker summary.
#' @export
plot_gnrh_conserved_markers <- function(
    files,
    dir=".",
    gene_col="gene",
    score_col=NULL,
    dataset_names=names(files),
    gene_case=c("upper","asis"),
    min_datasets=2L,
    top_n=50L,
    exclude_genes=NULL,
    scale_rows=TRUE,
    cluster_rows=TRUE,
    cluster_columns=FALSE,
    show_column_names=TRUE,
    fontsize_row=NULL,
    fontsize_col=NULL,
    fontsize_legend=10,
    filename=NULL,
    width=NULL,
    height=NULL
) {
  gene_case <- match.arg(gene_case)
  mat <- merge_gnrh_marker_scores(
    files,
    dataset_names=dataset_names,
    dir=dir,
    gene_col=gene_col,
    score_col=score_col,
    gene_case=gene_case,
    fill=NA_real_,
    require_same_score=TRUE
  )
  # ------------------------------------------------------------------------- #
  # Support statistics
  # ------------------------------------------------------------------------- #
  support <- rowSums(is.finite(mat))
  mean_score <- rowMeans(mat,na.rm=TRUE)
  median_score <- apply(mat,1,stats::median,na.rm=TRUE)
  min_score <- apply(mat,1,function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    min(x)
  })
  conserved <- data.frame(
    gene=rownames(mat),
    n_datasets=support,
    dataset_fraction=support/ncol(mat),
    mean_score=mean_score,
    median_score=median_score,
    min_score=min_score,
    stringsAsFactors=FALSE
  )
  # ------------------------------------------------------------------------- #
  # Filter
  # ------------------------------------------------------------------------- #
  conserved <- conserved[
    conserved$n_datasets>=min_datasets,
    ,
    drop=FALSE
  ]
  if (length(exclude_genes)) {
    exclude_genes <- if (gene_case=="upper") toupper(exclude_genes) else exclude_genes
    conserved <- conserved[
      !conserved$gene %in% exclude_genes,
      ,
      drop=FALSE
    ]
  }
  conserved <- conserved[
    order(
      -conserved$n_datasets,
      -conserved$median_score,
      -conserved$mean_score
    ),
    ,
    drop=FALSE
  ]
  if (!is.null(top_n)) conserved <- utils::head(conserved,as.integer(top_n))
  if (!nrow(conserved)) stop("No conserved GnRH markers meet `min_datasets`.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Heatmap
  # ------------------------------------------------------------------------- #
  selected <- mat[conserved$gene,,drop=FALSE]
  result <- .gnrh_marker_heatmap(
    selected,
    scale_rows=scale_rows,
    cluster_rows=cluster_rows,
    cluster_columns=cluster_columns,
    show_column_names=show_column_names,
    fontsize_row=fontsize_row,
    fontsize_col=fontsize_col,
    fontsize_legend=fontsize_legend,
    title="Conserved GnRH markers"
  )
  # ------------------------------------------------------------------------- #
  # Adaptive dimensions
  # ------------------------------------------------------------------------- #
  dims <- .gnrh_plot_dims(
    n_datasets=ncol(selected),
    n_genes=nrow(selected),
    type="heatmap"
  )
  if (is.null(width)) width <- unname(dims["width"])
  if (is.null(height)) height <- unname(dims["height"])
  .save_gnrh_heatmap(
    result$heatmap,
    filename=filename,
    width=width,
    height=height
  )
  c(
    result,
    list(
      matrix=mat,
      conserved=conserved,
      score_columns=attr(mat,"score_columns"),
      score_type=attr(mat,"score_type"),
      width=width,
      height=height
    )
  )
}
# ============================================================================= #
# Dataset-specific markers
# ============================================================================= #
#' Plot dataset-specific GnRH markers
#'
#' @inheritParams plot_gnrh_conserved_markers
#' @param top_n_per_dataset Maximum markers selected per dataset. If
#'   \code{NULL}, automatically adapted to dataset number.
#' @param max_datasets Maximum number of datasets supporting a marker.
#' @param min_target_score Optional minimum raw target-dataset score.
#'
#' @return List containing heatmap, score matrices, and specificity table.
#' @export
plot_gnrh_dataset_specific_markers <- function(
    files,
    dir=".",
    gene_col="gene",
    score_col=NULL,
    dataset_names=names(files),
    gene_case=c("upper","asis"),
    top_n_per_dataset=NULL,
    max_datasets=2L,
    min_target_score=NULL,
    exclude_genes=NULL,
    scale_rows=TRUE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=TRUE,
    fontsize_row=NULL,
    fontsize_col=NULL,
    fontsize_legend=10,
    filename=NULL,
    width=NULL,
    height=NULL
) {
  gene_case <- match.arg(gene_case)
  mat <- merge_gnrh_marker_scores(
    files,
    dataset_names=dataset_names,
    dir=dir,
    gene_col=gene_col,
    score_col=score_col,
    gene_case=gene_case,
    fill=NA_real_,
    require_same_score=TRUE
  )
  n_datasets <- ncol(mat)
  # ------------------------------------------------------------------------- #
  # Adaptive marker number
  # ------------------------------------------------------------------------- #
  if (is.null(top_n_per_dataset)) {
    top_n_per_dataset <- if (n_datasets<=4L) {
      25L
    } else if (n_datasets<=8L) {
      15L
    } else {
      10L
    }
  }
  # ------------------------------------------------------------------------- #
  # Eligibility
  # ------------------------------------------------------------------------- #
  support <- rowSums(is.finite(mat))
  eligible <- support<=max_datasets
  if (length(exclude_genes)) {
    exclude_genes <- if (gene_case=="upper") toupper(exclude_genes) else exclude_genes
    eligible[rownames(mat) %in% exclude_genes] <- FALSE
  }
  # ------------------------------------------------------------------------- #
  # Dataset specificity
  # ------------------------------------------------------------------------- #
  specificity <- matrix(
    -Inf,
    nrow=nrow(mat),
    ncol=ncol(mat),
    dimnames=dimnames(mat)
  )
  for (i in seq_len(ncol(mat))) {
    target <- mat[,i]
    for (g in seq_len(nrow(mat))) {
      if (!eligible[[g]] || !is.finite(target[[g]])) next
      other_values <- mat[g,-i,drop=TRUE]
      other_values <- other_values[is.finite(other_values)]
      contrast <- if (length(other_values)) {
        target[[g]]-mean(other_values)
      } else {
        target[[g]]
      }
      rarity_weight <- 1+(1-support[[g]]/ncol(mat))
      specificity[g,i] <- contrast*rarity_weight
    }
  }
  # ------------------------------------------------------------------------- #
  # Select genes
  # ------------------------------------------------------------------------- #
  selected_rows <- lapply(seq_len(ncol(specificity)),function(i) {
    value <- specificity[,i]
    valid <- is.finite(value)
    if (!is.null(min_target_score)) {
      valid <- valid & is.finite(mat[,i]) & mat[,i]>=min_target_score
    }
    value <- value[valid]
    if (!length(value)) return(character(0))
    value <- sort(value,decreasing=TRUE)
    names(utils::head(value,top_n_per_dataset))
  })
  selected <- unique(unlist(selected_rows,use.names=FALSE))
  if (!length(selected)) stop("No dataset-specific GnRH markers meet the requested criteria.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Best dataset
  # ------------------------------------------------------------------------- #
  selected_specificity <- specificity[selected,,drop=FALSE]
  best_index <- max.col(selected_specificity,ties.method="first")
  best_dataset <- colnames(mat)[best_index]
  best_score <- selected_specificity[
    cbind(seq_along(selected),best_index)
  ]
  raw_score <- mat[
    cbind(match(selected,rownames(mat)),best_index)
  ]
  order_index <- order(
    match(best_dataset,colnames(mat)),
    -best_score
  )
  selected <- selected[order_index]
  best_dataset <- best_dataset[order_index]
  best_score <- best_score[order_index]
  raw_score <- raw_score[order_index]
  selected_matrix <- mat[selected,,drop=FALSE]
  specific_table <- data.frame(
    gene=selected,
    best_dataset=best_dataset,
    n_datasets=support[selected],
    dataset_fraction=support[selected]/ncol(mat),
    specificity_score=best_score,
    raw_score=raw_score,
    stringsAsFactors=FALSE
  )
  # ------------------------------------------------------------------------- #
  # Heatmap
  # ------------------------------------------------------------------------- #
  row_split <- factor(
    best_dataset,
    levels=colnames(mat)
  )
  result <- .gnrh_marker_heatmap(
    selected_matrix,
    row_split=row_split,
    scale_rows=scale_rows,
    cluster_rows=cluster_rows,
    cluster_columns=cluster_columns,
    show_column_names=show_column_names,
    fontsize_row=fontsize_row,
    fontsize_col=fontsize_col,
    fontsize_legend=fontsize_legend,
    title="Dataset-specific GnRH markers"
  )
  # ------------------------------------------------------------------------- #
  # Adaptive dimensions
  # ------------------------------------------------------------------------- #
  dims <- .gnrh_plot_dims(
    n_datasets=ncol(selected_matrix),
    n_genes=nrow(selected_matrix),
    type="heatmap"
  )
  if (is.null(width)) width <- unname(dims["width"])
  if (is.null(height)) height <- unname(dims["height"])
  .save_gnrh_heatmap(
    result$heatmap,
    filename=filename,
    width=width,
    height=height
  )
  # ------------------------------------------------------------------------- #
  # Return
  # ------------------------------------------------------------------------- #
  c(
    result,
    list(
      matrix=mat,
      specificity=specificity,
      specific_table=specific_table,
      score_columns=attr(mat,"score_columns"),
      score_type=attr(mat,"score_type"),
      top_n_per_dataset=top_n_per_dataset,
      width=width,
      height=height
    )
  )
}
