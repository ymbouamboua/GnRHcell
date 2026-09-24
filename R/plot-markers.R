# ========================================================================= #
# Marker-program expression dot plot
# ========================================================================= #

#' Plot expression of GnRH marker programs
#'
#' Selects genes from a table returned by [find_gnrh_programs()] and displays
#' their expression across a Seurat metadata grouping. Genes are grouped by
#' their assigned biological program, ranked by integrated marker evidence,
#' and deduplicated across input datasets.
#'
#' @param object A Seurat object, or a collection returned by
#'   [run_gnrh_collection()]. For a collection, provide \code{selected_ids}.
#' @param programs Result returned by [find_gnrh_programs()]. When
#'   \code{object} is a collection, the default uses
#'   \code{object$comparisons$programs}.
#' @param selected_ids Dataset IDs or labels to combine when \code{object} is a
#'   collection.
#' @param gnrh_only Retain only GnRH-positive cells when combining a collection.
#' @param label_by Use dataset labels or IDs in the generated grouping column.
#' @param allow_saved Allow loading saved processed objects when collection
#'   objects are not retained in memory.
#' @param table Marker table to display: \code{"high_confidence"},
#'   \code{"context_specific"}, or \code{"candidate_table"}.
#' @param group.by Metadata column defining the displayed cell groups.
#' @param top_n Maximum number of genes displayed per biological program.
#' @param programs_keep Optional character vector of programs to retain.
#' @param program_order Display order for biological programs. Programs not
#'   listed are appended after the requested developmental sequence.
#' @param min_confidence Minimum confidence used with \code{"candidate_table"}.
#' @param include_unassigned Include genes without a program assignment.
#' @param scale Logical passed to [Seurat::DotPlot()]. The default
#'   \code{FALSE} preserves average-expression values and is recommended when
#'   comparing only a small number of groups.
#' @param title Optional plot title.
#' @param dot.scale Dot-size scaling passed to [gnrh_celldot()].
#' @param txtsize Base text size.
#' @param ... Additional arguments passed to [gnrh_celldot()].
#'
#' @return A \code{ggplot2} object. The selected gene groups are available in
#'   the \code{"gnrh_features"} attribute.
#'
#' @examples
#' \dontrun{
#' gnrh_program(
#'   combined_gnrh,
#'   programs = comparisons$programs,
#'   table = "high_confidence",
#'   group.by = "dataset"
#' )
#' }
#'
#' @export
gnrh_program <- function(
    object,
    programs = NULL,
    selected_ids = NULL,
    gnrh_only = TRUE,
    label_by = c("label", "id"),
    allow_saved = TRUE,
    table = c("high_confidence", "context_specific", "candidate_table"),
    group.by = "seurat_clusters",
    top_n = 8L,
    programs_keep = NULL,
    program_order = c("early", "migrating", "mature", "secretory"),
    min_confidence = c("medium", "high"),
    include_unassigned = FALSE,
    scale = FALSE,
    title = NULL,
    dot.scale = 6,
    txtsize = getOption("gnrhcell.base_size", 14),
    ...
) {
  if (is.list(object) && !inherits(object, "Seurat") && !is.null(object$results)) {
    collection <- object
    if (is.null(selected_ids)) {
      stop("`selected_ids` is required when `object` is a GnRHcell collection.", call. = FALSE)
    }
    if (is.null(programs)) programs <- collection$comparisons$programs %||% NULL
    object <- combine_gnrh_datasets(
      collection = collection,
      selected_ids = selected_ids,
      gnrh_only = gnrh_only,
      group.by = group.by,
      label_by = match.arg(label_by),
      allow_saved = allow_saved
    )
  }
  .check_seurat(object)
  table <- match.arg(table)

  if (!is.list(programs) || is.null(programs[[table]])) {
    stop("`programs` must be a result from `find_gnrh_programs()` containing `", table, "`.", call. = FALSE)
  }
  df <- programs[[table]]
  if (!is.data.frame(df) || !nrow(df)) {
    stop("The selected marker-program table is empty: ", table, ".", call. = FALSE)
  }
  required <- c("gene", "program")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Marker-program table is missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Grouping column not found: ", group.by, call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1L) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }

  df$gene <- toupper(trimws(as.character(df$gene)))
  df$program <- as.character(df$program)
  df <- df[!is.na(df$gene) & nzchar(df$gene) & df$gene %in% rownames(object), , drop = FALSE]

  if (identical(table, "candidate_table") && "confidence_level" %in% names(df)) {
    df <- df[as.character(df$confidence_level) %in% min_confidence, , drop = FALSE]
  }
  if (!isTRUE(include_unassigned)) {
    df <- df[!is.na(df$program) & nzchar(df$program) & df$program != "unassigned", , drop = FALSE]
  }
  if (!is.null(programs_keep)) {
    df <- df[df$program %in% programs_keep, , drop = FALSE]
  }
  if (!nrow(df)) {
    stop("No marker-program genes remain after filtering and matching object features.", call. = FALSE)
  }

  # Marker tables contain one row per gene and input dataset. Retain the
  # strongest evidence row before selecting the leading genes per program.
  rank_columns <- intersect(c("integrated_score", "marker_evidence_score", "n_datasets"), names(df))
  if (length(rank_columns)) {
    ordering <- lapply(rank_columns, function(column) -suppressWarnings(as.numeric(df[[column]])))
    ordering <- c(list(df$program), ordering, list(df$gene))
    df <- df[do.call(order, c(ordering, list(na.last = TRUE))), , drop = FALSE]
  } else {
    df <- df[order(df$program, df$gene), , drop = FALSE]
  }
  # A gene can be assigned to more than one program, but Seurat dot plots
  # require globally unique feature levels. Keep its highest-ranked assignment.
  df <- df[!duplicated(df$gene), , drop = FALSE]
  feature_groups <- split(df$gene, df$program)
  feature_groups <- lapply(feature_groups, utils::head, n = as.integer(top_n))
  feature_groups <- feature_groups[lengths(feature_groups) > 0L]
  program_rank <- match(
    tolower(names(feature_groups)),
    tolower(program_order),
    nomatch = length(program_order) + 1L
  )
  feature_groups <- feature_groups[order(program_rank, seq_along(feature_groups))]
  names(feature_groups) <- tools::toTitleCase(gsub("_", " ", names(feature_groups), fixed = TRUE))

  if (is.null(title)) {
    title <- switch(
      table,
      high_confidence = "High-confidence GnRH marker programs",
      context_specific = "Context-specific GnRH marker programs",
      candidate_table = "GnRH candidate marker programs"
    )
  }

  plot <- gnrh_celldot(
    object = object,
    features = feature_groups,
    group.by = group.by,
    scale = scale,
    dot.scale = dot.scale,
    flip = TRUE,
    txtsize = txtsize,
    title = title,
    x.ang = 60,
    ...
  )
  attr(plot, "gnrh_features") <- feature_groups
  plot
}

# ========================================================================= #
# GnRH association network
# ========================================================================= #
#' GnRH marker association network
#'
#' Build a GnRH marker network using co-expression, co-detection, or phenotype
#' association evidence. In \code{"auto"} mode, the most informative available
#' association metric is selected automatically. GNRH1 can be added as a
#' reference hub, and isolated genes are removed.
#'
#' @param df GnRH marker table.
#' @param top_n Maximum number of top-ranked marker genes considered.
#' @param threshold Minimum edge similarity.
#' @param mode Association mode. One of \code{"auto"}, \code{"coexpression"},
#'   \code{"codetection"}, or \code{"phenotype"}.
#' @param score_col Optional column used to scale node size.
#' @param association_col Optional association column overriding automatic
#'   column selection.
#' @param include_gnrh Logical. Add GNRH1 as a reference node.
#' @param gnrh_gene Name of the GnRH gene.
#' @param txtsize Base text size.
#' @param style Plot style passed to \code{gnrh_theme()}.
#' @param seed Random seed used by the graph layout.
#'
#' @return A ggraph object.
#'
#' @export
gnrh_network <- function(
    df,
    top_n=25L,
    threshold=0.10,
    mode=c("auto","coexpression","codetection","phenotype"),
    score_col=NULL,
    association_col=NULL,
    include_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize = getOption("gnrhcell.base_size", 14),
    style="void",
    seed=1234
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  if (!is.data.frame(df) || !nrow(df)) stop("`df` must be a non-empty data frame.",call.=FALSE)
  if (!"gene" %in% names(df)) stop("Missing column: gene",call.=FALSE)
  if (!is.numeric(top_n) || length(top_n)!=1L || !is.finite(top_n) || top_n<2L) stop("`top_n` must be >= 2.",call.=FALSE)
  if (!is.numeric(threshold) || length(threshold)!=1L || !is.finite(threshold)) stop("`threshold` must be a finite numeric value.",call.=FALSE)
  mode <- match.arg(mode)
  # ========================================================================= #
  # Association metric
  # ========================================================================= #
  if (is.null(association_col)) {
    if (mode=="coexpression") association_col <- intersect(c("coexpr_cor","coexpr"),names(df))[1L]
    if (mode=="codetection") association_col <- intersect(c("codetect_log2or","codetect_or"),names(df))[1L]
    if (mode=="phenotype") association_col <- intersect(c("phenotype_log2or","phenotype_or","specificity","spec"),names(df))[1L]
    if (mode=="auto") {
      candidates <- c(
        phenotype=intersect(c("phenotype_log2or","phenotype_or"),names(df))[1L],
        codetection=intersect(c("codetect_log2or","codetect_or"),names(df))[1L],
        coexpression=intersect(c("coexpr_cor","coexpr"),names(df))[1L]
      )
      candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
      if (!length(candidates)) stop("No association metric found.",call.=FALSE)
      quality <- vapply(candidates,function(x) {
        z <- suppressWarnings(as.numeric(df[[x]]))
        z <- z[is.finite(z)]
        if (length(z)<2L) return(-Inf)
        stats::sd(z,na.rm=TRUE)*sqrt(length(unique(z)))
      },numeric(1))
      resolved <- names(which.max(quality))
      association_col <- candidates[[resolved]]
      mode <- resolved
    }
  }
  if (!length(association_col) || is.na(association_col) || !association_col %in% names(df)) stop("No valid association column found.",call.=FALSE)
  # ========================================================================= #
  # Score
  # ========================================================================= #
  if (is.null(score_col)) score_col <- intersect(c("marker_score","final_score","association_score","specificity","avg_log2FC","avg_logFC"),names(df))[1L]
  if (!length(score_col) || is.na(score_col) || !score_col %in% names(df)) {
    score_col <- NULL
    df$.network_score <- 1
  } else {
    score <- suppressWarnings(as.numeric(df[[score_col]]))
    score[!is.finite(score)] <- 0
    df$.network_score <- abs(score)
    if (!any(df$.network_score>0)) df$.network_score <- 1
  }
  # ========================================================================= #
  # Select markers
  # ========================================================================= #
  df$gene <- as.character(df$gene)
  df$.association <- suppressWarnings(as.numeric(df[[association_col]]))
  if (grepl("_or$",association_col) && !grepl("log2",association_col)) df$.association <- log2(pmax(df$.association,.Machine$double.xmin))
  df <- df[!is.na(df$gene) & nzchar(df$gene) & is.finite(df$.association),,drop=FALSE]
  df <- df[!duplicated(df$gene),,drop=FALSE]
  if (isTRUE(include_gnrh)) df <- df[toupper(df$gene)!=toupper(gnrh_gene),,drop=FALSE]
  rank_col <- if (!is.null(score_col)) ".network_score" else ".association"
  df <- df[order(df[[rank_col]],df$.association,decreasing=TRUE,na.last=TRUE),,drop=FALSE]
  df <- utils::head(df,min(as.integer(top_n),nrow(df)))
  if (nrow(df)<2L) stop("At least two marker genes are required.",call.=FALSE)
  # ========================================================================= #
  # Normalize association
  # ========================================================================= #
  assoc <- pmax(df$.association,0)
  if (max(assoc,na.rm=TRUE)>0) assoc <- assoc/max(assoc,na.rm=TRUE)
  df$.association_scaled <- assoc
  genes <- df$gene
  # ========================================================================= #
  # Marker-marker similarity
  # ========================================================================= #
  mat <- outer(df$.association_scaled,df$.association_scaled,pmin)
  rownames(mat) <- genes
  colnames(mat) <- genes
  edges <- as.data.frame(as.table(mat),stringsAsFactors=FALSE)
  names(edges) <- c("from","to","weight")
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$weight <- suppressWarnings(as.numeric(edges$weight))
  edges$type <- "marker"
  edges <- edges[edges$from<edges$to & is.finite(edges$weight) & edges$weight>threshold,,drop=FALSE]
  # ========================================================================= #
  # GNRH1 reference edges
  # ========================================================================= #
  if (isTRUE(include_gnrh)) {
    gnrh_edges <- data.frame(from=gnrh_gene,to=df$gene,weight=df$.association_scaled,type="GNRH",stringsAsFactors=FALSE)
    gnrh_edges <- gnrh_edges[is.finite(gnrh_edges$weight) & gnrh_edges$weight>threshold,,drop=FALSE]
    edges <- rbind(edges,gnrh_edges)
  }
  if (!nrow(edges)) stop("No network edges remain above `threshold`.",call.=FALSE)
  # ========================================================================= #
  # Connected genes
  # ========================================================================= #
  connected_genes <- unique(c(edges$from,edges$to))
  df <- df[df$gene %in% connected_genes,,drop=FALSE]
  edges <- edges[edges$from %in% connected_genes & edges$to %in% connected_genes,,drop=FALSE]
  if (!nrow(df)) stop("No connected marker genes remain above `threshold`.",call.=FALSE)
  # ========================================================================= #
  # Nodes
  # ========================================================================= #
  nodes <- data.frame(
    name=df$gene,
    score=df$.network_score,
    association=df$.association_scaled,
    reference=FALSE,
    stringsAsFactors=FALSE
  )
  if (isTRUE(include_gnrh) && gnrh_gene %in% connected_genes) {
    gnrh_score <- if (length(nodes$score) && any(is.finite(nodes$score))) max(nodes$score,na.rm=TRUE)*1.25 else 1.25
    nodes <- rbind(
      data.frame(
        name=gnrh_gene,
        score=gnrh_score,
        association=1,
        reference=TRUE,
        stringsAsFactors=FALSE
      ),
      nodes
    )
  }
  # ========================================================================= #
  # Graph
  # ========================================================================= #
  g <- igraph::graph_from_data_frame(edges,vertices=nodes,directed=FALSE)
  isolated <- igraph::degree(g)==0L
  if (any(isolated)) g <- igraph::delete_vertices(g,igraph::V(g)[isolated])
  if (igraph::vcount(g)<2L || igraph::ecount(g)<1L) stop("The network contains too few connected genes.",call.=FALSE)
  # ========================================================================= #
  # Plot
  # ========================================================================= #
  set.seed(seed)
  p <- ggraph::ggraph(g,layout="fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width=.data$weight,alpha=.data$weight),
      colour="grey60",
      show.legend=FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size=.data$score,
        fill=.data$association,
        shape=.data$reference
      ),
      colour="black",
      stroke=0.3
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        label=.data$name,
        fontface=ifelse(.data$reference,"bold","italic")
      ),
      repel=TRUE,
      max.overlaps=Inf,
      size=txtsize/3
    ) +
    ggraph::scale_edge_width(range=c(0.25,1.6)) +
    ggraph::scale_edge_alpha(range=c(0.2,0.75)) +
    ggplot2::scale_size_continuous(
      range=c(3,10),
      name=if (is.null(score_col)) "Score" else score_col
    ) +
    ggplot2::scale_fill_viridis_c(
      name="Association",
      guide=ggplot2::guide_colourbar(
        frame.colour="black",
        frame.linewidth=0.35,
        ticks.colour="black",
        ticks.linewidth=0.35
      )
    ) +
    ggplot2::scale_shape_manual(
      values=c(`FALSE`=21,`TRUE`=23),
      labels=c(`FALSE`="Marker",`TRUE`=gnrh_gene)
    ) +
    ggplot2::guides(
      shape=ggplot2::guide_legend(title=NULL)
    ) +
    ggplot2::labs(
      title=paste0("GnRH ",mode," network"),
      subtitle=paste0("Association: ",association_col)
    ) +
    gnrh_theme(txtsize=txtsize,style=style) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",hjust=0.5),
      plot.subtitle=ggplot2::element_text(hjust=0.5)
    )
  # ========================================================================= #
  # Remove graph coordinates
  # ========================================================================= #
  if (identical(style,"void")) {
    p <- p +
      gnrh_theme(
        style="void", txtsize=txtsize, leg.pos="right",
        plot.margin=c(15,15,15,15)
      ) +
      ggplot2::theme(
        plot.title=ggplot2::element_text(face="bold",hjust=0.5),
        plot.subtitle=ggplot2::element_text(hjust=0.5),
        plot.margin=ggplot2::margin(t=15, r=15, b=15, l=15, unit="pt")
      ) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=0.10)) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0.10))
  }
  p
}



# ========================================================================= #
# Coexpression markers
# ========================================================================= #
#' Plot GNRH1-associated markers
#'
#' Visualize the strongest GNRH1-associated markers using co-expression,
#' differential-expression significance, and marker score.
#'
#' @param df GnRH marker table.
#' @param coexp_cutoff Minimum GNRH1 co-expression value.
#' @param top_n Maximum number of genes displayed. Use \code{NULL} to display
#'   all genes passing \code{coexp_cutoff}.
#' @param txtsize Base text size.
#' @param style Plot style.
#' @param coexpr_col Optional co-expression column.
#' @param score_col Optional marker-score column.
#' @param exclude_gnrh Logical. Exclude GNRH1 itself from the plot.
#' @param gnrh_gene GnRH gene name.
#' @param point.size Point-size range.
#' @param title Plot title.
#' @param subtitle Optional plot subtitle.
#' @param legend Logical. Display legends.
#'
#' @return A ggplot object.
#' @export
gnrh_coexpr <- function(
    df,
    coexp_cutoff=0.10,
    top_n=40L,
    txtsize = getOption("gnrhcell.base_size", 14),
    style="bw",
    coexpr_col=NULL,
    score_col=NULL,
    exclude_gnrh=TRUE,
    gnrh_gene="GNRH1",
    point.size=c(1.2,5),
    title="GNRH1-associated markers",
    subtitle=NULL,
    legend=TRUE
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  if (!is.data.frame(df) || !nrow(df)) stop("`df` must be a non-empty data frame.",call.=FALSE)
  required <- c("gene","p_val_adj")
  missing <- setdiff(required,names(df))
  if (length(missing)) stop("Missing columns: ",paste(missing,collapse=", "),call.=FALSE)
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n)!=1L || !is.finite(top_n) || top_n<1L)) stop("`top_n` must be NULL or a positive integer.",call.=FALSE)
  # ========================================================================= #
  # Columns
  # ========================================================================= #
  coexpr_col <- coexpr_col %||% intersect(c("coexpr_cor","coexpr","coexpression","coexpr_pct"),names(df))[1L]
  score_col <- score_col %||% intersect(c("final_score","marker_score","specificity_score","score","avg_log2FC","avg_logFC"),names(df))[1L]
  if (!length(coexpr_col) || is.na(coexpr_col)) stop("No co-expression column found.",call.=FALSE)
  if (!length(score_col) || is.na(score_col)) stop("No marker score column found.",call.=FALSE)
  # ========================================================================= #
  # Prepare
  # ========================================================================= #
  df$gene <- as.character(df$gene)
  df[[coexpr_col]] <- suppressWarnings(as.numeric(df[[coexpr_col]]))
  df[[score_col]] <- suppressWarnings(as.numeric(df[[score_col]]))
  df$p_val_adj <- suppressWarnings(as.numeric(df$p_val_adj))
  keep <- !is.na(df$gene) & nzchar(df$gene) & is.finite(df[[coexpr_col]]) & df[[coexpr_col]]>=coexp_cutoff & is.finite(df$p_val_adj)
  if (isTRUE(exclude_gnrh)) keep <- keep & toupper(df$gene)!=toupper(gnrh_gene)
  df <- df[keep,,drop=FALSE]
  if (!nrow(df)) stop("No genes pass `coexp_cutoff`.",call.=FALSE)
  # ========================================================================= #
  # Rank
  # ========================================================================= #
  df <- df[order(df[[coexpr_col]],df[[score_col]],decreasing=TRUE,na.last=TRUE),,drop=FALSE]
  df <- df[!duplicated(df$gene),,drop=FALSE]
  if (!is.null(top_n)) df <- utils::head(df,as.integer(top_n))
  df$log_padj <- -log10(pmax(df$p_val_adj,.Machine$double.xmin))
  df$marker_size <- abs(df[[score_col]])
  df$marker_size[!is.finite(df$marker_size)] <- 0
  df$gene <- factor(df$gene,levels=rev(df$gene))
  # ========================================================================= #
  # Plot
  # ========================================================================= #
  p <- ggplot2::ggplot(df,ggplot2::aes(x=.data[[coexpr_col]],y=.data$gene,size=.data$marker_size,fill=.data$log_padj)) +
    ggplot2::geom_point(shape=21,colour="black",stroke=0.18,alpha=0.92) +
    ggplot2::scale_size_continuous(name=score_col,range=point.size) +
    ggplot2::scale_fill_viridis_c(
      name=expression(-log[10]("adj. p")),
      guide=ggplot2::guide_colourbar(frame.colour="black",frame.linewidth=0.3,ticks.colour="black")
    ) +
    ggplot2::labs(x="Correlation with GNRH1",y=NULL,title=title,subtitle=subtitle) +
    gnrh_theme(style=style,txtsize=txtsize) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",hjust=0.5),
      axis.text.y=ggplot2::element_text(face="italic"),
      legend.position=if (isTRUE(legend)) "right" else "none",
      legend.key.height=grid::unit(0.38,"cm"),
      legend.spacing.y=grid::unit(0.08,"cm")
    )
  p
}



# ========================================================================= #
# GNRH1 co-detection
# ========================================================================= #
#' Plot GNRH1 co-detection marker associations
#'
#' Visualize binary co-detection enrichment between candidate markers and
#' GNRH1.
#'
#' @param markers Marker table returned by \code{gnrh_markers()}.
#' @param top_n Maximum number of genes to label.
#' @param min_or Minimum co-detection odds ratio.
#' @param max_fdr Maximum co-detection FDR.
#' @param min_specificity Optional minimum marker specificity.
#' @param exclude_gnrh Logical. Exclude GNRH1 itself.
#' @param gnrh_gene GnRH gene name.
#' @param txtsize Base text size.
#' @param label.size Label text size.
#' @param point.size Point-size range.
#' @param style Plot style.
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#' @param legend Logical. Display legends.
#' @param verbose Print plotting information.
#'
#' @return A ggplot object.
#' @export
gnrh_codetect <- function(
    markers,
    top_n=12L,
    min_or=2,
    max_fdr=0.05,
    min_specificity=0.05,
    exclude_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize = getOption("gnrhcell.base_size", 14),
    label.size=3,
    point.size=c(1.2,5),
    style="bw",
    title="GNRH1 co-detection",
    subtitle=NULL,
    legend=TRUE,
    verbose=TRUE
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  if (!is.data.frame(markers) || !nrow(markers)) stop("`markers` must be a non-empty data frame.",call.=FALSE)
  req <- c("gene","codetect_or","codetect_fdr")
  miss <- setdiff(req,names(markers))
  if (length(miss)) stop("Missing columns: ",paste(miss,collapse=", "),call.=FALSE)
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n)!=1L || !is.finite(top_n) || top_n<1L)) stop("`top_n` must be NULL or a positive integer.",call.=FALSE)
  # ========================================================================= #
  # Prepare
  # ========================================================================= #
  df <- markers
  df$gene <- as.character(df$gene)
  df$codetect_or <- suppressWarnings(as.numeric(df$codetect_or))
  df$codetect_fdr <- suppressWarnings(as.numeric(df$codetect_fdr))
  spec_col <- intersect(c("specificity","spec"),names(df))[1L]
  df$specificity <- if (!length(spec_col) || is.na(spec_col)) 1 else suppressWarnings(as.numeric(df[[spec_col]]))
  keep <- !is.na(df$gene) & nzchar(df$gene) & is.finite(df$codetect_or) & df$codetect_or>0 & is.finite(df$codetect_fdr) & is.finite(df$specificity)
  if (isTRUE(exclude_gnrh)) keep <- keep & toupper(df$gene)!=toupper(gnrh_gene)
  df <- df[keep,,drop=FALSE]
  if (!nrow(df)) stop("No valid genes available for co-detection plotting.",call.=FALSE)
  # ========================================================================= #
  # Variables
  # ========================================================================= #
  df$log2_or <- log2(df$codetect_or)
  df$log_fdr <- -log10(pmax(df$codetect_fdr,.Machine$double.xmin))
  df$significant <- df$codetect_or>=min_or & df$codetect_fdr<=max_fdr
  if (!is.null(min_specificity)) df$significant <- df$significant & df$specificity>=min_specificity
  df$label_score <- pmax(df$log2_or,0)*df$log_fdr*pmax(df$specificity,0)
  df$label_score[!is.finite(df$label_score)] <- 0
  # ========================================================================= #
  # Labels
  # ========================================================================= #
  lab <- df[df$significant,,drop=FALSE]
  if (nrow(lab)) {
    lab <- lab[order(-lab$label_score,-lab$specificity,lab$codetect_fdr),,drop=FALSE]
    lab <- lab[!duplicated(lab$gene),,drop=FALSE]
    if (!is.null(top_n)) lab <- utils::head(lab,as.integer(top_n))
  }
  # ========================================================================= #
  # Plot
  # ========================================================================= #
  label_col <- .gnrh_label_colour(style)
  p <- ggplot2::ggplot(df,ggplot2::aes(x=.data$log2_or,y=.data$log_fdr)) +
    ggplot2::geom_vline(xintercept=log2(min_or),linetype="dashed",linewidth=0.4) +
    ggplot2::geom_hline(yintercept=-log10(max_fdr),linetype="dashed",linewidth=0.4) +
    ggplot2::geom_point(
      ggplot2::aes(size=.data$specificity,fill=.data$log_fdr),
      shape=21,colour="black",stroke=0.15,alpha=0.85
    ) +
    ggrepel::geom_text_repel(
      data=lab,
      ggplot2::aes(label=.data$gene),
      size=label.size,
      fontface="italic",
      colour=label_col,
      box.padding=0.3,
      point.padding=0.15,
      min.segment.length=0,
      segment.colour=label_col,
      segment.alpha=0.5,
      max.overlaps=Inf,
      seed=1234,
      show.legend=FALSE
    ) +
    ggplot2::scale_size_continuous(name="Specificity",range=point.size) +
    ggplot2::scale_fill_viridis_c(
      name=expression(-log[10]("FDR")),
      guide=ggplot2::guide_colourbar(
        frame.colour=label_col,
        frame.linewidth=0.3,
        ticks.colour=label_col
      )
    ) +
    ggplot2::labs(
      x=expression(log[2]*"(co-detection odds ratio)"),
      y=expression(-log[10]*"(FDR)"),
      title=title,
      subtitle=subtitle
    ) +
    gnrh_theme(style=style,txtsize=txtsize) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",hjust=0.5),
      legend.position=if (isTRUE(legend)) "right" else "none"
    )
  if (isTRUE(verbose)) message("[INFO] Co-detection genes: ",nrow(df)," | enriched: ",sum(df$significant))
  p
}
