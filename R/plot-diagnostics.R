# ========================================================================= #
# GnRH phenotype association
# ========================================================================= #
#' Plot GnRH phenotype detection enrichment
#'
#' Visualize marker detection frequencies in GnRH target cells versus
#' reference cells.
#'
#' @param markers Marker table returned by \code{gnrh_markers()}.
#' @param top_n Maximum number of genes to label.
#' @param min_specificity Minimum target-versus-reference detection difference.
#' @param max_padj Maximum adjusted p-value.
#' @param exclude_gnrh Logical. Exclude GNRH1 itself.
#' @param gnrh_gene GnRH gene name.
#' @param txtsize Base text size.
#' @param label.size Label text size.
#' @param point.size Point-size range.
#' @param style Plot style.
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#' @param legend Logical. Display legends.
#'
#' @return A ggplot object.
#' @export
gnrh_pheno <- function(
    markers,
    top_n=12L,
    min_specificity=0.05,
    max_padj=0.05,
    exclude_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize = getOption("gnrhcell.base_size", 14),
    label.size=3,
    point.size=c(1.2,5),
    style="bw",
    title="GnRH phenotype association",
    subtitle=NULL,
    legend=TRUE
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  if (!is.data.frame(markers) || !nrow(markers)) stop("`markers` must be a non-empty data frame.",call.=FALSE)
  if (!"gene" %in% names(markers)) stop("Missing column: gene",call.=FALSE)
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n)!=1L || !is.finite(top_n) || top_n<1L)) stop("`top_n` must be NULL or a positive integer.",call.=FALSE)
  pct1_col <- intersect(c("pct1","pct.1"),names(markers))[1L]
  pct2_col <- intersect(c("pct2","pct.2"),names(markers))[1L]
  spec_col <- intersect(c("specificity","spec"),names(markers))[1L]
  if (!length(pct1_col) || is.na(pct1_col)) stop("No target detection column found (`pct1`/`pct.1`).",call.=FALSE)
  if (!length(pct2_col) || is.na(pct2_col)) stop("No reference detection column found (`pct2`/`pct.2`).",call.=FALSE)
  if (!length(spec_col) || is.na(spec_col)) stop("No specificity column found.",call.=FALSE)
  # ========================================================================= #
  # Prepare
  # ========================================================================= #
  df <- markers
  df$gene <- as.character(df$gene)
  df$pct1_plot <- suppressWarnings(as.numeric(df[[pct1_col]]))
  df$pct2_plot <- suppressWarnings(as.numeric(df[[pct2_col]]))
  df$specificity_plot <- suppressWarnings(as.numeric(df[[spec_col]]))
  padj_col <- intersect(c("p_val_adj","padj","FDR","fdr"),names(df))[1L]
  df$padj_plot <- if (length(padj_col) && !is.na(padj_col)) suppressWarnings(as.numeric(df[[padj_col]])) else NA_real_
  keep <- !is.na(df$gene) & nzchar(df$gene) & is.finite(df$pct1_plot) & is.finite(df$pct2_plot) & is.finite(df$specificity_plot)
  if (isTRUE(exclude_gnrh)) keep <- keep & toupper(df$gene)!=toupper(gnrh_gene)
  df <- df[keep,,drop=FALSE]
  if (!nrow(df)) stop("No valid genes available for detection plotting.",call.=FALSE)
  # ========================================================================= #
  # Variables
  # ========================================================================= #
  df$log_padj <- ifelse(is.finite(df$padj_plot),-log10(pmax(df$padj_plot,.Machine$double.xmin)),0)
  df$selected <- df$specificity_plot>=min_specificity
  if (any(is.finite(df$padj_plot))) df$selected <- df$selected & (!is.finite(df$padj_plot) | df$padj_plot<=max_padj)
  # ========================================================================= #
  # Labels
  # ========================================================================= #
  lab <- df[df$selected,,drop=FALSE]
  if (nrow(lab)) {
    lab <- lab[order(-lab$specificity_plot,-lab$log_padj,-lab$pct1_plot),,drop=FALSE]
    lab <- lab[!duplicated(lab$gene),,drop=FALSE]
    if (!is.null(top_n)) lab <- utils::head(lab,as.integer(top_n))
  }
  # ========================================================================= #
  # Plot
  # ========================================================================= #
  label_col <- .gnrh_label_colour(style)
  p <- ggplot2::ggplot(df,ggplot2::aes(x=.data$pct2_plot,y=.data$pct1_plot)) +
    ggplot2::geom_abline(intercept=0,slope=1,linetype="dashed",linewidth=0.4) +
    ggplot2::geom_point(
      ggplot2::aes(size=.data$specificity_plot,fill=.data$log_padj),
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
      name=expression(-log[10]("adj. p")),
      guide=ggplot2::guide_colourbar(
        frame.colour=label_col,
        frame.linewidth=0.3,
        ticks.colour=label_col
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x="Detection in reference cells",
      y="Detection in GnRH cells",
      title=title,
      subtitle=subtitle
    ) +
    gnrh_theme(style=style,txtsize=txtsize) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",hjust=0.5),
      legend.position=if (isTRUE(legend)) "right" else "none"
    )
  p
}

# ========================================================================= #
# Runtime
# ========================================================================= #

#' Plot gnrhcell runtime across datasets
#'
#' @param files Optional named vector of run-information files.
#' @param dir Directory searched when `files` is `NULL`.
#' @param pattern File-selection regular expression.
#' @param metric Runtime column to summarize.
#' @param x.ang X-axis label angle.
#' @param txtsize Base text size.
#' @param show_points Show dataset points.
#' @return A ggplot object.
#' @export
gnrh_runtime <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    metric = "total_sec",
    x.ang = 45,
    txtsize = getOption("gnrhcell.base_size", 14),
    show_points = TRUE
) {
  stats <- .load_gnrh_stats(
    files = files,
    dir = dir,
    pattern = pattern
  )

  if (!metric %in% names(stats))
    stop("Metric not found: ", metric, call. = FALSE)

  stats$pos <- as.numeric(stats$pos)
  stats$n_cells <- as.numeric(stats$n_cells)
  stats[[metric]] <- as.numeric(stats[[metric]])
  stats$pos[is.na(stats$pos)] <- 0

  stats <- stats::aggregate(
    cbind(pos, n_cells, runtime = stats[[metric]]) ~ dataset,
    data = stats,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  stats <- stats[order(stats$n_cells), , drop = FALSE]
  stats$label <- factor(
    paste0(
      stats$dataset, "\n",
      round(stats$pos), "/",
      round(stats$n_cells), " GnRH+"
    ),
    levels = paste0(
      stats$dataset, "\n",
      round(stats$pos), "/",
      round(stats$n_cells), " GnRH+"
    )
  )

  p <- ggplot2::ggplot(
    stats,
    ggplot2::aes(label, runtime, group = 1)
  ) +
    ggplot2::geom_line(
      linewidth = 0.8,
      colour = "#819DC7"
    )

  if (show_points)
    p <- p +
    ggplot2::geom_point(
      size = 2,
      colour = gnrh_palette("status")[["pos"]]
    )

  p +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(runtime, 1), "s")),
      vjust = -0.7,
      size = 3
    ) +
    ggplot2::labs(
      title = "gnrhcell runtime across datasets",
      x = "Dataset",
      y = metric
    ) +
    gnrh_theme(
      x.ang = x.ang,
      txtsize = txtsize
    )
}



# gnrh_cellmark <- function(
#     files = NULL,
#     dir = file.path("results", "tables"),
#     pattern = "_gnrh_run_info\\.tsv$",
#     x.ang = 45,
#     txtsize = getOption("gnrhcell.base_size", 14),
#     debug = FALSE
# ) {
#   stats <- .load_gnrh_stats(
#     files = files,
#     dir = dir,
#     pattern = pattern
#   )
#
#   if (debug)
#     print(stats[, intersect(
#       c("source_file", "dataset", "n_cells", "neg", "pos"),
#       names(stats)
#     )])
#
#   stats <- stats[, c("dataset", "pos")]
#   stats$pos <- as.numeric(stats$pos)
#
#   stats <- stats::aggregate(
#     pos ~ dataset,
#     stats,
#     function(x) mean(x, na.rm = TRUE)
#   )
#
#   stats$dataset <- factor(
#     stats$dataset,
#     levels = stats$dataset[order(stats$pos)]
#   )
#
#   ggplot2::ggplot(
#     stats,
#     ggplot2::aes(dataset, pos)
#   ) +
#     ggplot2::geom_col(
#       width = 0.7,
#       fill = gnrh_palette("status")[["pos"]]
#     ) +
#     ggplot2::geom_text(
#       ggplot2::aes(label = round(pos)),
#       vjust = -0.3,
#       size = 3
#     ) +
#     ggplot2::scale_y_continuous(
#       expand = ggplot2::expansion(mult = c(0, 0.08))
#     ) +
#     ggplot2::labs(
#       title = "Detected GnRH cells",
#       x = "Dataset",
#       y = "GnRH+ cells"
#     ) +
#     gnrh_theme(
#       x.ang = x.ang,
#       txtsize = txtsize
#     )
# }


# ========================================================================= #
# gnrh_cellmark
# ========================================================================= #
#' Plot detected GnRH-positive cells across datasets
#'
#' Plot the total number of cells analyzed in each dataset and annotate each
#' bar with the number and percentage of GnRH-positive cells detected.
#'
#' @param files Optional named vector of run-information files.
#' @param dir Directory searched when `files` is `NULL`.
#' @param pattern File-selection regular expression.
#' @param txtsize Base text size.
#' @param label.size Text size for bar annotations.
#' @param digits Number of decimal places shown for GnRH-positive percentages.
#' @param debug Print the imported summary columns.
#' @return A ggplot object.
#' @export
gnrh_cellmark <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    txtsize = getOption("gnrhcell.base_size", 14),
    label.size = 3.5,
    digits = 2,
    debug = FALSE
) {
  stats <- .load_gnrh_stats(files=files,dir=dir,pattern=pattern)
  if (debug) print(stats[,intersect(c("source_file","dataset","n_cells","neg","pos"),names(stats))])
  req <- c("dataset","n_cells","pos")
  miss <- setdiff(req,names(stats))
  if (length(miss)) stop("Missing required columns: ",paste(miss,collapse=", "),call.=FALSE)
  stats <- stats[,req]
  stats$n_cells <- as.numeric(stats$n_cells)
  stats$pos <- as.numeric(stats$pos)
  stats <- stats::aggregate(cbind(n_cells,pos)~dataset,stats,function(x) mean(x,na.rm=TRUE))
  stats$pct_gnrh <- 100*stats$pos/stats$n_cells
  #stats$dataset <- factor(stats$dataset,levels=rev(stats$dataset))
  stats$dataset <- factor(stats$dataset,levels=rev(stats$dataset[order(stats$pos,decreasing=TRUE)]))
  ggplot2::ggplot(stats,ggplot2::aes(.data$dataset,.data$n_cells)) +
    ggplot2::geom_col(width=0.7,fill=gnrh_palette("status")[["pos"]]) +
    ggplot2::geom_text(
      ggplot2::aes(label=paste0(
        scales::comma(round(.data$n_cells))," cells\n",
        scales::comma(round(.data$pos))," GnRH+ (",
        sprintf(paste0("%.",digits,"f"),.data$pct_gnrh),"%)"
      )),
      hjust=-0.08,size=label.size
    ) +
    ggplot2::coord_flip(clip="off") +
    ggplot2::scale_y_continuous(
      labels=scales::label_number(big.mark=","),
      expand=ggplot2::expansion(mult=c(0,0.30))
    ) +
    ggplot2::labs(
      title="Detected GnRH cells",
      x=NULL,
      y="Number of cells"
    ) +
    gnrh_theme(txtsize=txtsize) +
    ggplot2::theme(
      axis.text.y=ggplot2::element_text(face="bold"),
      plot.margin=ggplot2::margin(5.5,110,5.5,5.5)
    )
}



# ========================================================================= #
# Marker programs
# ========================================================================= #
#' Plot GnRH marker program results
#'
#' @param programs Result returned by \code{find_gnrh_programs()}.
#' @param table Result table to visualize.
#' @param type Plot type: bar, dot, or tile.
#' @param min_genes Minimum genes retained for summary plots.
#' @param require_coexpr Retain only markers with co-expression support.
#' @param top_n Maximum genes displayed per dataset for gene-level plots.
#'   If \code{NULL}, automatically adapted to the number of datasets.
#' @param facet_ncol Optional number of facet columns.
#' @param mode Light or dark display mode.
#' @param txtsize Base text size.
#' @param x.ang X-axis label angle.
#' @param style Plot theme.
#'
#' @return A ggplot object.
#' @export
gnrh_program_plot <- function(
    programs,
    table=c("summary","high_confidence","candidate_table","context_specific"),
    type=c("bar","dot","tile"),
    min_genes=1L,
    require_coexpr=FALSE,
    top_n=NULL,
    facet_ncol=NULL,
    mode=c("light","dark"),
    txtsize = getOption("gnrhcell.base_size", 14),
    x.ang=45,
    style=c("bw","test","classic","minimal")
) {
  # ========================================================================= #
  # Setup
  # ========================================================================= #
  table <- match.arg(table)
  type <- match.arg(type)
  style <- match.arg(style)
  mode <- match.arg(mode)
  df <- programs[[table]]
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) stop("Selected table is empty: ",table,call.=FALSE)
  # ========================================================================= #
  # Program summary
  # ========================================================================= #
  if (table=="summary") {
    required <- c("program_dimension","program","confidence_level","n_genes")
    missing <- setdiff(required,colnames(df))
    if (length(missing)) stop("Missing summary columns: ",paste(missing,collapse=", "),call.=FALSE)
    df$n_genes <- suppressWarnings(as.numeric(df$n_genes))
    df <- df[is.finite(df$n_genes) & df$n_genes>=min_genes,,drop=FALSE]
    if (!nrow(df)) stop("No rows remain after `min_genes`.",call.=FALSE)
    n_facets <- length(unique(df$program_dimension))
    if (is.null(facet_ncol)) facet_ncol <- min(3L,max(1L,ceiling(sqrt(n_facets))))
    # ----------------------------------------------------------------------- #
    # Bar
    # ----------------------------------------------------------------------- #
    if (type=="bar") {
      return(
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            x=.data$program,
            y=.data$n_genes,
            fill=.data$confidence_level
          )
        ) +
          ggplot2::geom_col(width=0.72) +
          ggplot2::facet_wrap(
            ~program_dimension,
            scales="free_x",
            ncol=facet_ncol
          ) +
          ggplot2::labs(
            title="GnRH marker programs",
            x=NULL,
            y="Genes",
            fill="Confidence"
          ) +
          gnrh_theme(
            mode=mode,
            txtsize=txtsize,
            x.ang=x.ang,
            style=style
          )
      )
    }
    # ----------------------------------------------------------------------- #
    # Dot
    # ----------------------------------------------------------------------- #
    if (type=="dot") {
      return(
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            x=.data$program,
            y=.data$confidence_level,
            size=.data$n_genes,
            colour=.data$program_dimension
          )
        ) +
          ggplot2::geom_point(alpha=0.9) +
          ggplot2::labs(
            title="GnRH marker programs",
            x=NULL,
            y="Confidence",
            size="Genes",
            colour="Dimension"
          ) +
          gnrh_theme(
            mode=mode,
            txtsize=txtsize,
            x.ang=x.ang,
            style=style
          )
      )
    }
    # ----------------------------------------------------------------------- #
    # Tile
    # ----------------------------------------------------------------------- #
    return(
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x=.data$program,
          y=.data$confidence_level,
          fill=.data$n_genes
        )
      ) +
        ggplot2::geom_tile(colour="white",linewidth=0.4) +
        ggplot2::geom_text(
          ggplot2::aes(label=.data$n_genes),
          size=max(2.5,min(3.5,txtsize/3.5))
        ) +
        ggplot2::facet_wrap(
          ~program_dimension,
          scales="free_x",
          ncol=facet_ncol
        ) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
          title="GnRH marker-program map",
          x=NULL,
          y="Confidence",
          fill="Genes"
        ) +
        gnrh_theme(
          mode=mode,
          txtsize=txtsize,
          x.ang=x.ang,
          style=style
        )
    )
  }
  # ========================================================================= #
  # Gene-level tables
  # ========================================================================= #
  required <- c("gene","program","confidence_level")
  missing <- setdiff(required,colnames(df))
  if (length(missing)) stop("Missing gene-level columns: ",paste(missing,collapse=", "),call.=FALSE)
  if (isTRUE(require_coexpr) && "coexpr_GNRH1" %in% colnames(df)) {
    df <- df[df$coexpr_GNRH1 %in% TRUE,,drop=FALSE]
  }
  if (!nrow(df)) stop("No genes remain after filtering.",call.=FALSE)
  score_col <- intersect(
    c("integrated_score","marker_evidence_score","marker_score","specificity_score"),
    colnames(df)
  )[1L]
  if (!length(score_col) || is.na(score_col)) {
    df$.score <- 1
    score_col <- ".score"
  }
  df[[score_col]] <- suppressWarnings(as.numeric(df[[score_col]]))
  df[[score_col]][!is.finite(df[[score_col]])] <- 0
  # ========================================================================= #
  # Adaptive selection
  # ========================================================================= #
  has_dataset <- "dataset" %in% colnames(df)
  n_datasets <- if (has_dataset) length(unique(df$dataset)) else length(unique(df$program_dimension))
  if (is.null(top_n)) {
    top_n <- if (n_datasets<=4L) {
      30L
    } else if (n_datasets<=8L) {
      20L
    } else {
      15L
    }
  }
  split_col <- if (has_dataset) "dataset" else "program_dimension"
  split_df <- split(df,df[[split_col]])
  df <- do.call(rbind,lapply(split_df,function(x) {
    x <- x[order(-x[[score_col]],x$gene),,drop=FALSE]
    x <- x[!duplicated(x$gene),,drop=FALSE]
    utils::head(x,top_n)
  }))
  rownames(df) <- NULL
  # ========================================================================= #
  # Adaptive facets / text
  # ========================================================================= #
  n_facets <- length(unique(df[[split_col]]))
  if (is.null(facet_ncol)) {
    facet_ncol <- if (n_facets<=4L) {
      2L
    } else if (n_facets<=9L) {
      3L
    } else {
      4L
    }
  }
  max_genes <- max(table(df[[split_col]]))
  gene_txt <- if (max_genes<=12L) {
    txtsize*0.90
  } else if (max_genes<=20L) {
    txtsize*0.72
  } else {
    txtsize*0.60
  }
  # ========================================================================= #
  # Plot
  # ========================================================================= #
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x=.data$program,
      y=.data$gene,
      colour=.data$confidence_level,
      size=.data[[score_col]]
    )
  ) +
    ggplot2::geom_point(alpha=0.9) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~",split_col)),
      scales="free_y",
      ncol=facet_ncol
    ) +
    ggplot2::labs(
      title="GnRH candidate markers by biological program",
      x="Program",
      y="Gene",
      colour="Confidence",
      size=score_col
    ) +
    gnrh_theme(
      mode=mode,
      txtsize=txtsize,
      x.ang=x.ang,
      style=style
    ) +
    ggplot2::theme(
      axis.text.y=ggplot2::element_text(
        face="italic",
        size=gene_txt
      ),
      strip.text=ggplot2::element_text(
        face="bold",
        size=min(txtsize,gene_txt+1)
      )
    )
  p
}


# ========================================================================= #
# Simple diagnostic plots
# ========================================================================= #

#' Plot GnRH classification counts
#'
#' @param object A Seurat object processed by `run_gnrh()`.
#' @param txtsize Base text size.
#' @return A ggplot object.
#' @export
gnrh_cellcount <- function(
    object,
    txtsize = getOption("gnrhcell.base_size", 14)
) {
  .validate_seurat(object)

  df <- object[[]] |>
    dplyr::count(.data$gnrh_class, name = "n") |>
    dplyr::mutate(pct = 100 * .data$n / sum(.data$n))

  cols <- gnrh_palette("class")

  ggplot2::ggplot(
    df,
    ggplot2::aes(.data$gnrh_class, .data$n, fill = .data$gnrh_class)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(
          "%s\n%.1f%%",
          scales::comma(.data$n),
          .data$pct
        )
      ),
      vjust = -0.25,
      size = 2.7
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = "GnRH classification",
      x = NULL,
      y = "Cells"
    ) +
    gnrh_theme(
      txtsize = txtsize,
      x.ang = 30,
      leg.pos = "none"
    )
}


#' Plot GnRH specificity landscape
#'
#' @param object A Seurat object processed by `run_gnrh()`.
#' @param txtsize Base text size.
#' @return A ggplot object.
#' @export
gnrh_specificity <- function(
    object,
    txtsize = 9
) {
  .validate_seurat(object)

  df <- object[[]]

  support_col <- if (
    "gnrh_support_score_raw" %in%
    colnames(df)
  ) {
    "gnrh_support_score_raw"
  } else {
    "gnrh_support_score"
  }

  required <- c(
    "gnrh_alternative_score",
    support_col,
    "gnrh_status"
  )

  missing <- setdiff(
    required,
    names(df)
  )

  if (length(missing))
    stop(
      "Missing metadata columns: ",
      paste(
        missing,
        collapse = ", "
      ),
      call. = FALSE
    )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x =
        .data$gnrh_alternative_score,

      y =
        .data[[support_col]],

      colour =
        .data$gnrh_status
    )
  ) +
    ggplot2::geom_point(
      size = 0.25,
      alpha = 0.35
    ) +
    ggplot2::scale_colour_manual(
      values =
        gnrh_palette("status")
    ) +
    ggplot2::labs(
      title =
        "GnRH specificity landscape",

      x =
        "Alternative neuronal program score",

      y =
        if (
          support_col ==
          "gnrh_support_score_raw"
        ) {
          "GnRH transcriptomic support"
        } else {
          "GnRH support score"
        },

      colour = "Status"
    ) +
    gnrh_theme(
      txtsize = txtsize
    )
}
