
# =========================================================
# Internal Helpers
# =========================================================
#' `%||%`
#' Return lhs unless NULL
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
#' Automatically determine rasterization
#' @keywords internal
#' @noRd
.auto_raster <- function(ncells, threshold = 1e5) {
  ncells > threshold
}

#' Automatically determine point size
#'
#' @keywords internal
#' @noRd
.auto_pt_size <- function(ncells, raster = FALSE) {

  if (!raster) {
    return(max(0.3, 1.5 / log10(ncells)))
  }

  if (ncells < 1e5) {
    0.1
  } else if (ncells < 2e5) {
    0.5
  } else {
    2
  }
}



#' Validate Seurat object
#' @keywords internal
#' @noRd
.validate_seurat <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
}
#' Get plotting defaults
#' @keywords internal
#' @noRd
.plot_defaults <- function(object, raster = NULL, pt.size = NULL) {
  ncells <- ncol(object)
  raster <- raster %||% .auto_raster(ncells)
  pt.size <- pt.size %||% .auto_pt_size(
    ncells = ncells,
    raster = raster
  )
  list(
    raster = raster,
    pt.size = pt.size,
    ncells = ncells
  )
}


# =========================================================
# Theme System
# =========================================================
#' Unified plotting theme
#'
#' Provides a consistent ggplot2 theme for GnRHcell visualizations.
#'
#' @param theme Theme preset. One of `"classic"`, `"minimal"`, or `"dark"`.
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param legend.position Position of legend.
#' @param x.angle Rotation angle for x-axis labels.
#' @param axis.text Logical; whether to display axis text.
#' @param axis.title Logical; whether to display axis titles.
#' @param grid Logical; whether to display major grid lines.
#'
#' @return A ggplot2 theme object.
#'
#' @export
celltheme <- function(
    theme = c("classic", "minimal", "dark"),
    base_size = 12,
    base_family = "Helvetica",
    legend.position = "right",
    x.angle = 0,
    axis.text = TRUE,
    axis.title = TRUE,
    grid = FALSE
) {
  theme <- match.arg(theme)
  base_theme <- switch(
    theme,
    classic = ggplot2::theme_classic(base_size = base_size),
    minimal = ggplot2::theme_minimal(base_size = base_size),
    dark = ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = "#1E1E1E"),
        panel.background = ggplot2::element_rect(fill = "#1E1E1E"),
        text             = ggplot2::element_text(color = "white"),
        axis.text        = ggplot2::element_text(color = "grey80"),
        axis.title       = ggplot2::element_text(color = "grey90"),
        legend.background = ggplot2::element_rect(fill = "#1E1E1E"),
        legend.key        = ggplot2::element_rect(fill = "#1E1E1E"),
        legend.text       = ggplot2::element_text(color = "white"),
        legend.title      = ggplot2::element_text(color = "white")
      )
  )
  base_theme +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family),
      axis.text.x = if (axis.text) {
        ggplot2::element_text(angle = x.angle, hjust = 1)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = if (axis.text) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      axis.title = if (axis.title) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      legend.position = legend.position,
      panel.grid.major = if (grid) {
        ggplot2::element_line(color = "grey85", linewidth = 0.2)
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )
}



# =========================================================
# Color Utilities
# =========================================================
#' Generate plotting palette
#'
#' @param n Number of colors.
#' @param preset Palette preset.
#' @param reverse Logical; whether to reverse color order.
#'
#' @return A character vector of colors.
#'
#' @export
cellpal <- function(
    n,
    preset = c("base", "bright", "pastel"),
    reverse = FALSE
) {
  preset <- match.arg(preset)
  palettes <- list(
    base = c(
      "#E41A1C", "#377EB8", "#4DAF4A",
      "#984EA3", "#FF7F00", "#FFFF33"
    ),
    bright = c(
      "#E6194B", "#3CB44B", "#FFE119",
      "#4363D8", "#F58231", "#911EB4"
    ),
    pastel = c(
      "#FFB3BA", "#BAE1FF", "#BAFFC9",
      "#FFFFBA", "#FFD6A5"
    )
  )
  cols <- grDevices::colorRampPalette(
    palettes[[preset]]
  )(n)
  if (reverse) {
    cols <- rev(cols)
  }
  cols
}



#' GnRH color palettes
#'
#' Returns predefined color palettes used throughout GnRHcell plots.
#'
#' @param type Palette type. One of `"status"`, `"class"`, or `"stage"`.
#'
#' @return A named character vector of colors.
#'
#' @export
gnrh_colors <- function(type = c("status", "class", "stage")) {
  type <- match.arg(type)
  switch(
    type,
    status = c(
      pos = "#E64B35",
      neg = "#D9D9D9"
    ),
    class = c(
      gnrh_neg  = "#D9D9D9",
      gnrh_low  = "#FDB863",
      gnrh_high = "#B2182B"
    ),
    stage = c(
      other         = "#D9D9D9",
      progenitor    = "#8dd3c7",
      migrating     = "#80b1d3",
      transitional  = "#fdb462",
      mature        = "#fb8072",
      secreting     = "#e31a1c"
    )
  )
}



# =========================================================
# Shared Scatter Plot Builder
# =========================================================
#' Internal scatter helper
#' @keywords internal
#' @noRd
.scatter_plot <- function(
    data,
    x,
    y,
    color = NULL,
    cols = NULL,
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    pt.size = 1,
    alpha = 0.8,
    theme = "classic"
) {

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  )
  if (!is.null(color)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[color]]),
        size = pt.size,
        alpha = alpha
      )
    if (!is.null(cols)) {
      p <- p +
        ggplot2::scale_color_manual(values = cols)
    }
  } else {
    p <- p +
      ggplot2::geom_point(
        size = pt.size,
        alpha = alpha
      )
  }
  p +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    celltheme(theme = theme)
}


# =========================================================
# DimPlot
# =========================================================
#' Enhanced Seurat DimPlot
#' Dimensional reduction plot for GnRH analysis
#'
#' Creates a scatter plot of cells in a reduced dimensional space
#' (e.g. UMAP or t-SNE), with flexible grouping and annotation options.
#'
#' This function is designed for Seurat objects and supports
#' standard dimensional reductions stored in the object.
#'
#' @param object A Seurat object containing dimensional reduction data.
#' @param group.by Metadata column used to color cells.
#' @param reduction Name of dimensional reduction to use (e.g. "umap", "tsne").
#' @param dims Numeric vector of dimensions to use (default: c(1, 2)).
#' @param cols Vector of colors for groups.
#' @param label Logical; whether to add cluster labels.
#' @param label.size Size of text labels.
#' @param show.n Logical; whether to display group sizes.
#' @param pt.size Point size for cells.
#' @param alpha Transparency of points.
#' @param raster Logical; whether to rasterize points for large datasets.
#' @param legend Logical; whether to show legend.
#' @param base_size Base font size for theme.
#' @param theme Optional ggplot theme override.
#' @param label.face Font face used for labels.
#' @param dark Logical; whether to use dark mode styling.
#' @param shuffle Logical; whether to randomly shuffle plotting order.
#'
#' @return A ggplot object.
#'
#' @export
dim_plot <- function(
    object,
    group.by = "ident",
    reduction = "umap",
    dims = c(1, 2),
    cols = NULL,
    label = FALSE,
    label.size = 4,
    label.face = "plain",
    show.n = TRUE,
    dark = FALSE,
    raster = NULL,
    pt.size = NULL,
    alpha = 1,
    shuffle = FALSE,
    legend = TRUE,
    base_size = 12,
    theme = c("classic", "minimal", "dark")
) {
  .validate_seurat(object)
  theme <- match.arg(theme)
  Seurat::Idents(object) <- group.by
  Seurat::Idents(object) <- droplevels(Seurat::Idents(object))
  groups <- levels(Seurat::Idents(object))
  if (is.null(cols)) {
    cols <- cellpal(length(groups))
    names(cols) <- groups
  }
  defaults <- .plot_defaults(
    object,
    raster = raster,
    pt.size = pt.size
  )
  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    dims = dims,
    group.by = group.by,
    cols = cols,
    raster = defaults$raster,
    pt.size = defaults$pt.size,
    alpha = alpha,
    shuffle = shuffle,
    label = FALSE
  )
  if (show.n) {
    counts <- table(Seurat::Idents(object))
    labels <- paste0(
      names(counts),
      " (",
      counts,
      ")"
    )
    p <- p +
      ggplot2::scale_color_manual(
        values = cols,
        labels = labels
      )
  }
  if (label) {
    emb <- Seurat::Embeddings(object, reduction)[, dims]
    centers <- aggregate(
      emb,
      by = list(cluster = Seurat::Idents(object)),
      FUN = median
    )
    p <- p +
      ggrepel::geom_text_repel(
        data = centers,
        ggplot2::aes(
          x = emb[,1],
          y = emb[,2],
          label = cluster
        ),
        inherit.aes = FALSE,
        size = label.size,
        fontface = label.face,
        seed = 42
      )
  }
  p <- p +
    celltheme(
      theme = if (dark) "dark" else theme,
      base_size = base_size
    )
  if (!legend) {
    p <- p + ggplot2::guides(color = "none")
  }
  p
}


# =========================================================
# Plot Embedding
# =========================================================
#' DimPlot Wrapper
#' @param obj Seurat object.
#' @param group.by Metadata column used for grouping.
#' @param reduction Dimensional reduction to use.
#' @param cols Optional colors.
#' @param label Logical; draw labels.
#' @param label.face Font face for labels.
#' @param dark Logical; use dark theme.
#' @export
plot_gnrh_embedding <- function(obj,
                                group.by = "all",
                                reduction = "umap",
                                cols = NULL,
                                label = FALSE,
                                label.face = "plain",
                                dark = FALSE) {

  allowed <- c("gnrh_status", "gnrh_class", "gnrh_stage")

  if (length(group.by) == 1 && group.by == "all") {
    group.by <- allowed
  }

  bad <- setdiff(group.by, allowed)
  if (length(bad) > 0) stop("Unsupported group.by")

  cmap <- list(
    gnrh_status = gnrh_colors("status"),
    gnrh_class  = gnrh_colors("class"),
    gnrh_stage  = gnrh_colors("stage")
  )

  if (length(group.by) == 1) {

    if (is.null(cols)) cols <- cmap[[group.by]]

    return(
      dim_plot(
        obj,
        reduction = reduction,
        group.by = group.by,
        label = label,
        label.face = label.face,
        cols = cols,
        theme = "classic",
        dark = dark
      )
    )
  }

  plist <- lapply(group.by, function(g) {
    dim_plot(
      obj,
      reduction = reduction,
      group.by = g,
      label = label,
      label.face = label.face,
      cols = cmap[[g]],
      theme = "classic",
      dark = dark
    )
  })

  patchwork::wrap_plots(plist)
}

# =========================================================
# Feature Plot
# =========================================================
#' Enhanced feature plot
#'
#' Wrapper around Seurat::FeaturePlot with GnRHcell styling.
#'
#' @param object A Seurat object.
#' @param features Features to visualize.
#' @param reduction Dimensional reduction to use.
#' @param cols Colors for low/high expression.
#' @param raster Logical; whether to rasterize points.
#' @param pt.size Point size.
#' @param order Logical; whether to plot high-expression cells on top.
#' @param dark Logical; whether to use dark theme.
#' @param base_size Base font size.
#' @param theme Theme preset.
#' @param ncol Number of columns for combined plots.
#'
#' @return A ggplot object or patchwork object.
#'
#' @export
plot_gnrh_feature <- function(
    object,
    features,
    reduction = "umap",
    cols = c("lightgrey", "red"),
    raster = NULL,
    pt.size = NULL,
    order = TRUE,
    dark = FALSE,
    base_size = 12,
    theme = c("classic", "minimal", "dark"),
    ncol = NULL
) {
  .validate_seurat(object)
  theme <- match.arg(theme)
  defaults <- .plot_defaults(
    object,
    raster = raster,
    pt.size = pt.size
  )
  Seurat::FeaturePlot(
    object = object,
    features = features,
    reduction = reduction,
    cols = cols,
    raster = defaults$raster,
    pt.size = defaults$pt.size,
    order = order,
    combine = TRUE,
    ncol = ncol
  ) &
    celltheme(
      theme = if (dark) "dark" else theme,
      base_size = base_size
    )
}



# =========================================================
# Violin Plot
# =========================================================
#' Enhanced violin plot
#'
#' Wrapper around Seurat::VlnPlot with GnRHcell styling.
#'
#' @param object A Seurat object.
#' @param features Features to visualize.
#' @param group.by Metadata column for grouping.
#' @param cols Colors for groups.
#' @param pt.size Point size for overlaid cells.
#' @param stack Logical; whether to stack violins.
#' @param flip Logical; whether to flip coordinates.
#' @param show.median Logical; whether to display median markers.
#' @param base_size Base font size.
#' @param theme Theme preset.
#' @param ncol Number of columns for combined plots.
#'
#' @return A ggplot object.
#'
#' @export
plot_gnrh_violin <- function(
    object,
    features,
    group.by = "seurat_clusters",
    cols = NULL,
    pt.size = 0,
    stack = FALSE,
    flip = FALSE,
    show.median = FALSE,
    base_size = 12,
    theme = c("classic", "minimal", "dark"),
    ncol = NULL
) {
  .validate_seurat(object)
  theme <- match.arg(theme)
  Seurat::Idents(object) <- group.by
  groups <- levels(Seurat::Idents(object))
  if (is.null(cols)) {
    cols <- cellpal(length(groups))
    names(cols) <- groups
  }
  p <- Seurat::VlnPlot(
    object = object,
    features = features,
    group.by = group.by,
    cols = cols,
    pt.size = pt.size,
    stack = stack,
    combine = TRUE,
    ncol = ncol
  )
  if (show.median) {
    p <- p +
      ggplot2::stat_summary(
        fun = median,
        geom = "point",
        shape = 95,
        size = 4
      )
  }
  if (flip) {
    p <- p + ggplot2::coord_flip()
  }
  p &
    Seurat::NoLegend() &
    celltheme(
      theme = theme,
      base_size = base_size,
      x.angle = 45
    )
}



# =========================================================
# ROC Plot
# =========================================================
#' Plot ROC curve
#'
#' Computes and visualizes ROC performance for GnRH classification.
#'
#' @param object A Seurat object containing `gnrh_truth_proxy`
#'   and `gnrh_expr` metadata fields.
#'
#' @return A ggplot object.
#'
#' @export
plot_gnrh_roc <- function(object) {
  truth <- object$gnrh_truth_proxy
  expr  <- object$gnrh_expr
  if (is.null(truth)) {
    stop("Missing truth proxy.")
  }
  roc_obj <- pROC::roc(truth, expr)
  roc_df <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities
  )
  auc_val <- round(
    pROC::auc(roc_obj),
    3
  )
  ggplot2::ggplot(
    roc_df,
    ggplot2::aes(fpr, tpr)
  ) +
    ggplot2::geom_line(
      color = "red",
      linewidth = 1
    ) +
    ggplot2::geom_abline(
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste0("ROC Curve (AUC = ", auc_val, ")"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    celltheme()
}



# =========================================================
# Threshold Diagnostics
# =========================================================
#' Plot threshold diagnostics
#'
#' Visualizes GnRH expression and composite scores relative
#' to threshold-based selection.
#'
#' @param object A Seurat object containing threshold diagnostics
#'   in `object@misc$gnrh_diag`.
#'
#' @return A ggplot object.
#'
#' @export
plot_gnrh_threshold <- function(object) {
  .validate_seurat(object)
  df <- object@misc$gnrh_diag
  if (is.null(df)) {
    stop("Missing gnrh_diag metadata.")
  }
  selected_pct <- round(
    mean(df$above_threshold) * 100,
    2
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(expr, score)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = class),
      alpha = 0.8,
      size = 1.5
    ) +
    ggplot2::geom_point(
      data = subset(df, above_threshold),
      color = "red",
      size = 2
    ) +
    ggplot2::scale_color_manual(
      values = gnrh_colors("class")
    ) +
    ggplot2::labs(
      title = paste0(
        "Threshold Diagnostics (",
        selected_pct,
        "% selected)"
      ),
      x = "GnRH expression",
      y = "Composite score"
    ) +
    celltheme()
}


# =========================================================
# Genes Network
# =========================================================

#' GnRH-centered coexpression star network
#'
#' Builds a gene coexpression network centered on a target gene
#' and visualizes top correlated partners.
#'
#' @param df Data frame containing gene, coexpr, and score columns.
#' @param target Character. Central gene (default: "GNRH1").
#' @param top_n Number of top edges to include.
#'
#' @return A ggraph object.
#'
#' @importFrom igraph graph_from_data_frame V
#' @import ggraph
#' @export
plot_star_network <- function(df, target = "GNRH1", top_n = 30) {

  # filter rows (base R)
  df2 <- df[df$gene != target, , drop = FALSE]

  # compute weight
  df2$weight <- df2$coexpr * df2$score

  # order and select top N
  df2 <- df2[order(-df2$weight), , drop = FALSE]
  df2 <- df2[seq_len(min(top_n, nrow(df2))), , drop = FALSE]

  # edge list
  edge_df <- data.frame(
    from = target,
    to = df2$gene,
    weight = df2$weight,
    stringsAsFactors = FALSE
  )

  # nodes MUST be passed into graph
  nodes <- data.frame(
    name = unique(c(target, df2$gene)),
    type = ifelse(unique(c(target, df2$gene)) == target,
                  "Target gene",
                  "Coexpressed genes"),
    stringsAsFactors = FALSE
  )

  # graph WITH vertices
  g <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = nodes
  )

  set.seed(123)

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      colour = "grey70",
      alpha = 0.5
    ) +
    ggraph::geom_node_point(
      aes(color = type, size = type)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      size = 3.5,
      fontface = "italic"
    ) +
    ggplot2::theme_void()
}



#' GnRH gene similarity network
#'
#' Constructs a similarity network based on coexpression structure
#' among top-ranked genes.
#'
#' @param df Data frame containing gene-level scores and coexpression.
#' @param top_n Number of genes to include.
#' @param threshold Minimum edge weight to retain.
#'
#' @return A ggraph object.
#'
#' @importFrom igraph graph_from_data_frame cluster_louvain V
#' @import ggraph
#' @export
plot_similarity_network <- function(df, top_n = 25, threshold = 0.4) {

  df <- df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
  genes <- df$gene

  mat <- outer(genes, genes, Vectorize(function(g1, g2) {
    a <- df$coexpr[df$gene == g1][1]
    b <- df$coexpr[df$gene == g2][1]
    min(a, b, na.rm = TRUE)
  }))

  rownames(mat) <- genes
  colnames(mat) <- genes

  edges <- as.data.frame(
    as.table(mat),
    stringsAsFactors = FALSE
  )

  colnames(edges) <- c("from", "to", "weight")

  edges <- edges[
    edges$from != edges$to & edges$weight > threshold,
    ,
    drop = FALSE
  ]

  nodes <- unique(data.frame(
    name = df$gene,
    score = df$score,
    stringsAsFactors = FALSE
  ))

  g <- igraph::graph_from_data_frame(
    edges,
    vertices = nodes,
    directed = FALSE
  )

  comm <- igraph::cluster_louvain(g)
  igraph::V(g)$community <- factor(comm$membership)

  set.seed(123)

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      colour = "grey75"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = community, size = score)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      max.overlaps = 30
    ) +
    ggplot2::theme_void()
}


