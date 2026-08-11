# ============================================================================= #
# GnRHcell visualization utilities
# ============================================================================= #


# ============================================================================= #
# Internal utilities
# ============================================================================= #

#' Return left-hand side unless NULL
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Validate a Seurat object
#' @keywords internal
#' @noRd
.validate_seurat <- function(object) {
  if (!inherits(object, "Seurat"))
    stop("`object` must be a Seurat object.", call. = FALSE)
  invisible(TRUE)
}

# Backward-compatible internal alias used by older plotting wrappers.
.check_seurat <- .validate_seurat


#' Automatically determine rasterization
#' @keywords internal
#' @noRd
.auto_raster <- function(ncells, threshold = 1e5)
  ncells > threshold


#' Automatically determine point size
#' @keywords internal
#' @noRd
.auto_pt_size <- function(ncells, raster = FALSE) {
  if (raster) {
    if (ncells < 1e5) return(0.15)
    if (ncells < 2e5) return(0.35)
    return(0.60)
  }

  if (ncells < 2500L) return(0.80)
  if (ncells < 10000L) return(0.50)
  if (ncells < 50000L) return(0.30)
  0.20
}


#' Internal plotting defaults
#' @keywords internal
#' @noRd
.plot_defaults <- function(object, raster = NULL, pt.size = NULL) {
  ncells <- ncol(object)
  raster <- raster %||% .auto_raster(ncells)
  pt.size <- pt.size %||% .auto_pt_size(ncells, raster)
  list(raster = raster, pt.size = pt.size)
}


# ============================================================================= #
# Theme
# ============================================================================= #

#' GnRHcell internal theme
#'
#' @keywords internal
#' @noRd
.gnrh_theme <- function(
    txtsize = 10,
    x.ang = 0,
    leg.pos = "right",
    title.position = c("left", "center", "right"),
    axes = TRUE,
    style = c("classic", "minimal", "bw", "test"),
    mode = c("light", "dark")
) {
  title.position <- match.arg(title.position)
  style <- match.arg(style)
  mode <- match.arg(mode)

  hjust <- switch(title.position, left = 0, center = 0.5, right = 1)

  base <- switch(
    style,
    classic = ggplot2::theme_classic(base_size = txtsize),
    minimal = ggplot2::theme_minimal(base_size = txtsize),
    bw = ggplot2::theme_bw(base_size = txtsize),
    test = ggplot2::theme_classic(base_size = txtsize)
  )

  fg <- if (mode == "dark") "white" else "#1A1A1A"
  bg <- if (mode == "dark") "#111111" else "white"

  base +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", colour = fg),
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),

      plot.title = ggplot2::element_text(
        face = "bold", size = txtsize + 2, hjust = hjust
      ),
      plot.subtitle = ggplot2::element_text(
        size = max(7, txtsize - 1), colour = if (mode == "dark") "grey80" else "grey35"
      ),

      axis.title = ggplot2::element_text(size = txtsize, colour = fg),
      axis.text.x = if (axes)
        ggplot2::element_text(
          size = txtsize - 1, angle = x.ang,
          hjust = if (x.ang == 0) 0.5 else 1,
          colour = fg
        )
      else ggplot2::element_blank(),
      axis.text.y = if (axes)
        ggplot2::element_text(size = txtsize - 1, colour = fg)
      else ggplot2::element_blank(),

      axis.ticks = if (axes)
        ggplot2::element_line(linewidth = 0.3, colour = fg)
      else ggplot2::element_blank(),
      axis.line = if (axes)
        ggplot2::element_line(linewidth = 0.3, colour = fg)
      else ggplot2::element_blank(),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.key.width = grid::unit(0.4, "cm"),
      legend.position = leg.pos,
      legend.title = ggplot2::element_text(face = "bold", size = txtsize),
      legend.text = ggplot2::element_text(size = txtsize - 1),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),

      panel.grid = if (style == "minimal")
        ggplot2::element_line(linewidth = 0.2, colour = if (mode == "dark") "grey30" else "grey90")
      else ggplot2::element_blank()
    )
}

# Backward-compatible internal wrapper used by feature and dot plots. Legacy
# arguments are accepted through `...`; visual labels remain controlled by the
# calling plot.
plot_theme <- function(
    style = "classic",
    txtsize = 10,
    x.ang = 0,
    leg.pos = "right",
    ...
) {
  .gnrh_theme(
    txtsize = txtsize,
    x.ang = x.ang,
    leg.pos = leg.pos,
    style = style
  )
}


# ============================================================================= #
# Palettes
# ============================================================================= #

#' GnRHcell color palettes
#'
#' @param type Palette type.
#'
#' @return Named character vector.
#' @export
gnrh_colors <- function(
    type = c("status", "confident", "class", "stage", "secretory")
) {
  switch(
    match.arg(type),

    status = c(
      neg = "#B0B0B0",
      pos = "#FF4D6D"
    ),

    confident = c(
      "FALSE" = "#B0B0B0",
      "TRUE" = "#E63946"
    ),

    class = c(
      neg = "#B0B0B0",
      dropout_rescue = "#F4A261",
      supported = "#4EA8DE",
      direct = "#E63946"
    ),

    stage = c(
      "non-gnrh" = "#B0B0B0",
      identity = "#1F78B4",
      migrating = "#6A3D9A",
      mature = "#E67E22",
      secreting = "#D81B60"
    ),

    secretory = c(
      "non-gnrh" = "#B0B0B0",
      limited = "#80B1D3",
      supported = "#D81B60"
    )
  )
}


#' Internal GnRH palette resolver
#' @keywords internal
#' @noRd
.gnrh_palette <- function(group_by) {
  type <- switch(
    group_by,
    gnrh_status = "status",
    gnrh_confident = "confident",
    gnrh_class = "class",
    gnrh_stage = "stage",
    gnrh_secretory = "secretory",
    NULL
  )

  if (is.null(type)) return(NULL)
  gnrh_colors(type)
}


#' Resolve colors for arbitrary levels
#' @keywords internal
#' @noRd
.resolve_colors <- function(levels, cols = NULL, group_by = NULL) {
  if (is.null(cols) && !is.null(group_by))
    cols <- .gnrh_palette(group_by)

  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(levels))
    names(cols) <- levels
  }

  if (is.null(names(cols)))
    names(cols) <- levels[seq_len(min(length(cols), length(levels)))]

  missing <- setdiff(levels, names(cols))
  if (length(missing)) {
    extra <- scales::hue_pal()(length(missing))
    names(extra) <- missing
    cols <- c(cols, extra)
  }

  cols[levels]
}


# ============================================================================= #
# Legend helpers
# ============================================================================= #

#' Create GnRHcell legend labels
#' @keywords internal
#' @noRd
.gnrh_legend_labels <- function(
    x,
    levels = NULL,
    percentage = FALSE,
    n_cells = TRUE
) {
  x <- as.character(x)
  x[is.na(x)] <- "Unknown"

  tab <- table(x)
  levels <- levels %||% names(tab)
  levels <- intersect(levels, names(tab))
  n <- as.integer(tab[levels])

  if (!n_cells)
    return(stats::setNames(levels, levels))

  labels <- if (percentage) {
    paste0(
      levels, " (",
      format(n, big.mark = ","),
      " | ",
      round(100 * n / sum(tab), 1),
      "%)"
    )
  } else {
    paste0(
      levels, " (",
      format(n, big.mark = ","),
      ")"
    )
  }

  stats::setNames(labels, levels)
}


#' Compact point legend
#' @keywords internal
#' @noRd
.point_legend <- function(size = 2.3) {
  list(
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = size, alpha = 1),
        keyheight = grid::unit(0.28, "cm")
      )
    ),
    ggplot2::theme(
      legend.spacing.y = grid::unit(0, "cm"),
      legend.key.width = grid::unit(0.35, "cm")
    )
  )
}


#' Compact line legend
#' @keywords internal
#' @noRd
.line_legend <- function(linewidth = 1.2) {
  list(
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(linewidth = linewidth),
        keyheight = grid::unit(0.30, "cm")
      )
    ),
    ggplot2::theme(
      legend.spacing.y = grid::unit(0, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )
  )
}


# ============================================================================= #
# Reduction resolver
# ============================================================================= #

#' Resolve dimensional reduction
#'
#' @keywords internal
#' @noRd
resolve_reduction <- function(
    object,
    reduction = NULL,
    fallback = c("umap", "umap_scvi", "umap.harmony", "umap.rpca", "pca")
) {
  .validate_seurat(object)

  available <- names(object@reductions)

  if (!length(available))
    stop("No dimensional reductions are available in `object`.", call. = FALSE)

  if (
    !is.null(reduction) &&
    length(reduction) == 1L &&
    !is.na(reduction) &&
    nzchar(reduction)
  ) {
    if (reduction %in% available)
      return(reduction)

    warning(
      sprintf(
        "Reduction `%s` not found. Using an available fallback.",
        reduction
      ),
      call. = FALSE
    )
  }

  candidate <- fallback[fallback %in% available]
  if (length(candidate)) return(candidate[[1L]])

  candidate <- available[grepl("umap", available, ignore.case = TRUE)]
  if (length(candidate)) return(candidate[[1L]])

  available[[1L]]
}


# ============================================================================= #
# Embedding plot
# ============================================================================= #

#' Plot GnRH embedding
#'
#' @param object A Seurat object.
#' @param group_by Metadata column used to colour cells.
#' @param reduction Dimensional reduction; an available reduction is selected when `NULL`.
#' @param dims Two reduction dimensions to display.
#' @param shuffle Randomize plotting order.
#' @param raster Use rasterized points; selected automatically when `NULL`.
#' @param raster.dpi Raster resolution passed to the raster geom.
#' @param alpha,background_alpha Opacity for highlighted and background cells.
#' @param n.cells,percentage Add cell counts or percentages to legend labels.
#' @param label,repel,label.size,label.face Cluster-label controls.
#' @param cols Optional named colour vector.
#' @param axes Show embedding axes.
#' @param plot.ttl Optional plot title.
#' @param legend Show the legend.
#' @param leg.ttl,leg.ttl.size Legend title and title size.
#' @param item.size,item.border Legend-key controls.
#' @param leg.pos,leg.dir,leg.size,leg.ncol Legend layout controls.
#' @param txtsize Base text size.
#' @param pt.size Point size; selected automatically when `NULL`.
#' @param dark Use dark display mode.
#' @param total.cells Include total-cell counts in the legend.
#' @param style Theme style.
#' @param ... Additional graphical arguments.
#' @return A ggplot object.
#' @export
plot_gnrh_embedding <- function(
    object,
    group_by = "gnrh_status",
    reduction = NULL,
    dims = c(1, 2),
    shuffle = FALSE,
    raster = NULL,
    raster.dpi = c(2048, 2048),
    alpha = 0.9,
    background_alpha = 0.18,
    n.cells = TRUE,
    percentage = FALSE,
    label = FALSE,
    repel = TRUE,
    label.size = 4,
    label.face = "plain",
    cols = NULL,
    axes = TRUE,
    plot.ttl = NULL,
    legend = TRUE,
    leg.ttl = NULL,
    leg.ttl.size = NULL,
    item.size = 3.5,
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = NULL,
    leg.ncol = NULL,
    item.border = TRUE,
    txtsize = 12,
    pt.size = NULL,
    dark = FALSE,
    total.cells = FALSE,
    style = "classic",
    ...
) {
  .validate_seurat(object)

  md <- object[[]]
  if (!group_by %in% names(md))
    stop("Metadata column `", group_by, "` not found.", call. = FALSE)

  if (length(dims) != 2L || any(!is.finite(dims)))
    stop("`dims` must contain exactly two valid dimensions.", call. = FALSE)

  reduction <- resolve_reduction(object, reduction)
  emb <- Seurat::Embeddings(object, reduction)

  if (max(dims) > ncol(emb))
    stop("Selected dimensions exceed those available in `", reduction, "`.", call. = FALSE)

  defaults <- .plot_defaults(object, raster, pt.size)
  raster <- defaults$raster
  pt.size <- defaults$pt.size

  orders <- list(
    gnrh_status = c("neg", "pos"),
    gnrh_class = c("neg", "dropout_rescue", "supported", "direct"),
    gnrh_confident = c("FALSE", "TRUE"),
    gnrh_stage = c("non-gnrh", "identity", "migrating", "mature", "secreting"),
    gnrh_secretory = c("non-gnrh", "limited", "supported")
  )

  values <- as.character(md[[group_by]])
  values[is.na(values)] <- "Unknown"

  lev <- orders[[group_by]] %||% sort(unique(values))
  lev <- c(lev[lev %in% values], sort(setdiff(unique(values), lev)))

  cols <- .resolve_colors(lev, cols, group_by)
  if ("Unknown" %in% names(cols)) cols["Unknown"] <- "grey70"

  df <- data.frame(
    x = emb[, dims[1]],
    y = emb[, dims[2]],
    group = factor(values, levels = lev)
  )

  background <- c("neg", "non-gnrh", "FALSE")
  df$.alpha <- ifelse(as.character(df$group) %in% background, background_alpha, alpha)

  if (shuffle) {
    set.seed(42)
    df <- df[sample.int(nrow(df)), , drop = FALSE]
  } else {
    df <- df[order(df$.alpha), , drop = FALSE]
  }

  present <- lev[lev %in% as.character(unique(df$group))]
  labels <- .gnrh_legend_labels(
    df$group,
    levels = present,
    percentage = percentage,
    n_cells = n.cells
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x, y, colour = group, alpha = .alpha)
  )

  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    p <- p + ggrastr::geom_point_rast(
      size = pt.size,
      raster.dpi = raster.dpi[1]
    )
  } else {
    p <- p + ggplot2::geom_point(size = pt.size)
  }

  p <- p +
    ggplot2::scale_colour_manual(
      values = cols[present],
      breaks = present,
      labels = labels,
      drop = FALSE
    ) +
    ggplot2::scale_alpha_identity()

  if (label) {
    lab <- stats::aggregate(cbind(x, y) ~ group, df, stats::median)

    if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab,
        ggplot2::aes(x, y, label = group),
        inherit.aes = FALSE,
        fontface = label.face,
        size = label.size,
        seed = 42
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x, y, label = group),
        inherit.aes = FALSE,
        fontface = label.face,
        size = label.size
      )
    }
  }

  if (total.cells)
    plot.ttl <- paste0(
      plot.ttl %||% group_by,
      " (n = ", format(ncol(object), big.mark = ","), ")"
    )

  leg.ttl <- leg.ttl %||% group_by
  leg.ttl.size <- leg.ttl.size %||% txtsize
  leg.size <- leg.size %||% max(8, txtsize - 2)
  leg.ncol <- leg.ncol %||% ifelse(length(present) > 30, 2, 1)

  p <- p +
    ggplot2::labs(
      title = plot.ttl,
      x = paste0(toupper(reduction), "_", dims[1]),
      y = paste0(toupper(reduction), "_", dims[2]),
      colour = leg.ttl
    ) +
    ggplot2::coord_equal() +
    .gnrh_theme(
      txtsize = txtsize,
      leg.pos = if (legend) leg.pos else "none",
      axes = axes,
      style = style,
      mode = if (dark) "dark" else "light"
    ) +
    ggplot2::theme(
      legend.direction = leg.dir,
      legend.title = ggplot2::element_text(
        size = leg.ttl.size,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(size = leg.size)
    )

  if (legend) {
    override <- if (item.border) {
      list(
        size = item.size,
        shape = 21,
        colour = if (dark) "white" else "black",
        fill = unname(cols[present]),
        stroke = 0.3,
        alpha = 1
      )
    } else {
      list(size = item.size, alpha = 1)
    }

    p <- p +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = override,
          ncol = leg.ncol,
          keyheight = grid::unit(0.38, "cm"),
          keywidth = grid::unit(0.38, "cm")
        )
      )
  }

  p
}


# ============================================================================= #
# Feature plots
# ============================================================================= #

#' Resolve GnRHcell feature presets
#' @keywords internal
#' @noRd
.gnrh_features <- function(
    type = c("core", "modules", "staging", "all")
) {
  type <- match.arg(type)

  presets <- list(
    core = c(
      "gnrh_raw",
      "gnrh_expr",
      "gnrh_score",
      "gnrh_knn"
    ),

    modules = c(
      "gnrh_core_hits",
      "gnrh_mig_hits",
      "gnrh_neuro_hits",
      "gnrh_knn"
    ),

    staging = c(
      "gnrh_identity_score",
      "gnrh_migrating_score",
      "gnrh_mature_score",
      "gnrh_secreting_score"
    )
  )

  if (type == "all") {
    return(unique(unlist(presets, use.names = FALSE)))
  }

  presets[[type]]
}


#' Cell Feature Plot for Seurat Objects
#'
#' A flexible wrapper around Seurat's `FeaturePlot` to visualize gene expression or metadata
#' features in a Seurat object. Supports custom color palettes, viridis, and hotspot/rainbow palettes.
#'
#' @param object A `Seurat` object.
#' @param features Character vector of features (genes or metadata columns) to plot.
#' @param preset Optional GnRHcell feature preset: `"core"`, `"modules"`,
#'   `"staging"`, or `"all"`.
#' @param cols Optional character vector of colors for plotting.
#' @param theme.cols Character. Predefined theme color palette (default: `"Reds"`). Options include `"Reds"`, `"Blues"`, etc., or custom list palettes `"hotspot"`,`"rainbow"`, etc.
#' @param rev.cols Logical. Reverse the color palette (default: `FALSE`).
#' @param na.col Color for NA or below-cutoff expression values (default: `"lightgray"`).
#' @param order Logical. If TRUE, plot high-expression cells on top (default: `FALSE`).
#' @param pt.size Numeric. Point size. If NULL, automatically calculated based on number of cells.
#' @param txtsize Numeric. Base font size for plot titles and axis labels (default: 10).
#' @param reduction Character. Dimensional reduction to use (default: first available in Seurat object).
#' @param na.cutoff Numeric. Minimum expression value for coloring; below this will be NA if palette requires (default: 1e-9).
#' @param raster Logical. If TRUE, rasterize points for faster plotting of large datasets.
#' @param raster.dpi Numeric vector of length 2. DPI for rasterization (default: c(512,512)).
#' @param split.by Character. Metadata column to split the plot.
#' @param ncol Numeric. Number of columns when combining multiple plots.
#' @param layer Character. Seurat assay slot to fetch data from (default: `"data"`).
#' @param label Logical. Whether to label clusters (default: FALSE).
#' @param axes Logical. Whether to show axes (default: TRUE).
#' @param combine Logical. Whether to return a single combined plot (default: TRUE).
#' @param blend Logical. Whether to blend exactly two features (default: FALSE).
#' @param merge.leg Logical. Whether to merge multiple legends into one (default: FALSE).
#' @param style Character. ggplot2 theme to apply (default: `"classic"`).
#' @param ... Additional arguments passed to `Seurat::FeaturePlot`.
#'
#' @return A `ggplot` object (or `patchwork` object if multiple features).
#' @export
#'
#' @examples
#' \dontrun{
#' path <- system.file("extdata", "hpsc.rds", package = "GnRHcell")
#' object <- readRDS(path)
#' plot_gnrh_feature(object, features = "GNRH1")
#' plot_gnrh_feature(object, preset = "core", ncol = 2)
#' }
#'
plot_gnrh_feature <- function(
    object,
    features = NULL,
    preset = NULL,
    cols = NULL,
    theme.cols = "gnrh",
    rev.cols = FALSE,
    na.col = "lightgray",
    order = FALSE,
    pt.size = NULL,
    txtsize = 10,
    reduction = NULL,
    na.cutoff = 1e-9,
    raster = NULL,
    raster.dpi = c(512, 512),
    split.by = NULL,
    ncol = NULL,
    layer = "data",
    label = FALSE,
    axes = TRUE,
    combine = TRUE,
    blend = FALSE,
    merge.leg = FALSE,
    style = "classic",
    ...
) {

  suppressWarnings(suppressMessages({

    # Checks
    if (!inherits(object, "Seurat"))
      stop("object must be a Seurat object")

    if (!is.null(features) && !is.null(preset)) {
      stop(
        "Supply either `features` or `preset`, not both.",
        call. = FALSE
      )
    }

    if (!is.null(preset)) {
      features <- .gnrh_features(preset)
    }

    if (is.null(features) || !length(features)) {
      stop(
        "Supply `features` or select a GnRH preset with `preset`.",
        call. = FALSE
      )
    }

    available <- c(
      rownames(object),
      colnames(object[[]])
    )

    missing_features <- setdiff(features, available)

    if (length(missing_features)) {
      warning(
        "Feature(s) not found and omitted: ",
        paste(missing_features, collapse = ", "),
        call. = FALSE
      )
    }

    features <- intersect(
      unique(as.character(features)),
      available
    )

    if (!length(features)) {
      stop(
        "None of the requested features were found in the object.",
        call. = FALSE
      )
    }

    features <- intersect(
      unique(as.character(features)),
      c(rownames(object), colnames(object@meta.data))
    )

    if (length(features) == 0)
      stop("No valid features found")

    reduction <- reduction %||% SeuratObject::DefaultDimReduc(object)

    # Colors
    if (!is.null(cols)) {
      if (length(cols) < 2)
        stop("colors_use must contain at least 2 colors")
      cols <- cols
      if (rev.cols) cols <- rev(cols)

    } else  {
      pal <- list(
        default = c("lightgrey", "blue"),
        gnrh = c("#F2F2F2", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),
        hotspot = c("navy","#2C7BB6","#27F5EB","green", "yellow", "orange", "red", "darkred"),
        rainbow = c("navy", "green", "yellow", "orange", "red"),
        viridis = c("#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF"),
        magma = c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF"),
        inferno = c("#000004", "#420A68", "#932667", "#DD513A", "#FCA50A", "#FCFFA4"),
        plasma = c("#0D0887", "#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
        blue_red = c("#2166AC", "#4393C3", "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
        teal_orange = c("#1B9E77", "#66C2A5", "#F7F7F7", "#FDD0A2", "#E6550D"),
        purple_yellow = c("#3F007D", "#6A51A3", "#9E9AC8", "#FEE391", "#FEC44F", "#D95F0E"),
        brain = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")
      )
      if (theme.cols %in% names(pal)) {
        cols <- pal[[theme.cols]]
      } else {
        cols <- RColorBrewer::brewer.pal(9, theme.cols)
      }
      if (rev.cols) cols <- rev(cols)
    }


    # Point size
    raster <- raster %||% (ncol(object) > 2e5)
    if (is.null(pt.size)) {
      pt.size <- if (raster) 1 else min(1583 / ncol(object), 1)
    }

    ## Global max for merged legend
    global_max <- NA
    if (merge.leg && is.null(split.by) && length(features) > 1 && !blend) {
      expr_data <- Seurat::FetchData(object, vars = features, layer = layer)
      global_max <- max(expr_data, na.rm = TRUE)
    }

    # Plot
    plt <- Seurat::FeaturePlot(
      object = object,
      features = features,
      reduction = reduction,
      order = order,
      pt.size = pt.size,
      raster = raster,
      raster.dpi = raster.dpi,
      split.by = split.by,
      combine = combine,
      blend = blend,
      label = label,
      ncol = ncol,
      ...
    )

    # Determine palette behavior
    like <- theme.cols %in% c("hotspot", "rainbow","magma","inferno","plasma", "viridis",
                              "blue_red","teal_orange","purple_yellow","brain","mult")

    # Brewer palettes should gray-out low expression
    if (!like && !blend) {

      plt <- plt & theme(legend.position = "none") & scale_color_gradientn(
        colors = cols,
        limits = c(na.cutoff, NA),
        oob = scales::censor,     # force < cutoff → NA
        na.value = na.col,        # show gray
        name = NULL
      )} else if (!blend) {
        # No gray background
        plt <- plt & scale_color_gradientn(
          colors = cols,
          limits = c(na.cutoff, NA),
          oob = scales::squish,     # clip instead of NA
          na.value = NA
        )
      }

    # Theme
    plt <- plt & plot_theme(style = style, txtsize = txtsize, ...) &
      theme(plot.title = element_text(
        hjust = 0.5, size = txtsize + 2, face = "bold.italic"
      ))

    # ONLY add colorbar guide for non-blend plots
    if (!blend) {
      plt <- plt &
        guides(
          color = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black"
          )
        )
    }

    # Manage axes
    if (!axes) {
      plt <- plt & Seurat::NoAxes()
    }

    ## Merge legend
    if (merge.leg && is.null(split.by) && length(features) > 1) {
      plt <- Seurat::CombinePlots(plots = plt, legend = "right", ncol = ncol)
    }

    return(plt)
  }))
}



# ============================================================================= #
# Dot plots
# ============================================================================= #

#' Create an enhanced Seurat dot plot
#'
#' Generates a customizable dot plot for gene expression patterns across
#' clusters or metadata-defined groups. This function wraps
#' \code{Seurat::DotPlot()} and adds improved color control, optional dot
#' outlines, flexible axis formatting, legend customization, and GnRHcell
#' theme support.
#'
#' @param object A Seurat object.
#' @param features Character vector or named list of features to plot.
#' @param group.by Metadata column used to group cells.
#' Default is \code{"seurat_clusters"}.
#' @param th.cols RColorBrewer palette name used for the expression gradient.
#' Default is \code{"RdYlBu"}.
#' @param rev.th.cols Logical; reverse the color gradient.
#' Default is \code{TRUE}.
#' @param dot.scale Numeric dot-size scaling factor passed to
#' \code{Seurat::DotPlot()}. Default is \code{4}.
#' @param x.ang Angle of x-axis labels in degrees. Default is \code{90}.
#' @param vjust.x,hjust.x Vertical and horizontal justification for x-axis labels.
#' @param flip Logical; flip x and y axes using \code{ggplot2::coord_flip()}.
#' Default is \code{FALSE}.
#' @param txtsize Base text size.
#' Default is \code{12}.
#' @param title Optional plot title.
#' @param leg.size Legend text size. Default is \code{10}.
#' @param leg.ttl.size Legend title size. Default is \code{10}.
#' @param leg.pos Legend position. One of \code{"right"}, \code{"left"},
#' \code{"top"}, \code{"bottom"}, or \code{NULL}.
#' Default is \code{"right"}.
#' @param leg.just Legend justification. Default is \code{"bottom"}.
#' @param leg.hjust Logical; reserved for legend layout customization.
#' Default is \code{FALSE}.
#' @param x.axis.pos Position of the x-axis. Default is \code{"bottom"}.
#' @param style Theme style.
#' Default is \code{"classic"}.
#' @param x.face,y.face Logical; italicize x- or y-axis labels.
#' Default is \code{FALSE}.
#' @param x.ttl,y.ttl Logical; show x- or y-axis titles.
#' Default is \code{FALSE}.
#' @param dot.outline Logical; draw outlines around dots.
#' Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{Seurat::DotPlot()} and
#' the internal GnRHcell theme.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Dot color represents average expression and dot size represents the
#' percentage of cells expressing each feature, following the standard
#' \code{Seurat::DotPlot()} convention.
#'
#' When \code{dot.outline = TRUE}, dots are drawn with a light outline using
#' shape 21. This can improve readability in publication figures.
#'
#' @examples
#' \dontrun{
#' celldot(pbmc, features = c("MS4A1", "CD3D"))
#'
#' celldot(
#'   pbmc,
#'   features = c("MS4A1", "CD14"),
#'   th.cols = "Blues",
#'   dot.outline = TRUE,
#'   flip = TRUE
#' )
#' }
#'
#' @export
plot_gnrh_dot <- function(
    object,
    features,
    group.by = "seurat_clusters",
    th.cols = "RdYlBu",
    rev.th.cols = TRUE,
    dot.scale = 4,
    x.ang = 90,
    vjust.x = NULL,
    hjust.x = NULL,
    flip = FALSE,
    txtsize = 12,
    title = NULL,
    leg.size = 10,
    leg.ttl.size = 10,
    leg.pos = "right",
    leg.just = "bottom",
    leg.hjust = FALSE,
    x.axis.pos = "bottom",
    style = "classic",
    x.face = FALSE,
    y.face = FALSE,
    x.ttl = FALSE,
    y.ttl = FALSE,
    dot.outline = FALSE,
    ...
) {

  .check_seurat(object)

  if (!group.by %in% colnames(object@meta.data)) {
    stop("Grouping column not found: ", group.by, call. = FALSE)
  }

  object <- Seurat::SetIdent(object, value = group.by)

  features <- unique(features)

  if (!length(features)) {
    stop("features must be provided.", call. = FALSE)
  }

  # --------------------------------------------------- # # #
  # Palette
  # --------------------------------------------------- # # #

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required.", call. = FALSE)
  }

  if (!th.cols %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Unknown RColorBrewer palette: ", th.cols, call. = FALSE)
  }

  n_pal <- min(
    9,
    RColorBrewer::brewer.pal.info[th.cols, "maxcolors"]
  )

  pal <- RColorBrewer::brewer.pal(n_pal, th.cols)

  if (rev.th.cols) {
    pal <- rev(pal)
  }

  outline_col <- if (dot.outline) "gray60" else NA
  outline_stroke <- if (dot.outline) 0.5 else 0

  # --------------------------------------------------- # # #
  # Plot
  # --------------------------------------------------- # # #

  plt <- suppressWarnings(
    suppressMessages(
      Seurat::DotPlot(
        object = object,
        features = features,
        dot.scale = dot.scale,
        ...
      )
    )
  )

  plt <- plt +
    ggplot2::scale_color_gradientn(
      colors = pal,
      oob = scales::squish
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[["pct.exp"]]),
      shape = 21,
      colour = outline_col,
      stroke = outline_stroke
    ) +
    ggplot2::labs(
      title = title,
      color = "Average\nExpression",
      size = "Percent\nExpressed"
    ) +
    plot_theme(
      style = style,
      txtsize = txtsize,
      x.ang = x.ang,
      leg.size = leg.size,
      leg.ttl.size = leg.ttl.size,
      leg.pos = leg.pos,
      x.hjust = hjust.x,
      x.vjust = vjust.x,
      xy.val = TRUE,
      xlab = TRUE,
      ylab = TRUE,
      x.ttl = x.ttl,
      y.ttl = y.ttl,
      ...
    ) +
    ggplot2::theme(
      axis.text.x = if (x.face || (flip && y.face)) {
        ggplot2::element_text(face = "italic")
      } else {
        ggplot2::element_text()
      },
      axis.text.y = if (y.face || (flip && x.face)) {
        ggplot2::element_text(face = "italic")
      } else {
        ggplot2::element_text()
      },
      axis.title = ggplot2::element_blank(),
      legend.spacing.y = grid::unit(0.05, "cm"),
      legend.spacing.x = grid::unit(0.05, "cm"),
      legend.box.spacing = grid::unit(0.05, "cm"),
      legend.margin = ggplot2::margin(2, 2, 2, 2)
    )

  if (flip) {
    plt <- plt + ggplot2::coord_flip()
  }

  if (is.list(features)) {
    plt <- plt +
      ggplot2::theme(
        strip.text.x = ggplot2::element_text(angle = 45)
      )
  }

  # --------------------------------------------------- # # #
  # Legend positioning
  # --------------------------------------------------- # # #

  if (!is.null(leg.pos)) {

    if (leg.pos == "right") {
      plt <- plt +
        ggplot2::theme(
          legend.position = "right",
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.box.margin = ggplot2::margin(0, 0, 0, 0)
        )
    }

    if (leg.pos == "left") {
      plt <- plt +
        ggplot2::theme(
          legend.position = "left",
          legend.justification = c("left", "bottom"),
          legend.box.just = "left",
          legend.box.margin = ggplot2::margin(0, 0, 0, 0)
        )
    }

    if (leg.pos == "top") {
      plt <- plt +
        ggplot2::theme(
          legend.position = "top",
          legend.justification = c("right", "top"),
          legend.box.just = "right"
        )
    }

    if (leg.pos == "bottom") {
      plt <- plt +
        ggplot2::theme(
          legend.position = "bottom",
          legend.justification = c("right", "bottom"),
          legend.box.just = "right"
        )
    }

    if (leg.pos == "none") {
      plt <- plt + Seurat::NoLegend()
    }
  }

  # --------------------------------------------------- # # #
  # Guides
  # --------------------------------------------------- # # #

  guide_color <- ggplot2::guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black"
  )

  guide_size <- ggplot2::guide_legend(
    override.aes = list(
      shape = 21,
      colour = outline_col,
      fill = "black"
    )
  )

  guide_color$order <- 1
  guide_size$order <- 2

  plt <- plt +
    ggplot2::guides(
      color = guide_color,
      size = guide_size
    )

  plt
}




# ============================================================================= #
# Distribution plots
# ============================================================================= #

#' Plot GnRHcell metadata distributions
#'
#' @param object A Seurat object.
#' @param group.by Metadata column defining categories.
#' @param split.by Optional metadata column defining bars.
#' @param cols Optional named colour vector.
#' @param proportion Display within-split proportions instead of counts.
#' @param position Bar position, such as `"stack"` or `"dodge"`.
#' @param label Add value labels.
#' @param label.size Label text size.
#' @param plot.ttl Optional plot title.
#' @param txtsize Base text size.
#' @param x.ang X-axis label angle.
#' @param flip Flip coordinates.
#' @return A ggplot object.
#' @export
plot_gnrh_distribution <- function(
    object,
    group.by,
    split.by = NULL,
    cols = NULL,
    proportion = FALSE,
    position = "stack",
    label = TRUE,
    label.size = 3,
    plot.ttl = NULL,
    txtsize = 10,
    x.ang = 45,
    flip = FALSE
) {
  .validate_seurat(object)

  md <- object[[]]

  if (!group.by %in% names(md))
    stop("`group.by` not found.", call. = FALSE)

  if (!is.null(split.by) && !split.by %in% names(md))
    stop("`split.by` not found.", call. = FALSE)

  if (is.null(split.by)) {
    df <- as.data.frame(
      table(md[[group.by]]),
      stringsAsFactors = FALSE
    )
    names(df) <- c("group", "n")
    df$split <- "all"
  } else {
    df <- as.data.frame(
      table(md[[split.by]], md[[group.by]]),
      stringsAsFactors = FALSE
    )
    names(df) <- c("split", "group", "n")
  }

  if (proportion) {
    df <- df |>
      dplyr::group_by(.data$split) |>
      dplyr::mutate(value = .data$n / sum(.data$n)) |>
      dplyr::ungroup()
  } else {
    df$value <- df$n
  }

  df$label <- if (proportion)
    scales::percent(df$value, accuracy = 0.1)
  else
    format(df$n, big.mark = ",")

  groups <- unique(as.character(df$group))
  cols <- .resolve_colors(groups, cols, group.by)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$split, .data$value, fill = .data$group)
  ) +
    ggplot2::geom_col(
      position = position,
      width = 0.75,
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      title = plot.ttl,
      x = if (is.null(split.by)) NULL else split.by,
      y = if (proportion) "Proportion" else "Cells",
      fill = group.by
    )

  if (label) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        position = if (position == "stack")
          ggplot2::position_stack(vjust = 0.5)
        else
          ggplot2::position_dodge(width = 0.75),
        size = label.size
      )
  }

  if (flip) p <- p + ggplot2::coord_flip()

  p +
    .gnrh_theme(
      txtsize = txtsize,
      x.ang = x.ang,
      leg.pos = if (is.null(split.by)) "none" else "right"
    )
}


# ============================================================================= #
# Module hits
# ============================================================================= #

#' Plot GnRH module hit distributions
#'
#' @param data Data frame containing plotting columns.
#' @param x Column mapped to the x axis.
#' @param fill Column mapped to fill colour.
#' @param palette Optional named colour vector.
#' @param type Display counts or fractions.
#' @param title Optional plot title.
#' @param txtsize Base text size.
#' @return A ggplot object.
#' @export
plot_gnrh_hits <- function(
    data,
    x,
    fill,
    palette = NULL,
    type = c("count", "fraction"),
    title = NULL,
    txtsize = 10
) {
  type <- match.arg(type)

  if (!is.data.frame(data) || !all(c(x, fill) %in% names(data)))
    stop("Invalid input data or metadata columns.", call. = FALSE)

  groups <- unique(as.character(data[[fill]]))

  if (is.null(palette))
    palette <- .resolve_colors(groups, group_by = fill)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[x]], fill = .data[[fill]])
  ) +
    ggplot2::geom_bar(
      position = if (type == "count") "stack" else "fill",
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(
      title = title,
      x = x,
      y = if (type == "count") "Cell count" else "Fraction",
      fill = fill
    )

  if (type == "fraction")
    p <- p +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format()
    )

  p + .gnrh_theme(txtsize = txtsize, x.ang = 45)
}


# ============================================================================= #
# Compact diagnostic report
# ============================================================================= #

#' Generate a publication-ready GnRHcell diagnostic report
#'
#' Creates a six-panel quality-control dashboard summarizing detection,
#' developmental staging, score separation, threshold behavior, detection
#' classes or ROC discrimination, and marker-program support.
#'
#' @param object A Seurat object processed with `run_gnrh()` and
#'   `gnrh_diagnostics()`.
#' @param style Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.
#' @param mode Display mode: `"light"` or `"dark"`.
#' @param truth Optional metadata column containing an independent binary
#'   reference classification for external ROC analysis.
#' @param roc_mode ROC behavior: `"auto"` uses external truth when supplied and
#'   otherwise restores internal status/stage discrimination curves;
#'   `"external"`, `"internal"`, and `"none"` force a specific behavior.
#' @param positive_truth Values in `truth` interpreted as positive.
#' @param txtsize Base text size.
#' @param max_points Maximum cells displayed in each scatter panel. All rare
#'   GnRH-positive cells are retained before negative cells are sampled.
#' @param seed Random seed used for display-only subsampling.
#' @param show_class_panel Show detection-class composition when ROC analysis is
#'   disabled or unavailable.
#' @param verbose Print progress messages.
#'
#' @return A patchwork object.
#' @export
#'
#' @examples
#' \dontrun{
#' report <- gnrh_report(wang, style = "bw")
#' report
#'
#' # External ROC using an independent manual/reference annotation
#' report <- gnrh_report(
#'   wang,
#'   truth = "manual_gnrh",
#'   roc_mode = "external",
#'   positive_truth = c("GnRH", "pos", "TRUE"),
#'   style = "minimal"
#' )
#'
#' # Explicit internal score-discrimination curves
#' report <- gnrh_report(wang, roc_mode = "internal", style = "bw")
#'
#' ggplot2::ggsave(
#'   "gnrh_report.pdf", report,
#'   width = 12, height = 7.5, units = "in",
#'   device = grDevices::cairo_pdf
#' )
#' }
gnrh_report <- function(
    object,
    style = c("test", "classic", "minimal", "bw"),
    mode = c("light", "dark"),
    truth = NULL,
    roc_mode = c("auto", "external", "internal", "none"),
    positive_truth = c("1", "TRUE", "true", "positive", "Positive", "pos", "Pos"),
    txtsize = 10,
    max_points = 100000L,
    seed = 1234L,
    show_class_panel = TRUE,
    verbose = TRUE) {

  style <- match.arg(style)
  mode <- match.arg(mode)
  roc_mode <- match.arg(roc_mode)
  resolved_roc_mode <- if (roc_mode == "auto") {
    if (is.null(truth)) "internal" else "external"
  } else {
    roc_mode
  }

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package `patchwork` is required.", call. = FALSE)
  }
  if (!is.numeric(max_points) || length(max_points) != 1L || max_points < 100L) {
    stop("`max_points` must be one number greater than or equal to 100.", call. = FALSE)
  }
  if (resolved_roc_mode == "external" && is.null(truth)) {
    stop("`roc_mode = \"external\"` requires an independent `truth` column.", call. = FALSE)
  }
  if (resolved_roc_mode == "external" && identical(truth, "gnrh_status")) {
    stop("Use `roc_mode = \"internal\"` to evaluate `gnrh_status` self-discrimination.", call. = FALSE)
  }

  if (verbose) message("[INFO] Generating publication-ready GnRH QC report")

  g <- object@misc$gnrh
  if (is.null(g) || is.null(g$diagnostics)) {
    stop(
      "GnRH diagnostics are missing. Run `run_gnrh()` or both ",
      "`detect_gnrh()` and `gnrh_diagnostics()` first.",
      call. = FALSE
    )
  }

  diagnostics <- as.data.frame(g$diagnostics)
  metadata <- object[[]]
  threshold_curve <- g$threshold_curve

  required <- c("expr", "score", "status", "core_hits", "mig_hits", "neuro_hits")
  missing <- setdiff(required, names(diagnostics))
  if (length(missing)) {
    stop("Missing diagnostic column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  # Add current metadata fields without assuming that diagnostics already carry
  # all outputs. Prefer barcode matching; fall back to row order only when safe.
  metadata_fields <- c(
    "gnrh_class", "gnrh_confident", "gnrh_stage", "gnrh_secretory",
    "gnrh_support_score", "gnrh_alternative_score",
    "gnrh_identity_score", "gnrh_migrating_score",
    "gnrh_mature_score", "gnrh_secreting_score", "gnrh_knn"
  )
  metadata_fields <- intersect(metadata_fields, names(metadata))
  diagnostic_cells <- rownames(diagnostics)
  can_match <- length(diagnostic_cells) == nrow(diagnostics) &&
    all(diagnostic_cells %in% rownames(metadata))

  for (field in metadata_fields) {
    if (!field %in% names(diagnostics)) {
      diagnostics[[field]] <- if (can_match) {
        metadata[diagnostic_cells, field, drop = TRUE]
      } else if (nrow(diagnostics) == nrow(metadata)) {
        metadata[[field]]
      } else {
        stop("Cannot align diagnostic rows with Seurat metadata.", call. = FALSE)
      }
    }
  }

  foreground <- if (mode == "dark") "#F5F5F5" else "#1A1A1A"
  background <- if (mode == "dark") "#111111" else "white"
  muted <- if (mode == "dark") "#BDBDBD" else "#555555"
  tile_low <- if (mode == "dark") "#252525" else "white"
  tile_border <- if (mode == "dark") "#111111" else "white"

  report_theme <- function(
    leg.pos = "right",
    x.ang = 0,
    axes = TRUE,
    title.position = "left") {
    .gnrh_theme(
      txtsize = txtsize,
      x.ang = x.ang,
      leg.pos = leg.pos,
      title.position = title.position,
      axes = axes,
      style = style,
      mode = mode
    ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = txtsize + 1,
          margin = ggplot2::margin(b = 4)
        ),
        plot.subtitle = ggplot2::element_text(
          size = max(7, txtsize - 1), color = muted,
          margin = ggplot2::margin(b = 5)
        ),
        legend.key.height = grid::unit(0.32, "cm"),
        legend.key.width = grid::unit(0.32, "cm"),
        plot.margin = ggplot2::margin(5, 6, 5, 5)
      )
  }

  empty_panel <- function(title, subtitle = NULL) {
    ggplot2::ggplot() +
      ggplot2::labs(title = title, subtitle = subtitle) +
      report_theme(leg.pos = "none", axes = FALSE) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank()
      )
  }

  legend_labels <- function(values, levels) {
    values <- as.character(values)
    values[is.na(values)] <- "Unknown"
    counts <- table(values)
    levels <- levels[levels %in% names(counts)]
    labels <- sprintf(
      "%s (n = %s; %.1f%%)",
      levels,
      format(as.integer(counts[levels]), big.mark = ","),
      100 * as.integer(counts[levels]) / sum(counts)
    )
    stats::setNames(labels, levels)
  }

  compact_point_guide <- function() {
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 2.3, alpha = 1),
        keyheight = grid::unit(0.32, "cm")
      )
    )
  }

  # Display-only sampling preserves all positive/rare cells and samples the
  # background. Statistical summaries always use the complete diagnostics.
  plot_data <- diagnostics
  if (nrow(plot_data) > max_points) {
    set.seed(seed)
    positive_index <- as.character(plot_data$status) == "pos"
    positive_rows <- which(positive_index)
    negative_rows <- which(!positive_index)
    effective_max <- max(as.integer(max_points), length(positive_rows))
    room <- max(0L, effective_max - length(positive_rows))
    sampled_negative <- if (length(negative_rows) > room) {
      sample(negative_rows, room)
    } else {
      negative_rows
    }
    selected <- c(sampled_negative, positive_rows)
    plot_data <- plot_data[selected, , drop = FALSE]
  }

  status_colors <- gnrh_colors("status")
  status_levels <- intersect(names(status_colors), unique(as.character(diagnostics$status)))
  status_labels <- legend_labels(diagnostics$status, status_levels)

  # A — Detection landscape
  p1 <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(.data$expr, .data$score, colour = .data$status)
  ) +
    ggplot2::geom_point(alpha = 0.42, size = 0.5, stroke = 0) +
    ggplot2::scale_colour_manual(
      values = status_colors,
      breaks = names(status_labels),
      labels = status_labels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "GnRH detection landscape",
      x = expression(italic(GNRH1)~"expression"),
      y = "Composite score",
      colour = "Status"
    ) +
    report_theme() +
    compact_point_guide()

  # B — Developmental stages
  if ("gnrh_stage" %in% names(diagnostics)) {
    stage_colors <- gnrh_colors("stage")
    stage_levels <- intersect(names(stage_colors), unique(as.character(diagnostics$gnrh_stage)))
    stage_labels <- legend_labels(diagnostics$gnrh_stage, stage_levels)

    p2 <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(.data$expr, .data$score, colour = .data$gnrh_stage)
    ) +
      ggplot2::geom_point(alpha = 0.42, size = 0.5, stroke = 0) +
      ggplot2::scale_colour_manual(
        values = stage_colors,
        breaks = names(stage_labels),
        labels = stage_labels,
        drop = FALSE
      ) +
      ggplot2::labs(
        title = "Developmental-stage landscape",
        x = expression(italic(GNRH1)~"expression"),
        y = "Composite score",
        colour = "Stage"
      ) +
      report_theme() +
      compact_point_guide()
  } else {
    p2 <- empty_panel("Developmental stages unavailable")
  }

  # C — Score separation
  density_data <- diagnostics[is.finite(diagnostics$score), , drop = FALSE]
  p3 <- ggplot2::ggplot(
    density_data,
    ggplot2::aes(.data$score, fill = .data$status, colour = .data$status)
  ) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 0.55, adjust = 1) +
    ggplot2::scale_fill_manual(values = status_colors, breaks = status_levels) +
    ggplot2::scale_colour_manual(values = status_colors, breaks = status_levels) +
    ggplot2::labs(
      title = "GnRH score separation",
      x = "Composite score", y = "Density",
      fill = "Status", colour = "Status"
    ) +
    report_theme()

  # D — Threshold performance
  curve_columns <- c("threshold", "sensitivity", "specificity", "F1")
  if (!is.null(threshold_curve) && all(curve_columns %in% names(threshold_curve))) {
    curve_data <- as.data.frame(threshold_curve)
    valid <- is.finite(curve_data$threshold) & is.finite(curve_data$F1)
    best <- if (any(valid)) {
      curve_data$threshold[which.max(replace(curve_data$F1, !valid, -Inf))]
    } else {
      NA_real_
    }

    long_curve <- rbind(
      data.frame(threshold = curve_data$threshold, metric = "Sensitivity", value = curve_data$sensitivity),
      data.frame(threshold = curve_data$threshold, metric = "Specificity", value = curve_data$specificity),
      data.frame(threshold = curve_data$threshold, metric = "F1", value = curve_data$F1)
    )
    metric_colors <- c(Sensitivity = "#0072B2", Specificity = "#009E73", F1 = "#D55E00")

    p4 <- ggplot2::ggplot(
      long_curve,
      ggplot2::aes(.data$threshold, .data$value, colour = .data$metric)
    ) +
      ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
      ggplot2::geom_vline(
        xintercept = best, linetype = "dashed",
        colour = muted, linewidth = 0.45, na.rm = TRUE
      ) +
      ggplot2::scale_colour_manual(values = metric_colors) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.03))) +
      ggplot2::labs(
        title = "Threshold performance",
        subtitle = if (is.finite(best)) sprintf("Optimal threshold: %.3f", best) else NULL,
        x = "Threshold", y = "Performance", colour = NULL
      ) +
      report_theme()
  } else {
    p4 <- empty_panel("Threshold curve unavailable")
  }

  # E — Restored ROC analysis or detection classes
  p5 <- NULL

  roc_add <- function(reference, predictor, label) {
    reference <- as.integer(reference)
    predictor <- as.numeric(predictor)
    valid <- !is.na(reference) & is.finite(predictor)
    if (sum(valid) < 2L || length(unique(reference[valid])) != 2L) return(NULL)

    roc <- pROC::roc(
      response = reference[valid], predictor = predictor[valid],
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
    data.frame(
      FPR = 1 - roc$specificities,
      TPR = roc$sensitivities,
      group = label,
      auc = as.numeric(pROC::auc(roc))
    )
  }

  if (resolved_roc_mode != "none") {
    if (!requireNamespace("pROC", quietly = TRUE)) {
      if (verbose) message("[INFO] Package `pROC` unavailable; using class-composition panel")
    } else {
      roc_list <- list()

      if (resolved_roc_mode == "external") {
        if (!truth %in% names(metadata)) {
          stop("Truth column `", truth, "` was not found in object metadata.", call. = FALSE)
        }
        truth_values <- if (can_match) {
          metadata[diagnostic_cells, truth, drop = TRUE]
        } else {
          metadata[[truth]]
        }
        reference <- as.character(truth_values) %in% as.character(positive_truth)
        roc_list[["External GnRH status"]] <- roc_add(
          reference, diagnostics$score, "External GnRH status"
        )
      } else {
        # These restored curves quantify internal score discrimination. The
        # classifications were derived from the same score system, so they are
        # intentionally labelled as diagnostic rather than external validation.
        roc_list[["GnRH status"]] <- roc_add(
          diagnostics$status == "pos", diagnostics$score, "GnRH status"
        )

        if ("gnrh_stage" %in% names(diagnostics)) {
          stage_scores <- c(
            identity = "gnrh_identity_score",
            migrating = "gnrh_migrating_score",
            mature = "gnrh_mature_score",
            secreting = "gnrh_secreting_score"
          )
          for (stage in names(stage_scores)) {
            score_column <- stage_scores[[stage]]
            if (
              stage %in% as.character(diagnostics$gnrh_stage) &&
              score_column %in% names(diagnostics)
            ) {
              roc_list[[stage]] <- roc_add(
                diagnostics$gnrh_stage == stage,
                diagnostics[[score_column]],
                stage
              )
            }
          }
        }
      }

      roc_list <- Filter(Negate(is.null), roc_list)
      if (length(roc_list)) {
        roc_data <- do.call(rbind, roc_list)
        auc_data <- unique(roc_data[c("group", "auc")])
        auc_data$label <- sprintf("%s (AUC = %.3f)", auc_data$group, auc_data$auc)
        roc_data$group <- factor(roc_data$group, levels = auc_data$group)

        roc_colors <- c(
          "External GnRH status" = status_colors[["pos"]],
          "GnRH status" = status_colors[["pos"]],
          gnrh_colors("stage")
        )
        missing_colors <- setdiff(auc_data$group, names(roc_colors))
        if (length(missing_colors)) {
          extra <- grDevices::hcl.colors(length(missing_colors), "Dark 3")
          names(extra) <- missing_colors
          roc_colors <- c(roc_colors, extra)
        }

        roc_title <- if (resolved_roc_mode == "external") {
          "External ROC performance"
        } else {
          "Internal score discrimination"
        }
        roc_subtitle <- if (resolved_roc_mode == "external") {
          paste0("Independent reference: ", truth)
        } else {
          "Diagnostic self-consistency; not independent validation"
        }

        p5 <- ggplot2::ggplot(
          roc_data,
          ggplot2::aes(.data$FPR, .data$TPR, colour = .data$group)
        ) +
          ggplot2::geom_abline(
            slope = 1, intercept = 0, linetype = "dashed",
            colour = muted, linewidth = 0.4
          ) +
          ggplot2::geom_line(linewidth = 0.8) +
          ggplot2::scale_colour_manual(
            values = roc_colors,
            breaks = auc_data$group,
            labels = stats::setNames(auc_data$label, auc_data$group)
          ) +
          ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
          ggplot2::labs(
            title = roc_title,
            subtitle = roc_subtitle,
            x = "False-positive rate",
            y = "True-positive rate",
            colour = NULL
          ) +
          report_theme() +
          ggplot2::guides(
            colour = ggplot2::guide_legend(
              override.aes = list(linewidth = 1.1),
              keyheight = grid::unit(0.32, "cm")
            )
          )
      }
    }
  }

  if (is.null(p5) && show_class_panel && "gnrh_class" %in% names(diagnostics)) {
    class_order <- c("neg", "supported", "direct")
    class_data <- as.data.frame(table(as.character(diagnostics$gnrh_class)), stringsAsFactors = FALSE)
    names(class_data) <- c("class", "n")
    class_data$class <- factor(class_data$class, levels = class_order)
    class_data <- class_data[!is.na(class_data$class), , drop = FALSE]
    class_data$pct <- 100 * class_data$n / sum(class_data$n)
    class_colors <- gnrh_colors("class")

    p5 <- ggplot2::ggplot(class_data, ggplot2::aes(.data$class, .data$n, fill = .data$class)) +
      ggplot2::geom_col(width = 0.68, colour = tile_border, linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%s\n%.1f%%", format(.data$n, big.mark = ","), .data$pct)),
        vjust = -0.2, size = txtsize / 3.2, colour = foreground
      ) +
      ggplot2::scale_fill_manual(values = class_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_comma(),
        expand = ggplot2::expansion(mult = c(0, 0.16))
      ) +
      ggplot2::labs(title = "GnRH detection classes", x = NULL, y = "Cells") +
      report_theme(leg.pos = "none", x.ang = 25)
  }
  if (is.null(p5)) p5 <- empty_panel("Detection-class panel unavailable")

  # F — Marker-program support
  module_data <- data.frame(
    status = rep(diagnostics$status, 3L),
    module = rep(c("Core", "Migration", "Neuroendocrine"), each = nrow(diagnostics)),
    hits = c(diagnostics$core_hits, diagnostics$mig_hits, diagnostics$neuro_hits)
  )
  module_summary <- stats::aggregate(
    hits ~ status + module,
    data = module_data,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  module_summary$module <- factor(
    module_summary$module,
    levels = c("Core", "Migration", "Neuroendocrine")
  )
  module_summary$status <- factor(module_summary$status, levels = c("neg", "pos"))

  p6 <- ggplot2::ggplot(
    module_summary,
    ggplot2::aes(.data$module, .data$status, fill = .data$hits)
  ) +
    ggplot2::geom_tile(colour = tile_border, linewidth = 0.45) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$hits)),
      size = txtsize / 3.2,
      colour = foreground
    ) +
    ggplot2::scale_fill_gradient(low = tile_low, high = status_colors[["pos"]]) +
    ggplot2::labs(
      title = "Marker-program support",
      x = NULL, y = NULL, fill = "Mean hits"
    ) +
    report_theme(x.ang = 25)

  # Assemble
  n_positive <- sum(as.character(diagnostics$status) == "pos", na.rm = TRUE)
  n_confident <- if ("gnrh_confident" %in% names(diagnostics)) {
    sum(as.character(diagnostics$gnrh_confident) %in% c("TRUE", "true", "1", "pos"), na.rm = TRUE)
  } else {
    NA_integer_
  }
  subtitle <- paste0(
    "Cells: ", format(ncol(object), big.mark = ","),
    "  |  Features: ", format(nrow(object), big.mark = ","),
    "  |  GnRH+: ", format(n_positive, big.mark = ","),
    if (!is.na(n_confident)) paste0("  |  Confident: ", format(n_confident, big.mark = ",")) else ""
  )

  annotation_theme <- ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = background, colour = NA),
    plot.title = ggplot2::element_text(
      family = "", face = "bold", size = txtsize + 5,
      colour = foreground, hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(
      family = "", size = txtsize, colour = muted, hjust = 0
    ),
    plot.tag = ggplot2::element_text(
      family = "", face = "bold", size = txtsize + 1,
      colour = foreground
    )
  )

  ((p1 | p2 | p3) / (p4 | p5 | p6)) +
    patchwork::plot_layout(guides = "keep", heights = c(1, 1)) +
    patchwork::plot_annotation(
      title = "GnRHcell diagnostic report",
      subtitle = subtitle,
      tag_levels = "A",
      theme = annotation_theme
    )
}


# ============================================================================= #
# Coexpression network
# ============================================================================= #

#' GnRH gene similarity network
#'
#' @param df Marker table containing `gene`, `coexpr`, and `score`.
#' @param top_n Maximum number of ranked genes included.
#' @param threshold Minimum co-expression edge weight.
#' @return A ggraph object.
#' @export
plot_network <- function(df, top_n = 25, threshold = 0.4) {
  req <- c("gene", "coexpr", "score")
  miss <- setdiff(req, names(df))

  if (length(miss))
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  df <- df[seq_len(min(top_n, nrow(df))), , drop = FALSE]

  if (nrow(df) < 2L)
    stop("At least two genes are required.", call. = FALSE)

  genes <- df$gene

  mat <- outer(
    df$coexpr,
    df$coexpr,
    pmin
  )

  rownames(mat) <- genes
  colnames(mat) <- genes

  edges <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(edges) <- c("from", "to", "weight")

  edges <- edges[
    edges$from < edges$to &
      is.finite(edges$weight) &
      edges$weight > threshold,
    ,
    drop = FALSE
  ]

  if (!nrow(edges))
    stop(
      "No network edges remain above `threshold`.",
      call. = FALSE
    )

  nodes <- data.frame(
    name = df$gene,
    score = df$score,
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(
    edges,
    vertices = nodes,
    directed = FALSE
  )

  igraph::V(g)$module <- factor(
    igraph::cluster_louvain(g)$membership
  )

  set.seed(123)

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      colour = "grey75",
      alpha = 0.7
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(colour = module, size = score)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      max.overlaps = 30
    ) +
    ggplot2::guides(edge_width = "none") +
    ggplot2::labs(
      title = "GnRH coexpression network",
      colour = "Module",
      size = "Score"
    ) +
    ggplot2::theme_void()
}


# ============================================================================= #
# Coexpression markers
# ============================================================================= #

#' Plot GnRH coexpression markers
#'
#' @param df Marker table containing `gene`, `coexpr`, `score`, and `p_val_adj`.
#' @param coexp_cutoff Minimum GNRH1 co-expression value.
#' @param txtsize Base text size.
#' @param style Theme style.
#' @return A ggplot object.
#' @export
plot_gnrh_coexpr <- function(
    df,
    coexp_cutoff = 0.25,
    txtsize = 10,
    style = "bw"
) {
  req <- c("gene", "coexpr", "score", "p_val_adj")
  miss <- setdiff(req, names(df))

  if (length(miss))
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  df <- df[
    df$gene != "GNRH1" &
      is.finite(df$coexpr) &
      df$coexpr >= coexp_cutoff,
    ,
    drop = FALSE
  ]

  if (!nrow(df))
    stop("No genes pass `coexp_cutoff`.", call. = FALSE)

  df$log_padj <- -log10(pmax(df$p_val_adj, 1e-300))
  df <- df[order(df$coexpr, decreasing = TRUE), , drop = FALSE]
  df$gene <- factor(df$gene, levels = rev(df$gene))

  ggplot2::ggplot(
    df,
    ggplot2::aes(coexpr, gene)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = score,
        colour = log_padj
      ),
      alpha = 0.9
    ) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(
      title = "GNRH1 coexpression",
      x = "Correlation coefficient",
      y = NULL,
      size = "Score",
      colour = expression(-log[10](adj.~p))
    ) +
    .gnrh_theme(txtsize = txtsize, style = style) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic")
    )
}


# ============================================================================= #
# Runtime
# ============================================================================= #

#' Plot GnRHcell runtime across datasets
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
plot_gnrh_runtime_curve <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    metric = "total_sec",
    x.ang = 45,
    txtsize = 10,
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
      colour = gnrh_colors("status")[["pos"]]
    )

  p +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(runtime, 1), "s")),
      vjust = -0.7,
      size = 3
    ) +
    ggplot2::labs(
      title = "GnRHcell runtime across datasets",
      x = "Dataset",
      y = metric
    ) +
    .gnrh_theme(
      x.ang = x.ang,
      txtsize = txtsize
    )
}


#' Plot detected GnRH-positive cells across datasets
#'
#' @param files Optional named vector of run-information files.
#' @param dir Directory searched when `files` is `NULL`.
#' @param pattern File-selection regular expression.
#' @param x.ang X-axis label angle.
#' @param txtsize Base text size.
#' @param debug Print the imported summary columns.
#' @return A ggplot object.
#' @export
plot_gnrh_detected <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    x.ang = 45,
    txtsize = 10,
    debug = FALSE
) {
  stats <- .load_gnrh_stats(
    files = files,
    dir = dir,
    pattern = pattern
  )

  if (debug)
    print(stats[, intersect(
      c("source_file", "dataset", "n_cells", "neg", "pos"),
      names(stats)
    )])

  stats <- stats[, c("dataset", "pos")]
  stats$pos <- as.numeric(stats$pos)

  stats <- stats::aggregate(
    pos ~ dataset,
    stats,
    function(x) mean(x, na.rm = TRUE)
  )

  stats$dataset <- factor(
    stats$dataset,
    levels = stats$dataset[order(stats$pos)]
  )

  ggplot2::ggplot(
    stats,
    ggplot2::aes(dataset, pos)
  ) +
    ggplot2::geom_col(
      width = 0.7,
      fill = gnrh_colors("status")[["pos"]]
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(pos)),
      vjust = -0.3,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(
      title = "Detected GnRH cells",
      x = "Dataset",
      y = "GnRH+ cells"
    ) +
    .gnrh_theme(
      x.ang = x.ang,
      txtsize = txtsize
    )
}


# ============================================================================= #
# Marker programs
# ============================================================================= #

#' Plot GnRH marker program results
#'
#' @param programs Result returned by `gnrh_marker_programs()`.
#' @param table Result table to visualize.
#' @param type Plot type: bar, dot, or tile.
#' @param min_genes Minimum genes retained for summary plots.
#' @param mode Light or dark display mode.
#' @param txtsize Base text size.
#' @param x.ang X-axis label angle.
#' @param style Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.
#' @return A ggplot object.
#' @export
plot_gnrh_marker_programs <- function(
    programs,
    table = c("summary", "high_confidence", "candidate_table"),
    type = c("bar", "dot", "tile"),
    min_genes = 1,
    mode = "light",
    txtsize = 12,
    x.ang = 45,
    style = c("bw","test", "classic", "minimal")
) {
  table <- match.arg(table)
  type <- match.arg(type)
  style <- match.arg(style)

  df <- programs[[table]]

  if (is.null(df) || !nrow(df))
    stop("Selected table is empty: ", table, call. = FALSE)

  if (table == "summary") {
    df$n_genes <- as.numeric(df$n_genes)
    df <- df[df$n_genes >= min_genes, , drop = FALSE]

    if (!nrow(df))
      stop("No rows after applying `min_genes`.", call. = FALSE)

    if (type == "bar") {
      return(
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            .data$program,
            .data$n_genes,
            fill = .data$confidence_level
          )
        ) +
          ggplot2::geom_col(width = 0.75) +
          ggplot2::facet_wrap(~dataset, scales = "free_y") +
          ggplot2::labs(
            title = "GnRH candidate marker programs",
            x = "Developmental program",
            y = "Number of genes",
            fill = "Confidence"
          ) +
          .gnrh_theme(
            mode = mode,
            txtsize = txtsize,
            x.ang = x.ang,
            style = style
          )
      )
    }

    if (type == "dot") {
      return(
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            .data$program,
            .data$dataset,
            size = .data$n_genes,
            colour = .data$confidence_level
          )
        ) +
          ggplot2::geom_point(alpha = 0.85) +
          ggplot2::labs(
            title = "GnRH marker program enrichment",
            x = "Developmental program",
            y = "Dataset",
            size = "Genes",
            colour = "Confidence"
          ) +
          .gnrh_theme(
            mode = mode,
            txtsize = txtsize,
            x.ang = x.ang,
            style = style
          )
      )
    }

    return(
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          .data$program,
          .data$dataset,
          fill = .data$n_genes
        )
      ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.4
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$n_genes),
          size = 3
        ) +
        ggplot2::scale_fill_gradient(
          low = "white",
          high = "#2563EB"
        ) +
        ggplot2::labs(
          title = "GnRH developmental marker program map",
          x = "Developmental program",
          y = "Dataset",
          fill = "Genes"
        ) +
        .gnrh_theme(
          mode = mode,
          txtsize = txtsize,
          x.ang = x.ang,
          style = style
        )
    )
  }

  if (!"coexpr_GNRH1" %in% names(df))
    stop("Column `coexpr_GNRH1` is required.", call. = FALSE)

  df <- df[df$coexpr_GNRH1 %in% TRUE, , drop = FALSE]

  if (!nrow(df))
    stop(
      "No GNRH1 co-expressed genes found in selected table.",
      call. = FALSE
    )

  if (!"specificity_score" %in% names(df))
    df$specificity_score <- 0

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      .data$program,
      .data$gene,
      colour = .data$confidence_level,
      size = .data$specificity_score
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~dataset, scales = "free_y") +
    ggplot2::labs(
      title = "GnRH candidate marker genes by developmental program",
      x = "Developmental program",
      y = "Gene",
      colour = "Confidence",
      size = "Specificity score"
    ) +
    .gnrh_theme(
      mode = mode,
      txtsize = txtsize,
      x.ang = x.ang,
      style = style
    )
}


# ============================================================================= #
# Simple diagnostic plots
# ============================================================================= #

#' Plot GnRH classification counts
#'
#' @param object A Seurat object processed by `run_gnrh()`.
#' @param txtsize Base text size.
#' @return A ggplot object.
#' @export
plot_class_counts <- function(object, txtsize = 9) {
  .validate_seurat(object)

  df <- object[[]] |>
    dplyr::count(.data$gnrh_class, name = "n") |>
    dplyr::mutate(pct = 100 * .data$n / sum(.data$n))

  cols <- gnrh_colors("class")

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
    .gnrh_theme(
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
plot_gnrh_specificity <- function(object, txtsize = 9) {
  .validate_seurat(object)

  df <- object[[]]

  req <- c(
    "gnrh_alternative_score",
    "gnrh_support_score",
    "gnrh_status"
  )

  miss <- setdiff(req, names(df))
  if (length(miss))
    stop(
      "Missing metadata columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      .data$gnrh_alternative_score,
      .data$gnrh_support_score,
      colour = .data$gnrh_status
    )
  ) +
    ggplot2::geom_point(
      size = 0.25,
      alpha = 0.35
    ) +
    ggplot2::scale_colour_manual(
      values = gnrh_colors("status")
    ) +
    ggplot2::labs(
      title = "GnRH specificity landscape",
      x = "Alternative neuronal program score",
      y = "GnRH support score",
      colour = "Status"
    ) +
    .gnrh_theme(txtsize = txtsize)
}
