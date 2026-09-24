# Embedding plot
# ========================================================================= #

#' Plot GnRH embedding
#'
#' Visualize GnRHcell categorical metadata on a dimensional reduction.
#'
#' @param object A Seurat object.
#' @param group_by One or more metadata columns used to colour cells.
#' @param split_by Optional metadata column used to facet the embedding.
#' @param reduction Dimensional reduction.
#' @param dims Two reduction dimensions to display.
#' @param shuffle Randomize plotting order.
#' @param raster Rasterize points; selected automatically when \code{NULL}.
#' @param raster.threshold Number of cells from which automatic rasterization
#'   is enabled.
#' @param raster.dpi Raster resolution.
#' @param alpha Opacity of highlighted cells.
#' @param background_alpha Opacity of background cells.
#' @param background.col Colour used for background groups such as
#'   \code{"neg"} and \code{"non-gnrh"}. When \code{NULL}, a light or dark
#'   mode default is selected automatically.
#' @param n.cells Add cell counts to legend labels.
#' @param percentage Add percentages to legend labels.
#' @param label Label groups on the embedding.
#' @param repel Use repelled labels.
#' @param label.size Label text size.
#' @param label.face Label font face.
#' @param cols Optional named colour vector.
#' @param axes Show embedding axis titles.
#' @param plot.ttl Optional plot title.
#' @param legend Show legend.
#' @param leg.ttl Legend title. Explicit \code{NULL} removes the title.
#' @param leg.ttl.size Legend-title size.
#' @param item.size Legend point size.
#' @param leg.pos Legend position.
#' @param leg.dir Legend direction.
#' @param leg.size Legend text size.
#' @param leg.ncol Number of legend columns.
#' @param item.border Draw bordered legend points.
#' @param txtsize Base text size.
#' @param pt.size Point size.
#' @param dark Use dark mode.
#' @param total.cells Add total cell number to title.
#' @param style Plot style.
#' @param facet.bg Draw facet-strip background.
#' @param ncol Number of columns for facets or multiple plots.
#' @param ... Additional graphical arguments.
#'
#' @return A ggplot or patchwork object.
#'
#' @export
gnrh_cellmap <- function(
    object,
    group_by = "gnrh_status",
    split_by = NULL,
    reduction = NULL,
    dims = c(1, 2),
    shuffle = FALSE,
    raster = NULL,
    raster.threshold = 3e5,
    raster.dpi = 600,
    alpha = 1,
    background_alpha = 0.5,
    background.col = NULL,
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
    item.size = 4,
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = NULL,
    leg.ncol = NULL,
    item.border = TRUE,
    txtsize = getOption("gnrhcell.base_size", 14),
    pt.size = NULL,
    dark = FALSE,
    total.cells = FALSE,
    style = c("test","classic", "minimal", "bw", "void", "dirty", "gray"),
    facet.bg = FALSE,
    ncol = NULL,
    ...
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  .validate_seurat(object)
  md <- object[[]]
  dots <- list(...)
  if (!is.null(dots$theme)) style <- dots$theme
  style <- match.arg(style)
  leg.ttl.missing <- missing(leg.ttl)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) stop("`alpha` must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(background_alpha) || length(background_alpha) != 1L || !is.finite(background_alpha) || background_alpha < 0 || background_alpha > 1) stop("`background_alpha` must be between 0 and 1.", call. = FALSE)
  if (is.null(background.col)) background.col <- if (dark) "grey55" else "grey75"
  # ========================================================================= #
  # Multiple metadata variables
  # ========================================================================= #
  if (length(group_by) > 1L) {
    plots <- lapply(group_by, function(g) {
      gnrh_cellmap(
        object = object,
        group_by = g,
        split_by = split_by,
        reduction = reduction,
        dims = dims,
        shuffle = shuffle,
        raster = raster,
        raster.threshold = raster.threshold,
        raster.dpi = raster.dpi,
        alpha = alpha,
        background_alpha = background_alpha,
        background.col = background.col,
        n.cells = n.cells,
        percentage = percentage,
        label = label,
        repel = repel,
        label.size = label.size,
        label.face = label.face,
        cols = cols,
        axes = axes,
        plot.ttl = g,
        legend = legend,
        leg.ttl = if (leg.ttl.missing) g else leg.ttl,
        leg.ttl.size = leg.ttl.size,
        item.size = item.size,
        leg.pos = leg.pos,
        leg.dir = leg.dir,
        leg.size = leg.size,
        leg.ncol = leg.ncol,
        item.border = item.border,
        txtsize = txtsize,
        pt.size = pt.size,
        dark = dark,
        total.cells = total.cells,
        style = style,
        facet.bg = facet.bg,
        ncol = ncol
      )
    })
    return(patchwork::wrap_plots(plots, ncol = ncol))
  }
  if (is.null(group_by) || length(group_by) != 1L || !group_by %in% colnames(md)) stop("Metadata column `", group_by, "` not found.", call. = FALSE)
  if (length(dims) != 2L || any(!is.finite(dims))) stop("`dims` must contain exactly two valid dimensions.", call. = FALSE)
  if (leg.ttl.missing) leg.ttl <- group_by
  # ========================================================================= #
  # Reduction
  # ========================================================================= #
  reduction <- gnrh_reduction(object, reduction)
  emb <- Seurat::Embeddings(object, reduction = reduction)
  if (max(dims) > ncol(emb)) stop("Selected dimensions exceed those available in `", reduction, "`.", call. = FALSE)
  # ========================================================================= #
  # Plot defaults
  # ========================================================================= #
  defaults <- .plot_defaults(
    object, raster = raster, pt.size = pt.size,
    raster.threshold = raster.threshold
  )
  raster <- defaults$raster
  pt.size <- defaults$pt.size
  # ========================================================================= #
  # Levels and colours
  # ========================================================================= #
  orders <- list(
    gnrh_status = c("neg", "pos"),
    gnrh_class = c("neg", "supported", "direct"),
    gnrh_confident = c("FALSE", "TRUE"),
    gnrh_stage = c(
      "non-gnrh", "early", "migrating", "post-migratory", "mature",
      "transitional", "undetermined"
    ),
    gnrh_stage_raw = c("non-gnrh", "early", "migrating", "mature"),
    gnrh_secretory = c("non-gnrh", "limited", "supported")
  )
  background_groups <- c("neg", "non-gnrh", "FALSE")
  values <- as.character(md[[group_by]])
  values[is.na(values)] <- "Unknown"
  lev <- orders[[group_by]] %||% sort(unique(values))
  lev <- c(lev[lev %in% values], sort(setdiff(unique(values), lev)))
  cols <- .resolve_colors(lev, cols = cols, group_by = group_by)
  background_present <- intersect(background_groups, names(cols))
  if (length(background_present)) cols[background_present] <- background.col
  if ("Unknown" %in% names(cols)) cols["Unknown"] <- if (dark) "grey50" else "grey70"
  # ========================================================================= #
  # Plotting data
  # ========================================================================= #
  df <- data.frame(
    x = emb[, dims[1]],
    y = emb[, dims[2]],
    group = factor(values, levels = lev),
    stringsAsFactors = FALSE
  )
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(md)) stop("Metadata column `", split_by, "` not found.", call. = FALSE)
    split_values <- md[[split_by]]
    df$split <- if (is.factor(split_values)) droplevels(split_values) else factor(split_values, levels = unique(split_values))
  }
  # ========================================================================= #
  # Background / foreground plotting order
  # ========================================================================= #
  df$.alpha <- ifelse(as.character(df$group) %in% background_groups, background_alpha, alpha)
  if (isTRUE(shuffle)) {
    set.seed(42)
    df <- df[sample.int(nrow(df)), , drop = FALSE]
  } else {
    df <- df[order(df$.alpha), , drop = FALSE]
  }
  # ========================================================================= #
  # Legend
  # ========================================================================= #
  present <- lev[lev %in% unique(as.character(df$group))]
  cols_use <- cols[present]
  labels <- .gnrh_legend_labels(
    x = values,
    levels = present,
    percentage = percentage,
    n_cells = n.cells
  )
  leg.ttl.size <- leg.ttl.size %||% txtsize
  leg.size <- leg.size %||% max(8, txtsize - 2)
  leg.ncol <- leg.ncol %||% ifelse(length(present) > 30L, 2L, 1L)
  # ========================================================================= #
  # Base plot
  # ========================================================================= #
  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      colour = .data$group,
      alpha = .data$.alpha
    )
  )
  if (isTRUE(raster)) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) stop("Package `ggrastr` is required when `raster = TRUE`.", call. = FALSE)
    plt <- plt + ggrastr::geom_point_rast(size = pt.size, raster.dpi = raster.dpi)
  } else {
    plt <- plt + ggplot2::geom_point(size = pt.size, shape = 16)
  }
  plt <- plt + ggplot2::scale_alpha_identity()
  # ========================================================================= #
  # Facets
  # ========================================================================= #
  if (!is.null(split_by)) {
    plt <- plt + ggplot2::facet_wrap(~split, ncol = ncol, scales = "fixed")
  }
  # ========================================================================= #
  # Colours
  # ========================================================================= #
  plt <- plt + ggplot2::scale_colour_manual(
    values = cols_use,
    breaks = present,
    labels = labels[present],
    drop = FALSE
  )
  # ========================================================================= #
  # Labels
  # ========================================================================= #
  if (isTRUE(label)) {
    if (is.null(split_by)) {
      lab_df <- stats::aggregate(cbind(x, y) ~ group, df, stats::median)
    } else {
      lab_df <- df |>
        dplyr::group_by(.data$split, .data$group) |>
        dplyr::summarise(
          x = stats::median(.data$x),
          y = stats::median(.data$y),
          .groups = "drop"
        )
    }
    if (isTRUE(repel) && requireNamespace("ggrepel", quietly = TRUE)) {
      plt <- plt + ggrepel::geom_text_repel(
        data = lab_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$group),
        inherit.aes = FALSE,
        colour = if (dark) "white" else "black",
        fontface = label.face,
        size = label.size,
        seed = 42
      )
    } else {
      plt <- plt + ggplot2::geom_text(
        data = lab_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$group),
        inherit.aes = FALSE,
        colour = if (dark) "white" else "black",
        fontface = label.face,
        size = label.size
      )
    }
  }
  # ========================================================================= #
  # Title
  # ========================================================================= #
  if (isTRUE(total.cells)) {
    plot.ttl <- paste0(
      plot.ttl %||% group_by,
      " (n = ",
      format(ncol(object), big.mark = ",", trim = TRUE),
      ")"
    )
  }
  # ========================================================================= #
  # Theme
  # ========================================================================= #
  plt <- plt +
    ggplot2::labs(
      title = plot.ttl,
      x = paste0(toupper(reduction), "_", dims[1]),
      y = paste0(toupper(reduction), "_", dims[2]),
      colour = leg.ttl
    ) +
    gnrh_theme(
      style = style,
      txtsize = txtsize,
      leg.pos = if (isTRUE(legend)) leg.pos else "none",
      mode = if (dark) "dark" else "light",
      xy.val = axes,
      xy.lab = axes,
      x.ttl = axes,
      y.ttl = axes,
      facet.bg = facet.bg
    ) +
    ggplot2::theme(
      legend.direction = leg.dir,
      legend.title = if (is.null(leg.ttl)) ggplot2::element_blank() else ggplot2::element_text(size = leg.ttl.size, face = "bold"),
      legend.text = ggplot2::element_text(size = leg.size),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.key.spacing.y = grid::unit(0, "pt"),
      strip.background = if (isTRUE(facet.bg)) ggplot2::element_rect(fill = if (dark) "grey25" else "grey90", colour = NA) else ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
  # ========================================================================= #
  # Compact legend
  # ========================================================================= #
  if (isTRUE(legend)) {
    override <- if (isTRUE(item.border)) {
      list(
        size = item.size,
        shape = 21,
        colour = if (dark) "white" else "black",
        fill = unname(cols_use),
        stroke = 0.25,
        alpha = 1
      )
    } else {
      list(
        size = item.size,
        alpha = 1
      )
    }
    plt <- plt + ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = override,
        ncol = leg.ncol,
        title = leg.ttl,
        keyheight = grid::unit(0.16, "cm"),
        keywidth = grid::unit(0.28, "cm")
      )
    )
  } else {
    plt <- plt + ggplot2::guides(colour = "none")
  }
  # ========================================================================= #
  # Axes
  # ========================================================================= #
  if (!isTRUE(axes)) plt <- plt + Seurat::NoAxes()
  plt
}


# ========================================================================= #
# Feature plots
# ========================================================================= #

#' Resolve GnRHcell feature presets
#'
#' @keywords internal
#' @noRd
.gnrh_features <- function(
    type = c(
      "core",
      "modules",
      "staging",
      "migration",
      "secretory",
      "hits",
      "all"
    )
) {

  type <- match.arg(type)

  presets <- list(

    # Core detection / confidence metrics
    core = c(
      "gnrh_expr",
      "gnrh_score",
      "gnrh_support_score",
      "gnrh_knn"
    ),

    # Biological support programs
    modules = c(
      "gnrh_identity_score",
      "gnrh_migration_score",
      "gnrh_neuro_score",
      "gnrh_alternative_score"
    ),

    # Developmental-stage scores
    staging = c(
      "gnrh_stage_early_score",
      "gnrh_stage_migrating_score",
      "gnrh_stage_mature_score"
    ),

    # Migration evidence
    migration = c(
      "gnrh_migration_score",
      "gnrh_migration_core_hits"
    ),

    # Secretory evidence
    secretory = c(
      "gnrh_secretory_core_hits",
      "gnrh_secretory_supportive_hits",
      "gnrh_secretory_hits"
    ),

    # Discrete program-hit metrics
    hits = c(
      "gnrh_migration_core_hits",
      "gnrh_secretory_core_hits",
      "gnrh_secretory_supportive_hits",
      "gnrh_secretory_hits"
    )
  )

  if (identical(type, "all")) {
    return(
      unique(
        unlist(
          presets,
          use.names = FALSE
        )
      )
    )
  }

  presets[[type]]
}



#' Plot GnRHcell features
#'
#' Visualizes GnRHcell scores, metadata, or genes on a dimensional reduction.
#'
#' @param object A Seurat object.
#' @param features Features or metadata columns.
#' @param preset Optional GnRHcell feature preset.
#' @param cols Optional colour gradient.
#' @param theme.cols Palette name.
#' @param rev.cols Reverse palette.
#' @param na.col Colour for missing or censored values.
#' @param order Plot high values last.
#' @param pt.size Point size.
#' @param txtsize Base text size.
#' @param reduction Dimensional reduction.
#' @param na.cutoff Lower colour-scale cutoff.
#' @param raster Rasterize large plots.
#' @param raster.dpi Raster resolution.
#' @param split.by Optional splitting metadata column.
#' @param ncol Number of plot columns.
#' @param layer Assay layer used when calculating a common legend range.
#' @param label Label clusters.
#' @param axes Show axes.
#' @param combine Combine plots.
#' @param blend Blend exactly two features.
#' @param merge.leg Use a common colour scale and collect legends.
#' @param style Plot theme.
#' @param ... Additional arguments passed to \code{Seurat::FeaturePlot()}.
#'
#' @return A ggplot, list of ggplots, or patchwork object.
#'
#' @export
gnrh_cellfeat <- function(
    object,
    features = NULL,
    preset = NULL,
    cols = NULL,
    theme.cols = "Reds",
    rev.cols = FALSE,
    na.col = "lightgray",
    order = FALSE,
    pt.size = NULL,
    txtsize = getOption("gnrhcell.base_size", 14),
    reduction = NULL,
    na.cutoff = 1e-9,
    raster = NULL,
    raster.dpi = c(2048, 2048),
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
    # ========================================================================= #
    # Checks
    # ========================================================================= #
    .validate_seurat(object)
    if (!is.null(features) && !is.null(preset)) stop("Supply either `features` or `preset`, not both.", call. = FALSE)
    if (!is.null(preset)) features <- .gnrh_features(preset)
    if (is.null(features) || !length(features)) stop("Supply `features` or `preset`.", call. = FALSE)
    features <- intersect(
      unique(as.character(features)),
      c(rownames(object), colnames(object[[]]))
    )
    if (!length(features)) stop("No valid features found.", call. = FALSE)
    reduction <- gnrh_reduction(object, reduction)
    if (isTRUE(blend) && length(features) != 2L) stop("`blend = TRUE` requires exactly two features.", call. = FALSE)
    # ========================================================================= #
    # Colors
    # ========================================================================= #
    if (!is.null(cols)) {
      if (length(cols) < 2L) stop("`cols` must contain at least 2 colors.", call. = FALSE)
      cols <- as.character(cols)
    } else {
      pal <- list(
        default = c("lightgrey", "blue"),
        gnrh = c("#F7F7F7", "#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#B30000", "darkred"),
        hotspot = c("navy", "#2C7BB6", "#27F5EB", "green", "yellow", "orange", "red", "darkred"),
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
        if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package `RColorBrewer` is required.", call. = FALSE)
        if (!theme.cols %in% rownames(RColorBrewer::brewer.pal.info)) stop("Unknown palette: ", theme.cols, call. = FALSE)
        n_pal <- min(9L, RColorBrewer::brewer.pal.info[theme.cols, "maxcolors"])
        cols <- RColorBrewer::brewer.pal(n_pal, theme.cols)
      }
    }
    if (isTRUE(rev.cols)) cols <- rev(cols)
    # ========================================================================= #
    # Point size
    # ========================================================================= #
    raster <- raster %||% (ncol(object) > 2e5)
    if (is.null(pt.size)) {
      pt.size <- if (raster) 1 else min(1583 / ncol(object), 1)
    }
    # ========================================================================= #
    # Global range for merged legend
    # ========================================================================= #
    global_max <- NULL
    if (isTRUE(merge.leg) && is.null(split.by) && length(features) > 1L && !isTRUE(blend)) {
      expr_data <- Seurat::FetchData(
        object,
        vars = features,
        layer = layer
      )
      global_max <- max(as.matrix(expr_data), na.rm = TRUE)
      if (!is.finite(global_max)) global_max <- NULL
    }
    # ========================================================================= #
    # FeaturePlot
    # ========================================================================= #
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
    if (!combine) return(plt)
    # Replace internal metadata names with concise, reader-facing panel titles.
    feature_title <- c(
      gnrh_secretory_core_hits = "Core secretory markers",
      gnrh_secretory_supportive_hits = "Supportive secretory markers",
      gnrh_secretory_hits = "Total secretory markers",
      gnrh_migration_core_hits = "Core migration markers"
    )
    for (i in seq_along(features)) {
      label_i <- unname(feature_title[features[[i]]])
      if (!length(label_i) || is.na(label_i)) {
        label_i <- if (features[[i]] %in% colnames(object[[]])) {
          tools::toTitleCase(gsub("_", " ", sub("^gnrh_", "", features[[i]])))
        } else {
          features[[i]]
        }
      }
      plt[[i]] <- plt[[i]] + ggplot2::labs(title = label_i)
    }
    # ========================================================================= #
    # Palette behavior
    # ========================================================================= #
    continuous_pal <- theme.cols %in% c(
      "gnrh",
      "hotspot",
      "rainbow",
      "magma",
      "inferno",
      "plasma",
      "viridis",
      "blue_red",
      "teal_orange",
      "purple_yellow",
      "brain",
      "mult"
    )
    if (!isTRUE(blend)) {
      limits <- if (!is.null(global_max)) c(na.cutoff, global_max) else c(na.cutoff, NA)
      if (!continuous_pal) {
        plt <- plt &
          ggplot2::scale_color_gradientn(
            colors = cols,
            limits = limits,
            oob = scales::censor,
            na.value = na.col,
            name = NULL
          )
      } else {
        plt <- plt &
          ggplot2::scale_color_gradientn(
            colors = cols,
            limits = limits,
            oob = scales::squish,
            na.value = na.col
          )
      }
    }
    # ========================================================================= #
    # Theme
    # ========================================================================= #
    plt <- plt &
      gnrh_theme(
        style = style,
        txtsize = txtsize
      ) &
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = txtsize + 2,
          face = "bold.italic"
        )
      )
    # ========================================================================= #
    # Colorbar
    # ========================================================================= #
    if (!isTRUE(blend)) {
      plt <- plt &
        ggplot2::guides(
          color = ggplot2::guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            barheight = grid::unit(32, "mm"),
            barwidth = grid::unit(4.5, "mm")
          )
        ) &
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = txtsize),
          legend.key.height = grid::unit(7, "mm")
        )
    }
    # ========================================================================= #
    # Axes
    # ========================================================================= #
    if (!isTRUE(axes)) plt <- plt & Seurat::NoAxes()
    # ========================================================================= #
    # Shared legend
    # ========================================================================= #
    if (isTRUE(merge.leg) && is.null(split.by) && length(features) > 1L && !isTRUE(blend)) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        plt <- plt +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = "right")
      }
    }
    plt
  }))
}



# ========================================================================= #
# Dot plots
# ========================================================================= #

#' Create an enhanced Seurat dot plot
#'
#' Generates a customizable dot plot for gene expression patterns across
#' clusters or metadata-defined groups. This function wraps
#' \code{Seurat::DotPlot()} and adds improved color control, optional dot
#' outlines, flexible axis formatting, legend customization, and gnrhcell
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
#' @param colors Optional custom gradient colours. Overrides `th.cols`.
#' @param n.colors Number of interpolated gradient colours.
#' @param col.min,col.max Limits for scaled average expression.
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
#' @param outline.col,outline.stroke Dot-outline colour and width.
#' @param ... Additional arguments passed to \code{Seurat::DotPlot()} and
#' the internal gnrhcell theme.
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
gnrh_celldot <- function(
    object,
    features,
    group.by = "seurat_clusters",
    th.cols = "Reds",
    rev.th.cols = FALSE,
    colors = NULL,
    n.colors = 100,
    col.min = -2.5,
    col.max = 2.5,
    dot.scale = 6,
    x.ang = 60,
    vjust.x = NULL,
    hjust.x = NULL,
    flip = FALSE,
    txtsize = getOption("gnrhcell.base_size", 14),
    title = NULL,
    leg.size = 10,
    leg.ttl.size = 10,
    leg.pos = "right",
    leg.just = "bottom",
    leg.hjust = FALSE,
    x.axis.pos = "bottom",
    style = "test",
    x.face = FALSE,
    y.face = FALSE,
    x.ttl = FALSE,
    y.ttl = FALSE,
    dot.outline = TRUE,
    outline.col = "gray40",
    outline.stroke = 0.35,
    ...
) {

  .check_seurat(object)

  if (!group.by %in% colnames(object@meta.data)) {
    stop("Grouping column not found: ", group.by, call. = FALSE)
  }

  object <- Seurat::SetIdent(object, value = group.by)

  if (is.list(features)) {
    # Seurat::DotPlot() converts the flattened feature list to a factor and
    # fails when one gene occurs in several named groups. Preserve the first
    # biological assignment and remove subsequent duplicates.
    seen_features <- character()
    features <- lapply(features, function(group_features) {
      group_features <- unique(as.character(group_features))
      group_features <- group_features[!group_features %in% seen_features]
      seen_features <<- c(seen_features, group_features)
      group_features
    })
    features <- features[lengths(features) > 0L]
  } else {
    features <- unique(
      as.character(features)
    )
  }

  if (!length(features)) {
    stop("features must be provided.", call. = FALSE)
  }

  # --------------------------------------------------- #
  # Palette
  # --------------------------------------------------- #

  pal <- .cellplot_gradient(
    palette = th.cols,
    colors = colors,
    n = n.colors,
    reverse = rev.th.cols
  )
  outline_col <- if (dot.outline) outline.col else NA
  outline_stroke <- if (dot.outline) outline.stroke else 0

  # --------------------------------------------------- #
  # Plot
  # --------------------------------------------------- #

  plt <- suppressWarnings(
    suppressMessages(
      Seurat::DotPlot(
        object = object,
        features = features,
        dot.scale = dot.scale,
        col.min = col.min,
        col.max = col.max,
        ...
      )
    )
  )

  # Seurat::DotPlot() installs its own colour scale. Remove it before applying
  # the GnRHcell palette so ggplot2 does not emit a replacement-scale message.
  plt$scales$scales <- Filter(
    function(scale) !any(c("colour", "color") %in% scale$aesthetics),
    plt$scales$scales
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
    gnrh_theme(
      style = style,
      txtsize = txtsize,
      x.ang = x.ang,
      leg.size = leg.size,
      leg.ttl.size = leg.ttl.size,
      leg.pos = leg.pos,
      hjust = hjust.x,
      vjust = vjust.x,
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

  # --------------------------------------------------- #
  # Legend positioning
  # --------------------------------------------------- #

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

  # --------------------------------------------------- #
  # Guides
  # --------------------------------------------------- #

  guide_color <- ggplot2::guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black"
  )

  guide_size <- ggplot2::guide_legend(
    override.aes = list(
      shape = 21,
      colour = outline_col,
      fill = "gray40"
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
