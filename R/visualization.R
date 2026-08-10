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
.gnrh_features <- function(type = c("core", "modules", "staging", "all")) {
  switch(
    match.arg(type),

    core = c(
      "gnrh_raw", "gnrh_expr",
      "gnrh_score", "gnrh_support_score"
    ),

    modules = c(
      "gnrh_core_hits", "gnrh_mig_hits",
      "gnrh_neuro_hits", "gnrh_knn"
    ),

    staging = c(
      "gnrh_identity_score",
      "gnrh_migration_score",
      "gnrh_neuro_score",
      "gnrh_hormone_score"
    ),

    all = c(
      "gnrh_raw", "gnrh_expr", "gnrh_score", "gnrh_support_score",
      "gnrh_core_hits", "gnrh_mig_hits", "gnrh_neuro_hits", "gnrh_knn",
      "gnrh_identity_score", "gnrh_migration_score",
      "gnrh_neuro_score", "gnrh_hormone_score"
    )
  )
}


#' Internal feature palette
#' @keywords internal
#' @noRd
.gnrh_feature_palette <- function(palette = "gnrh") {
  pals <- list(
    gnrh = c("#F2F2F2", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),
    red = c("#F7F7F7", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"),
    blue = c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"),
    purple = c("#FCFBFD", "#DADAEB", "#9E9AC8", "#6A51A3", "#3F007D"),
    viridis = c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725"),
    magma = c("#000004", "#51127C", "#B73779", "#FC8961", "#FCFDBF"),
    support = c("#F2F2F2", "#B2DF8A", "#33A02C", "#006D2C"),
    alternative = c("#F2F2F2", "#CAB2D6", "#984EA3", "#54278F")
  )

  if (!palette %in% names(pals))
    stop("Unknown palette: ", palette, call. = FALSE)

  pals[[palette]]
}


#' Plot GnRH features on an embedding
#' @export
plot_gnrh_feature <- function(
    object,
    features,
    cols = NULL,
    palette = "gnrh",
    rev.cols = FALSE,
    na.col = "grey90",
    order = FALSE,
    pt.size = NULL,
    reduction = NULL,
    dims = c(1, 2),
    min.cutoff = NULL,
    max.cutoff = NULL,
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
    legend = TRUE,
    plot.ttl = NULL,
    leg.ttl = NULL,
    leg.pos = "right",
    leg.size = 9,
    txtsize = 10,
    dark = FALSE,
    style = "classic",
    ...
) {
  .validate_seurat(object)

  if (!length(features))
    stop("`features` must contain at least one feature.", call. = FALSE)

  reduction <- resolve_reduction(object, reduction)
  emb <- Seurat::Embeddings(object, reduction)

  if (length(dims) != 2L || max(dims) > ncol(emb))
    stop("`dims` must contain two available dimensions.", call. = FALSE)

  defaults <- .plot_defaults(object, raster, pt.size)
  raster <- defaults$raster
  pt.size <- defaults$pt.size

  cols <- cols %||% .gnrh_feature_palette(palette)
  if (rev.cols) cols <- rev(cols)

  args <- list(
    object = object,
    features = features,
    reduction = reduction,
    dims = dims,
    order = order,
    pt.size = pt.size,
    raster = raster,
    raster.dpi = raster.dpi,
    split.by = split.by,
    combine = FALSE,
    blend = blend,
    label = label
  )

  if (!is.null(min.cutoff)) args$min.cutoff <- min.cutoff
  if (!is.null(max.cutoff)) args$max.cutoff <- max.cutoff

  plots <- suppressWarnings(
    suppressMessages(
      do.call(Seurat::FeaturePlot, c(args, list(...)))
    )
  )

  plots <- lapply(plots, function(p) {
    if (!blend)
      p <- p + ggplot2::scale_colour_gradientn(
        colours = cols,
        na.value = na.col,
        name = leg.ttl,
        oob = scales::squish
      )

    p +
      ggplot2::coord_equal() +
      .gnrh_theme(
        txtsize = txtsize,
        leg.pos = if (legend) leg.pos else "none",
        axes = axes,
        style = style,
        mode = if (dark) "dark" else "light"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold.italic",
          size = txtsize + 1
        ),
        legend.title = ggplot2::element_text(
          face = "bold",
          size = leg.size
        ),
        legend.text = ggplot2::element_text(size = leg.size)
      )
  })

  if (!combine) return(plots)

  p <- patchwork::wrap_plots(plots, ncol = ncol)

  if (merge.leg && is.null(split.by) && length(plots) > 1L)
    p <- p + patchwork::plot_layout(guides = "collect")

  if (!is.null(plot.ttl) && length(plots) > 1L)
    p <- p +
    patchwork::plot_annotation(
      title = plot.ttl,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = txtsize + 3
        )
      )
    )

  p
}


# ============================================================================= #
# Distribution plots
# ============================================================================= #

#' Plot GnRHcell metadata distributions
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

#' Generate a compact GnRH diagnostic report
#'
#' @param object Seurat object processed with GnRHcell.
#' @param style Plot theme style.
#' @param truth Optional metadata column containing independent binary truth.
#' @param verbose Print progress messages.
#'
#' @return A patchwork diagnostic dashboard.
#' @export
gnrh_report <- function(object, style = "test", truth = NULL, verbose = TRUE) {
  log <- .msg(verbose)
  log("Generating GnRH QC report")

  .validate_seurat(object)

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required.", call. = FALSE)

  g <- object@misc$gnrh
  if (is.null(g) || is.null(g$diagnostics))
    stop(
      "Missing GnRH diagnostics. Run detect_gnrh() and gnrh_diagnostics().",
      call. = FALSE
    )

  diag <- g$diagnostics
  md <- object[[]]
  curve <- g$threshold_curve %||% NULL

  add_meta <- c(
    "gnrh_class", "gnrh_confident", "gnrh_stage", "gnrh_secretory",
    "gnrh_support_score", "gnrh_alternative_score",
    "gnrh_identity_score", "gnrh_migration_score",
    "gnrh_neuro_score", "gnrh_hormone_score", "gnrh_knn"
  )

  for (x in intersect(add_meta, names(md)))
    if (!x %in% names(diag)) diag[[x]] <- md[[x]]

  req <- c("expr", "score", "status", "core_hits", "mig_hits", "neuro_hits")
  miss <- setdiff(req, names(diag))

  if (length(miss))
    stop(
      "Missing diagnostic columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )

  .empty <- function(title, subtitle = NULL)
    ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 8
      )
    )

  .roc_add <- function(ref, score, label) {
    ref <- as.integer(ref)
    ok <- !is.na(ref) & !is.na(score) & is.finite(score)

    if (sum(ok) < 2L || length(unique(ref[ok])) != 2L)
      return(NULL)

    roc <- pROC::roc(
      response = ref[ok],
      predictor = score[ok],
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    data.frame(
      FPR = 1 - roc$specificities,
      TPR = roc$sensitivities,
      group = label,
      auc = as.numeric(pROC::auc(roc))
    )
  }

  # ------------------------------------------------------------------------- #
  # 1. Detection landscape
  # ------------------------------------------------------------------------- #

  status_cols <- gnrh_colors("status")
  status_labels <- .gnrh_legend_labels(
    diag$status,
    levels = names(status_cols),
    percentage = TRUE
  )

  p1 <- ggplot2::ggplot(
    diag,
    ggplot2::aes(expr, score, colour = status)
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 0.65) +
    ggplot2::scale_colour_manual(
      values = status_cols,
      breaks = names(status_labels),
      labels = status_labels
    ) +
    ggplot2::labs(
      title = "GnRH detection landscape",
      x = "GNRH1 expression",
      y = "GnRH composite score",
      colour = "Status"
    ) +
    .gnrh_theme(leg.pos = "right") +
    .point_legend()

  # ------------------------------------------------------------------------- #
  # 2. Developmental stage landscape
  # ------------------------------------------------------------------------- #

  if ("gnrh_stage" %in% names(diag)) {
    stage_cols <- gnrh_colors("stage")
    stage_labels <- .gnrh_legend_labels(
      diag$gnrh_stage,
      levels = names(stage_cols),
      percentage = TRUE
    )

    p2 <- ggplot2::ggplot(
      diag,
      ggplot2::aes(expr, score, colour = gnrh_stage)
    ) +
      ggplot2::geom_point(alpha = 0.5, size = 0.65) +
      ggplot2::scale_colour_manual(
        values = stage_cols,
        breaks = names(stage_labels),
        labels = stage_labels
      ) +
      ggplot2::labs(
        title = "Developmental stage landscape",
        x = "GNRH1 expression",
        y = "GnRH composite score",
        colour = "Stage"
      ) +
      .gnrh_theme(leg.pos = "right") +
      .point_legend()
  } else {
    p2 <- .empty("Stage unavailable")
  }

  # ------------------------------------------------------------------------- #
  # 3. Score distribution
  # ------------------------------------------------------------------------- #

  p3 <- ggplot2::ggplot(
    diag,
    ggplot2::aes(score, fill = status)
  ) +
    ggplot2::geom_density(alpha = 0.45, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = status_cols) +
    ggplot2::labs(
      title = "GnRH score distribution",
      x = "GnRH composite score",
      y = "Density",
      fill = "Status"
    ) +
    .gnrh_theme(leg.pos = "right")

  # ------------------------------------------------------------------------- #
  # 4. Threshold performance
  # ------------------------------------------------------------------------- #

  if (
    !is.null(curve) &&
    all(c("threshold", "sensitivity", "specificity", "F1") %in% names(curve))
  ) {
    ok <- is.finite(curve$threshold) & is.finite(curve$F1)
    best <- if (any(ok))
      curve$threshold[which.max(replace(curve$F1, !ok, -Inf))]
    else NA_real_

    p4 <- ggplot2::ggplot(curve, ggplot2::aes(threshold)) +
      ggplot2::geom_line(
        ggplot2::aes(y = sensitivity, colour = "Sensitivity"),
        linewidth = 0.6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = specificity, colour = "Specificity"),
        linewidth = 0.6
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = F1, colour = "F1"),
        linewidth = 0.7
      ) +
      ggplot2::geom_vline(
        xintercept = best,
        linetype = "dashed",
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          F1 = "#EF476F",
          Sensitivity = "#3A86FF",
          Specificity = "#06D6A0"
        )
      ) +
      ggplot2::labs(
        title = "Threshold performance",
        subtitle = if (is.finite(best))
          paste0("Optimal threshold = ", round(best, 3))
        else NULL,
        x = "Threshold",
        y = "Performance",
        colour = NULL
      ) +
      .gnrh_theme(leg.pos = "right") +
      .line_legend()
  } else {
    p4 <- .empty("Threshold curve unavailable")
  }

  # ------------------------------------------------------------------------- #
  # 5. ROC performance
  # ------------------------------------------------------------------------- #

  p5 <- NULL

  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_list <- list()

    if (!is.null(truth)) {
      if (!truth %in% names(md))
        stop(
          "Truth column `", truth,
          "` not found in object metadata.",
          call. = FALSE
        )

      ref <- as.character(md[[truth]])

      ref[
        ref %in% c(
          "1", "TRUE", "true", "positive", "Positive", "pos", "Pos"
        )
      ] <- "pos"

      ref[
        ref %in% c(
          "0", "FALSE", "false", "negative", "Negative", "neg", "Neg"
        )
      ] <- "neg"

      roc_list$status <- .roc_add(
        ref == "pos",
        diag$score,
        "GnRH status"
      )
    } else {
      roc_list$status <- .roc_add(
        diag$status == "pos",
        diag$score,
        "GnRH status"
      )
    }

    if ("gnrh_stage" %in% names(diag)) {
      stage_scores <- c(
        identity = "gnrh_identity_score",
        migrating = "gnrh_migration_score",
        mature = "gnrh_neuro_score",
        secreting = "gnrh_hormone_score"
      )

      for (s in names(stage_scores)) {
        sc <- stage_scores[[s]]

        if (
          s %in% as.character(diag$gnrh_stage) &&
          sc %in% names(diag)
        ) {
          roc_list[[s]] <- .roc_add(
            diag$gnrh_stage == s,
            diag[[sc]],
            s
          )
        }
      }
    }

    roc_list <- Filter(Negate(is.null), roc_list)

    if (length(roc_list)) {
      roc_df <- do.call(rbind, roc_list)
      auc_df <- unique(roc_df[c("group", "auc")])

      auc_df$label <- paste0(
        auc_df$group,
        " (AUC = ",
        sprintf("%.3f", auc_df$auc),
        ")"
      )

      roc_cols <- c(
        "GnRH status" = gnrh_colors("status")[["pos"]],
        gnrh_colors("stage")
      )

      missing <- setdiff(auc_df$group, names(roc_cols))
      if (length(missing)) {
        extra <- grDevices::hcl.colors(length(missing), "Dark 3")
        names(extra) <- missing
        roc_cols <- c(roc_cols, extra)
      }

      roc_df$group <- factor(
        roc_df$group,
        levels = auc_df$group
      )

      p5 <- ggplot2::ggplot(
        roc_df,
        ggplot2::aes(FPR, TPR, colour = group)
      ) +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          colour = "grey65",
          linewidth = 0.4
        ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::scale_colour_manual(
          values = roc_cols,
          breaks = auc_df$group,
          labels = stats::setNames(
            auc_df$label,
            auc_df$group
          )
        ) +
        ggplot2::coord_equal(
          xlim = c(0, 1),
          ylim = c(0, 1)
        ) +
        ggplot2::labs(
          title = "ROC performance",
          x = "False positive rate",
          y = "True positive rate",
          colour = NULL
        ) +
        .gnrh_theme(leg.pos = "right") +
        .line_legend()
    }
  }

  if (is.null(p5))
    p5 <- .empty(
      "ROC unavailable",
      if (!requireNamespace("pROC", quietly = TRUE))
        "Install pROC to enable ROC analysis."
      else NULL
    )

  # ------------------------------------------------------------------------- #
  # 6. Marker program support
  # ------------------------------------------------------------------------- #

  module_df <- data.frame(
    status = rep(diag$status, 3),
    module = rep(
      c("Core", "Migration", "Neuroendocrine"),
      each = nrow(diag)
    ),
    hits = c(
      diag$core_hits,
      diag$mig_hits,
      diag$neuro_hits
    )
  )

  module_sum <- stats::aggregate(
    hits ~ status + module,
    module_df,
    mean,
    na.rm = TRUE
  )

  module_sum$module <- factor(
    module_sum$module,
    levels = c("Core", "Migration", "Neuroendocrine")
  )
  module_sum$status <- factor(
    module_sum$status,
    levels = c("neg", "pos")
  )

  p6 <- ggplot2::ggplot(
    module_sum,
    ggplot2::aes(module, status, fill = hits)
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(hits, 2)),
      size = 3
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = gnrh_colors("status")[["pos"]]
    ) +
    ggplot2::labs(
      title = "Marker program support",
      x = NULL,
      y = NULL,
      fill = "Mean hits"
    ) +
    .gnrh_theme()

  # ------------------------------------------------------------------------- #
  # Report summary
  # ------------------------------------------------------------------------- #

  n_pos <- sum(diag$status == "pos", na.rm = TRUE)

  n_conf <- if ("gnrh_confident" %in% names(diag))
    sum(diag$gnrh_confident %in% c(TRUE, "TRUE"), na.rm = TRUE)
  else NA_integer_

  subtitle <- paste0(
    "Cells: ", format(ncol(object), big.mark = ","),
    " | Genes: ", format(nrow(object), big.mark = ","),
    " | GnRH+: ", format(n_pos, big.mark = ","),
    if (!is.na(n_conf))
      paste0(" | Confident: ", format(n_conf, big.mark = ","))
    else ""
  )

  (p1 | p2 | p3) /
    (p4 | p5 | p6) +
    patchwork::plot_annotation(
      title = "GnRHcell Diagnostic Report",
      subtitle = subtitle
    )
}


# ============================================================================= #
# Coexpression network
# ============================================================================= #

#' GnRH gene similarity network
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
#' @export
plot_gnrh_coexpr <- function(
    df,
    coexp_cutoff = 0.25,
    txtsize = 10
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
    .gnrh_theme(txtsize = txtsize, style = "bw") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic")
    )
}


# ============================================================================= #
# Runtime
# ============================================================================= #

#' Plot GnRHcell runtime across datasets
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
#' @export
plot_gnrh_marker_programs <- function(
    programs,
    table = c("summary", "high_confidence", "candidate_table"),
    type = c("bar", "dot", "tile"),
    min_genes = 1,
    mode = "light",
    txtsize = 12,
    x.ang = 45
) {
  table <- match.arg(table)
  type <- match.arg(type)

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
            x.ang = x.ang
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
            x.ang = x.ang
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
          x.ang = x.ang
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
      x.ang = x.ang
    )
}


# ============================================================================= #
# Simple diagnostic plots
# ============================================================================= #

#' Plot GnRH classification counts
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

