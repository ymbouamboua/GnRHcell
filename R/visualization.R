# ========================================================================= #
# gnrhcell visualization utilities
# ========================================================================= #


# ========================================================================= #
# Internal utilities
# ========================================================================= #

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


# ========================================================================= #
# GnRHcell theme
# ========================================================================= #
#' GnRHcell plot theme
#' Apply a publication-ready ggplot2 theme
#'
#' Create a configurable ggplot2 theme with predefined visual styles, axis
#' controls, legend formatting, facet styling, grids, and light/dark modes.
#'
#' @param style Theme preset: `"classic"`, `"minimal"`, `"bw"`, `"test"`
#' @param txtsize Base text size.
#' @param xy.val Display axis values.
#' @param x.ang X-axis label angle.
#' @param hjust,vjust Optional x-axis label justification.
#' @param xlab,ylab Display x- and y-axis text/ticks.
#' @param xy.lab Display axis text/ticks on both axes.
#' @param facet.face Facet-label font face.
#' @param ttl.face Plot-title font face.
#' @param txt.face General text font face.
#' @param ttl.pos Plot-title position.
#' @param x.ttl,y.ttl Display x- and y-axis titles.
#' @param ticks Optional axis-tick control.
#' @param line Optional axis-line control.
#' @param border Optional panel-border control.
#' @param grid.major,grid.minor Optional grid controls.
#' @param panel.fill Panel background colour in light mode.
#' @param facet.bg Display facet-strip backgrounds.
#' @param mode Display mode: `"light"` or `"dark"`.
#' @param leg.pos Legend position.
#' @param leg.dir Legend direction.
#' @param leg.size Legend-text size.
#' @param leg.ttl Deprecated compatibility argument.
#' @param leg.ttl.size Legend-title size.
#' @param leg.just Legend justification.
#' @param leg.ttl.text Optional legend title.
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme object.
#'
#' @export
gnrh_theme <- function(
    style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
    txtsize = 12,
    xy.val = TRUE,
    x.ang = 0,
    hjust = NULL,
    vjust = NULL,
    xlab = TRUE,
    ylab = TRUE,
    xy.lab = TRUE,
    facet.face = "bold",
    ttl.face = "bold",
    txt.face = c("plain", "italic", "bold"),
    ttl.pos = c("center", "left", "right"),
    x.ttl = TRUE,
    y.ttl = TRUE,
    ticks = NULL,
    line = NULL,
    border = NULL,
    grid.major = NULL,
    grid.minor = NULL,
    panel.fill = "white",
    facet.bg = TRUE,
    mode = c("light", "dark"),
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = 10,
    leg.ttl = 10,
    leg.ttl.size = 10,
    leg.just = "center",
    leg.ttl.text = NULL,
    ...
) {
  # ========================================================================= #
  # Arguments
  # ========================================================================= #
  style <- match.arg(style)
  ttl.pos <- match.arg(ttl.pos)
  txt.face <- match.arg(txt.face)
  mode <- match.arg(mode)
  lw <- 0.3
  if (is.null(line)) line <- style == "classic"
  fg <- if (mode == "dark") "white" else "#1A1A1A"
  bg <- if (mode == "dark") "#111111" else panel.fill
  col.grid <- if (mode == "dark") "#444444" else "#D9D9D9"
  col.strip <- if (mode == "dark") "#383838" else "#EFEFEF"
  if (is.null(hjust) || is.null(vjust)) {
    pos <- switch(
      as.character(x.ang),
      `0` = c(0.5, 0.5),
      `45` = c(1, 1),
      `90` = c(1, 0.5),
      `270` = c(0, 0.5),
      c(1, 1)
    )
    if (is.null(hjust)) hjust <- pos[1]
    if (is.null(vjust)) vjust <- pos[2]
  }
  ttl.hjust <- switch(ttl.pos, left = 0, center = 0.5, right = 1)
  # ========================================================================= #
  # Preset
  # ========================================================================= #
  preset <- switch(
    style,
    classic = ggplot2::theme_classic(base_size = txtsize),
    minimal = ggplot2::theme_minimal(base_size = txtsize),
    bw = ggplot2::theme_bw(base_size = txtsize),
    test = ggplot2::theme_test(base_size = txtsize),
    void = ggplot2::theme_void(base_size = txtsize),
    dirty = ggplot2::theme_minimal(base_size = txtsize),
    gray = ggplot2::theme_gray(base_size = txtsize)
  )
  # ========================================================================= #
  # Base theme
  # ========================================================================= #
  th <- preset + ggplot2::theme(
    text = ggplot2::element_text(colour = fg, size = txtsize, face = txt.face, family = "Helvetica"),
    axis.text.x = ggplot2::element_text(colour = fg, size = txtsize, angle = x.ang, hjust = hjust, vjust = vjust),
    axis.text.y = ggplot2::element_text(colour = fg, size = txtsize),
    axis.title = ggplot2::element_text(colour = fg, size = txtsize),
    plot.title = ggplot2::element_text(hjust = ttl.hjust, face = ttl.face, size = txtsize + 2, colour = fg),
    plot.subtitle = ggplot2::element_text(colour = fg, size = max(7, txtsize - 1)),
    strip.text = ggplot2::element_text(face = facet.face, colour = fg),
    strip.background = ggplot2::element_rect(fill = col.strip, colour = NA),
    panel.background = ggplot2::element_rect(fill = bg, colour = NA),
    plot.background = ggplot2::element_rect(fill = bg, colour = NA),
    legend.title = ggplot2::element_text(size = leg.ttl.size + 2, face = "bold", colour = fg),
    legend.text = ggplot2::element_text(size = leg.size, colour = fg),
    legend.position = leg.pos,
    legend.direction = leg.dir,
    legend.justification = leg.just,
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.background = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    ...
  )
  # ========================================================================= #
  # Style-specific settings
  # ========================================================================= #
  if (style %in% c("bw", "test", "gray")) {
    th <- th + ggplot2::theme(
      panel.border = ggplot2::element_rect(linewidth = lw, colour = fg, fill = NA),
      axis.line = ggplot2::element_blank()
    )
  }
  if (isTRUE(line) && style == "classic") {
    th <- th + ggplot2::theme(
      axis.line.x = ggplot2::element_line(colour = fg, linewidth = lw),
      axis.line.y = ggplot2::element_line(colour = fg, linewidth = lw)
    )
  } else {
    th <- th + ggplot2::theme(axis.line = ggplot2::element_blank())
  }
  if (style == "gray") {
    th <- th + ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = if (mode == "dark") bg else "#EDEDED", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = col.grid, linewidth = lw),
      panel.grid.minor = ggplot2::element_line(colour = col.grid, linewidth = lw / 2)
    )
  }
  if (style == "dirty") {
    th <- th + ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  }
  if (style == "void") {
      th <- th + ggplot2::theme(
        axis.text=ggplot2::element_blank(),
        axis.title=ggplot2::element_blank(),
        axis.ticks=ggplot2::element_blank(),
        axis.ticks.length=grid::unit(0,"pt"),
        axis.line=ggplot2::element_blank(),
        panel.grid.major=ggplot2::element_blank(),
        panel.grid.minor=ggplot2::element_blank(),
        panel.border=ggplot2::element_blank(),
        strip.background=ggplot2::element_blank()
      )
    }
  # ========================================================================= #
  # Axis controls
  # ========================================================================= #
  if (!isTRUE(xy.val) || !isTRUE(xy.lab)) {
    th <- th + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  if (!isTRUE(xlab)) {
    th <- th + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  if (!isTRUE(ylab)) {
    th <- th + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  if (!isTRUE(x.ttl)) th <- th + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  if (!isTRUE(y.ttl)) th <- th + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  if (!is.null(ticks)) {
    th <- th + ggplot2::theme(
      axis.ticks = if (isTRUE(ticks)) ggplot2::element_line(colour = fg, linewidth = lw) else ggplot2::element_blank()
    )
  }
  # ========================================================================= #
  # Test style
  # ========================================================================= #
  if (style == "test") {
    th <- th + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = if (isTRUE(x.ttl)) ggplot2::element_text(colour = fg, size = txtsize) else ggplot2::element_blank(),
      axis.title.y = if (isTRUE(y.ttl)) ggplot2::element_text(colour = fg, size = txtsize) else ggplot2::element_blank()
    )
  }
  # ========================================================================= #
  # Panel controls
  # ========================================================================= #
  if (!is.null(border)) {
    th <- th + ggplot2::theme(
      panel.border = if (isTRUE(border)) ggplot2::element_rect(colour = col.grid, fill = NA, linewidth = lw) else ggplot2::element_blank()
    )
  }
  if (!is.null(grid.major)) {
    th <- th + ggplot2::theme(
      panel.grid.major = if (isTRUE(grid.major)) ggplot2::element_line(colour = col.grid, linewidth = lw) else ggplot2::element_blank()
    )
  }
  if (!is.null(grid.minor)) {
    th <- th + ggplot2::theme(
      panel.grid.minor = if (isTRUE(grid.minor)) ggplot2::element_line(colour = col.grid, linewidth = lw / 2) else ggplot2::element_blank()
    )
  }
  if (!isTRUE(facet.bg)) th <- th + ggplot2::theme(strip.background = ggplot2::element_blank())
  # ========================================================================= #
  # Legend title
  # ========================================================================= #
  if (!is.null(leg.ttl.text)) {
    th <- th + ggplot2::labs(
      colour = leg.ttl.text,
      color = leg.ttl.text,
      fill = leg.ttl.text
    )
  }
  th
}




# ========================================================================= #
# Palettes
# ========================================================================= #

#' Color palettes
#'
#' @param type Palette type.
#'
#' @return Named character vector.
#' @export
gnrh_colors <- function(
    type = c(
      "status",
      "confident",
      "class",
      "stage",
      "secretory"
    )
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
      supported = "#4EA8DE",
      direct = "#E63946"
    ),

    stage = c(
      "non-gnrh" = "#B0B0B0",
      identity = "#1F78B4",
      migrating = "#6A3D9A",
      mature = "#E67E22"
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


# ========================================================================= #
# Legend helpers
# ========================================================================= #

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


# ========================================================================= #
# Reduction resolver
# ========================================================================= #

#' Resolve an available dimensional reduction
#'
#' Selects a dimensional reduction from a Seurat object. If the requested
#' reduction is unavailable, the function searches a predefined fallback list,
#' then any reduction containing `"umap"`, and finally returns the first
#' available reduction.
#'
#' @param object A Seurat object containing dimensional reductions.
#' @param reduction Optional character string specifying the preferred
#'   dimensional reduction. If `NULL` or unavailable, a fallback is selected.
#' @param fallback Character vector defining the preferred fallback order.
#'   Defaults to `"umap"`, `"umap_scvi"`, `"umap.harmony"`, `"umap.rpca"`,
#'   and `"pca"`.
#'
#' @return A character string giving the name of the selected dimensional
#'   reduction.
#'
#' @details
#' Reduction selection follows this priority:
#' \enumerate{
#'   \item The explicitly requested `reduction`, when available.
#'   \item The first available reduction listed in `fallback`.
#'   \item The first available reduction whose name contains `"umap"`.
#'   \item The first dimensional reduction stored in the object.
#' }
#' An error is raised when the object contains no dimensional reductions.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' resolve_reduction(object)
#' resolve_reduction(object, "umap_scvi")
#' resolve_reduction(object, fallback = c("umap", "pca"))
#' }
resolve_reduction <- function(
    object,
    reduction = NULL,
    fallback = c("umap","umap_scvi","umap.harmony","umap.rpca","pca")
) {
  .validate_seurat(object)
  available <- names(object@reductions)
  if (!length(available)) {
    stop("No dimensional reductions are available in `object`.",call.=FALSE)
  }
  if (!is.null(reduction) && length(reduction) == 1L && !is.na(reduction) && nzchar(reduction)) {
    if (reduction %in% available) return(reduction)
    warning(sprintf("Reduction `%s` not found. Using an available fallback.",reduction),call.=FALSE)
  }
  candidate <- fallback[fallback %in% available]
  if (length(candidate)) return(candidate[[1L]])
  candidate <- available[grepl("umap",available,ignore.case=TRUE)]
  if (length(candidate)) return(candidate[[1L]])
  available[[1L]]
}

# ========================================================================= #
# Legend helpers
# ========================================================================= #
#' Legend label
#'
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
  if (!isTRUE(n_cells)) {
    return(stats::setNames(levels, levels))
  }
  n_fmt <- format(
    n,
    big.mark = ",",
    trim = TRUE,
    scientific = FALSE
  )
  labels <- if (isTRUE(percentage)) {
    paste0(
      levels,
      " (",
      n_fmt,
      " | ",
      round(100 * n / sum(tab), 1),
      "%)"
    )
  } else {
    paste0(
      levels,
      " (",
      n_fmt,
      ")"
    )
  }
  stats::setNames(labels, levels)
}
# ========================================================================= #
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
plot_gnrh_embedding <- function(
    object,
    group_by = "gnrh_status",
    split_by = NULL,
    reduction = NULL,
    dims = c(1, 2),
    shuffle = FALSE,
    raster = NULL,
    raster.dpi = 600,
    alpha = 1,
    background_alpha = 1,
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
    leg.size = 10,
    leg.ncol = NULL,
    item.border = TRUE,
    txtsize = 12,
    pt.size = NULL,
    dark = FALSE,
    total.cells = FALSE,
    style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
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
      plot_gnrh_embedding(
        object = object,
        group_by = g,
        split_by = split_by,
        reduction = reduction,
        dims = dims,
        shuffle = shuffle,
        raster = raster,
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
  reduction <- resolve_reduction(object, reduction)
  emb <- Seurat::Embeddings(object, reduction = reduction)
  if (max(dims) > ncol(emb)) stop("Selected dimensions exceed those available in `", reduction, "`.", call. = FALSE)
  # ========================================================================= #
  # Plot defaults
  # ========================================================================= #
  defaults <- .plot_defaults(object, raster = raster, pt.size = pt.size)
  raster <- defaults$raster
  pt.size <- defaults$pt.size
  # ========================================================================= #
  # Levels and colours
  # ========================================================================= #
  orders <- list(
    gnrh_status = c("neg", "pos"),
    gnrh_class = c("neg", "supported", "direct"),
    gnrh_confident = c("FALSE", "TRUE"),
    gnrh_stage = c("non-gnrh", "identity", "migrating", "mature"),
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
    ggplot2::coord_equal() +
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
      "gnrh_stage_identity_score",
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
plot_gnrh_feature <- function(
    object,
    features = NULL,
    preset = NULL,
    cols = NULL,
    theme.cols = "Reds",
    rev.cols = FALSE,
    na.col = "lightgray",
    order = FALSE,
    pt.size = NULL,
    txtsize = 10,
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
    reduction <- resolve_reduction(object, reduction)
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
            ticks.colour = "black"
          )
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
plot_gnrh_dot <- function(
    object,
    features,
    group.by = "seurat_clusters",
    th.cols = "RdYlBu",
    rev.th.cols = TRUE,
    dot.scale = 6,
    x.ang = 45,
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

  if (is.list(features)) {
    features <- lapply(
      features,
      unique
    )
  } else {
    features <- unique(
      as.character(features)
    )
  }

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
  outline_stroke <- if (dot.outline) 0.1 else 0

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
    gnrh_theme(
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




# ========================================================================= #
# Distribution plots
# ========================================================================= #

#' Plot gnrhcell metadata distributions
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
#' @param x.ang X-axis label angle. When `NULL`, the angle is selected from
#' @param style Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.
#'   the number and length of sample labels.
#' @param flip Flip coordinates.
#' @param adaptive Automatically use compact spacing, an economical legend
#'   layout, and recommended export dimensions based on the number of samples.
#' @param bar.width Bar width. When `NULL`, it is selected automatically.
#' @param bar.gap Gap between adjacent bar edges, in x-axis units. It is
#'   independent of `bar.width`; for example, `bar.width = 0.3` and
#'   `bar.gap = 0.2` produce centres 0.5 units apart. When `NULL`, an adaptive
#'   value is used.
#' @param legend.position Legend position. Use `"auto"` to place it according
#'   to the number of samples, or a standard ggplot2 legend position.
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
    x.ang = NULL,
    style = "classic",
    flip = FALSE,
    adaptive = TRUE,
    bar.width = NULL,
    bar.gap = NULL,
    legend.position = "auto"
) {
  .validate_seurat(object)

  style = match.arg(style)

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

  sample_labels <- unique(as.character(df$split))
  n_samples <- length(sample_labels)
  longest_label <- max(nchar(sample_labels), 1L)

  if (is.null(bar.width)) {
    bar.width <- if (isTRUE(adaptive)) {
      if (n_samples <= 4L) 0.58 else if (n_samples <= 10L) 0.68 else 0.78
    } else {
      0.75
    }
  }

  if (is.null(bar.gap)) bar.gap <- if (isTRUE(adaptive)) 0.12 else 0.25
  if (length(bar.gap) != 1L || !is.finite(bar.gap) || bar.gap < 0) {
    stop("`bar.gap` must be one non-negative number.", call. = FALSE)
  }

  split_levels <- unique(as.character(df$split))
  spacing <- bar.width + bar.gap
  split_positions <- stats::setNames(
    (seq_along(split_levels) - 1) * spacing + 1,
    split_levels
  )
  df$.split_position <- unname(split_positions[as.character(df$split)])
  outer_gap <- bar.gap / 2
  x_limits <- c(
    min(split_positions) - bar.width / 2 - outer_gap,
    max(split_positions) + bar.width / 2 + outer_gap
  )

  if (is.null(x.ang)) {
    x.ang <- if (isTRUE(flip)) {
      0
    } else if (n_samples <= 4L && longest_label <= 8L) {
      0
    } else if (n_samples <= 10L && longest_label <= 15L) {
      30
    } else {
      60
    }
  }

  if (identical(legend.position, "auto")) {
    legend.position <- if (is.null(split.by)) {
      "none"
    } else if (isTRUE(adaptive) && n_samples <= 8L) {
      "bottom"
    } else {
      "right"
    }
  }

  plot_width <- if (isTRUE(flip)) {
    max(5.5, min(9, 5 + longest_label / 10))
  } else {
    max(4.8, min(14, 3.4 + 0.62 * n_samples + longest_label / 30))
  }
  plot_height <- if (isTRUE(flip)) {
    max(4.2, min(12, 2.8 + 0.42 * n_samples))
  } else {
    if (identical(legend.position, "bottom")) 5.2 else 4.8
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$.split_position, .data$value, fill = .data$group)
  ) +
    ggplot2::geom_col(
      position = position,
      width = bar.width,
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_x_continuous(
      breaks = unname(split_positions),
      labels = names(split_positions),
      limits = x_limits,
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::labs(
      title = plot.ttl,
      x = if (is.null(split.by)) NULL else split.by,
      y = if (proportion) "Proportion" else "Cells",
      fill = group.by
    )

  if (isTRUE(proportion)) {
    p <- p + ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = function(x) ifelse(x == 1, "100", sprintf("%.2f", x)),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.015))
    )
  }

  if (label) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        position = if (position == "stack")
          ggplot2::position_stack(vjust = 0.5)
        else
          ggplot2::position_dodge(width = bar.width),
        size = label.size
      )
  }

  if (flip) p <- p + ggplot2::coord_flip()

  p <- p +
    gnrh_theme(
      txtsize = txtsize,
      x.ang = x.ang,
      leg.pos = legend.position,
      style = style
    ) +
    ggplot2::theme(
      legend.direction = if (identical(legend.position, "bottom")) "horizontal" else "vertical",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  attr(p, "recommended_size") <- c(width = plot_width, height = plot_height)
  attr(p, "n_samples") <- n_samples
  p
}


# ========================================================================= #
# Module hits
# ========================================================================= #

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

  p + gnrh_theme(txtsize = txtsize, x.ang = 45)
}


# ========================================================================= #
# Compact diagnostic report
# ========================================================================= #

#' Generate a publication-ready gnrhcell diagnostic report
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
    "gnrh_class",
    "gnrh_confident",

    "gnrh_stage",
    "gnrh_secretory",

    "gnrh_support_score",
    "gnrh_support_score_raw",

    "gnrh_alternative_score",

    "gnrh_identity_score",
    "gnrh_migration_score",
    "gnrh_neuro_score",

    "gnrh_stage_identity_score",
    "gnrh_stage_migrating_score",
    "gnrh_stage_mature_score",

    "gnrh_secretory_core_hits",
    "gnrh_secretory_supportive_hits",
    "gnrh_secretory_hits",

    "gnrh_knn",

    "gnrh_direct_signal",
    "gnrh_direct_isolated",
    "gnrh_transcriptomic_candidate"
  )

  if (
    !"gnrh_transcriptomic_candidate" %in%
    names(diagnostics) &&
    "gnrh_dropout_candidate" %in%
    names(metadata)
  ) {
    diagnostics$gnrh_transcriptomic_candidate <-
      metadata$gnrh_dropout_candidate
  }


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
    gnrh_theme(
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
    ggplot2::geom_point(alpha = 0.42, size = 1.5, stroke = 0) +
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
      ggplot2::geom_point(alpha = 0.42, size = 1.5, stroke = 0) +
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

  make_internal_threshold_curve <- function(score, reference, n_thresholds = 200L) {
    score <- suppressWarnings(as.numeric(score))
    reference <- as.logical(reference)
    ok <- is.finite(score) & !is.na(reference)
    score <- score[ok]; reference <- reference[ok]

    if (length(score) < 2L || length(unique(reference)) < 2L) return(NULL)

    rng <- range(score, finite = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) return(NULL)

    safe_div <- function(x, y) ifelse(y > 0, x / y, NA_real_)

    curve <- lapply(seq(rng[1], rng[2], length.out = n_thresholds), function(thr) {
      pred <- score >= thr
      TP <- sum(pred & reference); FP <- sum(pred & !reference)
      FN <- sum(!pred & reference); TN <- sum(!pred & !reference)

      sens <- safe_div(TP, TP + FN)
      spec <- safe_div(TN, TN + FP)
      prec <- safe_div(TP, TP + FP)
      F1 <- if (is.finite(prec) && is.finite(sens) && prec + sens > 0)
        2 * prec * sens / (prec + sens) else NA_real_

      data.frame(
        threshold = thr,
        sensitivity = sens,
        specificity = spec,
        precision = prec,
        F1 = F1
      )
    })

    do.call(rbind, curve)
  }


  # ----------------------------------------------------------------------- #-- #-- #
  # Resolve threshold source
  # ----------------------------------------------------------------------- #-- #-- #

  threshold_mode <- "none"
  report_threshold_curve <- NULL

  if (resolved_roc_mode == "external" && !is.null(threshold_curve)) {
    report_threshold_curve <- as.data.frame(threshold_curve)
    threshold_mode <- "external"

  } else if (resolved_roc_mode == "internal") {
    internal_predictor <- if ("gnrh_support_score_raw" %in% colnames(metadata)) {
      metadata$gnrh_support_score_raw
    } else if ("gnrh_support_score" %in% colnames(metadata)) {
      metadata$gnrh_support_score
    } else NULL

    if (!is.null(internal_predictor)) {
      report_threshold_curve <- make_internal_threshold_curve(
        internal_predictor,
        as.character(metadata$gnrh_status) == "pos"
      )
      if (!is.null(report_threshold_curve)) threshold_mode <- "internal"
    }

  } else if (!is.null(threshold_curve)) {
    report_threshold_curve <- as.data.frame(threshold_curve)
    threshold_mode <- "external"
  }


  # ----------------------------------------------------------------------- #-- #-- #
  # Plot
  # ----------------------------------------------------------------------- #-- #-- #

  req <- c("threshold", "sensitivity", "specificity", "F1")

  if (!is.null(report_threshold_curve) &&
      all(req %in% colnames(report_threshold_curve))) {

    d <- report_threshold_curve
    valid <- is.finite(d$threshold) & is.finite(d$F1)

    best_i <- if (any(valid))
      which.max(replace(d$F1, !valid, -Inf)) else NA_integer_

    best <- if (!is.na(best_i)) d$threshold[best_i] else NA_real_
    best_F1 <- if (!is.na(best_i)) d$F1[best_i] else NA_real_

    long <- rbind(
      data.frame(threshold = d$threshold, metric = "Sensitivity", value = d$sensitivity),
      data.frame(threshold = d$threshold, metric = "Specificity", value = d$specificity),
      data.frame(threshold = d$threshold, metric = "F1", value = d$F1)
    )

    cols <- c(
      Sensitivity = "#0072B2",
      Specificity = "#009E73",
      F1 = "#D55E00"
    )

    subtitle <- if (threshold_mode == "internal") {
      if (is.finite(best))
        sprintf("Internal consistency; optimal threshold %.3f", best)
      else
        "Internal consistency; not independent validation"
    } else {
      if (is.finite(best))
        sprintf("Optimal threshold %.3f", best)
      else
        "Independent validation"
    }

    metric_labels <- c(
      Sensitivity = "Sensitivity",
      Specificity = "Specificity",
      F1 = if (is.finite(best_F1))
        sprintf("F1 (max = %.3f)", best_F1)
      else
        "F1"
    )

    xlab <- if (
      threshold_mode == "internal" &&
      "gnrh_support_score_raw" %in% colnames(metadata)
    ) {
      "GnRH support-score threshold"
    } else if (threshold_mode == "internal") {
      "Standardized support-score threshold"
    } else {
      "Threshold"
    }

    p4 <- ggplot2::ggplot(
      long,
      ggplot2::aes(
        .data$threshold,
        .data$value,
        colour = .data$metric
      )
    ) +
      ggplot2::geom_line(
        linewidth = 0.7,
        na.rm = TRUE
      ) +
      ggplot2::geom_vline(
        xintercept = best,
        linetype = "dashed",
        colour = muted,
        linewidth = 0.45,
        na.rm = TRUE
      ) +
      ggplot2::scale_colour_manual(
        values = cols,
        breaks = c("Sensitivity", "Specificity", "F1"),
        labels = metric_labels
      )  +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        expand = ggplot2::expansion(mult = c(0, 0.03))
      ) +
      ggplot2::labs(
        title = "Threshold performance",
        subtitle = subtitle,
        x = xlab,
        y = "Performance",
        colour = NULL
      ) +
      report_theme() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(
          margin = ggplot2::margin(t = -20)
        )
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = list(linewidth = 1.1),
          keyheight = grid::unit(0.32, "cm")
        )
      )

  } else {
    p4 <- empty_panel(
      "Threshold curve unavailable",
      if (resolved_roc_mode == "internal")
        "GnRH support score unavailable for internal threshold analysis"
      else
        "Supply an independent truth annotation for threshold validation"
    )
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
            identity =
              "gnrh_stage_identity_score",

            migrating =
              "gnrh_stage_migrating_score",

            mature =
              "gnrh_stage_mature_score"
          )

          positive_stage <- as.character(
            diagnostics$gnrh_stage
          )

          for (stage in names(stage_scores)) {
            score_column <-
              stage_scores[[stage]]

            if (
              stage %in% positive_stage &&
              score_column %in%
              names(diagnostics)
            ) {
              roc_list[[stage]] <-
                roc_add(
                  positive_stage == stage,
                  diagnostics[[score_column]],
                  stage
                )
            }
          }
        }


        if (
          "gnrh_secretory" %in%
          names(diagnostics) &&
          "gnrh_secretory_hits" %in%
          names(diagnostics)
        ) {
          roc_list[["Secretory"]] <-
            roc_add(
              diagnostics$gnrh_secretory ==
                "supported",
              diagnostics$gnrh_secretory_hits,
              "Secretory"
            )
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

  if (is.null(p5) && show_class_panel) {
    status_order <- c("neg", "pos")
    status_data <- as.data.frame(
      table(factor(as.character(diagnostics$status), levels = status_order)),
      stringsAsFactors = FALSE
    )
    names(status_data) <- c("status", "n")
    status_data$pct <- 100 * status_data$n / sum(status_data$n)
    status_colors <- gnrh_colors("status")

    p5 <- ggplot2::ggplot(status_data, ggplot2::aes(.data$status, .data$n, fill = .data$status)) +
      ggplot2::geom_col(width = 0.68, colour = tile_border, linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%s\n%.1f%%", format(.data$n, big.mark = ","), .data$pct)),
        vjust = -0.2, size = txtsize / 3.2, colour = foreground
      ) +
      ggplot2::scale_fill_manual(values = status_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_comma(),
        expand = ggplot2::expansion(mult = c(0, 0.16))
      ) +
      ggplot2::labs(title = "GnRH detection status", x = NULL, y = "Cells") +
      report_theme(leg.pos = "none")
  }
  if (is.null(p5)) p5 <- empty_panel("Detection-status panel unavailable")

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
      title = "gnrhcell diagnostic report",
      subtitle = subtitle,
      tag_levels = "A",
      theme = annotation_theme
    )
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
plot_network <- function(
    df,
    top_n=25L,
    threshold=0.10,
    mode=c("auto","coexpression","codetection","phenotype"),
    score_col=NULL,
    association_col=NULL,
    include_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize=10,
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
      ggplot2::theme_void(base_size=txtsize) +
      ggplot2::theme(
        text=ggplot2::element_text(family="Helvetica",size=txtsize),
        plot.title=ggplot2::element_text(face="bold",hjust=0.5,size=txtsize+2),
        plot.subtitle=ggplot2::element_text(hjust=0.5,size=max(7,txtsize-1)),
        legend.position="right",
        legend.title=ggplot2::element_text(face="bold",size=txtsize+2),
        legend.text=ggplot2::element_text(size=txtsize)
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
plot_gnrh_coexpr <- function(
    df,
    coexp_cutoff=0.10,
    top_n=40L,
    txtsize=10,
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
plot_gnrh_codetect <- function(
    markers,
    top_n=12L,
    min_or=2,
    max_fdr=0.05,
    min_specificity=0.05,
    exclude_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize=10,
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
plot_gnrh_detection <- function(
    markers,
    top_n=12L,
    min_specificity=0.05,
    max_padj=0.05,
    exclude_gnrh=TRUE,
    gnrh_gene="GNRH1",
    txtsize=10,
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
      title = "gnrhcell runtime across datasets",
      x = "Dataset",
      y = metric
    ) +
    gnrh_theme(
      x.ang = x.ang,
      txtsize = txtsize
    )
}



# plot_gnrh_detected <- function(
#     files = NULL,
#     dir = file.path("results", "tables"),
#     pattern = "_gnrh_run_info\\.tsv$",
#     x.ang = 45,
#     txtsize = 10,
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
#       fill = gnrh_colors("status")[["pos"]]
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
# plot_gnrh_detected
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
plot_gnrh_detected <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    txtsize = 10,
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
  ggplot2::ggplot(stats,ggplot2::aes(dataset,n_cells)) +
    ggplot2::geom_col(width=0.7,fill=gnrh_colors("status")[["pos"]]) +
    ggplot2::geom_text(
      ggplot2::aes(label=paste0(
        scales::comma(round(n_cells))," cells\n",
        scales::comma(round(pos))," GnRH+ (",
        sprintf(paste0("%.",digits,"f"),pct_gnrh),"%)"
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
#' @param programs Result returned by \code{gnrh_marker_programs()}.
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
plot_gnrh_marker_programs <- function(
    programs,
    table=c("summary","high_confidence","candidate_table","context_specific"),
    type=c("bar","dot","tile"),
    min_genes=1L,
    require_coexpr=FALSE,
    top_n=NULL,
    facet_ncol=NULL,
    mode=c("light","dark"),
    txtsize=12,
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
plot_gnrh_specificity <- function(
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
        gnrh_colors("status")
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

