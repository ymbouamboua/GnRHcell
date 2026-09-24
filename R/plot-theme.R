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
.auto_raster <- function(ncells, threshold = 3e5)
  .cellplot_auto_raster(ncells, threshold)


#' Automatically determine point size
#' @keywords internal
#' @noRd
.auto_pt_size <- function(ncells, raster = FALSE) {
  .cellplot_auto_pt_size(ncells, raster)
}


#' Internal plotting defaults
#' @keywords internal
#' @noRd
.plot_defaults <- function(object, raster = NULL, pt.size = NULL, raster.threshold = 3e5) {
  params <- if (inherits(object, "Seurat")) {
    .cellplot_point_params(
      object, raster = raster, pt.size = pt.size,
      raster.threshold = raster.threshold
    )
  } else {
    ncells <- ncol(object)
    raster <- raster %||% .cellplot_auto_raster(ncells, raster.threshold)
    pt.size <- pt.size %||% .cellplot_auto_pt_size(ncells, raster)
    list(raster = raster, pt.size = pt.size)
  }
  list(raster = params$raster, pt.size = params$pt.size)
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
#' @param family Font family.
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
#' @param axis.size,axis.ttl.size Axis-text and axis-title sizes.
#' @param axis.face,axis.ttl.face Axis-text and axis-title faces.
#' @param ttl.size,subtitle.size,caption.size Title, subtitle, and caption sizes.
#' @param ttl.scope Align the title relative to the plotting panel (`"panel"`)
#'   or the complete figure including legends (`"plot"`).
#' @param facet.size Facet-label size.
#' @param linewidth Line width used by axes, borders, and grids.
#' @param tick.length Axis tick length in points.
#' @param grid.col Grid colour.
#' @param leg.key.size Legend-key size in centimetres.
#' @param leg.box Legend-box arrangement.
#' @param plot.margin Plot margins: top, right, bottom, left, in points.
#' @param panel.spacing Facet-panel spacing in centimetres.
#' @param aspect.ratio Optional panel aspect ratio.
#' @param facet.bg Display facet-strip backgrounds.
#' @param mode Display mode: `"light"` or `"dark"`.
#' @param leg.pos Legend position.
#' @param leg.dir Legend direction.
#' @param leg.ncol Optional number of legend columns; NULL uses the default.
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
    txtsize = getOption("gnrhcell.base_size", 14),
    family = getOption("gnrhcell.base_family", "Helvetica"),
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
    axis.size = NULL,
    axis.ttl.size = NULL,
    axis.face = "plain",
    axis.ttl.face = "plain",
    ttl.size = NULL,
    ttl.scope = c("panel", "plot"),
    subtitle.size = NULL,
    caption.size = NULL,
    facet.size = NULL,
    linewidth = 0.3,
    tick.length = 2.5,
    grid.col = NULL,
    facet.bg = TRUE,
    mode = c("light", "dark"),
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = NULL,
    leg.ttl = NULL,
    leg.ncol=NULL,
    leg.ttl.size = NULL,
    leg.just = "center",
    leg.ttl.text = NULL,
    leg.key.size = 0.4,
    leg.box = c("vertical", "horizontal"),
    plot.margin = c(5.5, 5.5, 5.5, 5.5),
    panel.spacing = 0.1,
    aspect.ratio = NULL,
    ...
) {
  # ========================================================================= #
  # Arguments
  # ========================================================================= #
  style <- match.arg(style)
  ttl.pos <- match.arg(ttl.pos)
  txt.face <- match.arg(txt.face)
  mode <- match.arg(mode)
  ttl.scope <- match.arg(ttl.scope)
  leg.box <- match.arg(leg.box)
  leg.size <- leg.size %||% max(10, txtsize - 1)
  leg.ttl.size <- leg.ttl.size %||% txtsize
  axis.size <- axis.size %||% txtsize
  axis.ttl.size <- axis.ttl.size %||% txtsize
  ttl.size <- ttl.size %||% (txtsize + 2)
  subtitle.size <- subtitle.size %||% max(7, txtsize - 1)
  caption.size <- caption.size %||% max(6, txtsize - 2)
  facet.size <- facet.size %||% txtsize
  lw <- linewidth
  if (is.null(line)) line <- style == "classic"
  fg <- if (mode == "dark") "white" else "#1A1A1A"
  bg <- if (mode == "dark") "#111111" else panel.fill
  col.grid <- grid.col %||% if (mode == "dark") "#444444" else "#D9D9D9"
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
    text = ggplot2::element_text(
      colour = fg,
      size = txtsize,
      face = txt.face,
      family = family
    ),
    axis.text.x = ggplot2::element_text(colour = fg, size = axis.size, face = axis.face, angle = x.ang, hjust = hjust, vjust = vjust),
    axis.text.y = ggplot2::element_text(colour = fg, size = axis.size, face = axis.face),
    axis.title = ggplot2::element_text(colour = fg, size = axis.ttl.size, face = axis.ttl.face),
    axis.ticks.length = grid::unit(tick.length, "pt"),
    plot.title = ggplot2::element_text(hjust = ttl.hjust, face = ttl.face, size = ttl.size, colour = fg),
    plot.subtitle = ggplot2::element_text(colour = fg, size = subtitle.size),
    plot.caption = ggplot2::element_text(colour = fg, size = caption.size),
    strip.text = ggplot2::element_text(face = facet.face, size = facet.size, colour = fg),
    strip.background = ggplot2::element_rect(fill = col.strip, colour = NA),
    panel.background = ggplot2::element_rect(fill = bg, colour = NA),
    plot.background = ggplot2::element_rect(fill = bg, colour = NA),
    legend.title = ggplot2::element_text(size = leg.ttl.size + 2, face = "bold", colour = fg),
    legend.text = ggplot2::element_text(size = leg.size, colour = fg),
    legend.position = leg.pos,
    legend.direction = leg.dir,
    legend.justification = leg.just,
    legend.key.height = grid::unit(leg.key.size, "cm"),
    legend.key.width = grid::unit(leg.key.size, "cm"),
    legend.background = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.box = leg.box,
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = ggplot2::margin(1, 1, 1, 1),
    panel.spacing = grid::unit(panel.spacing, "cm"),
    plot.margin = ggplot2::margin(
      plot.margin[1], plot.margin[2], plot.margin[3], plot.margin[4], unit = "pt"
    ),
    aspect.ratio = aspect.ratio,
    plot.title.position = ttl.scope,
    plot.caption.position = "plot",
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
  if (style=="test") {
    th <- th + ggplot2::theme(
      axis.title.x=if (isTRUE(x.ttl)) {
        ggplot2::element_text(
          colour=fg,
          size=axis.ttl.size,
          face=axis.ttl.face
        )
      } else {
        ggplot2::element_blank()
      },
      axis.title.y=if (isTRUE(y.ttl)) {
        ggplot2::element_text(
          colour=fg,
          size=axis.ttl.size,
          face=axis.ttl.face
        )
      } else {
        ggplot2::element_blank()
      }
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
  # Final void-style guarantee
  # ========================================================================= #
  # Keep this after all axis and panel controls. Thus optional arguments such
  # as `ticks`, `border`, or `grid.major` cannot partially restore axes when a
  # caller explicitly requests an axis-free plot.
  if (style == "void") {
    th <- th + ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0, "pt"),
      axis.line = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  }
  # --------------------------------------------------------------------------- #
  # Legend title
  # --------------------------------------------------------------------------- #
  if (!is.null(leg.ttl.text)) {
    warning(
      "`leg.ttl.text` cannot be applied through a theme and should be set with labs().",
      call.=FALSE
    )
  }

  # --------------------------------------------------------------------------- #
  # Legend columns
  # --------------------------------------------------------------------------- #
  if (!is.null(leg.ncol)) {
    if (!is.numeric(leg.ncol) || length(leg.ncol)!=1L || !is.finite(leg.ncol) || leg.ncol<1)
      stop("`leg.ncol` must be a positive integer or NULL.",call.=FALSE)

    return(list(
      th,
      ggplot2::guides(
        colour=ggplot2::guide_legend(ncol=as.integer(leg.ncol)),
        fill=ggplot2::guide_legend(ncol=as.integer(leg.ncol))
      )
    ))
  }

  th
}



# ========================================================================= #
# Palettes
# ========================================================================= #
#' Color palettes
#' @param type Palette type.
#' @return Named character vector.
#' @export
gnrh_palette <- function(type=c("status","confident","class","stage","developmental","secretory")) {
  switch(
    match.arg(type),
    status=c(neg="#B0B0B0",pos="#FF4D6D"),
    confident=c("FALSE"="#B0B0B0","TRUE"="#E63946"),
    class=c(neg="#B0B0B0",supported="#4EA8DE",direct="#E63946"),
    stage=c("non-gnrh"="#B0B0B0",early="#3B4CC0",migrating="#6A3D9A","post-migratory"="#7BC8A4",mature="#E67E22",transitional="#B8A9C9",undetermined="#D9D9D9"),
    developmental=c("non-gnrh"="#B0B0B0",early="#3B4CC0",migrating="#00A6CA","post-migratory"="#7BC8A4",mature="#F28E2B",transitional="#B8A9C9",undetermined="#D9D9D9"),
    secretory=c("non-gnrh"="#B0B0B0",limited="#80B1D3",supported="#D81B60")
  )
}
#' Internal GnRH palette resolver
#' @keywords internal
#' @noRd
.gnrh_palette <- function(group_by) {
  type <- switch(group_by,gnrh_status="status",gnrh_confident="confident",gnrh_class="class",gnrh_stage="stage",gnrh_stage_raw="stage",gnrh_secretory="secretory",NULL)
  if (is.null(type)) return(NULL)
  gnrh_palette(type)
}
#' Resolve colors for arbitrary levels
#' @keywords internal
#' @noRd
.resolve_colors <- function(levels,cols=NULL,group_by=NULL) {
  if (is.null(cols) && !is.null(group_by)) cols <- .gnrh_palette(group_by)
  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(levels))
    names(cols) <- levels
  }
  if (is.null(names(cols))) names(cols) <- levels[seq_len(min(length(cols),length(levels)))]
  missing <- setdiff(levels,names(cols))
  if (length(missing)) {
    extra <- scales::hue_pal()(length(missing))
    names(extra) <- missing
    cols <- c(cols,extra)
  }
  cols[levels]
}
# ========================================================================= #
# Legend helpers
# ========================================================================= #
#' Compact point legend
#' @keywords internal
#' @noRd
.point_legend <- function(size=2.3) {
  list(
    ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(size=size,alpha=1),keyheight=grid::unit(0.28,"cm"))),
    ggplot2::theme(legend.spacing.y=grid::unit(0,"cm"),legend.key.width=grid::unit(0.35,"cm"))
  )
}
#' Compact line legend
#' @keywords internal
#' @noRd
.line_legend <- function(linewidth=1.2) {
  list(
    ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(linewidth=linewidth),keyheight=grid::unit(0.30,"cm"))),
    ggplot2::theme(legend.spacing.y=grid::unit(0,"cm"),legend.key.width=grid::unit(0.45,"cm"))
  )
}
# ========================================================================= #
# Reduction resolver
# ========================================================================= #
#' Resolve an available dimensional reduction
#' @param object A Seurat object.
#' @param reduction Optional requested reduction name.
#' @param fallback Ordered character vector of fallback reductions.
#' @return The name of an available dimensional reduction.
#' @export
gnrh_reduction <- function(object,reduction=NULL,fallback=c("umap","umap_scvi","umap.harmony","umap.rpca","pca")) {
  .validate_seurat(object)
  available <- names(object@reductions)
  if (!length(available)) stop("No dimensional reductions are available in `object`.",call.=FALSE)
  if (!is.null(reduction) && length(reduction)==1L && !is.na(reduction) && nzchar(reduction)) {
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
# Legend labels
# ========================================================================= #
#' Legend label
#' @keywords internal
#' @noRd
.gnrh_legend_labels <- function(x,levels=NULL,percentage=FALSE,n_cells=TRUE) {
  x <- as.character(x)
  x[is.na(x)] <- "Unknown"
  tab <- table(x)
  levels <- intersect(levels %||% names(tab),names(tab))
  if (!isTRUE(n_cells)) return(stats::setNames(levels,levels))
  n <- as.integer(tab[levels])
  n_fmt <- format(n,big.mark=",",trim=TRUE,scientific=FALSE)
  labels <- if (isTRUE(percentage)) paste0(levels," (",n_fmt," | ",round(100*n/sum(tab),1),"%)") else paste0(levels," (",n_fmt,")")
  stats::setNames(labels,levels)
}




# ========================================================================= #
