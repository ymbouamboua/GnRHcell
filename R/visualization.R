#=========================================================#
# Internal plotting helpers
#=========================================================#

#' Return left-hand side unless NULL
#'
#' @param x Object to test.
#' @param y Fallback object.
#'
#' @return \code{x} if not \code{NULL}, otherwise \code{y}.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Automatically decide whether to rasterize
#'
#' @param ncells Number of cells.
#' @param threshold Cell number threshold above which rasterization is used.
#'
#' @return Logical.
#' @keywords internal
#' @noRd
.auto_raster <- function(ncells, threshold = 1e5) ncells > threshold


#' Automatically compute point size
#'
#' @param ncells Number of cells.
#' @param raster Logical; whether rasterization is used.
#'
#' @return Numeric point size.
#' @keywords internal
#' @noRd
.auto_pt_size <- function(ncells, raster = FALSE) {
  if (!raster) return(max(0.3, 1.5 / log10(ncells)))
  if (ncells < 1e5) 0.1 else if (ncells < 2e5) 0.5 else 2
}


#' Validate Seurat object
#'
#' @param obj Object to validate.
#'
#' @return Invisibly returns \code{TRUE}.
#' @keywords internal
#' @noRd
.validate_seurat <- function(obj) {
  if (!inherits(obj, "Seurat")) stop("Input must be Seurat object.", call. = FALSE)
  invisible(TRUE)
}



#' GnRHcell ggplot2 theme
#'
#' Flexible ggplot2 theme used across GnRHcell visualizations.
#'
#' @param style Theme style. One of \code{"classic"}, \code{"minimal"},
#' \code{"bw"}, \code{"test"}, \code{"void"}, \code{"dirty"}, or \code{"gray"}.
#' @param txtsize Base text size.
#' @param xy.val Show axis tick labels.
#' @param x.ang X-axis text angle.
#' @param hjust,vjust Horizontal and vertical justification for x-axis labels.
#' @param xlab,ylab Show x/y axis tick labels.
#' @param xy.lab Show all axis tick labels.
#' @param facet.face Facet label font face.
#' @param ttl.face Plot title font face.
#' @param txt.face Text font face.
#' @param ttl.pos Plot title position.
#' @param x.ttl,y.ttl Show x/y axis titles.
#' @param ticks,line,border Logical overrides for ticks, axis lines, and panel border.
#' @param grid.major,grid.minor Logical overrides for major/minor grid lines.
#' @param panel.fill Panel background fill.
#' @param facet.bg Show facet background.
#' @param mode Theme mode: \code{"light"} or \code{"dark"}.
#' @param leg.pos Legend position.
#' @param leg.dir Legend direction.
#' @param leg.size Legend text size.
#' @param leg.ttl Legend title size.
#' @param leg.ttl.size Legend title text size.
#' @param leg.just Legend justification.
#' @param leg.ttl.text Optional legend title text.
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @return A ggplot2 theme object.
#' @export
plot_theme <- function(
    style = c("classic","minimal","bw","test","void","dirty","gray"),
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
    txt.face = c("plain","italic","bold"),
    ttl.pos = c("center","left","right"),
    x.ttl = TRUE,
    y.ttl = TRUE,
    ticks = NULL,
    line = NULL,
    border = NULL,
    grid.major = NULL,
    grid.minor = NULL,
    panel.fill = "white",
    facet.bg = TRUE,
    mode = c("light","dark"),
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = 10,
    leg.ttl = 10,
    leg.ttl.size = 10,
    leg.just = "center",
    leg.ttl.text = NULL,
    ...
) {

  style <- match.arg(style)
  ttl.pos <- match.arg(ttl.pos)
  txt.face <- match.arg(txt.face)
  mode <- match.arg(mode)

  # Canonical line width (THIS FIXES YOUR PROBLEM)
  lw <- 0.3

  if (is.null(line)) {
    line <- style == "classic"
  }

  # Colors
  if (mode == "light") {
    col.txt   <- "#1A1A1A"
    col.grid  <- "#D9D9D9"
    col.panel <- panel.fill
    col.strip <- "#EFEFEF"
  } else {
    col.txt   <- "#DDDDDD"
    col.grid  <- "#444444"
    col.panel <- "#1E1E1E"
    col.strip <- "#383838"
  }

  # Automatic x-label alignment
  if (is.null(hjust) || is.null(vjust)) {
    if (x.ang == 0)   { hjust <- .5; vjust <- .5 }
    else if (x.ang == 45) { hjust <- 1; vjust <- 1 }
    else if (x.ang == 90) { hjust <- 1; vjust <- .5 }
    else if (x.ang == 270){ hjust <- 0; vjust <- .5 }
    else { hjust <- 1; vjust <- 1 }
  }


  ttl.pos <- switch(ttl.pos, left = 0, center = .5, right = 1)

  # Base theme components
  base <- theme(
    text = element_text(color = col.txt, size = txtsize, family = "Helvetica"),
    axis.text.x = element_text(color = col.txt, size = txtsize),
    axis.text.y = element_text(color = col.txt, size = txtsize),
    axis.title  = element_text(size = txtsize),
    plot.title  = element_text(
      hjust = ttl.pos, face = ttl.face,
      size = txtsize + 2, color = col.txt
    ),
    strip.text = element_text(face = facet.face, color = col.txt),
    legend.title = element_text(size = leg.ttl.size + 2, face = "bold"),
    legend.text  = element_text(size = leg.size),
    legend.position = leg.pos,
    legend.key.height = unit(.4, "cm"),
    legend.key.width  = unit(.4, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.box = "vertical",
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = ggplot2::margin(1,1,1,1),
    ...
  )

  # Preset themes
  preset <- switch(
    style,
    minimal = theme_minimal(base_size = txtsize),
    classic = theme_classic(base_size = txtsize),
    bw      = theme_bw(base_size = txtsize),
    test    = theme_test(base_size = txtsize),
    void    = ggplot2::theme_void(base_size = txtsize),
    dirty   = theme_minimal(base_size = txtsize) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()
      ),
    gray    = theme_gray(base_size = txtsize) +
      theme(
        panel.background = element_rect(fill = "#EDEDED", color = NA),
        panel.grid.major = element_line(color = "#CCCCCC", linewidth = lw),
        panel.grid.minor = element_line(color = "#DDDDDD", linewidth = lw/2)
      )
  )

  th <- preset + base

  # Normalize panel borders for border-based themes
  if (style %in% c("bw", "test", "gray")) {
    th <- th + theme(
      panel.border = element_rect(
        linewidth = lw,
        color = col.txt,
        fill = NA
      ),
      axis.line = element_blank()
    )
  }

  # Axis lines (classic-style)
  if (line && style == "classic") {
    th <- th + theme(
      axis.line.x = element_line(color = col.txt, linewidth = lw),
      axis.line.y = element_line(color = col.txt, linewidth = lw)
    )
  } else {
    th <- th + theme(axis.line = element_blank())
  }

  # X-axis angle
  if (xy.val && xlab) {
    th <- th + theme(
      axis.text.x = element_text(angle = x.ang, hjust = hjust, vjust = vjust)
    )
  }

  # Label / tick / grid overrides
  if (!xy.lab) th <- th + theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (!xlab)   th <- th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (!ylab)   th <- th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (!x.ttl)  th <- th + theme(axis.title.x = element_blank())
  if (!y.ttl)  th <- th + theme(axis.title.y = element_blank())

  if (!is.null(ticks)) {
    th <- th + if (ticks)
      theme(axis.ticks = element_line(color = col.txt, linewidth = lw))
    else
      theme(axis.ticks = element_blank())
  }

  if (!is.null(border)) {
    th <- th + if (border)
      theme(panel.border = element_rect(color = col.grid, fill = NA, linewidth = lw))
    else
      theme(panel.border = element_blank())
  }

  if (!is.null(grid.major)) {
    th <- th + if (grid.major)
      theme(panel.grid.major = element_line(color = col.grid, linewidth = lw))
    else
      theme(panel.grid.major = element_blank())
  }

  if (!is.null(grid.minor)) {
    th <- th + if (grid.minor)
      theme(panel.grid.minor = element_line(color = col.grid, linewidth = lw/2))
    else
      theme(panel.grid.minor = element_blank())
  }

  if (!facet.bg) {
    th <- th + theme(strip.background = element_blank())
  }

  if (!is.null(leg.ttl.text)) {
    th <- th + labs(color = leg.ttl.text)
  }

  # Void cleanup
  # Void cleanup
  if (style == "void") {
    th <- th + theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line    = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      strip.text   = element_blank(),
      strip.background = element_blank()
    )
    line <- FALSE
    ticks <- FALSE
    border <- FALSE
    grid.major <- FALSE
    grid.minor <- FALSE
    facet.bg <- FALSE
  }

  # Dirty theme
  if (style == "dirty") {
    line       <- TRUE
    ticks      <- FALSE
    border     <- FALSE
    grid.major <- FALSE
    grid.minor <- FALSE
    facet.bg   <- FALSE
  }

  th
}


#' Dark ggplot2 theme
#'
#' Dark theme presets for presentation-quality GnRHcell plots.
#'
#' @param style Dark theme style: \code{"soft"}, \code{"true"}, or \code{"paper"}.
#' @param txtsize Base text size.
#' @param family Font family.
#' @param axis.text Logical; show axis tick labels.
#' @param axis.title Logical; show axis titles.
#'
#' @return A ggplot2 theme object.
#' @export
dark_theme <- function(
    style = c("soft", "true", "paper"),
    txtsize = 12,
    family = "Helvetica",
    axis.text = TRUE,   # show/hide axis tick labels
    axis.title = TRUE   # show/hide axis titles
) {
  style <- match.arg(style)

  # Base common settings
  common <- theme(
    text = element_text(size = txtsize, family = family),
    legend.key.height  = unit(0.25, "cm"),
    legend.key.width   = unit(0.25, "cm"),
    legend.spacing.y   = unit(2, "pt"),
    legend.box.spacing = unit(2, "pt"),
    panel.grid.minor = element_blank()
  )

  # Helper to conditionally hide axis elements
  axis_theme <- theme(
    axis.text  = if(axis.text) element_text() else element_blank(),
    axis.title = if(axis.title) element_text() else element_blank()
  )

  # Theme switch
  th <- switch(
    style,

    ## SOFT DARK - general / Scanpy-like
    soft = common + theme(
      plot.background  = element_rect(fill = "#1E1E1E", colour = NA),
      panel.background = element_rect(fill = "#1E1E1E", colour = NA),
      panel.grid.major = element_blank(),
      legend.background = element_rect(fill = "#1E1E1E", colour = NA),
      legend.key        = element_rect(fill = "#1E1E1E", colour = NA),
      legend.text  = element_text(color = "grey90"),
      legend.title = element_text(color = "white", face = "bold"),
      plot.title   = element_text(color = "white", face = "bold"),
      plot.subtitle= element_text(color = "grey80"),
      axis.title = element_text(color = "grey85"),
      axis.text  = element_text(color = "grey70"),
      strip.text = element_text(color = "grey90")
    ),

    ## TRUE BLACK - slides
    true = common + theme(
      plot.background  = element_rect(fill = "black", colour = NA),
      panel.background = element_rect(fill = "black", colour = NA),
      panel.grid.major = element_blank(),
      legend.background = element_rect(fill = "black", colour = NA),
      legend.key        = element_rect(fill = "black", colour = NA),
      legend.text  = element_text(color = "white"),
      legend.title = element_text(color = "white", face = "bold"),
      plot.title   = element_text(color = "white", face = "bold"),
      plot.subtitle= element_text(color = "grey80"),
      axis.title = element_text(color = "white"),
      axis.text  = element_text(color = "grey80"),
      strip.text = element_text(color = "white")
    ),

    ## DARK PAPER - journal-friendly
    paper = common + theme(
      plot.background  = element_rect(fill = "#2A2A2A", colour = NA),
      panel.background = element_rect(fill = "#2A2A2A", colour = NA),
      panel.grid.major = element_line(color = "#3A3A3A", linewidth = 0.3),
      legend.background = element_rect(fill = "#2A2A2A", colour = NA),
      legend.key        = element_rect(fill = "#2A2A2A", colour = NA),
      legend.text  = element_text(color = "#E0E0E0"),
      legend.title = element_text(color = "white", face = "bold"),
      plot.title   = element_text(color = "white", face = "bold"),
      plot.subtitle= element_text(color = "#D0D0D0"),
      axis.title = element_text(color = "#E0E0E0"),
      axis.text  = element_text(color = "#BDBDBD"),
      strip.text = element_text(color = "#E0E0E0")
    )
  )

  th + axis_theme
}


#' Generate a discrete color palette
#'
#' This function returns a vector of colors from predefined palettes, thematic palettes (viridis, plasma, etc.),
#' or a user-provided set of base colors. It can oversample, remove grayish colors, adjust lightness and saturation,
#' and reverse the palette order.
#'
#' @param n Integer. Number of colors to return.
#' @param preset Character. Name of a predefined palette. Default NULL.
#' @param theme Character. Thematic palettes: "viridis", "magma", "plasma", "inferno", "cividis", or any RColorBrewer palette name.
#' @param base_colors Character vector. User-defined colors. Overrides preset/theme if provided.
#' @param space Character. Color interpolation space: "Lab", "rgb", or "HCL".
#' @param oversample_factor Numeric. How much to oversample before picking final n colors.
#' @param remove_gray Logical. Remove grayish colors when oversampling.
#' @param reverse Logical. Reverse the color order.
#' @param adjust_saturation Numeric. Factor to adjust saturation (1 = no change).
#' @param adjust_lightness Numeric. Factor to adjust lightness (1 = no change).
#'
#' @return Character vector of colors.
#' @export
cellpal <- function(
    n = 5,
    preset = "base",
    theme = NULL,
    base_colors = NULL,
    space = c("Lab", "rgb", "HCL"),
    oversample_factor = 1.3,
    remove_gray = TRUE,
    reverse = FALSE,
    adjust_saturation = 1,
    adjust_lightness = 1
) {

  space <- match.arg(space)

  # predefined palettes
  palettes <- list(
    base = c("#E41A1C", "#68618B", "#409388", "#57A156", "#8D5B96", "#CB7647", "#F38E38",
             "#F781BE", "#CC95C8", "#B27E85", "#9A6242", "#5FA3C9", "#3766A4", "#204E75",
             "#1B9D77", "#86CC84", "#D3EC90", "#FBF583", "#E7C715", "#F0A957", "#F57994",
             "#E7298A", "#A90D55", "#52587E", "#17CDD3", "#8ECDE0", "#BC6298", "#AE2373",
             "#5E4EA1", "#7E8D86", "#507C51", "#1F5917", "#BEC603", "#C5DD3B", "#A8DA83",
             "#8DD3C7")
  )

  # determine base colors
  if (!is.null(base_colors)) {
    cols <- base_colors
  } else if (!is.null(preset)) {
    if (!preset %in% names(palettes)) stop(
      "Unknown preset: ", preset, ". Available: ", paste(names(palettes), collapse = ", ")
    )
    cols <- palettes[[preset]]
  } else if (!is.null(theme)) {
    if (theme == "viridis") cols <- viridisLite::viridis(n)
    else if (theme == "magma") cols <- viridisLite::viridis(n, option = "magma")
    else if (theme == "plasma") cols <- viridisLite::viridis(n, option = "plasma")
    else if (theme == "inferno") cols <- viridisLite::viridis(n, option = "inferno")
    else if (theme == "cividis") cols <- viridisLite::viridis(n, option = "cividis")
    else if (theme %in% rownames(brewer.pal.info)) cols <- colorRampPalette(brewer.pal(8, theme))(n)
    else stop("Unknown theme: ", theme)
  } else {
    stop("Provide either preset, theme, or base_colors.")
  }

  # if enough colors, return first n
  if (n <= length(cols)) {
    final_colors <- cols[1:n]
    if (reverse) final_colors <- rev(final_colors)
    return(final_colors)
  }

  # oversample & interpolate
  extra_colors <- ceiling(n * oversample_factor)
  extended_colors <- grDevices::colorRampPalette(cols, space = space)(extra_colors)

  # remove grayish colors
  if (remove_gray) {
    is_grayish <- function(col) sd(grDevices::col2rgb(col)) < 15
    extended_colors <- Filter(function(c) !is_grayish(c), extended_colors)
  }

  # adjust saturation/lightness
  if (adjust_saturation != 1 || adjust_lightness != 1) {
    hcl_vals <- colorspace::coords(
      methods::as(colorspace::hex2RGB(extended_colors), "polarLUV")
    )
    hcl_vals[, "C"] <- pmax(0, hcl_vals[, "C"] * adjust_saturation)
    hcl_vals[, "L"] <- pmax(0, pmin(100, hcl_vals[, "L"] * adjust_lightness))
    extended_colors <- colorspace::hex(colorspace::polarLUV(hcl_vals))
  }

  # return first n colors
  final_colors <- extended_colors[1:n]
  if (reverse) final_colors <- rev(final_colors)
  return(final_colors)
}


#' GnRHcell color palettes
#'
#' Returns named color palettes for GnRHcell metadata variables.
#'
#' @param type Palette type. One of \code{"status"},
#' \code{"truth"}, or \code{"stage"}.
#'
#' @return Named character vector of colors.
#' @export
gnrh_colors <- function(type = c("status","truth","stage")) {

  switch(
    match.arg(type),

    # binary detection
    status = c(
      pos = "#FF4D6D",   # vivid GnRH positive (raspberry red)
      neg = "#B0B0B0"    # neutral cool gray
    ),

    # truth labels
    truth = c(
      pos = "#E63946",   # strong biologically "true positive"
      neg = "#B0B0B0"
    ),

    stage = c(
      "non-gnrh"    = "#B0B0B0",
      neurogenesis = "#1B9E77",
      identity     = "#1F78B4",
      migrating    = "#6A3D9A",
      mature       = "#E67E22",
      secreting    = "#D81B60"
    )
  )
}


#' @keywords internal
#' @noRd
.plot_defaults <- function(obj, raster = NULL, pt.size = NULL) {
  ncells <- ncol(obj)
  raster <- raster %||% .auto_raster(ncells)
  pt.size <- pt.size %||% .auto_pt_size(ncells, raster)
  list(raster = raster, pt.size = pt.size)
}



#' Internal scatter plot helper
#'
#' @param data Data frame.
#' @param x,y Column names for x and y axes.
#' @param color Optional column used for point color.
#' @param cols Optional named color vector.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param pt.size Point size.
#' @param alpha Point transparency.
#' @param style Theme style.
#'
#' @return A ggplot2 object.
#' @keywords internal
#' @noRd
.scatter_plot <- function(data, x, y, color = NULL, cols = NULL,
                          title = NULL, xlab = NULL, ylab = NULL,
                          pt.size = 1, alpha = 0.8, style = "classic") {

  p <- ggplot2::ggplot(data, ggplot2::aes(.data[[x]], .data[[y]]))

  if (!is.null(color)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color]]),
                          size = pt.size, alpha = alpha)

    if (!is.null(cols))
      p <- p + ggplot2::scale_color_manual(values = cols)
  } else {
    p <- p + ggplot2::geom_point(size = pt.size, alpha = alpha)
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    plot_theme(style = style)
}


#' Plot metadata distributions
#'
#' Publication-ready barplot utility for Seurat metadata.
#'
#' Supports:
#' \itemize{
#'   \item Counts or proportions
#'   \item Distribution by sample/group
#'   \item Stacked or dodged bars
#'   \item Horizontal plots
#' }
#'
#' @param object Seurat object.
#' @param group.by Metadata column to plot.
#' @param split.by Optional metadata column for grouped plots
#'   (e.g. sample, batch, condition).
#' @param cols Named color vector.
#' @param sort Sort bars.
#' @param decreasing Sort decreasing.
#' @param proportion Plot proportions instead of counts.
#' @param position Bar position:
#'   \code{"stack"} or \code{"dodge"}.
#' @param label Add labels.
#' @param label.size Label text size.
#' @param border Draw borders.
#' @param border.col Border color.
#' @param border.size Border linewidth.
#' @param width Bar width.
#' @param x.lab X axis label.
#' @param y.lab Y axis label.
#' @param plot.ttl Plot title.
#' @param txtsize Base text size.
#' @param x.ang X axis angle.
#' @param style Theme style.
#' @param flip Flip coordinates.
#' @param ... Additional arguments passed to \code{plot_theme()}.
#'
#' @return ggplot object.
#' @export
plot_gnrh_distribution <- function(
    object,
    group.by,
    split.by = NULL,
    cols = NULL,
    sort = FALSE,
    decreasing = TRUE,
    proportion = FALSE,
    position = "stack",
    label = TRUE,
    label.size = 3,
    border = TRUE,
    border.col = "black",
    border.size = 0.2,
    width = 0.7,
    x.lab = NULL,
    y.lab = NULL,
    plot.ttl = NULL,
    txtsize = 10,
    x.ang = 45,
    style = "test",
    flip = FALSE,
    ...
) {

  stopifnot(
    inherits(object, "Seurat"),
    group.by %in% colnames(object[[]])
  )

  if (!is.null(split.by)) {

    stopifnot(
      split.by %in% colnames(object[[]])
    )
  }

  md <- object[[]]

  # overall distribution
  if (is.null(split.by)) {

    vec <- md[[group.by]]

    tab <- table(vec)

    df <- data.frame(
      x = names(tab),
      fill = names(tab),
      n = as.numeric(tab),
      stringsAsFactors = FALSE
    )

  } else {

    # grouped distribution
    tab <- table(
      md[[split.by]],
      md[[group.by]]
    )

    df <- as.data.frame(tab)

    colnames(df) <- c(
      "x",
      "fill",
      "n"
    )
  }

  # proportions
  if (proportion) {

    if (is.null(split.by)) {

      df$value <- df$n / sum(df$n)

    } else {

      df$value <- ave(
        df$n,
        df$x,
        FUN = function(z) z / sum(z)
      )
    }

  } else {

    df$value <- df$n
  }

  # labels
  df$lab <- if (proportion) {

    scales::percent(
      df$value,
      accuracy = 0.1
    )

  } else {

    format(
      df$value,
      big.mark = ","
    )
  }

  # sorting
  if (sort && is.null(split.by)) {

    ord <- order(
      df$value,
      decreasing = decreasing
    )

    df <- df[ord, ]
  }

  df$x <- factor(
    df$x,
    levels = unique(df$x)
  )

  df$fill <- factor(
    df$fill,
    levels = unique(df$fill)
  )

  # colors
  if (is.null(cols)) {

    cols <- stats::setNames(
      cellpal(
        n = length(levels(df$fill)),
        preset = "base"
      ),
      levels(df$fill)
    )

  } else {

    if (is.null(names(cols))) {
      names(cols) <- levels(df$fill)
    }

    cols <- cols[levels(df$fill)]
  }

  # plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = x,
      y = value,
      fill = fill
    )
  ) +
    ggplot2::geom_col(
      position = position,
      width = width,
      color = if (border) border.col else NA,
      linewidth = border.size
    ) +
    ggplot2::scale_fill_manual(
      values = cols
    ) +
    ggplot2::labs(
      title = plot.ttl,
      x = x.lab %||% split.by %||% "",
      y = y.lab %||%
        ifelse(
          proportion,
          "Proportion",
          "Cells"
        ),
      fill = group.by
    )


  # labels (adaptive positioning)
  if (label) {

    if (proportion) {

      # proportions -> usually inside bars
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = lab),
          vjust = 0.5,
          size = label.size
        )

    } else {

      # counts -> above bars (your preferred behavior)
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = lab),
          vjust = -0.35,
          size = label.size
        )
    }
  }

  # flip
  if (flip) {

    p <- p +
      ggplot2::coord_flip()
  }

  # theme
  p +
    plot_theme(
      style = style,
      txtsize = txtsize,
      leg.pos = ifelse(
        is.null(split.by),
        "none",
        "right"
      ),
      x.ang = x.ang,
      ...
    )
}


#' Plot GnRH module hit distributions
#'
#' @description
#' Bar plots of GnRH module hits vs metadata variables
#'
#' @param data data.frame (usually object@meta.data)
#' @param x character. grouping variable (e.g. "total_hits_bin")
#' @param fill character. fill variable (e.g. "gnrh_status")
#' @param palette named vector of colors
#' @param type "count" or "fraction"
#' @param title plot title
#' @param rotate.x logical
#' @param txtsize theme text size
#' @param style theme style
#' @param ... Additional arguments passed to \code{plot_theme()}.
#'
#' @export
plot_gnrh_hits <- function(
    data,
    x,
    fill,
    palette = NULL,
    type = c("count", "fraction"),
    title = NULL,
    rotate.x = FALSE,
    txtsize = 12,
    style = "classic",
    ...
) {

  type <- match.arg(type)

  stopifnot(
    is.data.frame(data),
    x %in% colnames(data),
    fill %in% colnames(data)
  )

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[x]], fill = .data[[fill]])
  ) +
    ggplot2::geom_bar(
      position = if (type == "count") "stack" else "fill",
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::labs(
      x = x,
      y = if (type == "count") "Cell count" else "Fraction",
      title = title,
      fill = fill
    )

  if (type == "fraction") {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format()
    )
  }

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  if (isTRUE(rotate.x)) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  }

  p +
    plot_theme(
      style = style,
      txtsize = txtsize,
      leg.pos = "right",
      ...
    )
}



#' Visualize cells on a Seurat embedding
#'
#' Creates a customizable two-dimensional or three-dimensional embedding plot
#' for Seurat objects. This function wraps \code{Seurat::DimPlot()} and adds
#' improved color handling, automatic rasterization, custom legends, optional
#' cluster labels, figure-panel styling, dark mode, and support for plotting
#' multiple metadata variables.
#'
#' @param object A Seurat object.
#' @param group.by Metadata column used to color cells. If \code{NULL},
#' active identities are used. A character vector of multiple metadata columns
#' can be supplied to generate a patchwork of plots.
#' @param reduction Dimensional reduction to plot. Default is \code{"umap"}.
#' @param dims Numeric vector specifying dimensions to plot. Use two values for
#' 2D plots or three values for 3D plots. Default is \code{c(1, 2)}.
#' @param shuffle Logical; randomly shuffle plotting order of cells.
#' Default is \code{FALSE}.
#' @param raster Logical; rasterize points for large datasets. If \code{NULL},
#' rasterization is enabled automatically for objects with more than 100,000
#' cells. Default is \code{NULL}.
#' @param stroke.size Optional point stroke size passed to \code{Seurat::DimPlot()}.
#' @param raster.dpi Numeric vector of length 2 specifying rasterization
#' resolution. Default is \code{c(2048, 2048)}.
#' @param alpha Numeric point transparency. Default is \code{1}.
#' @param repel Logical; repel cluster labels when labels are drawn.
#' Default is \code{FALSE}.
#' @param n.cells Logical; append the number of cells per group to legend labels.
#' Default is \code{TRUE}.
#' @param label Logical; add cluster labels to the embedding.
#' Default is \code{FALSE}.
#' @param label.size Numeric label text size. Default is \code{4}.
#' @param label.face Font face for cluster labels. Default is \code{"plain"}.
#' @param cols Optional named color vector. Missing groups are assigned
#' \code{"gray70"}. If unnamed, colors are matched to group levels in order.
#' @param figplot Logical; generate a minimal figure-style plot with compact
#' arrow axes. Default is \code{FALSE}.
#' @param axes Logical; show embedding axes. Default is \code{TRUE}.
#' @param plot.ttl Optional plot title.
#' @param legend Logical; show legend. Default is \code{TRUE}.
#' @param leg.ttl Optional legend title. If \code{NULL}, \code{group.by} is used.
#' @param leg.ttl.size Numeric legend title size. Default is \code{txtsize}.
#' @param item.size Numeric legend item size. Default is \code{4}.
#' @param leg.pos Legend position passed to \code{\link{plot_theme}}.
#' Default is \code{"right"}.
#' @param leg.just Legend justification. Default is \code{"center"}.
#' @param leg.dir Legend direction. Default is \code{"vertical"}.
#' @param leg.size Numeric legend text size. Default is \code{10}.
#' @param leg.ncol Number of legend columns. If \code{NULL}, this is selected
#' automatically based on the number of groups.
#' @param item.border Logical; draw borders around legend keys.
#' Default is \code{TRUE}.
#' @param txtsize Base text size. Default is \code{12}.
#' @param pt.size Numeric point size. If \code{NULL}, a point size is chosen
#' automatically based on cell number and rasterization mode.
#' @param dark Logical; apply a dark theme. Default is \code{FALSE}.
#' @param total.cells Logical; append total cell number to the plot title.
#' Default is \code{FALSE}.
#' @param threeD Logical; generate an interactive 3D plot using \pkg{plotly}.
#' Default is \code{FALSE}.
#' @param style Theme style passed to \code{\link{plot_theme}}.
#' Default is \code{"test"}.
#' @param facet.bg Logical; show facet background in the applied theme.
#' Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{Seurat::DimPlot()} and
#' \code{\link{plot_theme}}.
#'
#' @return A \code{ggplot2} object for 2D plots, a patchwork object when
#' multiple \code{group.by} variables are supplied, or a plotly object for
#' 3D plots.
#'
#' @details
#' Cell identities are prepared internally to preserve factor-level order and
#' handle missing values as \code{"Unknown"}. Colors are strictly matched to
#' the groups present in the plot.
#'
#' If \code{figplot = TRUE}, a minimal publication-panel style is used and a
#' compact axis-arrow inset is added. If \code{dark = TRUE}, a dark theme is
#' applied after the main theme.
#'
#' @examples
#' \dontrun{
#' cellmap(obj, group.by = "seurat_clusters")
#'
#' cellmap(
#'   obj,
#'   group.by = c("gnrh_status", "gnrh_stage"),
#'   reduction = "umap",
#'   label = TRUE
#' )
#'
#' cellmap(
#'   obj,
#'   group.by = "gnrh_stage",
#'   dark = TRUE,
#'   legend = TRUE
#' )
#' }
#'
#' @export
cellmap <- function(
    object,
    group.by = NULL,
    reduction = "umap",
    dims = c(1,2),
    shuffle = FALSE,
    raster = NULL,
    stroke.size = NULL,
    #raster.dpi = c(512, 512),
    raster.dpi = c(2048, 2048),
    alpha = 1,
    repel = FALSE,
    n.cells = TRUE,
    label = FALSE,
    label.size = 4,
    label.face = "plain",
    cols = NULL,
    figplot = FALSE,
    axes = TRUE,
    plot.ttl = NULL,
    legend = TRUE,
    leg.ttl = NULL,
    leg.ttl.size = txtsize,
    item.size = 4,
    leg.pos = "right",
    leg.just = "center",
    leg.dir = "vertical",
    leg.size = 10,
    leg.ncol = NULL,
    item.border = TRUE,
    txtsize = 12,
    pt.size = NULL,
    dark = FALSE,
    total.cells = FALSE,
    threeD = FALSE,
    style = "test",
    facet.bg = FALSE,
    ...
) {

  if (!is.null(list(...)$theme)) {
    style <- list(...)$theme
  }

  # Helper: prepare object and colors
  .prepare_object <- function(obj, group, cols = NULL) {
    stopifnot(inherits(obj, "Seurat"))
    if (is.null(group)) group <- "ident"

    ## Extract values while preserving order
    if (group == "ident") {
      values <- as.character(Idents(obj))
      orig_levels <- levels(Idents(obj))
    } else {
      if (!group %in% colnames(obj@meta.data)) {
        stop(paste("Grouping column", group, "not found."))
      }
      values <- as.character(obj@meta.data[[group]])
      orig_levels <- if (is.factor(obj@meta.data[[group]]))
        levels(obj@meta.data[[group]])
      else
        unique(values)
    }

    ## NA handling
    values[is.na(values)] <- "Unknown"

    ## Preserve order, but only keep present levels
    levels_group <- intersect(orig_levels, unique(values))
    if ("Unknown" %in% values && !("Unknown" %in% levels_group)) {
      levels_group <- c(levels_group, "Unknown")
    }

    ## Re-factor with preserved order
    obj@meta.data[[group]] <- factor(values, levels = levels_group)
    Idents(obj) <- obj@meta.data[[group]]

    ## Colors (STRICTLY match present levels)
    if (is.null(cols)) {
      cols <- cellpal(n = length(levels_group), preset = "base")
      names(cols) <- levels_group
      if ("Unknown" %in% levels_group) cols["Unknown"] <- "gray70"
    } else {
      if (is.null(names(cols))) {
        cols <- setNames(cols[seq_along(levels_group)], levels_group)
      }
      missing <- setdiff(levels_group, names(cols))
      if (length(missing)) cols[missing] <- "gray70"
      cols <- cols[levels_group]
    }

    list(obj = obj, cols = cols, levels_group = levels_group)
  }

  # 3D plotting helper
  .plot_3D <- function(obj, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells){
    emb <- obj@reductions[[reduction]]@cell.embeddings
    df <- data.frame(x=emb[, dims[1]], y=emb[, dims[2]], z=emb[, dims[3]], cluster=Idents(obj))
    hover_labels <- if(n.cells){
      tbl <- table(df$cluster)
      paste0(df$cluster, " (", tbl[as.character(df$cluster)], ")")
    } else as.character(df$cluster)

    p3d <- plotly::plot_ly(df, x=~x, y=~y, z=~z, color=~cluster, colors=cols,
                           type="scatter3d", mode="markers",
                           marker=list(size=pt.size, opacity=alpha, line=list(width=0)),
                           text=hover_labels, hoverinfo="text")
    if(label){
      centers <- df %>% dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), z=median(z))
      p3d <- p3d %>% plotly::add_text(data=centers, x=~x, y=~y, z=~z, text=~cluster, textposition="top center")
    }
    p3d
  }

  if(length(group.by) > 1){
    plots <- lapply(group.by, function(g){
      cellmap(object = object, group.by = g, shuffle = shuffle,raster = raster,stroke.size = stroke.size,alpha = alpha,
              repel = repel,reduction = reduction,dims = dims,n.cells = n.cells,label = label,label.size = label.size,
              label.face = label.face,cols = cols,figplot = figplot,plot.ttl = g,legend = legend,leg.ttl = g,item.size = item.size,
              leg.pos = leg.pos,leg.just = leg.just,leg.dir = leg.dir,leg.ncol = leg.ncol,txtsize = txtsize,item.border = item.border,
              pt.size = pt.size,dark = dark,total.cells = total.cells,threeD = threeD,style = style,facet.bg = facet.bg,
              ...
      )
    })
    return(patchwork::wrap_plots(plots))
  }
  # Prepare object & colors
  prep <- .prepare_object(object, group.by, cols)
  object <- prep$obj
  cols <- prep$cols
  levels_group <- prep$levels_group
  if(is.null(leg.ncol)) leg.ncol <- if(length(levels_group) > 30) 2 else 1

  # Validate reduction
  if(!(reduction %in% names(object@reductions))){
    stop(paste0("Reduction '", reduction, "' not found. Available: ", paste(names(object@reductions), collapse=", ")))
  }
  emb <- object@reductions[[reduction]]@cell.embeddings
  if(max(dims) > ncol(emb)) stop("Selected dims exceed available dimensions in reduction.")

  # 3D plotting
  if(threeD || length(dims) == 3) return(.plot_3D(object, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells))

  # 2D plotting

  nb.cells <- ncol(object)

  if (is.null(raster)) {
    raster <- nb.cells > 1e5
  }

  if (is.null(pt.size)) {

    if (!raster) {

      # point mode (classic scaling)
      pt.size <- max(0.3, 1.5 / log10(nb.cells))

    } else {

      # raster mode (FIXED optimal density)
      pt.size <- dplyr::case_when(
        nb.cells < 1e5 ~ 0.1,
        nb.cells < 2e5 ~ 0.5,
        TRUE ~ 2
      )

    }
  }

  plt <- Seurat::DimPlot(object, group.by=group.by, shuffle=shuffle, raster=raster, pt.size=pt.size,
                         repel=repel, alpha=alpha, reduction=reduction, dims=dims, raster.dpi=raster.dpi,...)

  # reset Seurat's forced theme_classic
  plt <- plt + ggplot2::theme_void()

  present_levels <- levels(droplevels(object@active.ident))

  if (n.cells) {
    cell.nb <- table(object@active.ident)[present_levels]
    clust.lab <- paste0(present_levels, " (", cell.nb, ")")
  } else {
    clust.lab <- present_levels
  }

  cols_use <- cols[present_levels]

  leg.ttl <- if(is.null(leg.ttl)) group.by else leg.ttl

  plt <- plt + scale_color_manual(
    breaks = present_levels,
    labels = clust.lab,
    values = cols_use
  )

  if (legend) {
    plt <- plt +
      guides(
        color = guide_legend(
          override.aes = if (item.border)
            list(size = item.size, shape = 21,
                 color = if (dark) "white" else "black",
                 stroke = 0.3, fill = unname(cols)
            )
          else
            list(size = item.size),
          ncol = leg.ncol,
          title = leg.ttl,
          keyheight = unit(0.35, "cm"),
          keywidth  = unit(0.35, "cm")
        )
      )
  } else {
    plt <- plt + guides(color = "none")
  }


  # Plot title with total cells
  if(total.cells){
    plot.ttl <- paste0(plot.ttl, " (n=", format(ncol(object), big.mark=","), ")")
  }
  #if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl)
  if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl) else plt <- plt + labs(title = NULL)

  # Set default legend title size
  if (is.null(leg.ttl.size)) leg.ttl.size <- txtsize

  # Apply plot_theme using do.call
  theme_args <- list(
    style = style,
    txtsize = txtsize,
    leg.size  = leg.size,
    leg.pos   = leg.pos,
    leg.dir   = leg.dir,
    leg.ttl   = leg.ttl,
    leg.ttl.size = leg.ttl.size,
    facet.bg  = facet.bg
  )

  if(figplot){
    # Warn if the user specified a non-classic theme
    if(!missing(style) && style != "classic"){
      warning(sprintf(
        "figplot = TRUE ignores custom themes (style = '%s'). Use figplot = FALSE for full theming.",
        style
      ))
    }

    # Apply plot_theme for figplot figure
    plt <- plt &
      do.call(plot_theme, c(
        theme_args,
        list(
          x.ttl = FALSE,
          ticks = FALSE,
          line = FALSE,
          border = FALSE,
          grid.major = FALSE,
          grid.minor = FALSE,
          panel.fill = "white"
        ),
        list(...)
      ))
    text_col <- "black"
  } else {
    # Apply plot_theme normally
    plt <- plt & do.call(plot_theme, c(theme_args, list(...)))
  }


  # Add cluster labels if requested
  if(label){
    umap_data <- dplyr::tibble(x=emb[, dims[1]], y=emb[, dims[2]], cluster=as.character(object@active.ident)) %>%
      dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups="drop")
    plt <- plt + ggrepel::geom_text_repel(
      data=umap_data, aes(x, y, label=cluster),
      color = if(dark) "white" else "black",
      fontface = label.face,
      bg.color = if(dark) "#3A3A3A" else "grey95",
      bg.r = 0.1, size = label.size, seed = 42
    )
  }

  # figplot arrow axes (minimal figure)
  if(figplot){
    x.lab.reduc <- plt$labels$x %||% paste0(toupper(reduction), dims[1])
    y.lab.reduc <- plt$labels$y %||% paste0(toupper(reduction), dims[2])
    plt <- plt & Seurat::NoAxes()
    L <- 0.12
    axis.df <- data.frame(x0=c(0,0), y0=c(0,0), x1=c(L,0), y1=c(0,L))
    axis.plot <- ggplot(axis.df) +
      geom_segment(aes(x=x0, y=y0, xend=x1, yend=y1), linewidth=0.4, lineend="round") +
      xlab(x.lab.reduc) + ylab(y.lab.reduc) +
      coord_fixed() + theme_classic(base_size=txtsize) +
      theme(plot.background=element_rect(fill="transparent", colour=NA),
            panel.background=element_rect(fill="transparent", colour=NA),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            panel.border=element_blank(),
            axis.title=element_text(size=txtsize, face="plain"),
            plot.margin=ggplot2::margin(0,0,0,0))
    figure.layout <- c(patchwork::area(t=1,l=1,b=11,r=11), patchwork::area(t=10,l=1,b=11,r=2))
    return(plt + axis.plot + patchwork::plot_layout(design=figure.layout))
  }

  if (style == "test") {
    plt <- plt + plot_theme(style = "test", xy.val = F)
  }

  if(!legend) plt <- plt & Seurat::NoLegend()
  if(!axes) plt <- plt & Seurat::NoAxes()

  if (dark) {
    plt <- plt + dark_theme(style = "true", txtsize = txtsize, axis.text = FALSE, axis.title = FALSE)
  }

  if (style == "test" && axes && !figplot) {
    plt <- plt + theme_test(base_size = txtsize) +
      theme(axis.text  = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(color = "black", size = txtsize)
      )
  }

  # Compute axis titles from reduction and dims
  xlab <- paste0(toupper(reduction), dims[1])
  ylab <- paste0(toupper(reduction), dims[2])
  plt <- plt + labs(x = xlab, y = ylab)

  plt
}



#' Plot GnRHcell metadata on cell embeddings
#'
#' Visualizes GnRHcell metadata on a Seurat dimensional reduction.
#'
#' @param object A Seurat object processed by GnRHcell.
#' @param group.by Metadata variable to plot. Use \code{"all"} to plot
#' \code{gnrh_status}, \code{gnrh_truth}, and \code{gnrh_stage}.
#' @param reduction Dimensional reduction to use.
#' @param cols Optional color palette.
#' @param style Theme style.
#' @param dark Logical; use dark theme.
#' @param ncol Number of columns when plotting multiple panels.
#'
#' @return A ggplot2 or patchwork object.
#' @export
plot_gnrh_embedding <- function(object,
                                group.by = "all",
                                reduction = "umap",
                                cols = NULL,
                                style = "classic",
                                dark = FALSE,
                                ncol = NULL) {

  #allowed <- c("gnrh_status","gnrh_class","gnrh_truth","gnrh_stage")
  allowed <- c("gnrh_status","gnrh_truth","gnrh_stage")

  if (identical(group.by, "all")) group.by <- allowed
  stopifnot(all(group.by %in% allowed))

  cmap <- list(
    gnrh_status = gnrh_colors("status"),
    #gnrh_class  = gnrh_colors("class"),
    gnrh_truth  = gnrh_colors("truth"),
    gnrh_stage  = gnrh_colors("stage")
  )

  if (length(group.by) == 1) {
    return(cellmap(object, group.by = group.by,
                   reduction = reduction,
                   cols = cols %||% cmap[[group.by]],
                    style = style,
                   dark = dark)
    )
  }

  plots <- lapply(group.by, \(g)
                  cellmap(object, group.by = g,
                           reduction = reduction,
                           cols = cmap[[g]],
                          style = style,
                          dark = dark)
  )

  patchwork::wrap_plots(plots, ncol = ncol)
}


#' Cell Feature Plot for Seurat Objects
#'
#' A flexible wrapper around Seurat's `FeaturePlot` to visualize gene expression or metadata
#' features in a Seurat object. Supports custom color palettes, viridis, and hotspot/rainbow palettes.
#'
#' @param object A `Seurat` object.
#' @param features Character vector of features (genes or metadata columns) to plot.
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
#' @examples
#' \dontrun{
#' cellfeature(seurat_obj, features = "POMC")
#'
#' cellfeature(
#'   seurat_obj,
#'   features = c("POMC", "NPY"),
#'   merge.leg = TRUE
#' )
#'
#' cellfeature(
#'   seurat_obj,
#'   features = c("POMC", "NPY"),
#'   blend = TRUE
#' )
#' }
#' @export
#'
cellfeature <- function(
    object,
    features,
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
        oob = scales::censor,     # force < cutoff -> NA
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


#' Resolve GnRH feature presets
#'
#' @param type Feature preset: \code{"core"}, \code{"modules"}, or \code{"all"}.
#'
#' @return Character vector of feature names.
#' @keywords internal
#' @noRd
.gnrh_features <- function(type = "all") {

  switch(
    type,
    core = c("gnrh_raw", "gnrh_score"),
    modules = c("gnrh_core_hits", "gnrh_mig_hits", "gnrh_neuro_hits"),
    all = c("gnrh_raw", "gnrh_expr", "gnrh_score", "gnrh_knn"),
    stop("Unknown gnrh feature type")
  )
}


#' Plot GnRH feature maps
#'
#' Wrapper around \code{\link{cellfeature}} for GnRHcell diagnostic features.
#'
#' @param object A Seurat object processed by GnRHcell.
#' @param features Features to plot, or \code{"all"}.
#' @param feature_type Optional preset: \code{"core"}, \code{"modules"}, or \code{"all"}.
#' @param reduction Dimensional reduction to use.
#' @param ncol Number of columns.
#' @param style Theme style.
#' @param merge.leg Logical; merge legends across panels.
#' @param split.by Optional metadata column for splitting.
#' @param title Optional title.
#' @param ... Additional arguments passed to \code{\link{cellfeature}}.
#'
#' @return ggplot2 or patchwork object.
#' @export
plot_gnrh_feature <- function(
    object,
    features = "all",
    feature_type = NULL,
    reduction = "umap",
    ncol = 4,
    style = "classic",
    merge.leg = FALSE,
    split.by = NULL,
    title = NULL,
    ...
) {

  stopifnot(inherits(object, "Seurat"))

  # feature resolution logic (single source of truth)
  if (!is.null(feature_type)) {
    features <- .gnrh_features(feature_type)
  }

  if (identical(features, "all")) {
    features <- .gnrh_features("all")
  }

  features <- unique(as.character(features))

  # delegate rendering to cellfeature()
  p <- cellfeature(
    object = object,
    features = features,
    reduction = reduction,
    ncol = ncol,
    merge.leg = merge.leg,
    split.by = split.by,
    style = style,
    blend = FALSE,
    txtsize = 10,
    ...
  )

  return(p)
}



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
#' @param txtsize Base text size passed to \code{\link{plot_theme}}.
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
#' @param style Theme style passed to \code{\link{plot_theme}}.
#' Default is \code{"classic"}.
#' @param x.face,y.face Logical; italicize x- or y-axis labels.
#' Default is \code{FALSE}.
#' @param x.ttl,y.ttl Logical; show x- or y-axis titles.
#' Default is \code{FALSE}.
#' @param dot.outline Logical; draw outlines around dots.
#' Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{Seurat::DotPlot()} and
#' \code{\link{plot_theme}}.
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
celldot <- function(object,
                    features,
                    group.by="seurat_clusters",
                    #th.cols="Reds",
                    th.cols = "RdYlBu",
                    rev.th.cols=TRUE,
                    dot.scale=4,
                    x.ang=90,
                    vjust.x=NULL,
                    hjust.x=NULL,
                    flip=FALSE,
                    txtsize=12,
                    title=NULL,
                    leg.size=10,
                    leg.ttl.size = 10,
                    leg.pos="right",
                    leg.just="bottom",
                    leg.hjust=FALSE,
                    x.axis.pos="bottom",
                    style="classic",
                    x.face=FALSE,
                    y.face=FALSE,
                    x.ttl=FALSE,
                    y.ttl=FALSE,
                    dot.outline=FALSE,
                    ...) {

  stopifnot(inherits(object,"Seurat"))
  object <- Seurat::SetIdent(object, value=group.by)
  features <- unique(features)

  pal <- RColorBrewer::brewer.pal(9, th.cols)
  if (rev.th.cols) pal <- rev(pal)

  outline_col <- if (dot.outline) "gray60" else NA
  outline_stroke <- if (dot.outline) 0.5 else 0

  plt <- suppressWarnings({
    suppressMessages({
      Seurat::DotPlot(object, features=features, dot.scale=dot.scale, ...)
    })
  }) +
    scale_color_gradientn(colors=pal, oob=scales::squish) +
    geom_point(aes(size=pct.exp), shape=21, colour=outline_col, stroke=outline_stroke) +
    labs(title=title, color="Average\nExpression", size="Percent\nExpressed") +
    plot_theme(style=style, txtsize=txtsize, x.ang=x.ang, leg.size = leg.size,leg.ttl.size = leg.ttl.size,
               x.hjust=hjust.x, x.vjust=vjust.x, xy.val=TRUE, x.lab=TRUE, y.lab=TRUE, ...) +
    theme(
      axis.text.x = if (x.face || (flip && y.face)) element_text(face="italic") else element_text(),
      axis.text.y = if (y.face || (flip && x.face)) element_text(face="italic") else element_text(),
      axis.title = element_blank(),
      legend.spacing.y = unit(0.05, "cm"),
      legend.spacing.x = unit(0.05, "cm"),
      legend.box.spacing = unit(0.05, "cm"),
      legend.margin = margin(2,2,2,2)
    )

  if (flip) plt <- plt + coord_flip()
  if (is.list(features)) plt <- plt + theme(strip.text.x=element_text(angle=45))

  # Legend positioning
  # Legend positioning outside plot, bottom-right
  if (!is.null(leg.pos)) {
    if (leg.pos == "right") {
      plt <- plt + theme(
        legend.position = "right",
        legend.justification = c("right","bottom"),
        legend.box.just = "right",
        legend.box.margin = ggplot2::margin(0,0,0,0)
      )
    } else if (leg.pos == "left") {
      plt <- plt + theme(
        legend.position = "left",
        legend.justification = c("left","bottom"),
        legend.box.just = "left",
        legend.box.margin = ggplot2::margin(0,0,0,0)
      )
    } else if (leg.pos == "top") {
      plt <- plt + theme(
        legend.position = "top",
        legend.justification = c("right","top"),
        legend.box.just = "right"
      )
    } else if (leg.pos == "bottom") {
      plt <- plt + theme(
        legend.position = "bottom",
        legend.justification = c("right","bottom"),
        legend.box.just = "right"
      )
    }
  }

  # Keep original guides for color and size
  guide_color <- guide_colorbar(frame.colour="black", ticks.colour="black")
  guide_size <- guide_legend(override.aes=list(shape=21, colour=outline_col, fill="black"))
  guide_color$order <- 1
  guide_size$order  <- 2
  plt <- plt + guides(color=guide_color, size=guide_size)

  plt
}



#' Create violin plots for Seurat features
#'
#' Generates publication-ready violin plots for one or more genes or metadata
#' features from a Seurat object. This function wraps \code{Seurat::VlnPlot()}
#' and adds shared y-axis scaling, median markers, optional statistical tests,
#' flexible layout, global axis labels, and GnRHcell theme support.
#'
#' @param obj A Seurat object.
#' @param features Character vector of genes or metadata columns to plot.
#' @param ncol Number of columns in the patchwork layout. If \code{NULL},
#' selected automatically from the number of features.
#' @param shared.y Logical; use the same y-axis range for all features.
#' Default is \code{FALSE}.
#' @param ttl.pos Subplot title position. One of \code{"left"},
#' \code{"center"}, or \code{"right"}.
#' @param group.by Metadata column used to group cells.
#' Default is \code{"seurat_clusters"}.
#' @param split.by Optional metadata column used to split violin plots.
#' @param stack Logical; arrange plots in a single column.
#' Default is \code{FALSE}.
#' @param assay Assay used to retrieve gene expression.
#' Default is \code{"RNA"}.
#' @param slot Assay slot used for expression values. One of \code{"data"},
#' \code{"counts"}, or \code{"scale.data"}. Default is \code{"data"}.
#' @param log Logical; if \code{TRUE} and \code{slot = "counts"}, use
#' log-transformed counts. Default is \code{FALSE}.
#' @param cols Optional named color vector for groups.
#' @param med Logical; overlay median points on violin plots.
#' Default is \code{FALSE}.
#' @param med.size Size of median points. Default is \code{1}.
#' @param pt.size Size of jittered cells. Use \code{0} to hide points.
#' Default is \code{0}.
#' @param border.size Violin outline linewidth. Default is \code{0.1}.
#' @param style Theme style passed to \code{\link{plot_theme}}.
#' Default is \code{"classic"}.
#' @param leg.pos Legend position. Default is \code{"none"}.
#' @param x.ang X-axis text angle. Default is \code{45}.
#' @param title Optional global title for the combined plot.
#' @param rm.subttl Logical; remove individual feature titles.
#' Default is \code{FALSE}.
#' @param flip Logical; flip x and y axes. Default is \code{FALSE}.
#' @param auto.resize Logical; add dynamic width and height attributes to
#' the returned object. Default is \code{TRUE}.
#' @param ylab.global Optional global y-axis label. If \code{NULL}, inferred
#' from \code{slot}.
#' @param xlab.global Optional global x-axis label. Default is blank.
#' @param add.stats Logical; compute statistical comparisons.
#' Default is \code{FALSE}.
#' @param show.pval Logical; display p-values or significance labels.
#' Default is \code{FALSE}.
#' @param pairwise Logical; perform pairwise Wilcoxon tests between groups.
#' If \code{FALSE}, a global Kruskal-Wallis test is used.
#' Default is \code{FALSE}.
#' @param pval.label Label type passed to statistical annotation utilities.
#' Default is \code{"p.signif"}.
#' @param p.display How to display statistical results. One of
#' \code{"stars"}, \code{"value"}, or \code{"none"}.
#' Default is \code{"stars"}.
#' @param txtsize Base text size. Default is \code{10}.
#' @param subttl.size Feature subplot title size. Default is \code{10}.
#' @param ... Additional arguments passed to \code{Seurat::VlnPlot()} and
#' \code{\link{plot_theme}}.
#'
#' @return A patchwork/cowplot object containing the violin plots.
#'
#' @details
#' Features can be either assay genes or metadata columns. Missing features are
#' silently ignored; an error is raised if none of the requested features are
#' found.
#'
#' When \code{add.stats = TRUE} and \code{show.pval = TRUE}, statistical
#' annotations are added if at least two groups are available. Pairwise tests
#' use Wilcoxon rank-sum tests, while the global comparison uses a
#' Kruskal-Wallis test.
#'
#' @examples
#' \dontrun{
#' cellviolin(
#'   obj,
#'   features = c("GNRH1", "ISL1", "DLX5"),
#'   group.by = "gnrh_stage",
#'   pt.size = 0.1,
#'   ncol = 3
#' )
#'
#' cellviolin(
#'   obj,
#'   features = c("GNRH1", "ISL1"),
#'   group.by = "gnrh_status",
#'   add.stats = TRUE,
#'   show.pval = TRUE,
#'   pairwise = TRUE
#' )
#' }
#'
#' @export
cellviolin <- function(
    obj, features,
    ncol = NULL,
    shared.y = FALSE,
    ttl.pos = c("left","center","right"),
    group.by = "seurat_clusters",
    split.by = NULL,
    stack = FALSE,
    assay = "RNA",
    slot = "data",
    log = FALSE,
    cols = NULL,
    med = FALSE,
    med.size = 1,
    pt.size = 0,
    border.size = 0.1,
    style = "classic",
    leg.pos = "none",
    x.ang = 45,
    title = NULL,
    rm.subttl = FALSE,
    flip = FALSE,
    auto.resize = TRUE,
    ylab.global = NULL,
    xlab.global = NULL,
    add.stats = FALSE,
    show.pval = FALSE,
    pairwise = FALSE,
    pval.label = "p.signif",
    p.display = c("stars", "value", "none"),
    txtsize = 10,
    subttl.size = 10,
    ...
) {

  stopifnot(inherits(obj, "Seurat"))
  if (length(features) == 0) stop("features must be provided.")
  ttl.pos <- match.arg(ttl.pos)
  p.display <- match.arg(p.display)


  # determine ncol
  if (is.null(ncol)) {
    ncol <- if (stack) 1 else max(1, ceiling(sqrt(length(features) / 1.5)))
  }

  # Set identities if grouping
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(obj@meta.data))
      stop(paste(group.by, "not found in metadata."))
    Idents(obj) <- group.by
  }

  # Determine which features exist
  f.expr <- intersect(features, rownames(obj[[assay]]))
  f.meta <- intersect(features, colnames(obj@meta.data))
  features <- unique(c(f.expr, f.meta))
  if (!length(features)) stop("No features found in assay or metadata.")

  # Shared y-scale
  ymax <- NULL
  if (shared.y) {
    vals <- c(
      if (length(f.expr)) as.numeric(Seurat::GetAssayData(obj, assay, slot)[f.expr, ]),
      if (length(f.meta)) as.numeric(as.matrix(obj@meta.data[, f.meta, drop = FALSE]))
    )
    vals <- vals[is.finite(vals)]
    if (length(vals)) ymax <- max(vals)
  }

  # Colors
  if (is.null(cols)) {
    g <- tryCatch(unique(obj[[group.by]][, 1]), error = \(e) NULL)
    cols <- cellpal(n = length(g), preset = "base")
    #cols <- scales::hue_pal()(if (is.null(g)) 8 else length(g))
    #cols <- ggpubr::get_palette("npg", length(g))
    names(cols) <- g
  }

  # Helper to build a single violin with stats
  vln <- function(f) {

    if (log && slot == "counts") {
      obj[[assay]]@data[f, ] <- log1p(obj[[assay]]@counts[f, ])
      slot_use <- "data"
    } else {
      slot_use <- slot
    }

    df <- Seurat::FetchData(obj, vars = c(group.by, f))
    names(df)[2] <- "value"
    df[[group.by]] <- factor(df[[group.by]]) # ensure factor

    if (!is.null(split.by) && !is.null(cols)) {
      split_levels <- levels(obj@meta.data[[split.by]])
      cols <- cols[split_levels]
    }

    # Base plot
    # p <- ggplot(df, aes_string(group.by, "value", fill = group.by)) +
    #   geom_violin(scale = "width", color = "black", size = border.size) +
    #   scale_fill_manual(values = cols)

    p <- suppressWarnings({
      suppressMessages({Seurat::VlnPlot(
        obj,
        features = f,
        group.by = group.by,
        split.by = split.by,
        assay = assay,
        slot = slot,
        pt.size = pt.size,
        cols = cols,
        ...
      ) + scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.25))
      )
      })
    })

    # Points
    #if (pt.size > 0) p <- p + geom_jitter(width = 0.1, size = pt.size, alpha = 0.6)

    # Add statistics

    if (add.stats && show.pval) {

      if (nlevels(df[[group.by]]) > 1) {

        df$.grp <- df[[group.by]]

        y_max <- max(df$value, na.rm = TRUE)
        y_step <- (ymax %||% y_max) * 0.08

        if (pairwise && nlevels(df$.grp) > 1) {

          cmp <- utils::combn(levels(df$.grp), 2, simplify = FALSE)

          stat_df <- ggpubr::compare_means(
            value ~ .grp,
            data = df,
            method = "wilcox.test",
            comparisons = cmp
          )

          stat_df$y.position <- y_max + seq_len(nrow(stat_df)) * y_step

          # format labels
          if (p.display == "stars") {
            stat_df$label <- stat_df$p.signif
          }

          if (p.display == "value") {
            stat_df$label <- paste0("p = ", signif(stat_df$p, 3))
          }

          if (p.display == "none") {
            stat_df$label <- ""
          }

          p <- p + ggpubr::stat_pvalue_manual(
            stat_df,
            label = "label",
            y.position = "y.position",
            tip.length = 0.02,
            size = txtsize * 0.25
          )

        } else {

          # global test
          stat_df <- ggpubr::compare_means(
            value ~ .grp,
            data = df,
            method = "kruskal.test"
          )

          if (p.display == "stars") {
            label_text <- stat_df$p.signif
          }

          if (p.display == "value") {
            label_text <- paste0("p = ", signif(stat_df$p, 3))
          }

          if (p.display == "none") {
            label_text <- ""
          }

          p <- p + annotate(
            "text",
            x = 1,
            y = y_max + y_step,
            label = label_text,
            hjust = 0,
            size = txtsize * 0.25
          )
        }
      }
    }

    # Titles
    p <- if (!rm.subttl) p + labs(title = f) else p + labs(title = NULL)

    # Make feature titles bold + italic
    p <- p + theme(
      plot.title = element_text(face = "bold.italic", size = subttl.size)
    )

    .style_layers <- function(
    p,
    violin_lw = 0.15,
    point_size = NULL,
    jitter_width = NULL
    ) {
      for (i in seq_along(p$layers)) {

        layer <- p$layers[[i]]

        # Violin outline
        if (inherits(layer$geom, "GeomViolin")) {
          layer$aes_params$linewidth <- violin_lw
        }

        # Points (Seurat uses GeomPoint + position_jitterdodge)
        if (inherits(layer$geom, "GeomPoint")) {

          if (!is.null(point_size)) {
            layer$aes_params$size  <- point_size
            layer$aes_params$alpha <- 0.6
          }

          if (!is.null(jitter_width) &&
              inherits(layer$position, "PositionJitterdodge")) {
            layer$position$width <- jitter_width
          }
        }

        p$layers[[i]] <- layer
      }
      p
    }

    p <- .style_layers(
      p,
      violin_lw  = border.size,
      point_size = if (pt.size > 0) pt.size else NULL,
      jitter_width = 0.08
    )

    # Theme & formatting
    p <- p + plot_theme(style = style, txtsize = txtsize, x.ang = x.ang,
                        leg.pos = leg.pos, x.ttl = FALSE, ttl.pos = ttl.pos,...) +
      theme(
        plot.title = element_text(face = "bold.italic", size = subttl.size)
      )

    # Force final legend position
    if (!is.null(leg.pos)) {
      p <- p + theme(legend.position = leg.pos)
    }

    if (med) p <- p + stat_summary(fun = median, geom = "point", shape = 3, size = med.size)
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    if (flip) p <- p + coord_flip()
    p + ylab(NULL)
  }

  # number of cols
  if (is.null(ncol))
    ncol <- max(1, ceiling(sqrt(length(features) / 1.5)))

  plist <- lapply(features, vln)
  total <- length(plist)

  # Only show x-axis on bottom plots
  bottom <- sapply(1:ncol, \(i) max(seq(i, total, by = ncol)))
  for (i in seq_along(plist)) {
    if (!(i %in% bottom)) {
      plist[[i]] <- plist[[i]] +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }
  }

  # Layout
  combo <- patchwork::wrap_plots(plist, ncol = ncol) +
    patchwork::plot_layout(guides = "collect")

  if (!is.null(title)) {
    combo <- combo +
      patchwork::plot_annotation(
        title = title,
        theme = theme(title = element_text(face = "bold"))
      )
  }

  # Global labels
  auto_y <- switch(slot,
                   data = "Expression level",
                   counts = "Raw counts",
                   scale.data = "Scaled expression",
                   "Expression level")

  ylab <- if (is.null(ylab.global)) auto_y else ylab.global
  xlab <- if (is.null(xlab.global)) "" else xlab.global

  plt <- cowplot::ggdraw(combo) +
    cowplot::draw_label(ylab, x = -0.01, y = 0.55, angle = 90, size = txtsize) +
    cowplot::draw_label(xlab, x = 0.5, y = 0.02, size = txtsize) +
    theme(plot.margin = ggplot2::margin(15, 15, 15, 15))

  # Auto resize attributes
  if (auto.resize) {
    ng <- length(unique(obj[[group.by]][, 1]))
    attr(plt, "dynamic_width") <- 6 + ng * 0.3
    attr(plt, "dynamic_height") <- 4 + length(features) * 0.25
  }

  plt
}




#' Generate a GnRHcell diagnostic report
#'
#' Creates a multi-panel quality-control dashboard summarizing GnRH
#' detection results, diagnostic scores, threshold behavior, ROC
#' performance, module-hit enrichment, and class composition.
#'
#' @param object A Seurat object processed with \code{\link{detect_gnrh}}
#' and \code{\link{gnrh_diagnostics}}. The object must contain
#' \code{object@misc$gnrh} with diagnostic outputs.
#' @param style Plot theme style passed to \code{\link{plot_theme}}.
#' Default is \code{"test"}.
#' @param verbose Logical; print progress messages.
#' Default is \code{TRUE}.
#'
#' @return A patchwork object containing the GnRHcell diagnostic report.
#'
#' @details
#' The report includes:
#' \itemize{
#'   \item signal landscape: \code{GNRH1} expression versus composite score
#'   \item score density distribution by GnRH class
#'   \item threshold performance curve
#'   \item ROC curve when two truth classes are available
#'   \item module-hit signal versus score
#'   \item GnRH class distribution
#'   \item enrichment across total module hits
#' }
#'
#' The function requires diagnostic information generated by
#' \code{\link{gnrh_diagnostics}}. If ROC computation is not possible
#' because fewer than two truth classes are present, the ROC panel is
#' replaced by an empty placeholder.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{gnrh_diagnostics}},
#' \code{\link{plot_gnrh_distribution}}
#'
#' @examples
#' \dontrun{
#' obj <- detect_gnrh(obj)
#' obj <- gnrh_diagnostics(obj)
#' p <- gnrh_report(obj)
#' print(p)
#' }
#'
#' @export
gnrh_report <- function(object, style = "test", verbose = TRUE) {

  log <- .msg(verbose)
  log("Generating GnRH QC report")

  pkg_ok <- requireNamespace("patchwork", quietly = TRUE)
  if (!pkg_ok) stop("patchwork required")

  g <- object@misc$gnrh
  if (is.null(g)) stop("Missing @misc$gnrh - run detect_gnrh + gnrh_diagnostics first")

  diag <- g$diagnostics
  curve <- g$threshold_curve


  # SAFE ROC PREP
  # roc_plot <- NULL
  # if (!is.null(diag$truth) && length(unique(diag$truth)) == 2) {
  #
  #   roc_obj <- pROC::roc(diag$truth, diag$expr, quiet = TRUE)
  #
  #   df_roc <- data.frame(
  #     fpr = 1 - roc_obj$specificities,
  #     tpr = roc_obj$sensitivities
  #   )
  #
  #   roc_plot <-
  #     ggplot2::ggplot(df_roc, ggplot2::aes(fpr, tpr)) +
  #     ggplot2::geom_line(color = "red") +
  #     ggplot2::geom_abline(linetype = "dashed") +
  #     ggplot2::labs(
  #       title = paste0("ROC (AUC=", round(pROC::auc(roc_obj), 3), ")"),
  #       x = "FPR",
  #       y = "TPR"
  #     ) +
  #     plot_theme(style = style)
  # }

  roc_plot <- NULL

  if (!is.null(diag$truth) && length(unique(stats::na.omit(diag$truth))) == 2) {

    roc_obj <- pROC::roc(
      response = diag$truth,
      predictor = diag$score,
      levels = c("neg", "pos"),
      quiet = TRUE
    )

    df_roc <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    )

    roc_plot <-
      ggplot2::ggplot(df_roc, ggplot2::aes(fpr, tpr)) +
      ggplot2::geom_line(color = "#EF476F") +
      ggplot2::geom_abline(linetype = "dashed") +
      ggplot2::labs(
        title = paste0("ROC (AUC=", round(pROC::auc(roc_obj), 3), ")"),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      plot_theme(style = style)
  }

  # PANEL 1: expression vs score
  df <- object@misc$gnrh$diagnostics

  pct <- round(mean(df$status == "pos", na.rm = TRUE) * 100, 1)
  title_txt <- paste0("Signal landscape (", pct, "% GnRH+)")

  p1 <- ggplot2::ggplot(df, ggplot2::aes(expr, score)) +
    ggplot2::geom_point(ggplot2::aes(color = class),
      alpha = 0.8,size = 1.5) +
    ggplot2::labs(title = title_txt,
      x = "GNRH1 expression",
      y = "GnRH composite score"
    ) +
    ggplot2::scale_color_manual(values = gnrh_colors(type = "class")) +
    plot_theme(style = "test", leg.pos = c(0.2, 0.8))

  # PANEL 2: distribution
  p2 <- ggplot2::ggplot(diag, ggplot2::aes(score, fill = class)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::labs(title = "Score distribution", x = "Score", y = "Density") +
    ggplot2::scale_fill_manual(values = gnrh_colors(type = "class")) +
    plot_theme(style = "test", leg.pos = c(0.3, 0.8))

  # PANEL 3: threshold curve
  if (!is.null(curve)) {

    best <- curve$threshold[which.max(curve$F1)]

    p3 <- ggplot2::ggplot(curve, ggplot2::aes(threshold)) +
      ggplot2::geom_line(ggplot2::aes(y = sensitivity, color = "Sensitivity")) +
      ggplot2::geom_line(ggplot2::aes(y = specificity, color = "Specificity")) +
      ggplot2::geom_line(ggplot2::aes(y = F1, color = "F1")) +
      ggplot2::geom_vline(xintercept = best, linetype = "dashed") +
      ggplot2::labs(title = "Threshold curve", x = "Threshold", y = "Metric", color = NULL) +
      ggplot2::scale_color_manual(values = c(
        Sensitivity = "blue",
        Specificity = "darkgreen",
        F1 = "red"
      )) +
      plot_theme(style = "test")

  } else {
    p3 <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(title = "No threshold curve available")
  }

  # PANEL 4: ROC
  p4 <- if (!is.null(roc_plot)) roc_plot else
    ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::labs(title = "ROC unavailable")

  # PANEL 5: module breakdown (hits proxy)
  p5 <- ggplot2::ggplot(diag, ggplot2::aes(core_hits + mig_hits + neuro_hits, score)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::labs(title = "Module signal vs score", x = "Module hits", y = "Score") +
    plot_theme(style = "test")

  # PANEL 6: class composition

  p6 <- plot_gnrh_distribution(object, group.by = "gnrh_class",
                               style = "test",
                               plot.ttl = "Class distribution",
                               cols = gnrh_colors("class"), label = T)


  # PANEL 7: class  enrichment across module hits

  diag$hits_total <- diag$core_hits + diag$mig_hits + diag$neuro_hits

  p7 <- ggplot2::ggplot(diag, ggplot2::aes(x = factor(hits_total), fill = class)) +
    ggplot2::geom_bar(position = "fill", color = "black", linewidth = 0.2) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = gnrh_colors("class")) +
    ggplot2::labs(title = "Enrichment across module hits",
                  x = "Total hits", y = "Fraction") +
    plot_theme(style = "test", leg.pos = "none", x.ang = 90)


  # PANEL 8: Classification rule summary


  rule_df <- object@misc$gnrh$classify_summary

  if (!is.null(rule_df)) {

    rule_df$rule <- factor(rule_df$rule, levels = rev(rule_df$rule))

    p8 <- ggplot2::ggplot(rule_df, ggplot2::aes(rule, n_cells)) +
      ggplot2::geom_col(fill = "#7B61FF") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(
        ggplot2::aes(label = n_cells),
        hjust = -0.1,
        size = 3
      ) +
      ggplot2::labs(
        title = "Classification rule summary",
        x = NULL,
        y = "Number of cells"
      ) +
      plot_theme(style = style)
  }


  # COMBINE
  final_plot <-
    (p1 | p2 | p3 | p4) /
    (p5 | p6 | p7 ) +
    patchwork::plot_annotation(title = "GnRH Diagnostic Report")

  final_plot
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
plot_network <- function(df, top_n = 25, threshold = 0.4) {

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
    edges$from != edges$to & edges$weight > threshold, , drop = FALSE ]

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
  igraph::V(g)$module <- factor(comm$membership)

  set.seed(123)

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      colour = "grey75"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = module, size = score)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      max.overlaps = 30
    ) +
    ggplot2::theme_void()
}



#' Plot GnRH coexpression markers
#'
#' Visualizes genes coexpressed with \code{GNRH1} and ranks them by
#' coexpression strength.
#'
#' @param df A marker table returned by \code{\link{gnrh_markers}}.
#' Must contain the columns \code{gene}, \code{coexpr}, \code{score},
#' and \code{p_val_adj}.
#' @param coexp_cutoff Minimum coexpression value required for plotting.
#' Default is \code{0.25}.
#' @param txtsize Base text size passed to \code{\link{plot_theme}}.
#' Default is \code{12}.
#'
#' @return A \code{ggplot2} object showing coexpression strength,
#' marker score, and adjusted p-value significance.
#'
#' @details
#' \code{GNRH1} itself is removed before plotting. Genes are ordered
#' by decreasing coexpression strength.
#'
#' Point aesthetics:
#' \itemize{
#'   \item x-axis: correlation with \code{GNRH1}
#'   \item y-axis: coexpressed genes
#'   \item point size: composite marker score
#'   \item point color: \eqn{-log10(adjusted~p-value)}
#' }
#'
#' @examples
#' \dontrun{
#' p <- plot_gnrh_coexpr(
#'   markers,
#'   coexp_cutoff = 0.35
#' )
#'
#' print(p)
#' }
#'
#' @export
plot_gnrh_coexpr <- function(
    df,
    coexp_cutoff = 0.25,
    txtsize = 12
) {

  # Remove self-gene + apply cutoff
  df2 <- df[df$gene != "GNRH1" & df$coexpr >= coexp_cutoff, ,drop = FALSE]

  # Handle adjusted p-values
  df2$log_padj <- -log10(pmax(df2$p_val_adj, 1e-300))

  # Order by coexpression
  df2 <- df2[order(df2$coexpr, decreasing = TRUE), ,drop = FALSE ]
  df2$gene <- factor(df2$gene,levels = rev(df2$gene))

  # Plot
  ggplot2::ggplot(df2, ggplot2::aes(x = coexpr,y = gene)) +
    ggplot2::geom_point(
      ggplot2::aes(size = score, color = log_padj),alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::scale_color_viridis_c(option = "viridis", direction = 1) +
    ggplot2::labs(x = "Correlation coefficient", y = NULL, title = "GNRH1 correlation",
                  size = "Score",color = expression(-log[10](adj.~p))) +
    plot_theme(style = "bw", txtsize = txtsize) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
}


#' Plot GnRHcell runtime across datasets
#'
#' Plots runtime across datasets as a line plot. The x-axis displays
#' each dataset together with the number of detected GnRH-positive cells
#' over the total number of cells.
#'
#' @param files Optional character vector of GnRH run-info TSV files.
#' @param dir Directory containing run-info tables.
#' @param pattern Regex pattern used to find run-info files.
#' @param metric Runtime column to plot, e.g. \code{"total_sec"}.
#' @param style Plot theme style passed to \code{plot_theme()}.
#' @param x.ang X-axis label angle.
#' @param show_points Logical; show points on the curve.
#'
#' @return A \code{ggplot2} object.
#' @export
plot_gnrh_runtime_curve <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    metric = "total_sec",
    style = "classic",
    x.ang = 45,
    show_points = TRUE
) {
  stats <- .load_gnrh_stats(files = files, dir = dir, pattern = pattern)

  if (!metric %in% colnames(stats)) {
    stop("Metric not found: ", metric, call. = FALSE)
  }

  stats$pos <- as.numeric(stats$pos)
  stats$n_cells <- as.numeric(stats$n_cells)
  stats[[metric]] <- as.numeric(stats[[metric]])
  stats$pos[is.na(stats$pos)] <- 0

  stats <- stats::aggregate(
    cbind(pos, n_cells, runtime = stats[[metric]]) ~ dataset,
    data = stats,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  stats <- stats[order(stats$n_cells), ]
  stats$label <- factor(
    paste0(stats$dataset, "\n", round(stats$pos), "/", round(stats$n_cells), " GnRH+"),
    levels = paste0(stats$dataset, "\n", round(stats$pos), "/", round(stats$n_cells), " GnRH+")
  )

  p <- ggplot2::ggplot(stats, ggplot2::aes(label, runtime, group = 1)) +
    ggplot2::geom_line(linewidth = 0.8, color = "#819DC7")

  if (show_points) {
    p <- p + ggplot2::geom_point(size = 2, color = gnrh_colors("status")["pos"])
  }

  p +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(runtime, 1), "s")),
      vjust = -0.7,
      size = 3
    ) +
    plot_theme(style = style, x.ang = x.ang) +
    ggplot2::labs(
      title = "GnRHcell runtime across datasets",
      x = "Dataset",
      y = metric
    )
}


#' Plot detected GnRH-positive cells across datasets
#'
#' @param files Optional character vector of run-info TSV files.
#' @param dir Directory containing runtime summary tables.
#' @param pattern Regex pattern used to identify runtime files.
#' @param style Plot style.
#' @param x.ang X-axis text angle.
#' @param debug Print loaded table.
#'
#' @return A ggplot2 object.
#' @export
plot_gnrh_detected <- function(
    files = NULL,
    dir = file.path("results", "tables"),
    pattern = "_gnrh_run_info\\.tsv$",
    style = "classic",
    x.ang = 45,
    debug = TRUE
) {

  stats <- .load_gnrh_stats(
    files = files,
    dir = dir,
    pattern = pattern
  )

  if (debug) {
    print(stats[, c("source_file", "dataset", "n_cells", "neg", "pos")])
  }

  stats <- stats[, c("dataset", "pos")]

  stats <- aggregate(
    pos ~ dataset,
    data = stats,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  stats$dataset <- factor(
    stats$dataset,
    levels = stats$dataset[order(stats$pos)]
  )

  ggplot2::ggplot(stats, ggplot2::aes(x = dataset, y = pos)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = round(pos)),
      vjust = -0.3,
      size = 4
    ) +
    plot_theme(style = style, x.ang = x.ang) +
    ggplot2::labs(
      title = "Detected GnRH cells",
      x = "Dataset",
      y = "GnRH+ cells"
    )
}






