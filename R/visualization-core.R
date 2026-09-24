# ========================================================================= #
# Shared visualization engine
#
# Adapted from the author's ~/projects/pkg/cellplot.R.  GnRHcell keeps the
# biological plotting API, while these internal helpers provide one consistent
# implementation for sizing, rasterization, palettes, and plot metadata.
# ========================================================================= #

.cellplot_auto_raster <- function(ncells, threshold = 3e5) {
  if (!is.numeric(ncells) || length(ncells) != 1L || is.na(ncells) ||
      !is.finite(ncells) || ncells < 0) {
    stop("`ncells` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive number.", call. = FALSE)
  }
  ncells >= threshold
}

.cellplot_auto_pt_size <- function(ncells, raster = FALSE) {
  if (!is.numeric(ncells) || length(ncells) != 1L || is.na(ncells) ||
      !is.finite(ncells) || ncells < 1) {
    stop("`ncells` must be a single positive number.", call. = FALSE)
  }
  size <- if (ncells >= 5e5) {
    0.015
  } else if (ncells >= 2e5) {
    0.020
  } else if (ncells >= 1e5) {
    0.025
  } else if (ncells >= 5e4) {
    0.040
  } else if (ncells >= 2e4) {
    0.070
  } else if (ncells >= 1e4) {
    0.120
  } else if (ncells >= 5e3) {
    0.200
  } else if (ncells >= 1e3) {
    0.350
  } else {
    0.600
  }
  if (isTRUE(raster)) size <- size * 1.15
  size
}

.cellplot_point_params <- function(
    object,
    raster = NULL,
    pt.size = NULL,
    raster.threshold = 3e5
) {
  .validate_seurat(object)
  ncells <- ncol(object)
  if (is.null(raster)) raster <- .cellplot_auto_raster(ncells, raster.threshold)
  if (!is.logical(raster) || length(raster) != 1L || is.na(raster)) {
    stop("`raster` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  if (is.null(pt.size)) pt.size <- .cellplot_auto_pt_size(ncells, raster)
  if (!is.numeric(pt.size) || length(pt.size) != 1L ||
      !is.finite(pt.size) || pt.size <= 0) {
    stop("`pt.size` must be a single positive number.", call. = FALSE)
  }
  list(ncells = ncells, raster = isTRUE(raster), pt.size = pt.size)
}

.cellplot_gradient <- function(
    palette = "bwr",
    colors = NULL,
    n = 100L,
    reverse = FALSE
) {
  if (!is.null(colors)) {
    result <- grDevices::colorRampPalette(colors)(n)
  } else {
    key <- tolower(palette)
    custom <- list(
      bwr = c("#2166AC", "#F7F7F7", "#B2182B"),
      rwb = c("#B2182B", "#F7F7F7", "#2166AC"),
      coolwarm = c("#3B4CC0", "#F7F7F7", "#B40426"),
      tealred = c("#018571", "#F5F5F5", "#A6611A"),
      navyred = c("#053061", "#F7F7F7", "#67001F"),
      grayred = c("#F2F2F2", "#FB6A4A", "#A50F15")
    )
    if (key %in% names(custom)) {
      result <- grDevices::colorRampPalette(custom[[key]])(n)
    } else if (requireNamespace("RColorBrewer", quietly = TRUE) &&
               palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      base <- RColorBrewer::brewer.pal(min(max_n, 9L), palette)
      result <- grDevices::colorRampPalette(base)(n)
    } else {
      stop("Unknown palette: ", palette, ".", call. = FALSE)
    }
  }
  if (isTRUE(reverse)) rev(result) else result
}

.cellplot_attach_dimensions <- function(plot, width, height) {
  attr(plot, "gnrh_dimensions") <- c(width = width, height = height)
  attr(plot, "cellplot_dimensions") <- c(width = width, height = height)
  plot
}
