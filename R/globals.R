# CRAN notes for non-standard evaluation

#' @importFrom rlang .data

NULL

utils::globalVariables(
  c(
    "x", "y", "label", "type",
    "gene", "cluster", "avg_log2FC",
    "p_val_adj", "p_val",
    "score", "feature_rank",
    "coexpr", "pct1", "pct2", "spec",
    "from", "to", "weight",
    "name", "community",
    "sec", "threshold",
    "sensitivity", "specificity",
    "F1", "fpr", "tpr", "method_support",
    ".data", "expr", "median", "above_threshold",
    "pct.exp",
    "z",
    "x0", "y0", "x1", "y1",
    "core_hits", "mig_hits", "neuro_hits",
    "hits_total",
    "rule", "n_cells",
    "log_padj",
    "dataset", "pos", "runtime",
    "value", "fill", "lab",
    "module",
    "count"
  )
)

GNRH_SCALE_FACTOR <- 1e4
GNRH_DEFAULT_MIN_UMI <- 2
GNRH_DEFAULT_SCORE_Q <- 0.85
CELL_RASTER_THRESHOLD <- 1e5
CELL_DEFAULT_THEME <- "classic"
CELL_DEFAULT_BASE_SIZE <- 12


# Global package options

.onLoad <- function(libname, pkgname) {

  op <- options()

  defaults <- list(
    gnrhcell.raster_threshold = 1e5,
    gnrhcell.base_size = 12,
    gnrhcell.base_family = "Helvetica",
    gnrhcell.palette = "base"
  )

  toset <- !(names(defaults) %in% names(op))

  if (any(toset)) {
    options(defaults[toset])
  }

  invisible()
}
