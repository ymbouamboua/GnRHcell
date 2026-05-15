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
    ".data", "expr", "median", "above_threshold"
  )
)

GNRH_SCALE_FACTOR <- 1e4
GNRH_DEFAULT_MIN_UMI <- 2
GNRH_DEFAULT_SCORE_Q <- 0.85
CELL_RASTER_THRESHOLD <- 1e5
CELL_DEFAULT_THEME <- "classic"
CELL_DEFAULT_BASE_SIZE <- 12
