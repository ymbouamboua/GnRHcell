#' Format elapsed runtime into human-readable text
#'
#' Converts a duration in seconds into a compact human-readable string
#' suitable for console logging.
#'
#' Formatting rules:
#' \itemize{
#'   \item seconds < 60 → \code{"7.6s"}
#'   \item minutes < 60 → \code{"1m 0.6s"}
#'   \item hours ≥ 1 → \code{"1h 3m 0.4s"}
#' }
#'
#' @param seconds Numeric duration in seconds.
#'
#' @return A character string containing formatted elapsed time.
#'
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {

  if (is.null(seconds) || is.na(seconds)) {
    return(NA_character_)
  }

  seconds <- as.numeric(seconds)

  hrs <- floor(seconds / 3600)
  mins <- floor((seconds %% 3600) / 60)
  secs <- seconds %% 60

  if (hrs > 0) {
    sprintf("%dh %dm %.1fs", hrs, mins, secs)
  } else if (mins > 0) {
    sprintf("%dm %.1fs", mins, secs)
  } else {
    sprintf("%.1fs", secs)
  }
}



#' Styled GnRHcell logger
#'
#' CRAN-safe console logger using ASCII labels and optional ANSI colors.
#'
#' @param verbose Logical; print messages.
#' @param color Logical; use ANSI colors.
#'
#' @return A logging function.
#' @keywords internal
#' @noRd
.msg <- function(
    verbose = TRUE,
    color = interactive()
) {

  t0 <- Sys.time()

  col_reset  <- if (color) "\033[0m" else ""
  col_bold   <- if (color) "\033[1m" else ""

  # devtools/cli-like palette
  col_blue   <- if (color) "\033[34m" else ""
  col_green  <- if (color) "\033[32m" else ""
  col_yellow <- if (color) "\033[33m" else ""
  col_red    <- if (color) "\033[31m" else ""
  col_gray   <- if (color) "\033[90m" else ""
  col_cyan   <- if (color) "\033[36m" else ""
  col_purple <- if (color) "\033[35m" else ""

  function(
    ...,
    type = c("info", "step", "done", "warn", "error", "header"),
    duration = NULL
  ) {

    if (!isTRUE(verbose)) {
      return(invisible(NULL))
    }

    type <- match.arg(type)
    txt <- paste(..., collapse = " ")

    if (type == "done") {
      if (is.null(duration)) {
        duration <- as.numeric(
          difftime(Sys.time(), t0, units = "secs")
        )
      }

      txt <- sprintf(
        "%s Duration: %s",
        txt,
        .format_duration(duration)
      )
    }

    prefix <- switch(
      type,
      info   = "[INFO]",
      step   = "[STEP]",
      done   = "[DONE]",
      warn   = "[WARN]",
      error  = "[ERROR]",
      header = "[GNRH]"
    )

    prefix_col <- switch(
      type,
      info   = col_gray,
      step   = col_blue,
      done   = col_green,
      warn   = col_yellow,
      error  = col_red,
      header = paste0(col_purple, col_bold)
    )

    txt_col <- switch(
      type,
      info   = col_gray,
      step   = col_cyan,
      done   = col_green,
      warn   = col_yellow,
      error  = col_red,
      header = paste0(col_purple, col_bold)
    )

    cat(
      prefix_col,
      prefix,
      col_reset,
      " ",
      txt_col,
      txt,
      col_reset,
      "\n",
      sep = ""
    )

    invisible(NULL)
  }
}


#' Standardize numeric vector
#'
#' Internal helper that centers and scales a numeric vector.
#'
#' If the vector has zero or undefined standard deviation, a vector
#' of zeros is returned.
#'
#' @param x Numeric vector.
#'
#' @return Numeric standardized vector.
#'
#' @keywords internal
#' @noRd
.scale0 <- function(x) {

  s <- stats::sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) return(rep(0, length(x)))

  as.numeric(scale(x))
}

#' Match gene symbols case-insensitively
#'
#' Internal helper for matching requested features to available genes
#' regardless of letter case.
#'
#' @param features Character vector of requested gene symbols.
#' @param genes Character vector of available gene symbols.
#'
#' @return Character vector of matched genes using the original names
#' from \code{genes}.
#'
#' @keywords internal
#' @noRd
.match_genes <- function(features, genes) {

  idx <- match(toupper(features), toupper(genes))
  genes[stats::na.omit(idx)]
}


#' Extract expression matrix from a Seurat object
#'
#' Internal helper for retrieving assay data from a Seurat object.
#'
#' @param object A Seurat object.
#' @param assay Assay name. If \code{NULL}, the default assay is used.
#' @param layer Expression layer to extract. Default is \code{"counts"}.
#'
#' @return Expression matrix from the requested assay and layer.
#'
#' @keywords internal
#' @noRd
.get_expr <- function(object, assay = NULL, layer = "counts") {

  assay <- assay %||% Seurat::DefaultAssay(object)

  Seurat::GetAssayData(
    object,
    assay = assay,
    layer = layer
  )
}


#' Compute high-confidence GnRH truth labels
#'
#' Internal helper that derives high-confidence GnRH truth labels
#' from marker hits and raw GnRH expression.
#'
#' A cell is labeled \code{"pos"} when it has at least two core
#' marker hits, at least one migration or neuroendocrine marker hit,
#' and raw GnRH expression greater than or equal to \code{min_umi}.
#'
#' @param md Metadata data frame containing GnRH detection columns.
#' @param min_umi Minimum raw GnRH UMI count required.
#'
#' @return Factor vector with levels \code{"neg"} and \code{"pos"}.
#'
#' @keywords internal
#' @noRd
.compute_gnrh_truth <- function(md, min_umi) {

  required <- c(
    "gnrh_core_hits",
    "gnrh_mig_hits",
    "gnrh_neuro_hits",
    "gnrh_raw"
  )

  missing <- setdiff(required, colnames(md))

  if (length(missing) > 0) {
    stop(
      "Missing required columns for gnrh_truth: ",
      paste(missing, collapse = ", ")
    )
  }

  truth <- (
    md$gnrh_core_hits >= 2 &
      (md$gnrh_mig_hits >= 1 | md$gnrh_neuro_hits >= 1) &
      md$gnrh_raw >= min_umi
  )

  factor(
    ifelse(truth, "pos", "neg"),
    levels = c("neg", "pos")
  )
}


#' Add total GnRH marker hit summaries
#'
#' Internal helper that computes total marker hits and binned hit
#' categories from GnRH core and migration marker counts.
#'
#' @param object A Seurat object containing GnRH marker hit metadata.
#'
#' @return Seurat object updated with \code{total_hits} and
#' \code{total_hits_bin} metadata columns.
#'
#' @keywords internal
#' @noRd
.make_total_hits <- function(object) {

  md <- object@meta.data

  total <- md$gnrh_core_hits + md$gnrh_mig_hits

  bins <- cut(
    total,
    breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
    labels = c("0", "1", "2", "3", "4", "5+"),
    right = TRUE
  )

  object$total_hits <- total
  object$total_hits_bin <- factor(
    bins,
    levels = c("0", "1", "2", "3", "4", "5+")
  )

  object
}


#' Build GnRH ROC curve data
#'
#' Internal helper that computes ROC curve coordinates from GnRH
#' expression and GnRH classification labels.
#'
#' @param object A Seurat object containing \code{gnrh_expr} and
#' \code{gnrh_class} metadata columns.
#'
#' @return Data frame with false-positive rate (\code{fpr}) and
#' true-positive rate (\code{tpr}).
#'
#' @keywords internal
#' @noRd
.gnrh_build_roc <- function(object) {

  expr <- object$gnrh_expr

  #truth <- object$gnrh_class %in% c("high", "low")
  truth <- object$gnrh_class %in% c("pos")

  truth <- factor(truth, levels = c(FALSE, TRUE))

  if (length(unique(truth)) != 2) {
    stop("ROC requires 2 classes: found ",
         paste(unique(truth), collapse = ", "))
  }

  roc <- pROC::roc(
    response = truth,
    predictor = expr,
    quiet = TRUE
  )

  data.frame(
    fpr = 1 - roc$specificities,
    tpr = roc$sensitivities
  )
}


#' Initialize GnRHcell miscellaneous storage
#'
#' Internal helper that ensures \code{object@misc$gnrh} exists.
#'
#' @param object A Seurat object.
#'
#' @return Seurat object with initialized \code{object@misc$gnrh}.
#'
#' @keywords internal
#' @noRd
.init_gnrh_misc <- function(object) {
  if (is.null(object@misc$gnrh)) {
    object@misc$gnrh <- list()
  }
  object
}


#' Count factor levels
#'
#' Internal helper that counts values and returns them as a list.
#'
#' @param x Vector or factor to count.
#'
#' @return Named list containing counts per level or value.
#'
#' @keywords internal
#' @noRd
.count_factor <- function(x) {
  x <- as.character(x)
  as.list(table(x))
}


#' Format runtime duration
#'
#' Internal helper that formats elapsed time in seconds as a
#' human-readable string.
#'
#' @param seconds Numeric runtime in seconds.
#'
#' @return Character string formatted as seconds, minutes and seconds,
#' or hours, minutes and seconds.
#'
#' @keywords internal
#' @noRd
.format_runtime <- function(seconds) {

  if (is.null(seconds) || is.na(seconds)) {
    return(NA_character_)
  }

  seconds <- round(seconds)

  hrs <- seconds %/% 3600
  mins <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60

  if (hrs > 0) {
    sprintf("%dh %dm %ds", hrs, mins, secs)
  } else if (mins > 0) {
    sprintf("%dm %ds", mins, secs)
  } else {
    sprintf("%ds", secs)
  }
}

#' Add GnRHcell run information
#'
#' Internal helper that stores runtime, dataset size, classification
#' summaries, parameters, and session information in
#' \code{object@misc$gnrh$run_info}.
#'
#' @param object A Seurat object processed by GnRHcell.
#' @param step_times Named list of runtime values in seconds.
#' @param params Named list of analysis parameters.
#'
#' @return Seurat object with run information stored in
#' \code{object@misc$gnrh$run_info}.
#'
#' @keywords internal
#' @noRd
.add_run_info <- function(object, step_times = list(), params = list()) {

  object <- .init_gnrh_misc(object)
  md <- object[[]]

  timing_human <- lapply(step_times, .format_runtime)

  object@misc$gnrh$run_info <- list(
    date = as.character(Sys.time()),

    dataset = list(
      n_genes = nrow(object),
      n_cells = ncol(object)
    ),

    summary = list(
      status = if ("gnrh_status" %in% colnames(md)) .count_factor(md$gnrh_status) else NULL,
      truth  = if ("gnrh_truth" %in% colnames(md)) .count_factor(md$gnrh_truth) else NULL,
      stage  = if ("gnrh_stage" %in% colnames(md)) .count_factor(md$gnrh_stage) else NULL
    ),

    timing = list(
      seconds = step_times,
      human = timing_human
    ),

    params = params,

    session = list(
      r_version = R.version.string,
      platform = R.version$platform
    )
  )

  object
}


#' Load GnRHcell run statistics
#'
#' Internal helper that loads one or more GnRHcell run information
#' tables from tab-separated files.
#'
#' @param files Optional character vector of file paths. If \code{NULL},
#' files are discovered from \code{dir} using \code{pattern}.
#' @param dir Directory to search when \code{files = NULL}.
#' Default is \code{"."}.
#' @param pattern File name pattern used to identify run information
#' tables. Default is \code{"_gnrh_run_info\\\\.tsv$"}.
#'
#' @return Data frame containing combined GnRHcell run statistics.
#'
#' @keywords internal
#' @noRd
.load_gnrh_stats <- function(
    files = NULL,
    dir = ".",
    pattern = "_gnrh_run_info\\.tsv$"
) {

  if (is.null(files)) {
    files <- list.files(
      path = dir,
      pattern = pattern,
      full.names = TRUE
    )
  }

  if (length(files) == 0) {
    stop("No GnRH run info files found.", call. = FALSE)
  }

  missing_files <- files[!file.exists(files)]

  if (length(missing_files) > 0) {
    stop(
      "Missing files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  stats <- do.call(
    rbind,
    lapply(files, function(f) {

      x <- utils::read.delim(
        f,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      x$source_file <- basename(f)

      x
    })
  )

  numeric_cols <- c(
    "n_genes",
    "n_cells",
    "neg",
    "pos",
    "detect_sec",
    "stage_sec",
    "diagnostics_sec",
    "total_sec"
  )

  for (cc in intersect(numeric_cols, colnames(stats))) {
    stats[[cc]] <- suppressWarnings(as.numeric(stats[[cc]]))
  }

  if ("pos" %in% colnames(stats)) {
    stats$pos[is.na(stats$pos)] <- 0
  }

  if ("neg" %in% colnames(stats)) {
    stats$neg[is.na(stats$neg)] <- 0
  }

  stats
}



#' Extract GnRHcell run information
#'
#' Extracts runtime, dataset size, and detection summary statistics
#' from a Seurat object processed with \code{\link{run_gnrh}}.
#'
#' @param object A Seurat object containing
#' \code{object@misc$gnrh$run_info}.
#' @param dataset_name Optional dataset name. If \code{NULL},
#' \code{"dataset"} is used.
#'
#' @return A data frame containing dataset size, GnRH detection counts,
#' runtime values in seconds, and formatted runtime strings.
#'
#' @seealso
#' \code{\link{run_gnrh}},
#' \code{\link{detect_gnrh}},
#' \code{\link{stage_gnrh}}
#'
#' @export
extract_gnrh_run_info <- function(object,
                                  dataset_name = NULL) {

  `%||%` <- function(x, y) if (is.null(x)) y else x

  if (is.null(object@misc$gnrh$run_info)) {
    stop("Missing object@misc$gnrh$run_info. Run run_gnrh() first.", call. = FALSE)
  }

  info <- object@misc$gnrh$run_info
  sec <- info$timing$seconds %||% list()
  hum <- info$timing$human %||% list()

  md <- object[[]]

  pos_n <- if ("gnrh_status" %in% colnames(md)) {
    sum(as.character(md$gnrh_status) == "pos", na.rm = TRUE)
  } else {
    info$summary$status$pos %||% 0L
  }

  neg_n <- if ("gnrh_status" %in% colnames(md)) {
    sum(as.character(md$gnrh_status) == "neg", na.rm = TRUE)
  } else {
    info$summary$status$neg %||% 0L
  }

  data.frame(
    dataset = dataset_name %||% "dataset",
    n_genes = as.integer(info$dataset$n_genes),
    n_cells = as.integer(info$dataset$n_cells),

    neg = as.integer(neg_n),
    pos = as.integer(pos_n),

    detect_sec = as.numeric(sec$detect_sec %||% NA_real_),
    stage_sec = as.numeric(sec$stage_sec %||% NA_real_),
    diagnostics_sec = as.numeric(sec$diagnostics_sec %||% NA_real_),
    total_sec = as.numeric(sec$total_sec %||% NA_real_),

    detect_time = hum$detect_sec %||% NA_character_,
    stage_time = hum$stage_sec %||% NA_character_,
    diagnostics_time = hum$diagnostics_sec %||% NA_character_,
    total_time = hum$total_sec %||% NA_character_,

    stringsAsFactors = FALSE
  )
}

