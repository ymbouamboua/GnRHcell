# ============================================================================= #
# GnRHcell internal utilities
# ============================================================================= #


# ============================================================================= #
# Runtime / logging
# ============================================================================= #

#' Format elapsed runtime into human-readable text
#'
#' @param seconds Numeric duration in seconds.
#'
#' @return Character string containing formatted elapsed time.
#'
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {
  if (is.null(seconds) || !length(seconds) || !is.finite(seconds))
    return(NA_character_)

  seconds <- max(as.numeric(seconds), 0)

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


#' Format runtime duration
#'
#' Backward-compatible alias for \code{.format_duration()}.
#'
#' @keywords internal
#' @noRd
.format_runtime <- function(seconds) {
  .format_duration(seconds)
}


#' Styled GnRHcell logger
#'
#' @param verbose Logical; print messages.
#' @param color Logical; use ANSI colors.
#'
#' @return A logging function.
#'
#' @keywords internal
#' @noRd
.msg <- function(
    verbose = TRUE,
    color = interactive()
) {
  t0 <- Sys.time()

  col_reset  <- if (color) "\033[0m" else ""
  col_bold   <- if (color) "\033[1m" else ""
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
    if (!isTRUE(verbose))
      return(invisible(NULL))

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
      prefix_col, prefix, col_reset, " ",
      txt_col, txt, col_reset, "\n",
      sep = ""
    )

    invisible(NULL)
  }
}


# ============================================================================= #
# Numeric / gene helpers
# ============================================================================= #

#' Standardize numeric vector
#'
#' @param x Numeric vector.
#'
#' @return Numeric standardized vector.
#'
#' @keywords internal
#' @noRd
.scale0 <- function(x) {
  x <- as.numeric(x)

  finite <- is.finite(x)

  if (!any(finite))
    return(rep(0, length(x)))

  s <- stats::sd(x[finite])

  if (!is.finite(s) || s == 0)
    return(rep(0, length(x)))

  out <- rep(NA_real_, length(x))
  out[finite] <- (x[finite] - mean(x[finite])) / s
  out
}


#' Match gene symbols case-insensitively
#'
#' @param features Character vector of requested gene symbols.
#' @param genes Character vector of available gene symbols.
#'
#' @return Character vector of matched genes.
#'
#' @keywords internal
#' @noRd
.match_genes <- function(features, genes) {
  if (!length(features) || !length(genes))
    return(character(0))

  idx <- match(
    toupper(as.character(features)),
    toupper(as.character(genes))
  )

  unique(
    genes[idx[!is.na(idx)]]
  )
}


#' Extract expression matrix from a Seurat object
#'
#' @param object A Seurat object.
#' @param assay Assay name.
#' @param layer Expression layer.
#'
#' @return Expression matrix.
#'
#' @keywords internal
#' @noRd
.get_expr <- function(
    object,
    assay = NULL,
    layer = "counts"
) {
  if (!inherits(object, "Seurat"))
    stop("`object` must be a Seurat object.", call. = FALSE)

  assay <- assay %||% Seurat::DefaultAssay(object)

  if (!assay %in% names(object@assays))
    stop("Assay `", assay, "` not found.", call. = FALSE)

  Seurat::GetAssayData(
    object,
    assay = assay,
    layer = layer
  )
}


# ============================================================================= #
# GnRH confidence
# ============================================================================= #

#' Compute high-confidence GnRH classification
#'
#' High-confidence calls are intentionally stricter than
#' \code{gnrh_status}. Direct cells require strong independent GnRH identity.
#' Supported cells require strong identity together with independent
#' transcriptomic support.
#'
#' @param metadata Seurat metadata data frame.
#' @param min_umi Retained for backward compatibility.
#'
#' @return Logical vector.
#'
#' @keywords internal
#' @noRd
.compute_gnrh_confident <- function(
    metadata,
    min_umi = 2
) {
  required <- c(
    "gnrh_class",
    "gnrh_raw",
    "gnrh_identity_strong",
    "gnrh_independent_support"
  )

  missing <- setdiff(required, colnames(metadata))

  if (length(missing)) {
    stop(
      "Missing metadata columns required for confidence classification: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cls <- as.character(metadata$gnrh_class)

  raw <- as.numeric(metadata$gnrh_raw)

  identity_strong <-
    !is.na(metadata$gnrh_identity_strong) &
    metadata$gnrh_identity_strong

  independent_support <-
    !is.na(metadata$gnrh_independent_support) &
    metadata$gnrh_independent_support

  direct_confident <-
    cls == "direct" &
    raw >= min_umi &
    identity_strong

  supported_confident <-
    cls == "supported" &
    identity_strong &
    independent_support

  confident <-
    direct_confident |
    supported_confident

  confident[is.na(confident)] <- FALSE

  confident
}


# ============================================================================= #
# Marker-hit summaries
# ============================================================================= #

#' Add total GnRH marker hit summaries
#'
#' Computes identity, neuroendocrine, migration, and total marker evidence.
#'
#' @param object A Seurat object.
#'
#' @return Updated Seurat object.
#'
#' @keywords internal
#' @noRd
.make_total_hits <- function(object) {
  md <- object[[]]

  required <- c(
    "gnrh_core_hits",
    "gnrh_mig_hits",
    "gnrh_neuro_hits"
  )

  missing <- setdiff(required, colnames(md))

  if (length(missing)) {
    stop(
      "Missing GnRH marker-hit columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  identity_neuro <-
    md$gnrh_core_hits +
    md$gnrh_neuro_hits

  total <-
    identity_neuro +
    md$gnrh_mig_hits

  bin_hits <- function(x) {
    factor(
      cut(
        x,
        breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
        labels = c("0", "1", "2", "3", "4", "5+"),
        right = TRUE
      ),
      levels = c("0", "1", "2", "3", "4", "5+")
    )
  }

  object$gnrh_identity_neuro_hits <- identity_neuro
  object$gnrh_total_hits <- total
  object$gnrh_total_hits_bin <- bin_hits(total)

  # Backward compatibility
  object$total_hits <- total
  object$total_hits_bin <- object$gnrh_total_hits_bin

  object
}



# ============================================================================= #
# ROC validation
# ============================================================================= #

#' Build GnRH ROC curve data
#'
#' Computes ROC coordinates against an independent truth annotation.
#'
#' @param object A Seurat object.
#' @param truth Character name of the metadata truth column, or a vector
#'   containing binary truth values.
#' @param predictor Character name of the predictor metadata column.
#'   Default is \code{"gnrh_support_score_raw"}.
#' @param positive Optional value identifying the positive truth class.
#'
#' @return Data frame containing ROC coordinates, thresholds, and AUC.
#'
#' @keywords internal
#' @noRd
.gnrh_build_roc <- function(
    object,
    truth,
    predictor = "gnrh_support_score_raw",
    positive = NULL
) {
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required.", call. = FALSE)

  md <- object[[]]

  if (missing(truth) || is.null(truth))
    stop(
      "`truth` must contain an independent binary reference annotation.",
      call. = FALSE
    )

  y <- if (is.character(truth) && length(truth) == 1L) {
    if (!truth %in% colnames(md))
      stop("Truth column `", truth, "` not found.", call. = FALSE)

    md[[truth]]
  } else {
    truth
  }

  if (!predictor %in% colnames(md))
    stop(
      "Predictor column `", predictor, "` not found.",
      call. = FALSE
    )

  x <- as.numeric(md[[predictor]])

  if (length(y) != length(x))
    stop("`truth` must contain one value per cell.", call. = FALSE)

  keep <- !is.na(y) & is.finite(x)

  y <- y[keep]
  x <- x[keep]

  classes <- unique(as.character(y))

  if (length(classes) != 2L)
    stop(
      "ROC requires exactly two truth classes; found: ",
      paste(classes, collapse = ", "),
      call. = FALSE
    )

  if (!is.null(positive)) {
    if (!as.character(positive) %in% classes)
      stop("`positive` is not present in `truth`.", call. = FALSE)

    negative <- setdiff(classes, as.character(positive))

    y <- factor(
      as.character(y),
      levels = c(negative, as.character(positive))
    )
  } else {
    y <- factor(y)

    if (nlevels(y) != 2L)
      stop("Unable to determine binary truth levels.", call. = FALSE)
  }

  roc <- pROC::roc(
    response = y,
    predictor = x,
    levels = levels(y),
    direction = "<",
    quiet = TRUE
  )

  data.frame(
    threshold = roc$thresholds,
    fpr = 1 - roc$specificities,
    tpr = roc$sensitivities,
    auc = as.numeric(pROC::auc(roc)),
    predictor = predictor,
    stringsAsFactors = FALSE
  )
}



# ============================================================================= #
# Misc storage
# ============================================================================= #

#' Initialize GnRHcell miscellaneous storage
#'
#' @keywords internal
#' @noRd
.init_gnrh_misc <- function(object) {
  if (is.null(object@misc$gnrh) ||
      !is.list(object@misc$gnrh)) {
    object@misc$gnrh <- list()
  }

  object
}


#' Count values
#'
#' @keywords internal
#' @noRd
.count_factor <- function(x) {
  if (is.null(x))
    return(NULL)

  x <- as.character(x)
  x <- x[!is.na(x)]

  as.list(table(x))
}


#' Add GnRHcell run information
#'
#' @param object A Seurat object processed by GnRHcell.
#' @param step_times Named list of runtime values in seconds.
#' @param params Named list of analysis parameters.
#'
#' @return Updated Seurat object.
#'
#' @keywords internal
#' @noRd
.add_run_info <- function(
    object,
    step_times = list(),
    params = list()
) {
  object <- .init_gnrh_misc(object)
  md <- object[[]]

  count_if <- function(column) {
    if (column %in% colnames(md))
      .count_factor(md[[column]])
    else
      NULL
  }

  sum_if <- function(column) {
    if (column %in% colnames(md))
      sum(md[[column]] %in% TRUE, na.rm = TRUE)
    else
      NULL
  }

  timing_human <- lapply(
    step_times,
    .format_duration
  )

  object@misc$gnrh$run_info <- list(
    date = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z"
    ),

    dataset = list(
      n_genes = nrow(object),
      n_cells = ncol(object)
    ),

    summary = list(
      status = count_if("gnrh_status"),
      class = count_if("gnrh_class"),
      confident = count_if("gnrh_confident"),
      stage = count_if("gnrh_stage"),
      secretory = count_if("gnrh_secretory"),

      reference_positive =
        sum_if("gnrh_reference_positive"),

      transcriptomic_candidate =
        sum_if("gnrh_transcriptomic_candidate"),

      direct_isolated =
        sum_if("gnrh_direct_isolated"),

      direct_signal =
        sum_if("gnrh_direct_signal")
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



# ============================================================================= #
# Run statistics
# ============================================================================= #

#' Load GnRHcell run statistics
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

  if (!length(files))
    stop("No GnRH run info files found.", call. = FALSE)

  missing_files <- files[!file.exists(files)]

  if (length(missing_files))
    stop(
      "Missing files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )

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
    "direct",
    "supported",
    "confident",
    "reference_positive",
    "transcriptomic_candidate",
    "direct_isolated",
    "direct_signal",
    "detect_sec",
    "stage_sec",
    "diagnostics_sec",
    "total_sec"
  )

  for (cc in intersect(numeric_cols, colnames(stats))) {
    stats[[cc]] <- suppressWarnings(
      as.numeric(stats[[cc]])
    )
  }

  zero_cols <- intersect(
    c(
      "neg",
      "pos",
      "direct",
      "supported",
      "confident",
      "reference_positive",
      "transcriptomic_candidate",
      "direct_isolated",
      "direct_signal"
    ),
    colnames(stats)
  )

  for (cc in zero_cols)
    stats[[cc]][is.na(stats[[cc]])] <- 0

  stats
}



#' Extract GnRHcell run information
#'
#' @param object A Seurat object processed with \code{\link{run_gnrh}}.
#' @param dataset_name Optional dataset name.
#'
#' @return One-row data frame containing dataset and GnRHcell run statistics.
#'
#' @export
extract_gnrh_run_info <- function(
    object,
    dataset_name = NULL
) {
  if (is.null(object@misc$gnrh$run_info))
    stop(
      "Missing object@misc$gnrh$run_info. Run run_gnrh() first.",
      call. = FALSE
    )

  info <- object@misc$gnrh$run_info
  sec <- info$timing$seconds %||% list()
  hum <- info$timing$human %||% list()

  md <- object[[]]

  count_value <- function(column, value) {
    if (!column %in% colnames(md))
      return(0L)

    sum(
      as.character(md[[column]]) == value,
      na.rm = TRUE
    )
  }

  count_true <- function(column) {
    if (!column %in% colnames(md))
      return(0L)

    sum(
      md[[column]] %in% TRUE,
      na.rm = TRUE
    )
  }

  data.frame(
    dataset = dataset_name %||% "dataset",

    n_genes = as.integer(info$dataset$n_genes),
    n_cells = as.integer(info$dataset$n_cells),

    neg = count_value("gnrh_status", "neg"),
    pos = count_value("gnrh_status", "pos"),

    direct = count_value("gnrh_class", "direct"),
    supported = count_value("gnrh_class", "supported"),

    confident = count_true("gnrh_confident"),

    reference_positive =
      count_true("gnrh_reference_positive"),

    transcriptomic_candidate =
      count_true("gnrh_transcriptomic_candidate"),

    direct_signal =
      count_true("gnrh_direct_signal"),

    direct_isolated =
      count_true("gnrh_direct_isolated"),

    detect_sec =
      as.numeric(sec$detect_sec %||% NA_real_),

    stage_sec =
      as.numeric(sec$stage_sec %||% NA_real_),

    diagnostics_sec =
      as.numeric(sec$diagnostics_sec %||% NA_real_),

    total_sec =
      as.numeric(sec$total_sec %||% NA_real_),

    detect_time =
      hum$detect_sec %||%
      .format_duration(sec$detect_sec %||% NA_real_),

    stage_time =
      hum$stage_sec %||%
      .format_duration(sec$stage_sec %||% NA_real_),

    diagnostics_time =
      hum$diagnostics_sec %||%
      .format_duration(sec$diagnostics_sec %||% NA_real_),

    total_time =
      hum$total_sec %||%
      .format_duration(sec$total_sec %||% NA_real_),

    stringsAsFactors = FALSE
  )
}





#' Build gene sets from marker tables
#' @param files Character vector of marker-table filenames or paths.
#' @param dir Optional directory containing the marker files.
#' @param gene_col Column containing gene identifiers in each marker table.
#'
#' @export

# build_gene_sets <- function(
#     files,
#     dir,
#     gene_col = "gene"
# ) {
#   if (
#     !is.character(files) ||
#     !length(files) ||
#     is.null(names(files)) ||
#     any(names(files) == "") ||
#     anyDuplicated(names(files))
#   ) {
#     stop(
#       "`files` must be a named character vector with unique names.",
#       call. = FALSE
#     )
#   }
#
#   if (
#     length(dir) != 1L ||
#     is.na(dir) ||
#     !dir.exists(dir)
#   ) {
#     stop(
#       "`dir` must be an existing directory.",
#       call. = FALSE
#     )
#   }
#
#   if (
#     length(gene_col) != 1L ||
#     is.na(gene_col) ||
#     !nzchar(gene_col)
#   ) {
#     stop(
#       "`gene_col` must be a non-empty column name.",
#       call. = FALSE
#     )
#   }
#
#   stats::setNames(
#     lapply(
#       names(files),
#       function(dataset) {
#         file <- file.path(
#           dir,
#           files[[dataset]]
#         )
#
#         if (!file.exists(file)) {
#           warning(
#             "Missing file: ",
#             file,
#             call. = FALSE
#           )
#
#           return(character(0))
#         }
#
#         tab <- utils::read.delim(
#           file,
#           sep = "\t",
#           stringsAsFactors = FALSE,
#           check.names = FALSE
#         )
#
#         if (!gene_col %in% colnames(tab)) {
#           warning(
#             "Column `",
#             gene_col,
#             "` not found in ",
#             dataset,
#             ".",
#             call. = FALSE
#           )
#
#           return(character(0))
#         }
#
#         genes <- toupper(
#           trimws(
#             as.character(
#               tab[[gene_col]]
#             )
#           )
#         )
#
#         genes <- genes[
#           !is.na(genes) &
#             nzchar(genes)
#         ]
#
#         sort(unique(genes))
#       }
#     ),
#     names(files)
#   )
# }

# ========================================================================= #
# Build marker gene sets
# ========================================================================= #
build_gene_sets <- function(
    files,
    dir = ".",
    gene_col = "gene"
) {
  if (!is.character(files) || !length(files)) stop("`files` must be a non-empty character vector.", call. = FALSE)
  if (is.null(names(files)) || anyNA(names(files)) || any(!nzchar(names(files))) || anyDuplicated(names(files))) stop("`files` must have unique non-empty names.", call. = FALSE)
  out <- lapply(files, function(f) {
    path <- if (file.exists(f)) f else file.path(dir, f)
    if (!file.exists(path)) stop("Marker file not found: ", path, call. = FALSE)
    x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!gene_col %in% colnames(x)) stop("Column `", gene_col, "` not found in ", basename(path), ".", call. = FALSE)
    genes <- toupper(trimws(as.character(x[[gene_col]])))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    unique(genes)
  })
  names(out) <- names(files)
  out
}


# ========================================================================= #
# Gene-set overlap
# ========================================================================= #
#' Gene Set Overlap Analysis and Visualization
#'
#' @param gene_sets Named list of gene sets to compare.
#' @param min_size Minimum intersection size retained for visualization.
#' @param venn_title Plot title.
#' @param outdir Optional output directory.
#' @param save_plot Logical. Save PDF and PNG.
#' @param max_intersections Maximum intersections displayed. If \code{NULL},
#'   automatically determined from the number of gene sets.
#' @param plot_width Optional output width.
#' @param plot_height Optional output height.
#' @param dpi Raster resolution.
#'
#' @return Named overlap-analysis list.
#' @export
gnrh_gene_upset <- function(
    gene_sets,
    min_size=1,
    venn_title="Overlap of Gene Sets",
    outdir=".",
    save_plot=TRUE,
    max_intersections=Inf,
    plot_width=NULL,
    plot_height=NULL,
    dpi=600
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  for (pkg in c("ggplot2","ComplexUpset","patchwork")) {
    if (!requireNamespace(pkg,quietly=TRUE)) stop("Package '",pkg,"' is required.",call.=FALSE)
  }
  if (!is.list(gene_sets) || length(gene_sets)<2L) stop("`gene_sets` must be a named list containing at least two sets.",call.=FALSE)
  labels <- names(gene_sets)
  if (is.null(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) stop("`gene_sets` must have unique, non-empty names.",call.=FALSE)
  min_size <- as.integer(min_size)
  if (!is.finite(min_size) || min_size<1L) stop("`min_size` must be >= 1.",call.=FALSE)
  if (!is.null(max_intersections) && (length(max_intersections)!=1L || is.na(max_intersections) || max_intersections<1)) stop("`max_intersections` must be NULL, positive, or Inf.",call.=FALSE)
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  # ========================================================================= #
  # Clean sets
  # ========================================================================= #
  clean_genes <- function(x) {
    x <- toupper(trimws(as.character(x)))
    sort(unique(x[!is.na(x) & nzchar(x)]))
  }
  gene_sets <- lapply(gene_sets,clean_genes)
  if (!any(lengths(gene_sets)>0L)) stop("All gene sets are empty.",call.=FALSE)
  n_sets <- length(gene_sets)
  # ========================================================================= #
  # Set sizes
  # ========================================================================= #
  summary_df <- data.frame(
    dataset=labels,
    n_genes=lengths(gene_sets),
    stringsAsFactors=FALSE
  )
  utils::write.csv(summary_df,file.path(outdir,"gene_set_sizes.csv"),row.names=FALSE)
  # ========================================================================= #
  # Pairwise overlaps
  # ========================================================================= #
  combs <- utils::combn(labels,2,simplify=FALSE)
  pairwise_results <- list()
  pairwise_summary <- do.call(rbind,lapply(combs,function(cmb) {
    overlap <- intersect(gene_sets[[cmb[1]]],gene_sets[[cmb[2]]])
    pair_name <- paste(cmb,collapse="_vs_")
    pairwise_results[[pair_name]] <<- overlap
    safe_name <- gsub("[^A-Za-z0-9_\\-]+","_",pair_name)
    utils::write.csv(
      data.frame(gene=overlap),
      file.path(outdir,paste0("genes_",safe_name,".csv")),
      row.names=FALSE
    )
    data.frame(
      dataset1=cmb[1],
      dataset2=cmb[2],
      n_overlap=length(overlap),
      stringsAsFactors=FALSE
    )
  }))
  utils::write.csv(pairwise_summary,file.path(outdir,"pairwise_overlap_summary.csv"),row.names=FALSE)
  # ========================================================================= #
  # Common genes
  # ========================================================================= #
  common_all <- Reduce(intersect,gene_sets)
  utils::write.csv(
    data.frame(gene=common_all),
    file.path(outdir,"genes_common_all.csv"),
    row.names=FALSE
  )
  # ========================================================================= #
  # Unique genes
  # ========================================================================= #
  unique_results <- stats::setNames(vector("list",length(labels)),labels)
  unique_summary <- do.call(rbind,lapply(labels,function(lbl) {
    others <- gene_sets[names(gene_sets)!=lbl]
    other_genes <- if (length(others)) Reduce(union,others) else character(0)
    unique_genes <- setdiff(gene_sets[[lbl]],other_genes)
    unique_results[[lbl]] <<- unique_genes
    safe_lbl <- gsub("[^A-Za-z0-9_\\-]+","_",lbl)
    utils::write.csv(
      data.frame(gene=unique_genes),
      file.path(outdir,paste0("genes_unique_",safe_lbl,".csv")),
      row.names=FALSE
    )
    data.frame(dataset=lbl,n_unique=length(unique_genes),stringsAsFactors=FALSE)
  }))
  utils::write.csv(unique_summary,file.path(outdir,"unique_gene_summary.csv"),row.names=FALSE)
  # ========================================================================= #
  # Membership
  # ========================================================================= #
  all_genes <- sort(unique(unlist(gene_sets,use.names=FALSE)))
  membership <- data.frame(gene=all_genes,stringsAsFactors=FALSE,check.names=FALSE)
  for (lbl in labels) membership[[lbl]] <- all_genes %in% gene_sets[[lbl]]
  utils::write.csv(membership,file.path(outdir,"gene_set_membership.csv"),row.names=FALSE)
  signatures <- apply(
    membership[,labels,drop=FALSE],
    1L,
    function(x) paste(as.integer(x),collapse="")
  )
  intersection_sizes <- sort(table(signatures),decreasing=TRUE)
  intersection_sizes <- intersection_sizes[intersection_sizes>=min_size]
  n_available <- length(intersection_sizes)
  if (!n_available) stop("No intersections satisfy `min_size = ",min_size,"`.",call.=FALSE)
  # ========================================================================= #
  # Adaptive intersections
  # ========================================================================= #
  if (is.null(max_intersections)) {
    max_intersections <- if (n_sets<=4L) {
      25L
    } else if (n_sets<=6L) {
      35L
    } else if (n_sets<=10L) {
      50L
    } else {
      60L
    }
  }
  n_displayed <- if (is.infinite(max_intersections)) {
    n_available
  } else {
    min(n_available,as.integer(max_intersections))
  }
  n_intersections_arg <- n_displayed
  # ========================================================================= #
  # Adaptive layout
  # ========================================================================= #
  if (is.null(plot_width)) {
    plot_width <- max(9,min(18,7+0.16*n_displayed+0.25*n_sets))
  }
  if (is.null(plot_height)) {
    plot_height <- max(6,min(11,4.6+0.45*n_sets))
  }
  longest_label <- max(nchar(labels),1L)
  set_width_ratio <- min(0.32,max(0.19,0.17+longest_label/250))
  count_text_size <- if (n_displayed<=20L) {
    3.6
  } else if (n_displayed<=40L) {
    3.0
  } else {
    2.5
  }
  set_text_size <- if (n_sets<=6L) 9 else if (n_sets<=10L) 8 else 7
  # ========================================================================= #
  # UpSet
  # ========================================================================= #
  plot <- ComplexUpset::upset(
    data=membership,
    intersect=labels,
    min_size=min_size,
    n_intersections=n_intersections_arg,
    width_ratio=set_width_ratio,
    base_annotations=list(
      "Intersection size"=ComplexUpset::intersection_size(
        counts=TRUE,
        bar_number_threshold=0.82,
        text=list(
          size=count_text_size,
          fontface="bold"
        )
      )
    ),
    set_sizes=ComplexUpset::upset_set_size(
      geom=ggplot2::geom_bar(
        width=0.68,
        fill="#4D4D4D"
      )
    ),
    sort_sets="descending",
    sort_intersections_by="cardinality"
  ) +
    patchwork::plot_annotation(
      title=venn_title,
      theme=ggplot2::theme(
        plot.title=ggplot2::element_text(
          hjust=0.5,
          face="bold",
          size=max(12,min(17,19-0.4*n_sets)),
          margin=ggplot2::margin(b=8)
        )
      )
    ) &
    ggplot2::theme(
      axis.text=ggplot2::element_text(size=set_text_size)
    )
  # ========================================================================= #
  # Save
  # ========================================================================= #
  if (isTRUE(save_plot)) {
    for (ext in c("pdf","png")) {
      ggplot2::ggsave(
        filename=file.path(outdir,paste0("gene_overlap_plot.",ext)),
        plot=plot,
        width=plot_width,
        height=plot_height,
        dpi=dpi,
        bg="white"
      )
    }
  }
  # ========================================================================= #
  # Return
  # ========================================================================= #
  list(
    gene_sets=gene_sets,
    summary=summary_df,
    pairwise=pairwise_results,
    pairwise_summary=pairwise_summary,
    common_all=common_all,
    unique=unique_results,
    unique_summary=unique_summary,
    membership=membership,
    intersections_available=n_available,
    intersections_displayed=n_displayed,
    plot_width=plot_width,
    plot_height=plot_height,
    plot=plot,
    upset_plot=plot,
    plot_type="upset"
  )
}


#' GnRH module weights
#'
#' Internal weights used to combine curated GnRH marker modules.
#'
#' @return Named numeric vector of module weights.
#'
#' @keywords internal
#' @noRd
.gnrh_module_weights <- function() {
  c(
    identity_primary          = 3.00,
    identity_supportive       = 1.50,
    migration_primary         = 1.00,
    migration_supportive      = 0.50,
    neuroendocrine_primary    = 1.50,
    neuroendocrine_supportive = 0.75,
    hormone_supportive        = 0.25,
    guidance_environment      = 0.00
  )
}


#' Compute a weighted GnRH module score
#'
#' @param scores Named list or data frame containing module scores.
#' @param weights Named numeric vector of module weights.
#'
#' @return Numeric vector containing the weighted GnRH score.
#'
#' @keywords internal
#' @noRd
.weight_gnrh_scores <- function(
    scores,
    weights = .gnrh_module_weights(),
    normalize = FALSE
) {
  if (!is.list(scores) && !is.data.frame(scores)) {
    stop(
      "`scores` must be a named list or data frame.",
      call. = FALSE
    )
  }

  if (is.null(names(scores))) {
    stop(
      "`scores` must contain named module scores.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(weights) ||
    is.null(names(weights))
  ) {
    stop(
      "`weights` must be a named numeric vector.",
      call. = FALSE
    )
  }

  common <- intersect(
    names(weights),
    names(scores)
  )

  if (!length(common)) {
    stop(
      "No GnRH module scores matched the supplied weights.",
      call. = FALSE
    )
  }

  w <- weights[common]

  keep <- is.finite(w) & w != 0
  common <- common[keep]
  w <- w[keep]

  if (!length(common)) {
    stop(
      "All matched GnRH module weights are zero or non-finite.",
      call. = FALSE
    )
  }

  score_matrix <- do.call(
    cbind,
    lapply(scores[common], as.numeric)
  )

  if (is.null(dim(score_matrix))) {
    score_matrix <- matrix(
      score_matrix,
      ncol = 1L
    )
  }

  weighted <- sweep(
    score_matrix,
    MARGIN = 2L,
    STATS = w,
    FUN = "*"
  )

  score <- rowSums(
    weighted,
    na.rm = TRUE
  )

  if (isTRUE(normalize)) {
    score <- score / sum(abs(w))
  }

  score
}


