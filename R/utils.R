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


#' Compute high-confidence GnRH classification
#'
#' High-confidence GnRH cells must contain detected GNRH1.
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
    "gnrh_identity_moderate"
  )

  missing <- setdiff(
    required,
    colnames(metadata)
  )

  if (length(missing) > 0L) {
    stop(
      "Missing metadata columns required for confidence classification: ",
      paste(
        missing,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  gnrh_class <- as.character(
    metadata$gnrh_class
  )

  gnrh_raw <- metadata$gnrh_raw

  identity_moderate <-
    !is.na(metadata$gnrh_identity_moderate) &
    metadata$gnrh_identity_moderate

  # --------------------------------------------------------------------------- #
  # Direct candidates
  # --------------------------------------------------------------------------- #

  direct_confident <-
    gnrh_class == "direct" &
    (
      identity_moderate |
        gnrh_raw >= (min_umi + 1L)
    )

  # --------------------------------------------------------------------------- #
  # Supported candidates
  #
  # Supported cells already passed independent identity and support gates
  # during classification.
  # --------------------------------------------------------------------------- #

  supported_confident <-
    gnrh_class == "supported" &
    identity_moderate

  # --------------------------------------------------------------------------- #
  # Final confidence
  # --------------------------------------------------------------------------- #

  confident <-
    direct_confident |
    supported_confident

  confident[
    is.na(confident)
  ] <- FALSE

  confident
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
#' \code{gnrh_status} metadata columns.
#'
#' @return Data frame with false-positive rate (\code{fpr}) and
#' true-positive rate (\code{tpr}).
#'
#' @keywords internal
#' @noRd
.gnrh_build_roc <- function(object) {

  expr <- object$gnrh_expr

  true_gnrh <- object$gnrh_status %in% c("pos")

  true_gnrh <- factor(true_gnrh, levels = c(FALSE, TRUE))

  if (length(unique(true_gnrh)) != 2) {
    stop("ROC requires 2 true_gnrh: found ",
         paste(unique(true_gnrh), collapse = ", "))
  }

  roc <- pROC::roc(
    response = true_gnrh,
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
      confident  = if ("gnrh_confident" %in% colnames(md)) .count_factor(md$gnrh_confident) else NULL,
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



#' Gene Set Overlap Analysis and Visualization
#'
#' Perform overlap analysis between multiple gene sets, export overlap
#' tables, compute unique/common genes, and always generate a
#' publication-quality UpSet plot.
#'
#' @param gene_sets A named list of gene vectors. Each element should contain
#'   a character vector of gene symbols.
#' @param min_size Integer. Minimum intersection size to display in the UpSet
#'   plot. Default is \code{1}.
#' @param venn_title Character string specifying the UpSet plot title. The
#'   argument name is retained for backward compatibility.
#'   Default is \code{"Overlap of Gene Sets"}.
#' @param outdir Output directory where CSV tables and figures will be saved.
#'   Default is current working directory.
#' @param save_plot Logical indicating whether plots should be exported.
#'   Default is \code{TRUE}.
#' @param plot_width Numeric width of exported figures in inches.
#'   Default is \code{10}.
#' @param plot_height Numeric height of exported figures in inches.
#'   Default is \code{8}.
#' @param dpi Numeric resolution for PNG export. Default is \code{600}.
#'
#' @return A list containing cleaned gene sets, overlap tables, unique genes,
#'   common genes, the membership table, and the generated UpSet plot.
#'
#' @details
#' Prior to overlap analysis, gene symbols are standardized by:
#' \itemize{
#'   \item removing duplicated entries,
#'   \item removing missing values,
#'   \item removing empty strings,
#'   \item converting all gene symbols to uppercase.
#' }
#'
#' This ensures robust overlap comparisons across datasets originating
#' from different species or annotation conventions.
#'
#' @examples
#' \dontrun{
#' gene_sets <- list(
#'   Dataset_A = c("GNRH1", "KISS1", "TAC3"),
#'   Dataset_B = c("GNRH1", "TAC3", "PAX6"),
#'   Dataset_C = c("GNRH1", "DLX1", "DLX2")
#' )
#'
#' results <- gene_upset(
#'   gene_sets = gene_sets,
#'   venn_title = "GnRH Marker Overlap",
#'   outdir = "results/gene_overlap"
#' )
#'
#' results$upset_plot
#' }
#'
#' @export
gnrh_gene_upset <- function(
    gene_sets,
    min_size = 1,
    venn_title = "Overlap of Gene Sets",
    outdir = ".",
    save_plot = TRUE,
    plot_width = 10,
    plot_height = 8,
    dpi = 600
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' is required.", call. = FALSE)
  }

  if (!is.list(gene_sets)) {
    stop("'gene_sets' must be a named list.", call. = FALSE)
  }

  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("'gene_sets' must be a named list.", call. = FALSE)
  }

  if (length(gene_sets) < 2) {
    stop("'gene_sets' must contain at least two gene sets.", call. = FALSE)
  }

  min_size <- as.integer(min_size)
  if (is.na(min_size) || min_size < 1) {
    min_size <- 1L
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  labels <- names(gene_sets)

  gene_sets <- lapply(gene_sets, function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x)]
    x <- x[x != ""]
    # standardisation
    x <- toupper(x)
    sort(unique(x))
  })

  summary_df <- data.frame(
    dataset = labels,
    n_genes = vapply(gene_sets, length, integer(1)),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    summary_df,
    file.path(outdir, "gene_set_sizes.csv"),
    row.names = FALSE
  )

  combs <- utils::combn(labels, 2, simplify = FALSE)

  pairwise_results <- list()
  pairwise_summary <- data.frame(
    dataset1 = character(0),
    dataset2 = character(0),
    n_overlap = integer(0),
    stringsAsFactors = FALSE
  )

  for (cmb in combs) {

    overlap <- intersect(gene_sets[[cmb[1]]], gene_sets[[cmb[2]]])
    pair_name <- paste(cmb, collapse = "_vs_")

    pairwise_results[[pair_name]] <- overlap

    safe_pair_name <- gsub("[^A-Za-z0-9_\\-]+", "_", pair_name)

    utils::write.csv(
      data.frame(gene = overlap),
      file.path(outdir, paste0("genes_", safe_pair_name, ".csv")),
      row.names = FALSE
    )

    pairwise_summary <- rbind(
      pairwise_summary,
      data.frame(
        dataset1 = cmb[1],
        dataset2 = cmb[2],
        n_overlap = length(overlap),
        stringsAsFactors = FALSE
      )
    )
  }

  utils::write.csv(
    pairwise_summary,
    file.path(outdir, "pairwise_overlap_summary.csv"),
    row.names = FALSE
  )

  common_all <- Reduce(intersect, gene_sets)

  utils::write.csv(
    data.frame(gene = common_all),
    file.path(outdir, "genes_common_all.csv"),
    row.names = FALSE
  )

  unique_results <- list()
  unique_summary <- data.frame(
    dataset = character(0),
    n_unique = integer(0),
    stringsAsFactors = FALSE
  )

  for (lbl in labels) {

    others <- gene_sets[names(gene_sets) != lbl]

    unique_genes <- setdiff(
      gene_sets[[lbl]],
      Reduce(union, others)
    )

    unique_results[[lbl]] <- unique_genes

    safe_lbl <- gsub("[^A-Za-z0-9_\\-]+", "_", lbl)

    utils::write.csv(
      data.frame(gene = unique_genes),
      file.path(outdir, paste0("genes_unique_", safe_lbl, ".csv")),
      row.names = FALSE
    )

    unique_summary <- rbind(
      unique_summary,
      data.frame(
        dataset = lbl,
        n_unique = length(unique_genes),
        stringsAsFactors = FALSE
      )
    )
  }

  utils::write.csv(
    unique_summary,
    file.path(outdir, "unique_gene_summary.csv"),
    row.names = FALSE
  )

  all_genes <- sort(unique(unlist(gene_sets, use.names = FALSE)))

  membership <- data.frame(
    gene = all_genes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (lbl in labels) {
    membership[[lbl]] <- all_genes %in% gene_sets[[lbl]]
  }

  utils::write.csv(
    membership,
    file.path(outdir, "gene_set_membership.csv"),
    row.names = FALSE
  )

  plot <- ComplexUpset::upset(
    data = membership,
    intersect = labels,
    min_size = min_size,
    width_ratio = 0.18,
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        counts = TRUE,
        text = list(size = 3.5)
      )
    ),
    set_sizes = ComplexUpset::upset_set_size(),
    sort_sets = "descending",
    sort_intersections_by = "cardinality"
  ) +
    patchwork::plot_annotation(
      title = venn_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        )
      )
    )

  if (save_plot) {

    ggplot2::ggsave(
      filename = file.path(outdir, "gene_overlap_plot.pdf"),
      plot = plot,
      width = plot_width,
      height = plot_height,
      dpi = dpi,
      bg = "white"
    )

    ggplot2::ggsave(
      filename = file.path(outdir, "gene_overlap_plot.png"),
      plot = plot,
      width = plot_width,
      height = plot_height,
      dpi = dpi,
      bg = "white"
    )
  }

  list(
    gene_sets = gene_sets,
    summary = summary_df,
    pairwise = pairwise_results,
    pairwise_summary = pairwise_summary,
    common_all = common_all,
    unique = unique_results,
    unique_summary = unique_summary,
    membership = membership,
    plot = plot,
    upset_plot = plot,
    plot_type = "upset"
  )
}


#' Build gene sets from marker tables
#'
#' Reads marker tables from multiple datasets and extracts unique gene
#' symbols into a named list suitable for overlap analysis, UpSet plots, or
#' marker comparison workflows.
#'
#' Gene names are automatically standardized to uppercase to ensure
#' consistent comparisons across species and datasets.
#'
#' @param files Named character vector containing marker table filenames.
#' Names correspond to dataset identifiers and values correspond to file
#' names.
#' @param dir Character. Directory containing marker tables.
#' @param gene_col Character. Name of the column containing gene symbols.
#' Default is \code{"gene"}.
#'
#' @details
#' For each dataset:
#' \itemize{
#'   \item Marker tables are imported using
#'   \code{\link[utils]{read.delim}}.
#'   \item Missing values and empty gene names are removed.
#'   \item Duplicate genes are removed.
#'   \item Gene symbols are converted to uppercase.
#'   \item Gene symbols are sorted alphabetically.
#' }
#'
#' This standardization ensures robust overlap analysis between datasets
#' originating from different species or annotation conventions.
#'
#' @return
#' A named list where each element contains a character vector of unique
#' gene symbols for a dataset.
#'
#' @examples
#' \dontrun{
#'
#' files <- c(
#'   "HuDeCa Nose" = "gnrh_nose_markers.tsv",
#'   "HPSC Wang 2022" = "gnrh_wang_markers.tsv",
#'   "Human HypoMap" = "gnrh_human_hypomap_markers.tsv"
#' )
#'
#' gene_sets <- build_gene_sets(
#'   files = files,
#'   dir = file.path(outdir, "tables")
#' )
#'
#' names(gene_sets)
#' lengths(gene_sets)
#'
#' }
#'
#' @seealso
#' \code{\link{gene_upset}},
#' \code{\link{gnrh_marker_programs}}
#'
#' @export
build_gene_sets <- function(
    files,
    dir,
    gene_col = "gene"
) {

  if (is.null(names(files))) {
    stop("'files' must be a named vector.", call. = FALSE)
  }

  gene_sets <- stats::setNames(

    lapply(names(files), function(dataset) {

      file <- file.path(dir, files[[dataset]])

      if (!file.exists(file)) {
        warning("Missing file: ", file, call. = FALSE)
        return(character(0))
      }

      tab <- utils::read.delim(
        file,
        sep = "\t",
        stringsAsFactors = FALSE
      )

      if (!gene_col %in% colnames(tab)) {
        warning(
          "Column '", gene_col,
          "' not found in ", dataset,
          call. = FALSE
        )
        return(character(0))
      }

      genes <- unique(as.character(tab[[gene_col]]))

      genes <- genes[!is.na(genes)]
      genes <- genes[genes != ""]

      genes <- toupper(genes)

      sort(unique(genes))
    }),

    names(files)
  )

  gene_sets
}
