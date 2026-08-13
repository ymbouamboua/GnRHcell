#' Merge GnRH marker scores across datasets
#'
#' Builds a gene-by-dataset score matrix from marker tables returned by
#' [gnrh_markers()]. Marker tables may be supplied as data frames, as one list
#' of data frames, or as a named vector of TSV files.
#'
#' @param ... Marker data frames, one list of marker data frames, or a character
#'   vector of marker-table paths.
#' @param dataset_names Optional dataset labels.
#' @param dir Directory prepended to relative file paths.
#' @param gene_col Gene-symbol column.
#' @param score_col Score column. `NULL` selects the first available column from
#'   `score_candidates` independently for each table.
#' @param score_candidates Preferred GnRH marker-score columns.
#' @param gene_case Standardize genes to uppercase or preserve their case.
#' @param fill Value used for genes absent from a dataset. `NA_real_` is safest.
#'
#' @return A numeric gene-by-dataset matrix. Attribute `score_columns` records
#'   the score column selected for each dataset.
#' @export
merge_gnrh_marker_scores <- function(
    ...,
    dataset_names = NULL,
    dir = ".",
    gene_col = "gene",
    score_col = NULL,
    score_candidates = c("score", "avg_log2FC", "coexpr", "avg_logFC"),
    gene_case = c("upper", "asis"),
    fill = NA_real_
) {
  gene_case <- match.arg(gene_case)
  inputs <- list(...)
  if (length(inputs) == 1L && is.character(inputs[[1L]])) {
    inputs <- as.list(inputs[[1L]])
  }
  if (length(inputs) == 1L && is.list(inputs[[1L]]) &&
      !is.data.frame(inputs[[1L]])) inputs <- inputs[[1L]]
  if (!length(inputs)) stop("Supply at least one marker table.", call. = FALSE)

  if (is.character(inputs)) inputs <- as.list(inputs)
  if (all(vapply(inputs, function(x) is.character(x) && length(x) == 1L, logical(1)))) {
    paths <- vapply(inputs, function(x) file.path(dir, x), character(1))
    missing <- paths[!file.exists(paths)]
    if (length(missing)) stop("Missing marker files: ", paste(missing, collapse = ", "), call. = FALSE)
    inputs <- lapply(paths, utils::read.delim, check.names = FALSE, stringsAsFactors = FALSE)
  }

  dataset_names <- dataset_names %||% names(inputs)
  if (is.null(dataset_names) || any(!nzchar(dataset_names))) {
    dataset_names <- paste0("dataset_", seq_along(inputs))
  }
  if (length(dataset_names) != length(inputs)) {
    stop("`dataset_names` must match the number of marker tables.", call. = FALSE)
  }

  selected_columns <- character(length(inputs))
  scores <- lapply(seq_along(inputs), function(i) {
    tab <- inputs[[i]]
    if (!is.data.frame(tab)) stop("Every marker input must be a data frame or TSV file.", call. = FALSE)
    if (!gene_col %in% names(tab)) stop("Column `", gene_col, "` is missing from ", dataset_names[[i]], ".", call. = FALSE)
    selected <- score_col
    if (is.null(selected)) selected <- score_candidates[score_candidates %in% names(tab)][1L]
    if (!length(selected) || is.na(selected) || !selected %in% names(tab)) {
      stop("No supported score column found in ", dataset_names[[i]], ". Available columns: ", paste(names(tab), collapse = ", "), call. = FALSE)
    }
    selected_columns[[i]] <<- selected
    gene <- as.character(tab[[gene_col]])
    if (gene_case == "upper") gene <- toupper(gene)
    value <- suppressWarnings(as.numeric(tab[[selected]]))
    keep <- !is.na(gene) & nzchar(gene) & is.finite(value)
    out <- data.frame(gene = gene[keep], value = value[keep], stringsAsFactors = FALSE)
    out <- out[order(out$value, decreasing = TRUE), , drop = FALSE]
    out <- out[!duplicated(out$gene), , drop = FALSE]
    stats::setNames(out$value, out$gene)
  })

  genes <- sort(unique(unlist(lapply(scores, names), use.names = FALSE)))
  matrix <- do.call(cbind, lapply(scores, function(x) unname(x[genes])))
  if (is.null(dim(matrix))) matrix <- matrix(matrix, ncol = 1L)
  matrix[is.na(matrix)] <- fill
  rownames(matrix) <- genes
  colnames(matrix) <- dataset_names
  attr(matrix, "score_columns") <- stats::setNames(selected_columns, dataset_names)
  matrix
}

.gnrh_marker_heatmap <- function(matrix, row_split = NULL, scale_rows = TRUE,
                                 cluster_rows = TRUE, cluster_columns = FALSE,
                                 show_column_names = TRUE, fontsize_row = 8,
                                 fontsize_col = 9, fontsize_legend = 10,
                                 title = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Packages `ComplexHeatmap` and `circlize` are required for marker heatmaps.", call. = FALSE)
  }
  plot_matrix <- matrix
  if (isTRUE(scale_rows)) {
    plot_matrix <- t(vapply(seq_len(nrow(matrix)), function(i) {
      values <- matrix[i, ]
      observed <- is.finite(values)
      scaled <- rep(NA_real_, length(values))
      if (sum(observed) == 1L) {
        # A singleton is maximally dataset-specific but has no estimable SD.
        scaled[observed] <- 2
      } else if (sum(observed) > 1L) {
        scaled[observed] <- as.numeric(scale(values[observed]))
      }
      scaled
    }, numeric(ncol(matrix))))
    dimnames(plot_matrix) <- dimnames(matrix)
    colour <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
    legend_title <- "z-score"
  } else {
    rng <- range(plot_matrix, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0) rng <- c(-1, 1)
    colour <- circlize::colorRamp2(c(rng[1], mean(rng), rng[2]), c("#F7FBFF", "#6BAED6", "#08306B"))
    legend_title <- "Marker score"
  }

  # Preserve NA values for display, but use row-mean imputation only for
  # dendrogram calculation because stats::dist()/hclust() reject missing data.
  clustering_matrix <- plot_matrix
  for (i in seq_len(nrow(clustering_matrix))) {
    missing <- !is.finite(clustering_matrix[i, ])
    if (any(missing)) {
      observed_mean <- mean(clustering_matrix[i, !missing], na.rm = TRUE)
      if (!is.finite(observed_mean)) observed_mean <- 0
      clustering_matrix[i, missing] <- observed_mean
    }
  }
  row_clustering <- cluster_rows
  if (isTRUE(cluster_rows) && nrow(clustering_matrix) > 1L) {
    row_clustering <- stats::hclust(stats::dist(clustering_matrix))
  } else if (isTRUE(cluster_rows)) {
    row_clustering <- FALSE
  }
  column_clustering <- cluster_columns
  if (isTRUE(cluster_columns) && ncol(clustering_matrix) > 1L) {
    column_clustering <- stats::hclust(stats::dist(t(clustering_matrix)))
  } else if (isTRUE(cluster_columns)) {
    column_clustering <- FALSE
  }
  heatmap <- ComplexHeatmap::Heatmap(
    plot_matrix, name = legend_title, col = colour, row_split = row_split,
    cluster_rows = row_clustering, cluster_columns = column_clustering,
    show_column_names = show_column_names, show_row_names = TRUE,
    column_title = title, na_col = "#F2F2F2", border = TRUE,
    row_names_gp = grid::gpar(fontsize = fontsize_row),
    column_names_gp = grid::gpar(fontsize = fontsize_col),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = fontsize_legend, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = fontsize_legend - 1)
    )
  )
  list(
    heatmap = heatmap,
    plot_matrix = plot_matrix,
    clustering_matrix = clustering_matrix
  )
}

#' Plot conserved GnRH markers across datasets
#'
#' @inheritParams merge_gnrh_marker_scores
#' @param files Named character vector of GnRH marker TSV filenames or paths.
#'   Names are used as dataset labels unless `dataset_names` is supplied.
#' @param min_datasets Minimum datasets supporting a gene.
#' @param top_n Maximum genes displayed.
#' @param scale_rows Row-standardize scores before plotting.
#' @param cluster_rows,cluster_columns Cluster heatmap rows or columns.
#' @param show_column_names Show dataset labels.
#' @param fontsize_row,fontsize_col,fontsize_legend Text sizes.
#' @param filename Optional PDF filename. Nothing is written when `NULL`.
#' @param width,height PDF dimensions.
#'
#' @return A list with `heatmap`, `matrix`, `plot_matrix`, and `conserved`.
#' @export
plot_gnrh_conserved_markers <- function(
    files, dir = ".", gene_col = "gene", score_col = NULL,
    dataset_names = names(files), gene_case = c("upper", "asis"),
    min_datasets = 2, top_n = 50, scale_rows = TRUE,
    cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
    fontsize_row = 8, fontsize_col = 9, fontsize_legend = 10,
    filename = NULL, width = 7, height = 9
) {
  matrix <- merge_gnrh_marker_scores(files, dataset_names = dataset_names,
    dir = dir, gene_col = gene_col, score_col = score_col,
    gene_case = match.arg(gene_case), fill = NA_real_)
  support <- rowSums(!is.na(matrix))
  mean_score <- rowMeans(matrix, na.rm = TRUE)
  table <- data.frame(gene = rownames(matrix), mean_score, n_datasets = support,
    stringsAsFactors = FALSE)
  table <- table[table$n_datasets >= min_datasets, , drop = FALSE]
  table <- table[order(-table$n_datasets, -table$mean_score), , drop = FALSE]
  table <- utils::head(table, top_n)
  if (!nrow(table)) stop("No conserved GnRH markers meet `min_datasets`.", call. = FALSE)
  selected <- matrix[table$gene, , drop = FALSE]
  result <- .gnrh_marker_heatmap(selected, scale_rows = scale_rows,
    cluster_rows = cluster_rows, cluster_columns = cluster_columns,
    show_column_names = show_column_names, fontsize_row = fontsize_row,
    fontsize_col = fontsize_col, fontsize_legend = fontsize_legend,
    title = "Conserved GnRH markers")
  if (!is.null(filename)) {
    grDevices::pdf(filename, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(result$heatmap)
  }
  c(result, list(matrix = matrix, conserved = table))
}

#' Plot dataset-specific GnRH markers
#'
#' @inheritParams plot_gnrh_conserved_markers
#' @param top_n_per_dataset Maximum markers selected for each dataset.
#' @param max_datasets Maximum datasets supporting a dataset-specific marker.
#'
#' @return A list with `heatmap`, matrices, specificity scores, and
#'   `specific_table`.
#' @export
plot_gnrh_dataset_specific_markers <- function(
    files, dir = ".", gene_col = "gene", score_col = NULL,
    dataset_names = names(files), gene_case = c("upper", "asis"),
    top_n_per_dataset = 20, max_datasets = 2, scale_rows = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
    fontsize_row = 7, fontsize_col = 9, fontsize_legend = 10,
    filename = NULL, width = 8, height = 10
) {
  matrix <- merge_gnrh_marker_scores(files, dataset_names = dataset_names,
    dir = dir, gene_col = gene_col, score_col = score_col,
    gene_case = match.arg(gene_case), fill = NA_real_)
  support <- rowSums(!is.na(matrix))
  eligible <- support <= max_datasets
  specificity <- vapply(seq_len(ncol(matrix)), function(i) {
    target <- matrix[, i]
    others <- matrix[, -i, drop = FALSE]
    other_mean <- rowMeans(others, na.rm = TRUE)
    other_mean[!is.finite(other_mean)] <- 0
    value <- target - other_mean
    value[is.na(target) | !eligible] <- -Inf
    value
  }, numeric(nrow(matrix)))
  dimnames(specificity) <- dimnames(matrix)
  selected <- unique(unlist(lapply(seq_len(ncol(specificity)), function(i) {
    value <- sort(specificity[, i], decreasing = TRUE)
    names(utils::head(value[is.finite(value)], top_n_per_dataset))
  }), use.names = FALSE))
  if (!length(selected)) stop("No dataset-specific GnRH markers meet `max_datasets`.", call. = FALSE)
  best_index <- max.col(specificity[selected, , drop = FALSE], ties.method = "first")
  best_dataset <- colnames(matrix)[best_index]
  best_score <- specificity[cbind(match(selected, rownames(specificity)), best_index)]
  order_index <- order(match(best_dataset, colnames(matrix)), -best_score)
  selected <- selected[order_index]
  best_dataset <- best_dataset[order_index]
  best_score <- best_score[order_index]
  selected_matrix <- matrix[selected, , drop = FALSE]
  specific_table <- data.frame(
    gene = selected, best_dataset = best_dataset,
    n_datasets = support[selected], specificity_score = best_score,
    raw_score = selected_matrix[cbind(seq_along(selected), match(best_dataset, colnames(matrix)))],
    stringsAsFactors = FALSE
  )
  split <- factor(best_dataset, levels = colnames(matrix))
  result <- .gnrh_marker_heatmap(selected_matrix, row_split = split,
    scale_rows = scale_rows, cluster_rows = cluster_rows,
    cluster_columns = cluster_columns, show_column_names = show_column_names,
    fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    fontsize_legend = fontsize_legend,
    title = "Dataset-specific GnRH markers")
  if (!is.null(filename)) {
    grDevices::pdf(filename, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(result$heatmap)
  }
  c(result, list(matrix = matrix, specificity = specificity,
    specific_table = specific_table))
}
