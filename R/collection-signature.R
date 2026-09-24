# ========================================================================= #
# External-signature validation across GnRHcell datasets
# ========================================================================= #

#' Plot an external GnRH signature across datasets
#'
#' Quantifies a supplied gene signature within cells classified as GnRH-positive
#' in each selected dataset. Point size represents the percentage of detected
#' cells expressing a gene. Point colour represents either the mean normalized
#' expression or its gene-wise standardized value across datasets. Rendering is
#' delegated to [gnrh_celldot()] so signature plots share the package-wide
#' dotplot style and behaviour.
#'
#' @param collection Result returned by [run_gnrh_collection()].
#' @param genes Character vector of signature genes in the desired display
#'   order.
#' @param selected_ids Dataset IDs or unique labels. When `NULL`, all datasets
#'   whose `species` value is `"Human"` are used.
#' @param signature_name Optional figure title.
#' @param status_field Metadata column defining GnRH detection.
#' @param positive_value Value in `status_field` identifying detected cells.
#' @param layer Seurat expression layer passed to [Seurat::FetchData()].
#' @param scale_expression Standardize mean expression separately for each gene
#'   across datasets. Recommended when datasets were normalized independently.
#' @param show_title Display `signature_name` as the main plot title. The
#'   default is `FALSE`, matching [gnrh_celldot()].
#' @param txtsize Base text size.
#' @param x.ang Rotation of dataset labels.
#' @param dot.scale Maximum point size.
#' @param compact_labels Use concise, multi-line dataset labels. Disabled by
#'   default to preserve the standard [gnrh_celldot()] display.
#' @param legend_position Legend position, either `"bottom"` or `"right"`.
#' @param th.cols,colors,style Palette and style passed to [gnrh_celldot()].
#' @param filename Optional output filename. When supplied, the plot is saved
#'   with [gnrh_save()].
#' @param width,height Optional output dimensions. When `NULL`, they are adapted
#'   to the number and length of dataset and gene labels.
#' @param units Output units passed to [gnrh_save()].
#' @param dpi Resolution for raster output.
#' @param ... Additional arguments passed to [gnrh_celldot()].
#'
#' @return A ggplot object. The plotting table, missing genes, selected dataset
#'   table, and recommended dimensions are stored as attributes.
#' @export
gnrh_signature <- function(
    collection,
    genes,
    selected_ids = NULL,
    signature_name = "External GnRH signature",
    status_field = "gnrh_status",
    positive_value = "pos",
    layer = "data",
    scale_expression = TRUE,
    show_title = FALSE,
    txtsize = getOption("gnrhcell.base_size", 14),
    x.ang = 60,
    dot.scale = 6,
    compact_labels = FALSE,
    legend_position = c("right", "bottom"),
    th.cols = "Reds",
    colors = NULL,
    style = "test",
    filename = NULL,
    width = NULL,
    height = NULL,
    units = "in",
    dpi = 600,
    ...
) {
  if (!is.list(collection) || is.null(collection$datasets)) {
    stop("`collection` must be returned by `run_gnrh_collection()`.", call. = FALSE)
  }
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) stop("`genes` must contain at least one gene.", call. = FALSE)
  if (!is.logical(scale_expression) || length(scale_expression) != 1L || is.na(scale_expression)) {
    stop("`scale_expression` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_title) || length(show_title) != 1L || is.na(show_title)) {
    stop("`show_title` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(compact_labels) || length(compact_labels) != 1L || is.na(compact_labels)) {
    stop("`compact_labels` must be TRUE or FALSE.", call. = FALSE)
  }
  legend_position <- match.arg(legend_position)

  datasets <- as.data.frame(collection$datasets, stringsAsFactors = FALSE)
  if (!all(c("id", "label") %in% names(datasets))) {
    stop("`collection$datasets` must contain `id` and `label`.", call. = FALSE)
  }
  if (is.null(selected_ids)) {
    if (!"species" %in% names(datasets)) {
      stop("Supply `selected_ids` because `collection$datasets$species` is unavailable.", call. = FALSE)
    }
    selected_ids <- datasets$id[tolower(trimws(datasets$species)) == "human"]
  }
  selected_ids <- as.character(selected_ids)
  resolve_row <- function(value) {
    hit <- unique(c(which(datasets$id == value), which(datasets$label == value)))
    if (length(hit) == 1L) hit[[1L]] else NA_integer_
  }
  rows <- vapply(selected_ids, resolve_row, integer(1))
  if (anyNA(rows)) {
    stop(
      "Dataset selection not found or ambiguous: ",
      paste(selected_ids[is.na(rows)], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  datasets_sel <- datasets[rows, , drop = FALSE]
  if (anyDuplicated(datasets_sel$id)) stop("Duplicated datasets in `selected_ids`.", call. = FALSE)

  all_objects <- .get_collection_objects(collection, allow_saved = TRUE)
  objects <- all_objects[datasets_sel$id]
  unavailable <- datasets_sel$id[vapply(objects, is.null, logical(1))]
  if (length(unavailable)) {
    stop("Processed object(s) unavailable: ", paste(unavailable, collapse = ", "), ".", call. = FALSE)
  }

  missing_by_dataset <- stats::setNames(vector("list", nrow(datasets_sel)), datasets_sel$id)
  signature_data <- dplyr::bind_rows(lapply(seq_len(nrow(datasets_sel)), function(i) {
    id <- datasets_sel$id[[i]]
    label <- datasets_sel$label[[i]]
    object <- objects[[id]]
    metadata <- object[[]]
    if (!status_field %in% colnames(metadata)) {
      warning("Metadata field `", status_field, "` is absent from ", id, ".", call. = FALSE)
      return(NULL)
    }
    positive <- as.character(metadata[[status_field]]) == positive_value
    cells <- rownames(metadata)[!is.na(positive) & positive]
    if (!length(cells)) {
      warning("No `", positive_value, "` cells were found in ", id, ".", call. = FALSE)
      return(NULL)
    }
    available <- intersect(genes, rownames(object))
    missing_by_dataset[[id]] <<- setdiff(genes, available)
    if (!length(available)) {
      warning("None of the signature genes are present in ", id, ".", call. = FALSE)
      return(NULL)
    }
    expression <- Seurat::FetchData(
      object = object,
      vars = available,
      cells = cells,
      layer = layer
    )
    tibble::tibble(
      dataset = id,
      label = label,
      gene = available,
      n_gnrh = length(cells),
      pct_detected = 100 * colMeans(expression > 0, na.rm = TRUE),
      mean_expr = colMeans(expression, na.rm = TRUE)
    )
  }))
  if (!nrow(signature_data)) stop("No signature values could be calculated.", call. = FALSE)

  signature_data <- signature_data |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(
      scaled_mean_expr = {
        value <- .data$mean_expr
        if (sum(is.finite(value)) > 1L && stats::sd(value, na.rm = TRUE) > 0) {
          as.numeric(scale(value))
        } else {
          rep(0, length(value))
        }
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      label = factor(.data$label, levels = datasets_sel$label),
      gene = factor(.data$gene, levels = rev(genes)),
      plot_expression = if (isTRUE(scale_expression)) .data$scaled_mean_expr else .data$mean_expr
    )

  compact_dataset_label <- function(value) {
    value <- sub("^Human\\s+", "", value, ignore.case = TRUE)
    value <- sub("^hPSC-derived\\s+GnRH", "hPSC", value, ignore.case = TRUE)
    value <- sub("^hPSC-derived", "hPSC", value, ignore.case = TRUE)
    value <- sub("^fetal nose", "Fetal nose", value, ignore.case = TRUE)
    value <- sub("^fetal hypothalamus", "Fetal hypothalamus", value, ignore.case = TRUE)
    value <- sub("^adult hypothalamus", "Adult hypothalamus", value, ignore.case = TRUE)
    value <- sub("^adult ME", "Adult ME", value, ignore.case = TRUE)
    sub(" \\(", "\n(", value)
  }
  display_labels <- if (isTRUE(compact_labels)) {
    vapply(datasets_sel$label, compact_dataset_label, character(1))
  } else {
    datasets_sel$label
  }
  names(display_labels) <- datasets_sel$label

  n_datasets <- nrow(datasets_sel)
  n_genes <- length(unique(signature_data$gene))
  if (is.null(width)) width <- max(8, min(14, 5.5 + 1.25 * n_datasets))
  if (is.null(height)) height <- max(7, min(14, 3.2 + 0.205 * n_genes))

  group_field <- ".gnrh_signature_dataset"
  combined <- combine_gnrh_datasets(
    collection = collection,
    selected_ids = datasets_sel$id,
    gnrh_only = FALSE,
    status_col = status_field,
    positive_values = positive_value,
    group.by = group_field,
    label_by = "label",
    allow_saved = TRUE
  )
  group_values <- as.character(combined[[group_field]][, 1L])
  positive_values <- as.character(combined[[status_field]][, 1L]) == positive_value
  group_values[is.na(positive_values) | !positive_values] <- ".__non_gnrh__"
  combined[[group_field]] <- factor(
    group_values,
    levels = c(datasets_sel$label, ".__non_gnrh__")
  )
  if (!identical(layer, "data")) {
    assay <- Seurat::DefaultAssay(combined)
    combined <- SeuratObject::SetAssayData(
      combined,
      assay = assay,
      layer = "data",
      new.data = SeuratObject::LayerData(combined, assay = assay, layer = layer)
    )
  }
  plot_features <- genes[genes %in% unique(signature_data$gene)]
  if (length(plot_features) == 1L) {
    auxiliary <- setdiff(rownames(combined), plot_features)
    if (length(auxiliary)) plot_features <- c(plot_features, auxiliary[[1L]])
  }
  plot <- gnrh_celldot(
    object = combined,
    features = plot_features,
    group.by = group_field,
    th.cols = th.cols,
    colors = colors,
    dot.scale = dot.scale,
    x.ang = x.ang,
    txtsize = txtsize,
    title = if (isTRUE(show_title)) signature_name else NULL,
    #leg.size = max(9, txtsize - 3),
    #leg.ttl.size = max(10, txtsize - 2),
    leg.pos = legend_position,
    style = style,
    scale = scale_expression,
    ...
  )

  # Merged Seurat objects contain the union of features. Remove combinations
  # where a gene was absent from the original dataset rather than displaying
  # an artificial zero-expression dot.
  valid_keys <- paste(as.character(signature_data$label), as.character(signature_data$gene), sep = "\r")
  plot_keys <- paste(as.character(plot$data$id), as.character(plot$data$features.plot), sep = "\r")
  plot$data <- plot$data[plot_keys %in% valid_keys, , drop = FALSE]
  plot$data$features.plot <- droplevels(plot$data$features.plot)

  if (isTRUE(compact_labels)) {
    plot <- plot + ggplot2::scale_x_discrete(labels = display_labels, drop = FALSE)
  }

  attr(plot, "gnrh_signature_data") <- signature_data
  attr(plot, "gnrh_missing_genes") <- missing_by_dataset
  attr(plot, "gnrh_datasets") <- datasets_sel
  attr(plot, "gnrh_dimensions") <- c(width = width, height = height)
  attr(plot, "gnrh_renderer") <- "gnrh_celldot"
  if (!is.null(filename)) {
    gnrh_save(
      plot,
      filename,
      width = width,
      height = height,
      units = units,
      dpi = dpi
    )
  }
  plot
}
