# ========================================================================= #== #
# GnRHcell multi-dataset workflows
# ========================================================================= #== #


# ========================================================================= #== #
# Internal helpers
# ========================================================================= #== #

#' Resolve a metadata column for dataset splitting
#'
#' @keywords internal
#' @noRd
resolve_split_column <- function(
    object,
    split_by = NULL
) {
  md <- object[[]]

  candidates <- unique(
    c(
      split_by,
      "orig.ident",
      "sample",
      "library_id"
    )
  )

  candidates <- candidates[
    !is.na(candidates) &
      nzchar(candidates)
  ]

  found <- candidates[
    candidates %in% colnames(md)
  ]

  if (!length(found))
    return(NULL)

  found[[1L]]
}


#' Save GnRHcell table
#'
#' @keywords internal
#' @noRd
.save_gnrh_table <- function(
    x,
    filename
) {
  utils::write.table(
    x = x,
    file = paste0(filename, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  invisible(filename)
}


#' Save a figure, table, or R object
#'
#' Chooses the writer from the filename extension. Supported figure formats
#' include PDF, SVG, EPS, PS, PNG, JPEG, TIFF, BMP, and WebP. Tables can be
#' written as CSV, TSV, TXT, or XLSX (when `writexl` is installed). Any R
#' object can be serialized as RDS, RDA, or RData. A figure filename without
#' an extension retains the GnRHcell convention of writing both PDF and PNG;
#' tables default to TSV and other objects default to RDS.
#'
#' @param x Object to save.
#' @param filename Output filename, with or without an extension.
#' @param format Optional output format. By default it is inferred from
#'   `filename` or from the object class when no extension is present.
#' @param width,height Figure dimensions.
#' @param units Figure dimension units passed to [ggplot2::ggsave()].
#' @param dpi Resolution for raster figures. Vector output ignores this value.
#' @param bg Figure background colour.
#' @param row.names Include row names in text tables.
#' @param object_name Name assigned to an object stored in RDA/RData output.
#' @param compress Compression passed to [saveRDS()] or [save()].
#' @param ... Additional arguments passed to [ggplot2::ggsave()] for figures.
#' @return The written filename or filenames, invisibly.
#' @export
gnrh_save <- function(x, filename, format = NULL, width = 7, height = 5,
                      units = "in", dpi = 600, bg = "white",
                      row.names = FALSE, object_name = "object",
                      compress = TRUE, ...) {
  if (is.null(x)) return(invisible(NULL))
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) || !nzchar(filename)) {
    stop("`filename` must be one non-empty path.", call. = FALSE)
  }
  is_plot <- inherits(x, c("ggplot", "patchwork", "grob", "gtable"))
  is_table <- is.data.frame(x) || is.matrix(x)
  ext <- tolower(tools::file_ext(filename))
  requested <- if (is.null(format)) "" else tolower(sub("^\\.", "", as.character(format)[1L]))
  aliases_match <- (requested %in% c("jpg", "jpeg") && ext %in% c("jpg", "jpeg")) ||
    (requested %in% c("tif", "tiff") && ext %in% c("tif", "tiff"))
  if (nzchar(requested) && nzchar(ext) && requested != ext && !aliases_match) {
    stop("`format` does not match the extension of `filename`.", call. = FALSE)
  }
  if (nzchar(requested) && !nzchar(ext)) {
    filename <- paste0(filename, ".", requested)
    ext <- requested
  }
  if (!nzchar(ext)) {
    if (is_plot) return(.gnrh_save_plot(x, filename, width, height, units, dpi, bg, ...))
    ext <- if (is_table) "tsv" else "rds"
    filename <- paste0(filename, ".", ext)
  }
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  figure_formats <- c("pdf", "svg", "eps", "ps", "png", "jpg", "jpeg",
                      "tif", "tiff", "bmp", "webp")
  if (ext %in% figure_formats) {
    if (!is_plot) stop("Figure formats require a ggplot, patchwork, grob, or gtable object.", call. = FALSE)
    return(.gnrh_save_plot(x, filename, width, height, units, dpi, bg, ...))
  }
  if (ext %in% c("csv", "tsv", "tab", "txt")) {
    if (!is_table) stop("Text-table formats require a data frame or matrix.", call. = FALSE)
    utils::write.table(
      x, filename, sep = if (ext == "csv") "," else "\t",
      quote = ext == "csv", row.names = row.names, col.names = TRUE
    )
  } else if (ext == "xlsx") {
    if (!is_table) stop("XLSX output requires a data frame or matrix.", call. = FALSE)
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Install the optional `writexl` package to save XLSX files.", call. = FALSE)
    }
    writexl::write_xlsx(as.data.frame(x), filename)
  } else if (ext == "rds") {
    saveRDS(x, filename, compress = compress)
  } else if (ext %in% c("rda", "rdata")) {
    object_name <- make.names(as.character(object_name)[1L])
    environment <- new.env(parent = emptyenv())
    assign(object_name, x, envir = environment)
    save(list = object_name, file = filename, envir = environment, compress = compress)
  } else {
    stop("Unsupported output format: `", ext, "`.", call. = FALSE)
  }
  invisible(filename)
}


# ========================================================================= #
# Adaptive plot dimensions
# ========================================================================= #
#' Plot dimensions
#'
#' @keywords internal
#' @noRd
.gnrh_plot_dims <- function(
    n_datasets=1L,
    n_genes=1L,
    type=c("heatmap","upset","program")
) {
  type <- match.arg(type)
  n_datasets <- max(1L,as.integer(n_datasets))
  n_genes <- max(1L,as.integer(n_genes))
  if (type=="heatmap") {
    width <- max(7,min(16,4.5+0.85*n_datasets))
    height <- max(5,min(20,3.5+0.13*n_genes))
  } else if (type=="upset") {
    width <- max(7,min(18,5.5+0.8*n_datasets))
    height <- max(5.5,min(11,4.5+0.4*n_datasets))
  } else {
    width <- max(8,min(20,5+2.4*ceiling(sqrt(n_datasets))))
    height <- max(6,min(20,4+2.2*ceiling(n_datasets/3)))
  }
  c(width=width,height=height)
}

#' Adaptive dot-plot dimensions
#' @keywords internal
#' @noRd
.gnrh_dot_dims <- function(labels, n_genes, n_facets = 1L) {
  labels <- as.character(labels)
  n_idents <- max(1L, length(labels))
  n_genes <- max(1L, as.integer(n_genes))
  n_facets <- max(1L, as.integer(n_facets))
  max_chars <- if (length(labels)) max(nchar(labels), na.rm = TRUE) else 0
  if (!is.finite(max_chars)) max_chars <- 0
  label_space <- min(1.8, max(0, (max_chars - 12) * 0.04))

  width <- if (n_facets == 1L) {
    4.5 + 0.85 * n_idents + label_space
  } else {
    4.0 + 0.70 * n_idents * n_facets + label_space
  }
  width <- max(6.5, min(18, width))
  height <- 3.6 + 0.23 * n_genes + min(2.2, max_chars * 0.035)
  height <- max(6, min(16, height))
  x_angle <- if (max_chars > 28) 55 else if (max_chars > 18) 40 else 20

  c(width = width, height = height, x_angle = x_angle)
}
# ========================================================================= #
# Adaptive text
# ========================================================================= #
.gnrh_text_size <- function(n,base=10,min_size=5,max_size=12) {
  size <- base*sqrt(20/max(20,n))
  max(min_size,min(max_size,size))
}
# ========================================================================= #
# Contrast-aware label colour
# ========================================================================= #
.gnrh_label_colour <- function(style="bw") {
  dark <- tolower(style) %in% c("dark","dirty")
  if (dark) "#F2F2F2" else "#111111"
}


# ========================================================================= #
# Save GnRHcell plot
# ========================================================================= #
#' Save GnRHcell ggplot
#'
#' @keywords internal
#' @noRd
.gnrh_save_plot <- function(
    plot,
    filename,
    width=7,
    height=5,
    units="in",
    dpi=600,
    bg="white",
    ...
) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE)
  ext <- tolower(tools::file_ext(filename))
  if (!nzchar(ext)) {
    files <- paste0(filename, c(".pdf", ".png"))
    ggplot2::ggsave(
      filename=files[[1L]],
      plot=plot,
      width=width,
      height=height,
      units=units,
      bg=bg,
      ...
    )
    ggplot2::ggsave(
      filename=files[[2L]],
      plot=plot,
      width=width,
      height=height,
      units=units,
      dpi=dpi,
      bg=bg,
      ...
    )
    return(invisible(files))
  } else {
    if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "bmp", "webp")) {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        bg = bg,
        ...
      )
    } else {
      # Vector devices do not use raster resolution. Recent ggplot2 versions
      # reject `dpi = NULL`, so omit the argument entirely for PDF/SVG output.
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = units,
        bg = bg,
        ...
      )
    }
  }
  invisible(filename)
}


#' Recover processed collection objects
#'
#' Uses in-memory objects when available and optionally falls back to
#' saved objects under output_dir/objects.
#'
#' @keywords internal
#' @noRd
.get_collection_objects <- function(
    collection,
    allow_saved = TRUE
) {
  ids <- as.character(
    collection$datasets$id
  )

  out <- stats::setNames(
    vector("list", length(ids)),
    ids
  )

  for (id in ids) {
    x <- collection$results[[id]]$object %||% NULL

    if (
      is.null(x) &&
      isTRUE(allow_saved)
    ) {
      file <- file.path(
        collection$output_dir,
        "objects",
        paste0(id, "_gnrh.rds")
      )

      if (file.exists(file))
        x <- readRDS(file)
    }

    out[[id]] <- x
  }

  out
}


#' Combine selected datasets from a GnRHcell collection
#'
#' Extracts processed Seurat objects from a result returned by
#' [run_gnrh_collection()], optionally retains only detected GnRH cells, adds a
#' dataset-origin metadata column, and merges the selected objects.
#'
#' @param collection Result returned by [run_gnrh_collection()].
#' @param selected_ids Character vector of dataset IDs or unique dataset labels.
#' @param gnrh_only Logical; retain only cells classified as GnRH-positive.
#' @param status_col Metadata column containing the detection status.
#' @param positive_values Values in \code{status_col} treated as positive.
#' @param group.by Name of the metadata column recording dataset origin.
#' @param label_by Use dataset \code{"label"} or \code{"id"} as origin values.
#' @param allow_saved Logical; if in-memory objects were cleaned, load objects
#'   saved under \code{collection$output_dir/objects} when available.
#'
#' @return A Seurat object containing the selected datasets.
#' @export
combine_gnrh_datasets <- function(
    collection,
    selected_ids,
    gnrh_only = TRUE,
    status_col = "gnrh_status",
    positive_values = "pos",
    group.by = "dataset",
    label_by = c("label", "id"),
    allow_saved = TRUE
) {
  label_by <- match.arg(label_by)
  if (!is.list(collection) || is.null(collection$results) || is.null(collection$datasets)) {
    stop("`collection` must be a result returned by `run_gnrh_collection()`.", call. = FALSE)
  }
  if (!is.character(selected_ids) || !length(selected_ids) || anyNA(selected_ids)) {
    stop("`selected_ids` must be a non-empty character vector.", call. = FALSE)
  }
  datasets <- as.data.frame(collection$datasets, stringsAsFactors = FALSE)
  if (!all(c("id", "label") %in% names(datasets))) {
    stop("`collection$datasets` must contain `id` and `label` columns.", call. = FALSE)
  }

  resolve_one <- function(value) {
    found <- unique(c(which(datasets$id == value), which(datasets$label == value)))
    if (!length(found)) return(NA_integer_)
    if (length(found) > 1L) {
      stop("Dataset selection is ambiguous: ", value, ". Use its unique `id`.", call. = FALSE)
    }
    found[[1L]]
  }
  rows <- vapply(selected_ids, resolve_one, integer(1))
  if (anyNA(rows)) {
    stop("Dataset(s) not found: ", paste(selected_ids[is.na(rows)], collapse = ", "), ".", call. = FALSE)
  }
  selected <- datasets[rows, , drop = FALSE]
  if (anyDuplicated(selected$id)) {
    stop("`selected_ids` resolves to duplicated datasets.", call. = FALSE)
  }

  all_objects <- .get_collection_objects(collection, allow_saved = allow_saved)
  objects <- all_objects[selected$id]
  missing_objects <- selected$id[vapply(objects, is.null, logical(1))]
  if (length(missing_objects)) {
    stop(
      "Processed object(s) unavailable: ", paste(missing_objects, collapse = ", "),
      ". Re-run with `clean_objects = FALSE` or enable saved objects.",
      call. = FALSE
    )
  }

  for (i in seq_along(objects)) {
    object <- objects[[i]]
    .check_seurat(object)
    if (isTRUE(gnrh_only)) {
      md <- object[[]]
      if (!status_col %in% names(md)) {
        stop("Metadata column not found in ", selected$id[[i]], ": ", status_col, ".", call. = FALSE)
      }
      keep <- rownames(md)[
        !is.na(md[[status_col]]) & as.character(md[[status_col]]) %in% positive_values
      ]
      if (!length(keep)) {
        stop("No GnRH-positive cells found in dataset: ", selected$id[[i]], ".", call. = FALSE)
      }
      object <- subset(object, cells = keep)
    }
    object[[group.by]] <- selected[[label_by]][[i]]
    objects[[i]] <- object
  }

  if (length(objects) == 1L) return(objects[[1L]])
  combined <- merge(
    x = objects[[1L]],
    y = objects[-1L],
    add.cell.ids = paste0("dataset", seq_along(objects)),
    merge.data = TRUE
  )
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    combined <- SeuratObject::JoinLayers(combined)
  }
  combined
}



#' Prepare a GnRHcell dataset configuration table
#' @param datasets Data frame describing datasets to be processed by the
#'   multi-dataset GnRHcell workflow.
#' @param check_files Logical. Check whether dataset files exist.
#' @param remove_missing Logical. Remove datasets whose input files cannot be
#'   found.
#'
#' @export
prepare_gnrh_datasets <- function(
    datasets,
    check_files = TRUE,
    remove_missing = FALSE
) {
  if (!is.data.frame(datasets))
    stop(
      "`datasets` must be a data frame or tibble.",
      call. = FALSE
    )

  required_columns <- c(
    "id",
    "label",
    "species",
    "file",
    "split_by",
    "reduction"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(datasets)
  )

  if (length(missing_columns))
    stop(
      "Missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )

  datasets <- tibble::as_tibble(
    datasets
  )

  # ----------------------------------------------------------------------- # #-- #
  # Standardize
  # ----------------------------------------------------------------------- # #-- #

  datasets$id <- trimws(
    as.character(datasets$id)
  )

  datasets$label <- trimws(
    as.character(datasets$label)
  )

  datasets$species <- stringr::str_to_title(
    trimws(
      as.character(
        datasets$species
      )
    )
  )

  datasets$file <- path.expand(
    trimws(
      as.character(
        datasets$file
      )
    )
  )

  datasets$split_by <- trimws(
    as.character(
      datasets$split_by
    )
  )

  datasets$reduction <- trimws(
    as.character(
      datasets$reduction
    )
  )

  # ----------------------------------------------------------------------- # #-- #
  # Required values
  # ----------------------------------------------------------------------- # #-- #

  for (column in c(
    "id",
    "label",
    "file"
  )) {
    bad <-
      is.na(datasets[[column]]) |
      !nzchar(datasets[[column]])

    if (any(bad))
      stop(
        "Column `",
        column,
        "` contains missing or empty values.",
        call. = FALSE
      )
  }

  if (anyDuplicated(datasets$id))
    stop(
      "Duplicated dataset IDs: ",
      paste(
        unique(
          datasets$id[
            duplicated(datasets$id)
          ]
        ),
        collapse = ", "
      ),
      call. = FALSE
    )

  # ----------------------------------------------------------------------- # #-- #
  # Species
  # ----------------------------------------------------------------------- # #-- #

  invalid_species <- setdiff(
    unique(
      datasets$species
    ),
    c(
      "Human",
      "Mouse"
    )
  )

  invalid_species <- invalid_species[
    !is.na(invalid_species) &
      nzchar(invalid_species)
  ]

  if (length(invalid_species))
    warning(
      "Unrecognized species: ",
      paste(
        invalid_species,
        collapse = ", "
      ),
      call. = FALSE
    )

  # ----------------------------------------------------------------------- # #-- #
  # File availability
  # ----------------------------------------------------------------------- # #-- #

  datasets$exists <- file.exists(
    datasets$file
  )

  if (
    isTRUE(check_files) &&
    any(!datasets$exists)
  ) {
    warning(
      "Missing dataset files: ",
      paste(
        datasets$id[
          !datasets$exists
        ],
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (isTRUE(remove_missing)) {
    datasets <- datasets[
      datasets$exists,
      ,
      drop = FALSE
    ]
  }

  datasets
}
