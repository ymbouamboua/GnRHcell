# ============================================================================= #
# GnRHcell multi-dataset workflows
# ============================================================================= #


#' Resolve a dimensional reduction
#'
#' Selects the requested dimensional reduction from a Seurat object. If the
#' requested reduction is unavailable, the function falls back to `"umap"`
#' when present.
#'
#' @param object A Seurat object.
#' @param reduction Character scalar giving the requested dimensional
#'   reduction. Default is `"umap"`.
#'
#' @return A character scalar containing the name of the reduction to use.
#'
#' @keywords internal
#' @noRd
resolve_reduction <- function(object, reduction = "umap") {
  available <- names(object@reductions)
  
  if (reduction %in% available) {
    return(reduction)
  }
  
  if ("umap" %in% available) {
    warning(
      "Reduction '", reduction,
      "' not found; using 'umap'.",
      call. = FALSE
    )
    return("umap")
  }
  
  stop(
    "No usable UMAP reduction found. Available reductions: ",
    paste(available, collapse = ", "),
    call. = FALSE
  )
}


#' Resolve a metadata column for dataset splitting
#'
#' Identifies an available metadata column that can be used to split or group
#' cells in downstream GnRHcell plots. The requested column is preferred,
#' followed by `"orig.ident"`, `"sample"`, and `"library_id"`.
#'
#' @param object A Seurat object.
#' @param split_by Character scalar giving the preferred metadata column.
#'
#' @return A character scalar containing the selected metadata column, or
#'   `NULL` if none of the candidate columns are present.
#'
#' @keywords internal
#' @noRd
resolve_split_column <- function(object, split_by) {
  candidates <- unique(c(
    split_by,
    "orig.ident",
    "sample",
    "library_id"
  ))
  
  found <- candidates[
    candidates %in% colnames(object@meta.data)
  ]
  
  if (length(found) == 0L) {
    return(NULL)
  }
  
  found[[1]]
}


#' Run GnRHcell analysis on a single dataset
#'
#' Runs the GnRHcell workflow on a Seurat object and generates dataset-level
#' diagnostic plots, embeddings, feature plots, distribution plots, marker
#' tables, co-expression plots, and marker networks.
#'
#' The function is primarily used internally by [run_gnrh_collection()] but
#' can also be called directly for individual datasets.
#'
#' @param object A Seurat object containing the dataset to analyze.
#' @param dataset_id Character scalar giving a short unique identifier for the
#'   dataset. This identifier is used in output file names and directories.
#' @param dataset_label Character scalar giving a human-readable dataset name.
#' @param split_by Character scalar giving the metadata column used to split
#'   GnRH distributions. If unavailable, common alternatives such as
#'   `"orig.ident"`, `"sample"`, and `"library_id"` are considered.
#' @param reduction Character scalar giving the dimensional reduction used for
#'   embedding plots. Default is `"umap"`.
#' @param output_dir Character scalar giving the root output directory.
#' @param run_markers Logical. Whether to identify GnRH-associated markers and
#'   generate marker-based plots. Default is `TRUE`.
#' @param clean_object Logical. Whether to trigger garbage collection after
#'   processing the dataset. Default is `TRUE`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{object}{The processed Seurat object.}
#'   \item{run_info}{Dataset-level GnRHcell run information.}
#'   \item{markers}{GnRH marker results, or `NULL` when marker analysis was
#'   disabled.}
#'   \item{reduction}{The dimensional reduction used for plotting.}
#'   \item{split_by}{The metadata column used for distribution plots, or
#'   `NULL`.}
#' }
#'
#' @seealso
#' [GnRHcell::run_gnrh_collection()], [GnRHcell::run_gnrh()]
#'
#' @export
run_gnrh_dataset <- function(
    object,
    dataset_id,
    dataset_label,
    split_by = "orig.ident",
    reduction = "umap",
    output_dir,
    run_markers = TRUE,
    clean_object = TRUE
) {
  message("Running GnRHcell: ", dataset_label)
  dataset_dir <- file.path(
    output_dir,
    "figures",
    dataset_id
  )
  table_dir <- file.path(
    output_dir,
    "tables"
  )
  marker_dir <- file.path(
    output_dir,
    "markers"
  )
  invisible(lapply(
    c(
      dataset_dir,
      table_dir,
      marker_dir
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  # --------------------------------------------------------------------------- # #
  # Local helpers
  # --------------------------------------------------------------------------- # #
  save_table <- function(x, filename) {
    utils::write.table(
      x = x,
      file = paste0(filename, ".tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }
  save_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 300
  ) {
    ggplot2::ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf
    )
    ggplot2::ggsave(
      filename = paste0(filename, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi
    )
    invisible(NULL)
  }
  reduction <- resolve_reduction(
    object,
    reduction
  )
  split_by <- resolve_split_column(
    object,
    split_by
  )
  # --------------------------------------------------------------------------- # #
  # Run pipeline
  # --------------------------------------------------------------------------- # #
  object <- GnRHcell::run_gnrh(
    object
  )
  run_info <- GnRHcell::extract_gnrh_run_info(
    object,
    dataset_name = dataset_label
  )
  save_table(
    run_info,
    file.path(
      table_dir,
      paste0(
        dataset_id,
        "_gnrh_run_info"
      )
    )
  )
  # --------------------------------------------------------------------------- # #
  # QC report
  # --------------------------------------------------------------------------- # #
  report_plot <- GnRHcell::gnrh_report(
    object
  )
  save_plot(
    plot = report_plot,
    filename = file.path(
      dataset_dir,
      paste0(
        dataset_id,
        "_gnrh_report"
      )
    ),
    width = 16,
    height = 8
  )
  # --------------------------------------------------------------------------- # #
  # Embeddings
  # --------------------------------------------------------------------------- # #
  embedding_plot <- GnRHcell::plot_gnrh_embedding(
    object,
    group.by = c(
      "gnrh_status",
      "gnrh_confident",
      "gnrh_stage"
    ),
    reduction = reduction
  )
  save_plot(
    plot = embedding_plot,
    filename = file.path(
      dataset_dir,
      paste0(
        dataset_id,
        "_gnrh_embedding"
      )
    ),
    width = 15,
    height = 5
  )
  feature_plot <- GnRHcell::plot_gnrh_feature(
    object,
    feature_type = "all",
    reduction = reduction
  )
  save_plot(
    plot = feature_plot,
    filename = file.path(
      dataset_dir,
      paste0(
        dataset_id,
        "_gnrh_features"
      )
    ),
    width = 16,
    height = 4
  )
  # --------------------------------------------------------------------------- # #
  # Distribution plots
  # --------------------------------------------------------------------------- # #
  if (!is.null(split_by)) {
    status_plot <- GnRHcell::plot_gnrh_distribution(
      object,
      label = FALSE,
      proportion = TRUE,
      group.by = "gnrh_status",
      split.by = split_by,
      cols = GnRHcell::gnrh_colors("status")
    )
    save_plot(
      plot = status_plot,
      filename = file.path(
        dataset_dir,
        paste0(
          dataset_id,
          "_status_distribution"
        )
      ),
      width = 5,
      height = 4
    )
    stage_plot <- GnRHcell::plot_gnrh_distribution(
      object,
      label = FALSE,
      proportion = TRUE,
      group.by = "gnrh_stage",
      split.by = split_by,
      cols = GnRHcell::gnrh_colors("stage")
    )
    save_plot(
      plot = stage_plot,
      filename = file.path(
        dataset_dir,
        paste0(
          dataset_id,
          "_stage_distribution"
        )
      ),
      width = 5,
      height = 4
    )
  }
  # --------------------------------------------------------------------------- # #
  # Markers
  # --------------------------------------------------------------------------- # #
  markers <- NULL
  if (isTRUE(run_markers)) {
    markers <- GnRHcell::gnrh_markers(
      object
    )
    invisible(gc())
    save_table(
      markers,
      file.path(
        marker_dir,
        paste0(
          "gnrh_",
          dataset_id,
          "_markers"
        )
      )
    )
    if (
      "coexpr_flag" %in% colnames(markers) &&
      any(
        markers$coexpr_flag %in% TRUE,
        na.rm = TRUE
      )
    ) {
      coexpr <- markers[
        !is.na(markers$coexpr_flag) &
          markers$coexpr_flag,
        ,
        drop = FALSE
      ]
      coexpr_plot <- GnRHcell::plot_gnrh_coexpr(
        coexpr,
        coexp_cutoff = 0.3
      )
      save_plot(
        plot = coexpr_plot,
        filename = file.path(
          dataset_dir,
          paste0(
            dataset_id,
            "_coexpression"
          )
        ),
        width = 5,
        height = 6
      )
      network_plot <- GnRHcell::plot_network(
        markers,
        top_n = 50,
        threshold = 0.1
      )
      save_plot(
        plot = network_plot,
        filename = file.path(
          dataset_dir,
          paste0(
            dataset_id,
            "_network"
          )
        ),
        width = 10,
        height = 10
      )
    }
  }
  if (isTRUE(clean_object)) {
    invisible(gc())
  }
  list(
    object = object,
    run_info = run_info,
    markers = markers,
    reduction = reduction,
    split_by = split_by
  )
}



#' Run GnRHcell across multiple datasets
#'
#' Applies the GnRHcell workflow to a collection of Seurat datasets and
#' optionally performs cross-dataset marker comparisons and conserved marker
#' program analysis.
#'
#' Each row of `datasets` represents one dataset. Dataset files are loaded
#' sequentially to limit memory usage.
#'
#' @param datasets A data frame or tibble containing the columns `id`, `label`,
#'   `species`, `file`, `split_by`, and `reduction`.
#' @param output_dir Character scalar giving the root output directory.
#' @param run_markers Logical. Whether to identify GnRH markers for each
#'   dataset. Default is `TRUE`.
#' @param run_comparisons Logical. Whether to perform cross-dataset
#'   comparisons after all datasets have been processed. Default is `TRUE`.
#' @param run_programs Logical. Whether to identify conserved marker programs
#'   during cross-dataset comparisons. Default is `TRUE`.
#' @param clean_objects Logical. Whether processed Seurat objects should be
#'   removed from the returned dataset-level results to reduce memory usage.
#'   Default is `TRUE`.
#' @param save_objects Logical. Whether processed Seurat objects should be
#'   saved to disk. Default is `FALSE`.
#' @param verbose Logical. Whether to print progress messages. Default is
#'   `TRUE`.
#'
#' @return An object of class `"gnrh_collection"` containing:
#' \describe{
#'   \item{datasets}{The dataset configuration table used for the analysis.}
#'   \item{results}{Named list of dataset-level GnRHcell results.}
#'   \item{comparisons}{Cross-dataset comparison results, or `NULL`.}
#'   \item{output_dir}{Normalized output directory.}
#' }
#'
#' @seealso
#' [GnRHcell::prepare_gnrh_datasets()], [GnRHcell::run_gnrh_dataset()],
#' [GnRHcell::compare_gnrh_datasets()]
#'
#' @examples
#' \dontrun{
#' datasets <- data.frame(
#'   id = c("human_hpsc", "mouse_hypomap"),
#'   label = c("Human hPSC GnRH", "Mouse HypoMap"),
#'   species = c("Human", "Mouse"),
#'   file = c("human_hpsc.rds", "mouse_hypomap.rds"),
#'   split_by = c("orig.ident", "orig.ident"),
#'   reduction = c("umap", "umap")
#' )
#'
#' datasets <- prepare_gnrh_datasets(datasets)
#'
#' results <- run_gnrh_collection(
#'   datasets = datasets,
#'   output_dir = "gnrh_results"
#' )
#' }
#'
#' @export
run_gnrh_collection <- function(
    datasets,
    output_dir,
    run_markers = TRUE,
    run_comparisons = TRUE,
    run_programs = TRUE,
    clean_objects = TRUE,
    save_objects = FALSE,
    verbose = TRUE
) {
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
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing dataset columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  datasets$exists <- file.exists(
    path.expand(datasets$file)
  )
  
  unavailable <- datasets$id[
    !datasets$exists
  ]
  
  if (
    length(unavailable) > 0L &&
    isTRUE(verbose)
  ) {
    warning(
      "Skipping unavailable datasets: ",
      paste(unavailable, collapse = ", "),
      call. = FALSE
    )
  }
  
  datasets <- datasets[
    datasets$exists,
    ,
    drop = FALSE
  ]
  
  if (nrow(datasets) == 0L) {
    stop(
      "No dataset files were found.",
      call. = FALSE
    )
  }
  
  directories <- c(
    output_dir,
    file.path(output_dir, "figures"),
    file.path(output_dir, "tables"),
    file.path(output_dir, "markers"),
    file.path(output_dir, "comparisons"),
    file.path(output_dir, "objects")
  )
  
  invisible(lapply(
    directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  
  dataset_results <- setNames(
    vector("list", nrow(datasets)),
    datasets$id
  )
  
  for (i in seq_len(nrow(datasets))) {
    info <- datasets[
      i,
      ,
      drop = FALSE
    ]
    
    if (isTRUE(verbose)) {
      message(
        "[",
        i,
        "/",
        nrow(datasets),
        "] ",
        info$label
      )
    }
    
    object <- readRDS(
      path.expand(info$file)
    )
    
    result <- run_gnrh_dataset(
      object = object,
      dataset_id = info$id,
      dataset_label = info$label,
      split_by = info$split_by,
      reduction = info$reduction,
      output_dir = output_dir,
      run_markers = run_markers,
      clean_object = clean_objects
    )
    
    if (isTRUE(save_objects)) {
      saveRDS(
        result$object,
        file.path(
          output_dir,
          "objects",
          paste0(
            info$id,
            "_gnrh.rds"
          )
        )
      )
    }
    
    if (isTRUE(clean_objects)) {
      result$object <- NULL
    }
    
    dataset_results[[info$id]] <- result
    
    rm(object)
    invisible(gc())
  }
  
  comparison_results <- NULL
  
  if (isTRUE(run_comparisons)) {
    comparison_results <- compare_gnrh_datasets(
      datasets = datasets,
      output_dir = output_dir,
      run_programs = run_programs
    )
  }
  
  structure(
    list(
      datasets = datasets,
      results = dataset_results,
      comparisons = comparison_results,
      output_dir = normalizePath(
        output_dir,
        mustWork = FALSE
      )
    ),
    class = "gnrh_collection"
  )
}


#' Compare GnRHcell results across datasets
#'
#' Combines dataset-level GnRHcell outputs to compare runtime, detected GnRH
#' cells, marker overlap, and optionally conserved marker programs.
#'
#' The function expects outputs previously generated by
#' [GnRHcell::run_gnrh_collection()] or [GnRHcell::run_gnrh_dataset()].
#'
#' @param datasets A dataset configuration data frame containing at least
#'   `id` and `label`.
#' @param output_dir Character scalar giving the GnRHcell output directory.
#' @param run_programs Logical. Whether to compute conserved marker programs.
#'   Default is `TRUE`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{run_info_files}{Run-information files included in the comparison.}
#'   \item{marker_files}{Marker files included in the comparison.}
#'   \item{runtime_plot}{Runtime comparison plot, if available.}
#'   \item{detected_plot}{Detected GnRH cell comparison plot, if available.}
#'   \item{overlap}{Marker overlap results, if at least two marker files are
#'   available.}
#'   \item{programs}{Conserved marker program results, if requested and
#'   available.}
#' }
#'
#' @seealso
#' [GnRHcell::run_gnrh_collection()], [GnRHcell::gnrh_marker_programs()]
#'
#' @export
compare_gnrh_datasets <- function(
    datasets,
    output_dir,
    run_programs = TRUE
) {
  table_dir <- file.path(
    output_dir,
    "tables"
  )
  
  marker_dir <- file.path(
    output_dir,
    "markers"
  )
  
  comparison_dir <- file.path(
    output_dir,
    "comparisons"
  )
  
  dir.create(
    comparison_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  run_info_files <- file.path(
    table_dir,
    paste0(
      datasets$id,
      "_gnrh_run_info.tsv"
    )
  )
  
  run_info_files <- run_info_files[
    file.exists(run_info_files)
  ]
  
  marker_files <- setNames(
    paste0(
      "gnrh_",
      datasets$id,
      "_markers.tsv"
    ),
    datasets$label
  )
  
  marker_files <- marker_files[
    file.exists(
      file.path(
        marker_dir,
        marker_files
      )
    )
  ]
  
  runtime_plot <- NULL
  detected_plot <- NULL
  overlap <- NULL
  programs <- NULL
  
  if (length(run_info_files) > 0L) {
    runtime_plot <- plot_gnrh_runtime_curve(
      files = run_info_files,
      show_points = TRUE,
      txtsize = 8,
      x.ang = 60
    )
    
    detected_plot <- plot_gnrh_detected(
      files = run_info_files,
      txtsize = 8,
      x.ang = 60
    )
  }
  
  if (length(marker_files) >= 2L) {
    gene_sets <- build_gene_sets(
      files = marker_files,
      dir = marker_dir
    )
    
    overlap <- gene_upset(
      gene_sets = gene_sets,
      venn_title = "Overlap of GnRH co-expressed markers",
      outdir = file.path(
        comparison_dir,
        "marker_overlap"
      )
    )
  }
  
  if (
    isTRUE(run_programs) &&
    !is.null(overlap) &&
    length(marker_files) >= 2L
  ) {
    programs <- gnrh_marker_programs(
      files = marker_files,
      results = overlap,
      outdir = file.path(
        comparison_dir,
        "marker_programs"
      ),
      min_high_score = 0,
      min_medium_score = 0,
      write_output = TRUE
    )
  }
  
  list(
    run_info_files = run_info_files,
    marker_files = marker_files,
    runtime_plot = runtime_plot,
    detected_plot = detected_plot,
    overlap = overlap,
    programs = programs
  )
}


#' Prepare a GnRHcell dataset configuration table
#'
#' Validates and standardizes a dataset configuration table before running
#' multi-dataset GnRHcell analyses.
#'
#' @param datasets A data frame or tibble containing dataset configuration
#'   information.
#' @param check_files Logical. Whether to warn when dataset files are missing.
#'   Default is \code{TRUE}.
#' @param remove_missing Logical. Whether rows corresponding to missing files
#'   should be removed. Default is \code{FALSE}.
#'
#' @return A tibble containing the standardized dataset configuration and an
#'   additional logical column named \code{exists}.
#'
#' @export
prepare_gnrh_datasets <- function(
    datasets,
    check_files = TRUE,
    remove_missing = FALSE
) {
  
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
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  datasets <- tibble::as_tibble(
    datasets
  )
  
  # --------------------------------------------------------------------------- #
  # Standardize columns
  # --------------------------------------------------------------------------- #
  
  datasets$id <- as.character(
    datasets$id
  )
  
  datasets$label <- as.character(
    datasets$label
  )
  
  datasets$species <- stringr::str_to_title(
    as.character(
      datasets$species
    )
  )
  
  datasets$file <- path.expand(
    as.character(
      datasets$file
    )
  )
  
  datasets$split_by <- as.character(
    datasets$split_by
  )
  
  datasets$reduction <- as.character(
    datasets$reduction
  )
  
  datasets$exists <- file.exists(
    datasets$file
  )
  
  # --------------------------------------------------------------------------- #
  # Validate IDs
  # --------------------------------------------------------------------------- #
  
  duplicated_ids <- unique(
    datasets$id[
      duplicated(
        datasets$id
      )
    ]
  )
  
  if (length(duplicated_ids) > 0L) {
    stop(
      "Duplicated dataset IDs: ",
      paste(
        duplicated_ids,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # --------------------------------------------------------------------------- #
  # Validate species
  # --------------------------------------------------------------------------- #
  
  invalid_species <- setdiff(
    unique(
      datasets$species
    ),
    c(
      "Human",
      "Mouse"
    )
  )
  
  if (length(invalid_species) > 0L) {
    warning(
      "Unrecognized species: ",
      paste(
        invalid_species,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # --------------------------------------------------------------------------- #
  # Check dataset files
  # --------------------------------------------------------------------------- #
  
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
  
  # --------------------------------------------------------------------------- #
  # Optionally remove missing datasets
  # --------------------------------------------------------------------------- #
  
  if (isTRUE(remove_missing)) {
    datasets <- datasets[
      datasets$exists,
      ,
      drop = FALSE
    ]
  }
  
  datasets
}
