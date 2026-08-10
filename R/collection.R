# =========================================================================== #
# GnRHcell multi-dataset workflows
# =========================================================================== #


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
  # --------------------------------------------------------------------------- # # #
  # Local helpers
  # --------------------------------------------------------------------------- # # #
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
  # --------------------------------------------------------------------------- # # #
  # Run pipeline
  # --------------------------------------------------------------------------- # # #
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
  # --------------------------------------------------------------------------- # # #
  # QC report
  # --------------------------------------------------------------------------- # # #
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
  # --------------------------------------------------------------------------- # # #
  # Embeddings
  # --------------------------------------------------------------------------- # # #
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
  # --------------------------------------------------------------------------- # # #
  # Distribution plots
  # --------------------------------------------------------------------------- # # #
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
  # --------------------------------------------------------------------------- # # #
  # Markers
  # --------------------------------------------------------------------------- # # #
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

  # --------------------------------------------------------------------------- # #
  # Standardize columns
  # --------------------------------------------------------------------------- # #

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

  # --------------------------------------------------------------------------- # #
  # Validate IDs
  # --------------------------------------------------------------------------- # #

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

  # --------------------------------------------------------------------------- # #
  # Validate species
  # --------------------------------------------------------------------------- # #

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

  # --------------------------------------------------------------------------- # #
  # Check dataset files
  # --------------------------------------------------------------------------- # #

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

  # --------------------------------------------------------------------------- # #
  # Optionally remove missing datasets
  # --------------------------------------------------------------------------- # #

  if (isTRUE(remove_missing)) {
    datasets <- datasets[
      datasets$exists,
      ,
      drop = FALSE
    ]
  }

  datasets
}




#' Validate a GnRHcell multi-dataset collection
#'
#' Performs cross-dataset validation of GnRH detection and developmental
#' staging results generated by \code{\link{run_gnrh_collection}}.
#'
#' The function summarizes GnRH detection classes, consistency between
#' detection status and class, evidence scores, developmental stages,
#' stage refinement, migration-specific evidence, and selected biological
#' marker expression across datasets.
#'
#' Processed Seurat objects must be present in the collection object.
#' Therefore, \code{run_gnrh_collection()} should normally be called with
#' \code{clean_objects = FALSE} before running this function.
#'
#' @param collection An object of class \code{"gnrh_collection"} returned by
#'   \code{\link{run_gnrh_collection}}.
#' @param output_dir Optional directory in which validation tables are written.
#'   Default is \code{file.path(collection$output_dir, "validation")}.
#' @param positive_classes Character vector defining positive GnRH detection
#'   classes. Default is
#'   \code{c("direct", "supported")}.
#' @param validation_markers Character vector of genes used for independent
#'   biological marker validation. If \code{NULL}, a predefined marker panel
#'   covering GnRH identity, migration, and neuroendocrine maturation is used.
#' @param assay Assay used for biological marker expression summaries.
#'   Default is \code{"RNA"}.
#' @param layer Expression layer used for biological marker expression
#'   summaries. Default is \code{"data"}.
#' @param write_output Logical. Whether validation tables should be written
#'   to disk. Default is \code{TRUE}.
#' @param verbose Logical. Whether to print progress messages.
#'   Default is \code{TRUE}.
#'
#' @return An object of class \code{"gnrh_validation"} containing:
#' \describe{
#'   \item{\code{datasets}}{Dataset metadata used for validation.}
#'   \item{\code{input_summary}}{Number of cells and features per dataset.}
#'   \item{\code{detection}}{GnRH detection summary.}
#'   \item{\code{dropout_candidates}}{
#'   Summary of GNRH1-negative cells showing strong GnRH-like transcriptomic
#'   evidence. These cells are diagnostic candidates only and are not counted
#'   as GnRH-positive.}
#'   \item{\code{status_class}}{Consistency between detection status and class.}
#'   \item{\code{scores}}{Summary of GnRH evidence scores by detection class.}
#'   \item{\code{stages}}{Developmental-stage distribution among GnRH-positive
#'   cells.}
#'   \item{\code{stage_class}}{Developmental-stage distribution within each
#'   GnRH detection class.}
#'   \item{\code{stage_refinement}}{Summary of raw-to-final stage refinement.}
#'   \item{\code{stage_reassignment}}{Detailed raw-to-final stage transitions.}
#'   \item{\code{migration_core}}{Migration-core evidence by raw stage.}
#'   \item{\code{migration_refinement}}{
#'   Comparison of migration-core evidence between retained and reassigned
#'   migrating cells.}
#'   \item{\code{migration_refinement_summary}}{
#'   Dataset-level summary of migration-stage refinement.}
#'   \item{\code{biological_markers}}{
#'   Expression of selected biological validation markers by positive
#'   GnRH detection class.}
#'   \item{\code{output_dir}}{Validation output directory, or \code{NULL}.}
#' }
#'
#' @seealso
#' \code{\link{run_gnrh_collection}},
#' \code{\link{run_gnrh_dataset}},
#' \code{\link{stage_gnrh}}
#'
#' @examples
#' \dontrun{
#' results <- run_gnrh_collection(
#'   datasets = datasets,
#'   output_dir = "gnrh_results",
#'   clean_objects = FALSE
#' )
#'
#' validation <- validate_gnrh_collection(
#'   results
#' )
#'
#' validation$detection
#' validation$stage_refinement
#' validation$migration_core
#' validation$migration_refinement
#' validation$migration_refinement_summary
#' }
#'
#' @export
validate_gnrh_collection <- function(
    collection,
    output_dir = file.path(
      collection$output_dir,
      "validation"
    ),
    positive_classes = c(
      "direct",
      "supported"
    ),
    validation_markers = NULL,
    assay = "RNA",
    layer = "data",
    write_output = TRUE,
    verbose = TRUE
) {

  # --------------------------------------------------------------------------- #
  # Input validation
  # --------------------------------------------------------------------------- #

  if (!inherits(collection, "gnrh_collection")) {
    stop(
      "`collection` must be an object of class `gnrh_collection`.",
      call. = FALSE
    )
  }

  if (
    is.null(collection$results) ||
    length(collection$results) == 0L
  ) {
    stop(
      "`collection$results` is empty.",
      call. = FALSE
    )
  }

  if (
    is.null(collection$datasets) ||
    !is.data.frame(collection$datasets)
  ) {
    stop(
      "`collection$datasets` is missing or invalid.",
      call. = FALSE
    )
  }

  if (
    length(positive_classes) == 0L ||
    anyNA(positive_classes)
  ) {
    stop(
      "`positive_classes` must contain at least one non-missing class.",
      call. = FALSE
    )
  }

  positive_classes <- unique(
    as.character(positive_classes)
  )

  allowed_positive_classes <- c(
    "direct",
    "supported"
  )

  invalid_positive_classes <- setdiff(
    positive_classes,
    allowed_positive_classes
  )

  if (length(invalid_positive_classes) > 0L) {
    stop(
      "Unsupported positive GnRH class(es): ",
      paste(
        invalid_positive_classes,
        collapse = ", "
      ),
      ". Allowed classes are: direct, supported.",
      call. = FALSE
    )
  }

  required_dataset_columns <- c(
    "id",
    "label",
    "species"
  )

  missing_dataset_columns <- setdiff(
    required_dataset_columns,
    colnames(collection$datasets)
  )

  if (length(missing_dataset_columns) > 0L) {
    stop(
      "Missing dataset metadata columns: ",
      paste(
        missing_dataset_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # --------------------------------------------------------------------------- #
  # Local logging
  # --------------------------------------------------------------------------- #

  log <- function(...) {
    if (isTRUE(verbose)) {
      message(...)
    }
  }

  log("==== GNRH COLLECTION VALIDATION START ====")


  # --------------------------------------------------------------------------- #
  # Extract processed Seurat objects
  # --------------------------------------------------------------------------- #

  gnrh_list <- lapply(
    collection$results,
    function(x) {
      x$object
    }
  )

  missing_objects <- names(gnrh_list)[
    vapply(
      gnrh_list,
      is.null,
      logical(1)
    )
  ]

  if (length(missing_objects) > 0L) {
    stop(
      "Processed Seurat objects are unavailable for: ",
      paste(
        missing_objects,
        collapse = ", "
      ),
      ". Re-run `run_gnrh_collection()` with ",
      "`clean_objects = FALSE`.",
      call. = FALSE
    )
  }

  valid_seurat <- vapply(
    gnrh_list,
    inherits,
    logical(1),
    what = "Seurat"
  )

  if (!all(valid_seurat)) {
    stop(
      "Invalid processed object(s): ",
      paste(
        names(gnrh_list)[!valid_seurat],
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # --------------------------------------------------------------------------- #
  # Dataset metadata
  # --------------------------------------------------------------------------- #

  datasets <- tibble::as_tibble(
    collection$datasets
  )

  datasets$id <- as.character(
    datasets$id
  )

  dataset_metadata <- datasets |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species
    )

  if (is.null(names(gnrh_list))) {
    stop(
      "`collection$results` must be named using dataset IDs.",
      call. = FALSE
    )
  }

  expected_ids <- datasets$id[
    datasets$id %in% names(gnrh_list)
  ]

  if (length(expected_ids) == 0L) {
    stop(
      "No dataset IDs in `collection$datasets` match ",
      "`collection$results`.",
      call. = FALSE
    )
  }

  gnrh_list <- gnrh_list[
    expected_ids
  ]

  dataset_metadata <- dataset_metadata |>
    dplyr::filter(
      .data$id %in% expected_ids
    )


  # --------------------------------------------------------------------------- #
  # Required GnRH metadata
  # --------------------------------------------------------------------------- #

  required_columns <- c(
    "gnrh_status",
    "gnrh_class",
    "gnrh_stage"
  )

  for (id in names(gnrh_list)) {

    md <- gnrh_list[[id]][[]]

    missing_columns <- setdiff(
      required_columns,
      colnames(md)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Dataset `",
        id,
        "` is missing GnRHcell metadata column(s): ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }


  # --------------------------------------------------------------------------- #
  # Default biological validation markers
  # --------------------------------------------------------------------------- #

  if (is.null(validation_markers)) {

    validation_markers <- c(
      "GNRH1",
      "GNRHR",
      "FEZF1",
      "ISL1",
      "OTX2",
      "SIX3",
      "SIX6",
      "CHGA",
      "CHGB",
      "SCG2",
      "PCSK1",
      "PCSK2",
      "CPE",
      "VGF",
      "SYP",
      "ANOS1",
      "PROKR2",
      "PROK2",
      "NRP1",
      "NRP2",
      "ROBO1",
      "ROBO2",
      "DCX"
    )
  }

  validation_markers <- unique(
    as.character(validation_markers)
  )


  # =========================================================================== #
  # Input summary
  # =========================================================================== #

  log("Summarizing datasets")

  input_summary <- purrr::imap_dfr(
    gnrh_list,
    function(object, id) {

      tibble::tibble(
        id = id,
        n_cells = ncol(object),
        n_features = nrow(object)
      )
    }
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      .data$n_cells,
      .data$n_features
    )


  # =========================================================================== #
  # Detection summary
  # =========================================================================== #

  log("Validating GnRH detection")

  summarise_detection <- function(
    object,
    id
  ) {

    md <- object[[]]

    classes <- table(
      as.character(
        md$gnrh_class
      ),
      useNA = "no"
    )

    count_class <- function(x) {

      if (x %in% names(classes)) {
        return(
          unname(
            as.integer(
              classes[[x]]
            )
          )
        )
      }

      0L
    }

    direct <- count_class(
      "direct"
    )

    supported <- count_class(
      "supported"
    )

    positive <- sum(
      as.character(
        md$gnrh_class
      ) %in%
        positive_classes,
      na.rm = TRUE
    )

    total <- nrow(md)

    dropout_candidates <-
      if (
        "gnrh_dropout_candidate" %in%
        colnames(md)
      ) {

        sum(
          md$gnrh_dropout_candidate %in% TRUE,
          na.rm = TRUE
        )

      } else {

        0L
      }

    direct_supported <-
      if (
        "gnrh_direct_supported" %in%
        colnames(md)
      ) {

        sum(
          md$gnrh_direct_supported %in% TRUE,
          na.rm = TRUE
        )

      } else {

        NA_integer_
      }

    direct_isolated <-
      if (
        "gnrh_direct_isolated" %in%
        colnames(md)
      ) {

        sum(
          md$gnrh_direct_isolated %in% TRUE,
          na.rm = TRUE
        )

      } else {

        NA_integer_
      }

    confident <-
      if (
        "gnrh_confident" %in%
        colnames(md)
      ) {

        sum(
          md$gnrh_confident %in% TRUE,
          na.rm = TRUE
        )

      } else {

        NA_integer_
      }

    tibble::tibble(
      id = id,

      n_cells = total,

      direct = direct,
      supported = supported,

      direct_supported =
        direct_supported,

      direct_isolated =
        direct_isolated,

      dropout_candidates =
        dropout_candidates,

      gnrh_pos = positive,

      gnrh_confident =
        confident,

      pct_gnrh =
        if (total > 0L) {
          100 * positive / total
        } else {
          NA_real_
        },

      pct_direct =
        if (positive > 0L) {
          100 * direct / positive
        } else {
          NA_real_
        },

      pct_supported =
        if (positive > 0L) {
          100 * supported / positive
        } else {
          NA_real_
        },

      pct_direct_supported =
        if (
          direct > 0L &&
          !is.na(direct_supported)
        ) {
          100 *
            direct_supported /
            direct
        } else {
          NA_real_
        },

      pct_direct_isolated =
        if (
          direct > 0L &&
          !is.na(direct_isolated)
        ) {
          100 *
            direct_isolated /
            direct
        } else {
          NA_real_
        },

      pct_confident =
        if (
          positive > 0L &&
          !is.na(confident)
        ) {
          100 *
            confident /
            positive
        } else {
          NA_real_
        },

      pct_direct_all =
        if (total > 0L) {
          100 * direct / total
        } else {
          NA_real_
        },

      pct_supported_all =
        if (total > 0L) {
          100 * supported / total
        } else {
          NA_real_
        },

      pct_dropout_candidates =
        if (total > 0L) {
          100 *
            dropout_candidates /
            total
        } else {
          NA_real_
        }
    )
  }

  detection <- purrr::imap_dfr(
    gnrh_list,
    summarise_detection
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      dplyr::everything()
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$pct_gnrh
      )
    )



  # =========================================================================== #
  # GNRH1-negative transcriptomic candidates
  #
  # Diagnostic only.
  # These cells remain gnrh_status == "neg".
  # =========================================================================== #

  log(
    "Summarizing GNRH1-negative transcriptomic candidates"
  )

  has_dropout_candidate <- vapply(
    gnrh_list,
    function(object) {

      "gnrh_dropout_candidate" %in%
        colnames(
          object[[]]
        )
    },
    logical(1)
  )

  if (any(has_dropout_candidate)) {

    dropout_candidates <- purrr::imap_dfr(
      gnrh_list[has_dropout_candidate],
      function(object, id) {

        md <- object[[]] |>
          tibble::as_tibble()

        x <- md |>
          dplyr::filter(
            .data$gnrh_dropout_candidate %in% TRUE
          )

        if (nrow(x) == 0L) {

          return(
            tibble::tibble(
              id = id,
              n_cells = 0L,
              pct_dataset = 0,
              pct_GNRH1_detected = 0,
              median_GNRH1 = 0,
              median_support = NA_real_,
              median_identity_primary = NA_real_,
              median_core = NA_real_,
              median_neuro = NA_real_,
              median_knn = NA_real_,
              median_alternative = NA_real_
            )
          )
        }

        get_median <- function(column) {

          if (!column %in% colnames(x)) {
            return(NA_real_)
          }

          stats::median(
            x[[column]],
            na.rm = TRUE
          )
        }

        tibble::tibble(
          id = id,

          n_cells =
            nrow(x),

          pct_dataset =
            100 *
            nrow(x) /
            nrow(md),

          pct_GNRH1_detected =
            if (
              "gnrh_raw" %in%
              colnames(x)
            ) {

              100 *
                mean(
                  x$gnrh_raw > 0,
                  na.rm = TRUE
                )

            } else {

              NA_real_
            },

          median_GNRH1 =
            get_median(
              "gnrh_raw"
            ),

          median_support =
            get_median(
              "gnrh_support_score"
            ),

          median_identity_primary =
            get_median(
              "gnrh_identity_primary_hits"
            ),

          median_core =
            get_median(
              "gnrh_core_hits"
            ),

          median_neuro =
            get_median(
              "gnrh_neuro_hits"
            ),

          median_knn =
            get_median(
              "gnrh_knn"
            ),

          median_alternative =
            get_median(
              "gnrh_alternative_score"
            )
        )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        dplyr::everything()
      )

  } else {

    dropout_candidates <- tibble::tibble()
  }


  # =========================================================================== #
  # Status-class consistency
  # =========================================================================== #

  log("Checking status/class consistency")

  status_class <- purrr::imap_dfr(
    gnrh_list,
    function(object, id) {

      object[[]] |>
        tibble::as_tibble() |>
        dplyr::mutate(
          gnrh_status = as.character(
            .data$gnrh_status
          ),
          gnrh_class = as.character(
            .data$gnrh_class
          )
        ) |>
        dplyr::count(
          .data$gnrh_status,
          .data$gnrh_class,
          name = "n_cells"
        ) |>
        dplyr::mutate(
          id = id,
          pct_cells =
            100 *
            .data$n_cells /
            sum(.data$n_cells)
        )
    }
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      .data$gnrh_status,
      .data$gnrh_class,
      .data$n_cells,
      .data$pct_cells
    ) |>
    dplyr::arrange(
      .data$id,
      dplyr::desc(
        .data$n_cells
      )
    )


  # =========================================================================== #
  # Classification consistency checks
  # =========================================================================== #

  classification_consistency <- purrr::imap_dfr(
    gnrh_list,
    function(object, id) {

      md <- object[[]]

      class <- as.character(
        md$gnrh_class
      )

      status <- as.character(
        md$gnrh_status
      )

      raw <-
        if (
          "gnrh_raw" %in%
          colnames(md)
        ) {
          md$gnrh_raw
        } else {
          rep(NA_real_, nrow(md))
        }

      dropout <-
        if (
          "gnrh_dropout_candidate" %in%
          colnames(md)
        ) {
          md$gnrh_dropout_candidate %in% TRUE
        } else {
          rep(FALSE, nrow(md))
        }

      tibble::tibble(
        id = id,

        n_positive_without_GNRH1 =
          sum(
            status == "pos" &
              raw <= 0,
            na.rm = TRUE
          ),

        n_negative_positive_class =
          sum(
            status == "neg" &
              class %in%
              positive_classes,
            na.rm = TRUE
          ),

        n_positive_negative_class =
          sum(
            status == "pos" &
              !class %in%
              positive_classes,
            na.rm = TRUE
          ),

        n_dropout_called_positive =
          sum(
            dropout &
              status == "pos",
            na.rm = TRUE
          )
      )
    }
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      dplyr::everything()
    )

  # =========================================================================== #
  # Score validation
  # =========================================================================== #

  log("Summarizing GnRH evidence scores")

  candidate_score_columns <- c(
    "gnrh_score",
    "gnrh_support_score",

    "gnrh_identity_score",
    "gnrh_migration_score",
    "gnrh_neuro_score",
    "gnrh_hormone_score",
    "gnrh_guidance_score",
    "gnrh_alternative_score",

    "gnrh_identity_primary_hits",
    "gnrh_identity_supportive_hits",
    "gnrh_core_hits",

    "gnrh_migration_primary_hits",
    "gnrh_migration_supportive_hits",
    "gnrh_mig_hits",

    "gnrh_neuro_primary_hits",
    "gnrh_neuro_supportive_hits",
    "gnrh_neuro_hits",

    "gnrh_migration_core_hits",

    "gnrh_migrating_score",
    "gnrh_mature_score",
    "gnrh_secreting_score",

    "gnrh_knn"
  )

  available_score_columns <- Reduce(
    intersect,
    lapply(
      gnrh_list,
      function(object) {
        colnames(
          object[[]]
        )
      }
    )
  )

  score_columns <- intersect(
    candidate_score_columns,
    available_score_columns
  )

  if (length(score_columns) > 0L) {

    scores <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {

        object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_class = as.character(
              .data$gnrh_class
            )
          ) |>
          dplyr::filter(
            .data$gnrh_class %in%
              positive_classes
          ) |>
          dplyr::group_by(
            .data$gnrh_class
          ) |>
          dplyr::summarise(
            n_cells = dplyr::n(),

            dplyr::across(
              dplyr::all_of(
                score_columns
              ),
              list(
                median = ~ stats::median(
                  .x,
                  na.rm = TRUE
                ),
                q25 = ~ stats::quantile(
                  .x,
                  0.25,
                  na.rm = TRUE,
                  names = FALSE
                ),
                q75 = ~ stats::quantile(
                  .x,
                  0.75,
                  na.rm = TRUE,
                  names = FALSE
                )
              )
            ),

            .groups = "drop"
          ) |>
          dplyr::mutate(
            id = id
          )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$gnrh_class,
        dplyr::everything()
      )

  } else {

    scores <- tibble::tibble()
  }


  # =========================================================================== #
  # Stage distribution
  # =========================================================================== #

  log("Validating developmental stages")

  stages <- purrr::imap_dfr(
    gnrh_list,
    function(object, id) {

      x <- object[[]] |>
        tibble::as_tibble() |>
        dplyr::mutate(
          gnrh_class = as.character(
            .data$gnrh_class
          ),
          gnrh_stage = as.character(
            .data$gnrh_stage
          )
        ) |>
        dplyr::filter(
          .data$gnrh_class %in%
            positive_classes,
          !is.na(.data$gnrh_stage),
          .data$gnrh_stage != "non-gnrh"
        )

      if (nrow(x) == 0L) {
        return(
          tibble::tibble()
        )
      }

      x |>
        dplyr::count(
          .data$gnrh_stage,
          name = "n_cells"
        ) |>
        dplyr::mutate(
          id = id,
          pct_positive =
            100 *
            .data$n_cells /
            sum(.data$n_cells)
        )
    }
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      .data$gnrh_stage,
      .data$n_cells,
      .data$pct_positive
    ) |>
    dplyr::arrange(
      .data$id,
      dplyr::desc(
        .data$n_cells
      )
    )


  # =========================================================================== #
  # Stage by GnRH detection class
  # =========================================================================== #

  stage_class <- purrr::imap_dfr(
    gnrh_list,
    function(object, id) {

      x <- object[[]] |>
        tibble::as_tibble() |>
        dplyr::mutate(
          gnrh_class = as.character(
            .data$gnrh_class
          ),
          gnrh_stage = as.character(
            .data$gnrh_stage
          )
        ) |>
        dplyr::filter(
          .data$gnrh_class %in%
            positive_classes,
          !is.na(.data$gnrh_stage),
          .data$gnrh_stage != "non-gnrh"
        )

      if (nrow(x) == 0L) {
        return(
          tibble::tibble()
        )
      }

      x |>
        dplyr::count(
          .data$gnrh_class,
          .data$gnrh_stage,
          name = "n_cells"
        ) |>
        dplyr::group_by(
          .data$gnrh_class
        ) |>
        dplyr::mutate(
          pct_class =
            100 *
            .data$n_cells /
            sum(.data$n_cells)
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          id = id
        )
    }
  ) |>
    dplyr::left_join(
      dataset_metadata,
      by = "id"
    ) |>
    dplyr::select(
      .data$id,
      .data$label,
      .data$species,
      .data$gnrh_class,
      .data$gnrh_stage,
      .data$n_cells,
      .data$pct_class
    )


  # =========================================================================== #
  # Stage refinement summary
  # =========================================================================== #

  refinement_required <- c(
    "gnrh_stage_raw",
    "gnrh_stage_reassigned",
    "gnrh_stage_reason"
  )

  has_refinement <- vapply(
    gnrh_list,
    function(object) {

      all(
        refinement_required %in%
          colnames(
            object[[]]
          )
      )
    },
    logical(1)
  )

  if (any(has_refinement)) {

    stage_refinement <- purrr::imap_dfr(
      gnrh_list[has_refinement],
      function(object, id) {

        md <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::filter(
            .data$gnrh_status == "pos"
          )

        n_positive <- nrow(md)

        n_reassigned <- sum(
          md$gnrh_stage_reassigned,
          na.rm = TRUE
        )

        n_migration_filtered <- sum(
          as.character(
            md$gnrh_stage_reason
          ) == "migration_core_filter",
          na.rm = TRUE
        )

        tibble::tibble(
          id = id,

          n_positive =
            n_positive,

          n_reassigned =
            n_reassigned,

          pct_reassigned =
            if (n_positive > 0L) {
              100 *
                n_reassigned /
                n_positive
            } else {
              NA_real_
            },

          n_migration_filtered =
            n_migration_filtered
        )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        dplyr::everything()
      )

    stage_reassignment <- purrr::imap_dfr(
      gnrh_list[has_refinement],
      function(object, id) {

        object[[]] |>
          tibble::as_tibble() |>
          dplyr::filter(
            .data$gnrh_status == "pos"
          ) |>
          dplyr::mutate(
            gnrh_stage_raw = as.character(
              .data$gnrh_stage_raw
            ),
            gnrh_stage = as.character(
              .data$gnrh_stage
            ),
            gnrh_stage_reason = as.character(
              .data$gnrh_stage_reason
            )
          ) |>
          dplyr::count(
            .data$gnrh_stage_raw,
            .data$gnrh_stage,
            .data$gnrh_stage_reason,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id
          )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$gnrh_stage_raw,
        .data$gnrh_stage,
        .data$gnrh_stage_reason,
        .data$n_cells
      )

  } else {

    stage_refinement <- tibble::tibble()
    stage_reassignment <- tibble::tibble()
  }


  # =========================================================================== #
  # Migration-core validation
  # =========================================================================== #

  has_migration_hits <- vapply(
    gnrh_list,
    function(object) {

      all(
        c(
          "gnrh_migration_core_hits",
          "gnrh_stage_raw"
        ) %in%
          colnames(
            object[[]]
          )
      )
    },
    logical(1)
  )

  if (any(has_migration_hits)) {

    migration_core <- purrr::imap_dfr(
      gnrh_list[has_migration_hits],
      function(object, id) {

        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_stage_raw = as.character(
              .data$gnrh_stage_raw
            )
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(
              .data$gnrh_stage_raw
            )
          )

        if (nrow(x) == 0L) {
          return(
            tibble::tibble()
          )
        }

        x |>
          dplyr::group_by(
            .data$gnrh_stage_raw
          ) |>
          dplyr::summarise(
            n = dplyr::n(),

            median_hits =
              stats::median(
                .data$gnrh_migration_core_hits,
                na.rm = TRUE
              ),

            q25_hits =
              stats::quantile(
                .data$gnrh_migration_core_hits,
                0.25,
                na.rm = TRUE,
                names = FALSE
              ),

            q75_hits =
              stats::quantile(
                .data$gnrh_migration_core_hits,
                0.75,
                na.rm = TRUE,
                names = FALSE
              ),

            pct_ge1 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 1,
                na.rm = TRUE
              ),

            pct_ge2 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 2,
                na.rm = TRUE
              ),

            pct_ge3 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 3,
                na.rm = TRUE
              ),

            .groups = "drop"
          ) |>
          dplyr::mutate(
            id = id
          )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$gnrh_stage_raw,
        .data$n,
        .data$median_hits,
        .data$q25_hits,
        .data$q75_hits,
        .data$pct_ge1,
        .data$pct_ge2,
        .data$pct_ge3
      )

  } else {

    migration_core <- tibble::tibble()
  }


  # =========================================================================== #
  # Migration refinement validation
  # =========================================================================== #

  has_migration_refinement <- vapply(
    gnrh_list,
    function(object) {

      all(
        c(
          "gnrh_stage_raw",
          "gnrh_stage",
          "gnrh_stage_reason",
          "gnrh_migration_core_hits"
        ) %in%
          colnames(
            object[[]]
          )
      )
    },
    logical(1)
  )

  if (any(has_migration_refinement)) {

    migration_refinement <- purrr::imap_dfr(
      gnrh_list[has_migration_refinement],
      function(object, id) {

        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_stage_raw = as.character(
              .data$gnrh_stage_raw
            ),
            gnrh_stage = as.character(
              .data$gnrh_stage
            ),
            gnrh_stage_reason = as.character(
              .data$gnrh_stage_reason
            )
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            .data$gnrh_stage_raw == "migrating"
          ) |>
          dplyr::mutate(
            migration_outcome = dplyr::case_when(
              .data$gnrh_stage == "migrating" ~
                "retained_migrating",

              .data$gnrh_stage != "migrating" ~
                "reassigned",

              TRUE ~
                NA_character_
            )
          )

        if (nrow(x) == 0L) {
          return(
            tibble::tibble()
          )
        }

        x |>
          dplyr::group_by(
            .data$migration_outcome,
            .data$gnrh_stage
          ) |>
          dplyr::summarise(
            n_cells = dplyr::n(),

            median_hits =
              stats::median(
                .data$gnrh_migration_core_hits,
                na.rm = TRUE
              ),

            q25_hits =
              stats::quantile(
                .data$gnrh_migration_core_hits,
                0.25,
                na.rm = TRUE,
                names = FALSE
              ),

            q75_hits =
              stats::quantile(
                .data$gnrh_migration_core_hits,
                0.75,
                na.rm = TRUE,
                names = FALSE
              ),

            mean_hits =
              mean(
                .data$gnrh_migration_core_hits,
                na.rm = TRUE
              ),

            pct_zero =
              100 *
              mean(
                .data$gnrh_migration_core_hits == 0,
                na.rm = TRUE
              ),

            pct_ge1 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 1,
                na.rm = TRUE
              ),

            pct_ge2 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 2,
                na.rm = TRUE
              ),

            pct_ge3 =
              100 *
              mean(
                .data$gnrh_migration_core_hits >= 3,
                na.rm = TRUE
              ),

            .groups = "drop"
          ) |>
          dplyr::mutate(
            id = id
          )
      }
    ) |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$migration_outcome,
        .data$gnrh_stage,
        .data$n_cells,
        .data$median_hits,
        .data$q25_hits,
        .data$q75_hits,
        .data$mean_hits,
        .data$pct_zero,
        .data$pct_ge1,
        .data$pct_ge2,
        .data$pct_ge3
      ) |>
      dplyr::arrange(
        .data$id,
        .data$migration_outcome,
        .data$gnrh_stage
      )


    # ------------------------------------------------------------------------- #
    # Dataset-level migration refinement summary
    # ------------------------------------------------------------------------- #

    migration_refinement_summary <- migration_refinement |>
      dplyr::group_by(
        .data$id,
        .data$label,
        .data$species
      ) |>
      dplyr::summarise(
        mean_hits_retained = {
          idx <- .data$migration_outcome == "retained_migrating"

          if (any(idx)) {
            stats::weighted.mean(
              .data$mean_hits[idx],
              .data$n_cells[idx],
              na.rm = TRUE
            )
          } else {
            NA_real_
          }
        },

        n_raw_migrating = sum(
          .data$n_cells,
          na.rm = TRUE
        ),

        n_retained = sum(
          .data$n_cells[
            .data$migration_outcome == "retained_migrating"
          ],
          na.rm = TRUE
        ),

        n_reassigned = sum(
          .data$n_cells[
            .data$migration_outcome == "reassigned"
          ],
          na.rm = TRUE
        ),

        n_to_identity = sum(
          .data$n_cells[
            .data$migration_outcome == "reassigned" &
              .data$gnrh_stage == "identity"
          ],
          na.rm = TRUE
        ),

        n_to_mature = sum(
          .data$n_cells[
            .data$migration_outcome == "reassigned" &
              .data$gnrh_stage == "mature"
          ],
          na.rm = TRUE
        ),

        n_to_secreting = sum(
          .data$n_cells[
            .data$migration_outcome == "reassigned" &
              .data$gnrh_stage == "secreting"
          ],
          na.rm = TRUE
        ),

        pct_retained =
          if (n_raw_migrating > 0L) {
            100 * n_retained / n_raw_migrating
          } else {
            NA_real_
          },

        pct_reassigned =
          if (n_raw_migrating > 0L) {
            100 * n_reassigned / n_raw_migrating
          } else {
            NA_real_
          },

        pct_to_identity =
          if (n_reassigned > 0L) {
            100 * n_to_identity / n_reassigned
          } else {
            NA_real_
          },

        pct_to_mature =
          if (n_reassigned > 0L) {
            100 * n_to_mature / n_reassigned
          } else {
            NA_real_
          },

        pct_to_secreting =
          if (n_reassigned > 0L) {
            100 * n_to_secreting / n_reassigned
          } else {
            NA_real_
          },

        .groups = "drop"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$n_raw_migrating,
        .data$n_retained,
        .data$n_reassigned,
        .data$pct_retained,
        .data$pct_reassigned,
        .data$n_to_identity,
        .data$n_to_mature,
        .data$n_to_secreting,
        .data$pct_to_identity,
        .data$pct_to_mature,
        .data$pct_to_secreting,
        .data$mean_hits_retained
      ) |>
      dplyr::arrange(
        dplyr::desc(
          .data$pct_reassigned
        )
      )

  } else {

    migration_refinement <- tibble::tibble()
    migration_refinement_summary <- tibble::tibble()
  }


  # =========================================================================== #
  # Independent biological marker validation
  # =========================================================================== #

  log("Summarizing biological marker expression")

  summarise_marker_expression <- function(
    object,
    id
  ) {

    md <- object[[]]

    cells <- rownames(md)[
      as.character(
        md$gnrh_class
      ) %in%
        positive_classes
    ]

    if (length(cells) == 0L) {
      return(
        tibble::tibble()
      )
    }

    if (!assay %in% names(object@assays)) {

      warning(
        "Assay `",
        assay,
        "` not found in dataset `",
        id,
        "`; skipping biological marker validation.",
        call. = FALSE
      )

      return(
        tibble::tibble()
      )
    }

    genes <- intersect(
      validation_markers,
      rownames(
        object[[assay]]
      )
    )

    if (length(genes) == 0L) {
      return(
        tibble::tibble()
      )
    }

    mat <- tryCatch(
      SeuratObject::LayerData(
        object = object,
        assay = assay,
        layer = layer
      ),
      error = function(e) {

        warning(
          "Unable to retrieve layer `",
          layer,
          "` from assay `",
          assay,
          "` in dataset `",
          id,
          "`: ",
          conditionMessage(e),
          call. = FALSE
        )

        NULL
      }
    )

    if (is.null(mat)) {
      return(
        tibble::tibble()
      )
    }

    available_cells <- intersect(
      cells,
      colnames(mat)
    )

    if (length(available_cells) == 0L) {
      return(
        tibble::tibble()
      )
    }

    available_genes <- intersect(
      genes,
      rownames(mat)
    )

    if (length(available_genes) == 0L) {
      return(
        tibble::tibble()
      )
    }

    mat <- mat[
      available_genes,
      available_cells,
      drop = FALSE
    ]

    classes <- as.character(
      md[
        available_cells,
        "gnrh_class",
        drop = TRUE
      ]
    )

    purrr::map_dfr(
      positive_classes,
      function(cl) {

        idx <- which(
          classes == cl
        )

        if (length(idx) == 0L) {
          return(
            tibble::tibble()
          )
        }

        x <- mat[
          ,
          idx,
          drop = FALSE
        ]

        tibble::tibble(
          id = id,

          gnrh_class = cl,

          gene =
            rownames(x),

          avg_expression =
            Matrix::rowMeans(x),

          pct_expressing =
            100 *
            Matrix::rowSums(
              x > 0
            ) /
            ncol(x)
        )
      }
    )
  }

  biological_markers <- purrr::imap_dfr(
    gnrh_list,
    summarise_marker_expression
  )

  if (nrow(biological_markers) > 0L) {

    biological_markers <- biological_markers |>
      dplyr::left_join(
        dataset_metadata,
        by = "id"
      ) |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$gnrh_class,
        .data$gene,
        .data$avg_expression,
        .data$pct_expressing
      )

  } else {

    biological_markers <- tibble::tibble()
  }


  # =========================================================================== #
  # Export
  # =========================================================================== #

  validation_tables <- list(
    input_summary = input_summary,
    detection = detection,
    dropout_candidates = dropout_candidates,
    classification_consistency = classification_consistency,
    status_class = status_class,
    scores = scores,
    stages = stages,
    stage_class = stage_class,
    stage_refinement = stage_refinement,
    stage_reassignment = stage_reassignment,
    migration_core = migration_core,
    migration_refinement = migration_refinement,
    migration_refinement_summary = migration_refinement_summary,
    biological_markers = biological_markers
  )

  if (isTRUE(write_output)) {

    if (
      is.null(output_dir) ||
      length(output_dir) != 1L ||
      is.na(output_dir) ||
      !nzchar(output_dir)
    ) {
      stop(
        "`output_dir` must be a valid directory when ",
        "`write_output = TRUE`.",
        call. = FALSE
      )
    }

    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    purrr::iwalk(
      validation_tables,
      function(x, name) {

        if (
          is.data.frame(x) &&
          ncol(x) > 0L
        ) {

          readr::write_csv(
            x,
            file.path(
              output_dir,
              paste0(
                name,
                ".csv"
              )
            )
          )
        }
      }
    )
  }


  # =========================================================================== #
  # Final summary
  # =========================================================================== #

  total_cells <- sum(
    detection$n_cells,
    na.rm = TRUE
  )

  total_positive <- sum(
    detection$gnrh_pos,
    na.rm = TRUE
  )

  total_dropout_candidates <- sum(
    detection$dropout_candidates,
    na.rm = TRUE
  )

  log(
    "Validated ",
    length(gnrh_list),
    " datasets; ",
    format(
      total_cells,
      big.mark = ","
    ),
    " cells; ",
    format(
      total_positive,
      big.mark = ","
    ),
    " GnRH-positive cells."
  )

  log(
    "GNRH1-negative transcriptomic candidates: ",
    format(
      total_dropout_candidates,
      big.mark = ","
    ),
    " (diagnostic only; not counted as GnRH-positive)."
  )

  if (nrow(stage_refinement) > 0L) {

    total_reassigned <- sum(
      stage_refinement$n_reassigned,
      na.rm = TRUE
    )

    log(
      "Developmental-stage refinements: ",
      format(
        total_reassigned,
        big.mark = ","
      ),
      " cells."
    )
  }

  if (nrow(migration_refinement_summary) > 0L) {

    total_raw_migrating <- sum(
      migration_refinement_summary$n_raw_migrating,
      na.rm = TRUE
    )

    total_migration_reassigned <- sum(
      migration_refinement_summary$n_reassigned,
      na.rm = TRUE
    )

    pct_migration_reassigned <-
      if (total_raw_migrating > 0L) {
        100 *
          total_migration_reassigned /
          total_raw_migrating
      } else {
        NA_real_
      }

    log(
      "Migration refinement: ",
      format(
        total_migration_reassigned,
        big.mark = ","
      ),
      " / ",
      format(
        total_raw_migrating,
        big.mark = ","
      ),
      " raw migrating cells reassigned (",
      sprintf(
        "%.2f",
        pct_migration_reassigned
      ),
      "%)."
    )
  }

  log(
    "==== GNRH COLLECTION VALIDATION DONE ===="
  )


  # =========================================================================== #
  # Return
  # =========================================================================== #

  result <- c(
    list(
      datasets =
        dataset_metadata
    ),
    validation_tables,
    list(
      output_dir =
        if (isTRUE(write_output)) {
          normalizePath(
            output_dir,
            mustWork = FALSE
          )
        } else {
          NULL
        }
    )
  )

  structure(
    result,
    class = "gnrh_validation"
  )
}
