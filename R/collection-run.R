  #' Run GnRHcell across multiple datasets
  #'
  #' @param datasets Dataset configuration table.
  #' @param output_dir Root output directory.
  #' @param run_markers Run marker discovery.
  #' @param run_comparisons Run cross-dataset comparisons.
  #' @param run_programs Run marker-program analysis.
  #' @param run_figures Generate minimal status/stage figures.
  #' @param run_report Generate dataset-level diagnostic reports.
  #' @param clean_objects Remove processed Seurat objects from returned results.
  #' @param save_objects Save processed Seurat objects.
  #' @param detect_args Arguments passed to \code{detect_gnrh()}.
  #' @param stage_args Arguments passed to \code{stage_gnrh()}.
  #' @param diagnostic_args Arguments passed to \code{gnrh_diagnostics()}.
  #' @param marker_args Arguments passed to \code{gnrh_markers()}.
  #' @param comparison_args Arguments passed to \code{compare_gnrh_datasets()}.
  #' @param verbose Print progress.
  #'
  #' @return Object of class \code{"gnrh_collection"}.
  #'
  #' @export
  run_gnrh_collection <- function(
    datasets,
    output_dir,
    run_markers = TRUE,
    run_comparisons = TRUE,
    run_programs = TRUE,
    run_figures = TRUE,
    run_report = TRUE,
    clean_objects = TRUE,
    save_objects = FALSE,
    detect_args = list(),
    stage_args = list(),
    diagnostic_args = list(),
    marker_args = list(),
    comparison_args = list(),
    verbose = TRUE
  ) {
    # ========================================================================= #
    # Dataset configuration
    # ========================================================================= #
    datasets <- prepare_gnrh_datasets(
      datasets,
      check_files = TRUE,
      remove_missing = TRUE
    )
    if (!nrow(datasets)) {
      stop("No available dataset files.", call. = FALSE)
    }
    datasets$id <- as.character(datasets$id)
    datasets$file <- as.character(datasets$file)
    if (anyNA(datasets$id) || any(!nzchar(datasets$id))) {
      stop("`datasets$id` cannot contain missing or empty values.", call. = FALSE)
    }
    if (anyDuplicated(datasets$id)) {
      stop("`datasets$id` must contain unique dataset identifiers.", call. = FALSE)
    }
    # ========================================================================= #
    # Output directories
    # ========================================================================= #
    output_dir <- path.expand(output_dir)
    dirs <- c(
      output_dir,
      file.path(output_dir, "figures"),
      file.path(output_dir, "tables"),
      file.path(output_dir, "markers"),
      file.path(output_dir, "comparisons"),
      file.path(output_dir, "objects")
    )
    invisible(lapply(
      dirs,
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    ))
    # ========================================================================= #
    # Initialize results
    # ========================================================================= #
    results <- stats::setNames(
      vector("list", nrow(datasets)),
      datasets$id
    )
    log <- .msg(verbose)
    # ========================================================================= #
    # Process datasets
    # ========================================================================= #
    for (i in seq_len(nrow(datasets))) {
      info <- datasets[i, , drop = FALSE]
      id <- info$id[[1L]]
      label <- info$label[[1L]]
      file <- info$file[[1L]]
      split_by <- info$split_by[[1L]]
      reduction <- info$reduction[[1L]]
      log(
        sprintf("[%d/%d] %s", i, nrow(datasets), label),
        type = "step"
      )
      object <- readRDS(file)
      if (!inherits(object, "Seurat")) {
        stop(
          "Dataset `", id, "` is not a Seurat object.",
          call. = FALSE
        )
      }
      result <- run_gnrh_dataset(
        object = object,
        dataset_id = id,
        dataset_label = label,
        split_by = split_by,
        reduction = reduction,
        output_dir = output_dir,
        run_markers = run_markers,
        run_figures = run_figures,
        run_report = run_report,
        detect_args = detect_args,
        stage_args = stage_args,
        diagnostic_args = diagnostic_args,
        marker_args = marker_args,
        clean_object = FALSE,
        verbose = verbose
      )
      # ----------------------------------------------------------------------- #
      # Save processed object
      # ----------------------------------------------------------------------- #
      if (isTRUE(save_objects)) {
        saveRDS(
          result$object,
          file.path(
            output_dir,
            "objects",
            paste0(id, "_gnrh.rds")
          )
        )
      }
      # ----------------------------------------------------------------------- #
      # Store result
      # ----------------------------------------------------------------------- #
      results[[id]] <- result
      rm(object)
      invisible(gc())
    }
    # ========================================================================= #
    # Cross-dataset comparisons
    # ========================================================================= #
    comparisons <- NULL
    if (isTRUE(run_comparisons)) {
      # ----------------------------------------------------------------------- #
      # Collect processed objects
      # ----------------------------------------------------------------------- #
      comparison_objects <- lapply(results, `[[`, "object")
      keep <- !vapply(comparison_objects, is.null, logical(1))
      comparison_objects <- comparison_objects[keep]
      if (length(comparison_objects)) {
        names(comparison_objects) <- names(results)[keep]
      }
      if (anyDuplicated(names(comparison_objects))) {
        stop("Dataset IDs must be unique for cross-dataset comparisons.", call. = FALSE)
      }
      # ----------------------------------------------------------------------- #
      # Comparison arguments
      # ----------------------------------------------------------------------- #
      default_comparison_args <- list(
        datasets = datasets,
        objects = comparison_objects,
        output_dir = output_dir,
        run_programs = run_programs,
        run_gallery = isTRUE(run_figures)
      )
      comparison_args <- utils::modifyList(
        default_comparison_args,
        comparison_args
      )
      # ----------------------------------------------------------------------- #
      # Run comparisons
      # ----------------------------------------------------------------------- #
      comparisons <- do.call(
        compare_gnrh_datasets,
        comparison_args
      )
    }
    # ========================================================================= #
    # Remove heavy Seurat objects
    # ========================================================================= #
    if (isTRUE(clean_objects)) {
      results <- lapply(
        results,
        function(x) {
          x$object <- NULL
          x
        }
      )
    }
    # ========================================================================= #
    # Return
    # ========================================================================= #
    structure(
      list(
        datasets = datasets,
        results = results,
        comparisons = comparisons,
        parameters = list(
          detect_args = detect_args,
          stage_args = stage_args,
          diagnostic_args = diagnostic_args,
          marker_args = marker_args,
          run_markers = run_markers,
          run_comparisons = run_comparisons,
          run_programs = run_programs,
          run_figures = run_figures,
          run_report = run_report
        ),
        output_dir = normalizePath(
          output_dir,
          mustWork = FALSE
        )
      ),
      class = "gnrh_collection"
    )
  }



  # ========================================================================= #
