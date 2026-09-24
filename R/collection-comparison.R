  # Named file-vector validation
  # ========================================================================= #
  .validate_named_files <- function(x, name = "files") {
    if (!is.character(x)) {
      stop("`", name, "` must be a character vector.", call. = FALSE)
    }
    if (!length(x)) {
      return(invisible(TRUE))
    }
    if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
      stop("`", name, "` must be a named character vector with unique names.", call. = FALSE)
    }
    invisible(TRUE)
  }



  # ========================================================================= #
  # Compare GnRHcell results across datasets
  # ========================================================================= #
  #' Compare GnRHcell results across datasets
  #'
  #' @param datasets Dataset configuration table.
  #' @param output_dir GnRHcell output directory.
  #' @param run_programs Run marker-program integration.
  #' @param run_marker_heatmaps Generate conserved and dataset-specific heatmaps.
  #' @param marker_score_col Score used for quantitative cross-dataset marker comparisons.
  #' @param specific_max_datasets Maximum number of marker lists containing a
  #'   gene for it to be called dataset-specific. Default \code{1L} matches the
  #'   strict unique sets reported by the UpSet analysis.
  #' @param run_gallery Build cross-dataset gallery.
  #' @param objects Optional processed Seurat objects.
  #' @param comparison_dir Optional directory for comparison outputs. When
  #'   `NULL`, defaults to `file.path(output_dir, "comparisons")`.
  #' @param show_titles Display main titles on generated comparison figures.
  #'
  #' @return Named cross-dataset comparison list.
  #'
  #' @export
  compare_gnrh_datasets <- function(
    datasets,
    output_dir,
    run_programs = TRUE,
    run_marker_heatmaps = TRUE,
    marker_score_col = NULL,
    specific_max_datasets = 1L,
    run_gallery = FALSE,
    objects = NULL,
    comparison_dir = NULL,
    show_titles = TRUE
  ) {
    # ========================================================================= #
    # Validation
    # ========================================================================= #
    if (!is.data.frame(datasets)) {
      stop("`datasets` must be a data frame or tibble.", call. = FALSE)
    }
    required <- c("id", "label")
    missing <- setdiff(required, colnames(datasets))
    if (length(missing)) {
      stop("Missing dataset columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    datasets$id <- trimws(as.character(datasets$id))
    datasets$label <- trimws(as.character(datasets$label))
    if (anyNA(datasets$id) || any(!nzchar(datasets$id))) {
      stop("`datasets$id` cannot contain missing or empty values.", call. = FALSE)
    }
    if (anyDuplicated(datasets$id)) {
      stop("`datasets$id` must contain unique values.", call. = FALSE)
    }
    if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles)) {
      stop("`show_titles` must be TRUE or FALSE.", call. = FALSE)
    }
    specific_max_datasets <- as.integer(specific_max_datasets)
    if (length(specific_max_datasets) != 1L || is.na(specific_max_datasets) || specific_max_datasets < 1L) {
      stop("`specific_max_datasets` must be a positive integer.", call. = FALSE)
    }
    # ========================================================================= #
    # Directories
    # ========================================================================= #
    output_dir <- path.expand(output_dir)
    table_dir <- file.path(output_dir, "tables")
    marker_dir <- file.path(output_dir, "markers")
    if (is.null(comparison_dir)) comparison_dir <- file.path(output_dir,"comparisons")
    comparison_dir <- path.expand(comparison_dir)
    dir.create(comparison_dir,recursive=TRUE,showWarnings=FALSE)
    # ========================================================================= #
    # Run-info files
    # ========================================================================= #
    run_info_files <- stats::setNames(
      file.path(table_dir, paste0(datasets$id, "_gnrh_run_info.tsv")),
      datasets$id
    )
    run_info_files <- run_info_files[file.exists(unname(run_info_files))]
    .validate_named_files(run_info_files, "run_info_files")
    # ========================================================================= #
    # Marker files
    # ========================================================================= #
    marker_files <- stats::setNames(
      file.path(marker_dir, paste0("gnrh_", datasets$id, "_markers.tsv")),
      datasets$id
    )
    marker_files <- marker_files[file.exists(unname(marker_files))]
    .validate_named_files(marker_files, "marker_files")
    # ========================================================================= #
    # Marker score
    # ========================================================================= #
    if (length(marker_files)) {
      marker_cols <- lapply(marker_files,function(x) names(utils::read.delim(x,nrows=1,check.names=FALSE)))
      common_cols <- Reduce(intersect,marker_cols)
      if (is.null(marker_score_col)) {
        score_candidates <- c("marker_score","final_score","avg_log2FC","avg_logFC","specificity_score","specificity")
        marker_score_col <- intersect(score_candidates,common_cols)[1L]
      }
      if (!length(marker_score_col) || is.na(marker_score_col)) {
        stop(
          "No common marker score column found. Common columns: ",
          paste(common_cols,collapse=", "),
          call.=FALSE
        )
      }
      if (!marker_score_col %in% common_cols) {
        stop(
          "`marker_score_col = ",marker_score_col,"` is not present in all marker tables. ",
          "Common columns: ",paste(common_cols,collapse=", "),
          call.=FALSE
        )
      }
      message("[INFO] Cross-dataset marker score: ",marker_score_col)
    }
    # ========================================================================= #
    # Initialize
    # ========================================================================= #
    runtime_plot <- NULL
    detected_plot <- NULL
    overlap <- NULL
    programs <- NULL
    conserved_markers <- NULL
    dataset_specific_markers <- NULL
    gallery <- NULL
    # ========================================================================= #
    # Runtime / detection
    # ========================================================================= #
    if (length(run_info_files)) {
      runtime_plot <- gnrh_runtime(
        files = run_info_files,
        show_points = TRUE,
        txtsize = getOption("gnrhcell.base_size", 14),
        x.ang = 60
      )
      detected_plot <- gnrh_cellmark(
        files = run_info_files,
        txtsize = 8
      )
      if (!isTRUE(show_titles)) {
        runtime_plot <- runtime_plot + ggplot2::labs(title = NULL, subtitle = NULL)
        detected_plot <- detected_plot + ggplot2::labs(title = NULL, subtitle = NULL)
      }
      .gnrh_save_plot(
        runtime_plot,
        file.path(comparison_dir, "gnrh_runtime_across_datasets"),
        width = 7,
        height = 6
      )
      .gnrh_save_plot(
        detected_plot,
        file.path(comparison_dir, "gnrh_detected_across_datasets"),
        width = 7,
        height = 6
      )
    } else {
      warning("No GnRHcell run-info files were found.", call. = FALSE)
    }
    # ========================================================================= #
    # Marker overlap
    # ========================================================================= #
    if (length(marker_files) >= 2L) {
      marker_basenames <- stats::setNames(
        basename(unname(marker_files)),
        names(marker_files)
      )
      gene_sets <- build_gene_sets(
        files = marker_basenames,
        dir = marker_dir
      )
      if (is.null(names(gene_sets)) || anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
        names(gene_sets) <- names(marker_files)
      }
      overlap_dir <- file.path(
        comparison_dir,
        "marker_overlap"
      )
      dir.create(
        overlap_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      overlap <- gnrh_upset(
        gene_sets = gene_sets,
        venn_title = if (isTRUE(show_titles)) "Overlap of GnRH-associated markers" else NULL,
        outdir = overlap_dir
      )
    } else {
      warning(
        "Fewer than two marker tables are available; cross-dataset marker analyses were skipped.",
        call. = FALSE
      )
    }
    # ========================================================================= #
    # Marker programs
    # ========================================================================= #
    program_plot <- NULL
    program_high_plot <- NULL
    if (isTRUE(run_programs) && !is.null(overlap) && length(marker_files)>=2L) {
      program_dir <- file.path(comparison_dir,"marker_programs")
      dir.create(program_dir,recursive=TRUE,showWarnings=FALSE)
      programs <- find_gnrh_programs(
        files=marker_files,
        results=overlap,
        outdir=program_dir,
        write_output=TRUE
      )
      # ----------------------------------------------------------------------- #
      # Adaptive dimensions
      # ----------------------------------------------------------------------- #
      n_program_datasets <- length(marker_files)
      n_program_genes <- if (
        !is.null(programs$high_confidence) &&
        is.data.frame(programs$high_confidence) &&
        nrow(programs$high_confidence) &&
        "dataset" %in% colnames(programs$high_confidence)
      ) {
        max(table(programs$high_confidence$dataset))
      } else {
        15L
      }
      program_width <- max(
        8,
        min(
          18,
          5 + 2.2*min(4L,ceiling(sqrt(n_program_datasets)))
        )
      )
      program_height <- max(
        6,
        min(
          18,
          4 + 0.22*n_program_genes*ceiling(n_program_datasets/3)
        )
      )
      summary_height <- max(
        5.5,
        min(
          10,
          4.5 + 0.35*n_program_datasets
        )
      )
      # ----------------------------------------------------------------------- #
      # Summary plot
      # ----------------------------------------------------------------------- #
      if (!is.null(programs$summary) && is.data.frame(programs$summary) && nrow(programs$summary)) {
        program_plot <- tryCatch(
          gnrh_program_plot(
            programs=programs,
            table="summary",
            type="bar",
            min_genes=1L,
            txtsize = getOption("gnrhcell.base_size", 14),
            style="bw"
          ),
          error=function(e) {
            warning(
              "Marker-program summary plot skipped: ",
              conditionMessage(e),
              call.=FALSE
            )
            NULL
          }
        )
        if (!isTRUE(show_titles) && !is.null(program_plot)) {
          program_plot <- program_plot + ggplot2::labs(title = NULL, subtitle = NULL)
        }
        if (!is.null(program_plot)) {
          .gnrh_save_plot(
            program_plot,
            file.path(program_dir,"gnrh_marker_programs"),
            width=program_width,
            height=summary_height
          )
        }
      }
      # ----------------------------------------------------------------------- #
      # High-confidence marker plot
      # ----------------------------------------------------------------------- #
      if (!is.null(programs$high_confidence) && is.data.frame(programs$high_confidence) && nrow(programs$high_confidence)) {
        program_high_plot <- tryCatch(
          gnrh_program_plot(
            programs=programs,
            table="high_confidence",
            type="dot",
            require_coexpr=FALSE,
            top_n=NULL,
            txtsize = getOption("gnrhcell.base_size", 14),
            style="bw"
          ),
          error=function(e) {
            warning(
              "High-confidence marker-program plot skipped: ",
              conditionMessage(e),
              call.=FALSE
            )
            NULL
          }
        )
        if (!isTRUE(show_titles) && !is.null(program_high_plot)) {
          program_high_plot <- program_high_plot + ggplot2::labs(title = NULL, subtitle = NULL)
        }
        if (!is.null(program_high_plot)) {
          .gnrh_save_plot(
            program_high_plot,
            file.path(program_dir,"gnrh_high_confidence_marker_programs"),
            width=program_width,
            height=program_height
          )
        }
      }
    }
    # ========================================================================= #
    # Conserved / dataset-specific markers
    # ========================================================================= #
    if (isTRUE(run_marker_heatmaps) && length(marker_files)>=2L) {
      conserved_pdf <- file.path(comparison_dir,"conserved_gnrh_markers.pdf")
      conserved_png <- file.path(comparison_dir,"conserved_gnrh_markers.png")
      dataset_pdf <- file.path(comparison_dir,"dataset_specific_gnrh_markers.pdf")
      dataset_png <- file.path(comparison_dir,"dataset_specific_gnrh_markers.png")
      message("[INFO] Generating conserved marker heatmap...")
      conserved_markers <- gnrh_conserved(
        files=marker_files,
        score_col=marker_score_col,
        min_datasets=2L,
        top_n=50L,
        exclude_genes=NULL,
        filename=conserved_pdf,
        title=if (isTRUE(show_titles)) "Conserved GnRH markers" else NULL
      )
      gnrh_conserved(
        files=marker_files,
        score_col=marker_score_col,
        min_datasets=2L,
        top_n=50L,
        exclude_genes=NULL,
        filename=conserved_png,
        title=if (isTRUE(show_titles)) "Conserved GnRH markers" else NULL
      )
      message("[INFO] Conserved markers: ",conserved_pdf)
      message("[INFO] Generating dataset-specific marker heatmap...")
      dataset_specific_markers <- gnrh_specific(
        files=marker_files,
        score_col=marker_score_col,
        top_n_per_dataset=20L,
        max_datasets=specific_max_datasets,
        exclude_genes=NULL,
        filename=dataset_pdf,
        title=if (isTRUE(show_titles)) "Dataset-specific GnRH markers" else NULL
      )
      gnrh_specific(
        files=marker_files,
        score_col=marker_score_col,
        top_n_per_dataset=20L,
        max_datasets=specific_max_datasets,
        exclude_genes=NULL,
        filename=dataset_png,
        title=if (isTRUE(show_titles)) "Dataset-specific GnRH markers" else NULL
      )
      message("[INFO] Dataset-specific markers: ",dataset_pdf)
    }
    # ========================================================================= #
    # Gallery
    # ========================================================================= #
    if (isTRUE(run_gallery)) {
      if (is.null(objects)) {
        objects <- .load_collection_gallery_objects(
          datasets,
          output_dir
        )
      }
      if (is.null(objects) || !length(objects)) {
        warning("Gallery generation requires processed objects.", call. = FALSE)
      } else {
        gallery <- .build_collection_gallery(
          objects = objects,
          datasets = datasets,
          output_dir = file.path(comparison_dir, "gallery"),
          show_titles = show_titles
        )
      }
    }
    # ========================================================================= #
    # Return
    # ========================================================================= #
    list(
      run_info_files = run_info_files,
      marker_files = marker_files,
      runtime_plot = runtime_plot,
      detected_plot = detected_plot,
      overlap = overlap,
      programs = programs,
      conserved_markers = conserved_markers,
      dataset_specific_markers = dataset_specific_markers,
      gallery = gallery
    )
  }
