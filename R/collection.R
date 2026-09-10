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


# ========================================================================= #
# Adaptive plot dimensions
# ========================================================================= #
#' Plot dimmensions
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
    width <- max(9,min(24,7+1.1*n_datasets))
    height <- max(5.5,min(12,5+0.35*n_datasets))
  } else {
    width <- max(8,min(20,5+2.4*ceiling(sqrt(n_datasets))))
    height <- max(6,min(20,4+2.2*ceiling(n_datasets/3)))
  }
  c(width=width,height=height)
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
.save_gnrh_plot <- function(
    plot,
    filename,
    width=7,
    height=5,
    dpi=300,
    bg="white"
) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE)
  ext <- tolower(tools::file_ext(filename))
  if (!nzchar(ext)) {
    ggplot2::ggsave(
      filename=paste0(filename,".pdf"),
      plot=plot,
      width=width,
      height=height,
      units="in",
      bg=bg
    )
    ggplot2::ggsave(
      filename=paste0(filename,".png"),
      plot=plot,
      width=width,
      height=height,
      units="in",
      dpi=dpi,
      bg=bg
    )
  } else {
    ggplot2::ggsave(
      filename=filename,
      plot=plot,
      width=width,
      height=height,
      units="in",
      dpi=if (ext %in% c("png","jpg","jpeg","tiff")) dpi else NULL,
      bg=bg
    )
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



#' Run GnRHcell analysis on a single dataset
#'
#' Runs GnRH detection, developmental staging, diagnostics, essential
#' visualization, marker discovery, marker co-expression, and marker-network
#' analysis for one Seurat dataset.
#'
#' @param object A Seurat object.
#' @param dataset_id Short unique dataset identifier.
#' @param dataset_label Human-readable dataset label.
#' @param split_by Metadata column used for grouped summaries.
#' @param reduction Dimensional reduction used for embedding plots.
#' @param output_dir Root output directory.
#' @param run_markers Run marker discovery.
#' @param run_figures Generate essential GnRH figures.
#' @param run_report Generate GnRHcell diagnostic report.
#' @param detect_args Named list passed to \code{\link{detect_gnrh}}.
#' @param stage_args Named list passed to \code{\link{stage_gnrh}}.
#' @param diagnostic_args Named list passed to \code{\link{gnrh_diagnostics}}.
#' @param marker_args Named list passed to \code{\link{gnrh_markers}}.
#' @param clean_object Trigger garbage collection after analysis.
#' @param verbose Print progress messages.
#'
#' @return Named dataset-level result list.
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
    run_figures = TRUE,
    run_report = TRUE,
    detect_args = list(),
    stage_args = list(),
    diagnostic_args = list(),
    marker_args = list(),
    clean_object = TRUE,
    verbose = TRUE
) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  args <- list(
    detect_args = detect_args,
    stage_args = stage_args,
    diagnostic_args = diagnostic_args,
    marker_args = marker_args
  )
  if (any(!vapply(args, is.list, logical(1)))) {
    stop(
      "`detect_args`, `stage_args`, `diagnostic_args`, and `marker_args` must be lists.",
      call. = FALSE
    )
  }
  log <- .msg(verbose)
  log("Running GnRHcell:", dataset_label, type = "header")
  # ========================================================================= #
  # Directories
  # ========================================================================= #
  dataset_dir <- file.path(output_dir, "figures", dataset_id)
  table_dir <- file.path(output_dir, "tables")
  marker_dir <- file.path(output_dir, "markers")
  invisible(lapply(
    c(dataset_dir, table_dir, marker_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  # ========================================================================= #
  # Plotting metadata
  # ========================================================================= #
  reduction <- resolve_reduction(object, reduction)
  split_by <- resolve_split_column(object, split_by)
  # ========================================================================= #
  # Run GnRHcell
  # ========================================================================= #
  object <- GnRHcell::run_gnrh(
    object = object,
    detect_args = detect_args,
    stage_args = stage_args,
    diagnostic_args = diagnostic_args,
    verbose = verbose
  )
  md <- object[[]]
  # ========================================================================= #
  # Run information
  # ========================================================================= #
  run_info <- GnRHcell::extract_gnrh_run_info(
    object,
    dataset_name = dataset_label
  )
  .save_gnrh_table(
    run_info,
    file.path(table_dir, paste0(dataset_id, "_gnrh_run_info"))
  )
  # ========================================================================= #
  # Figures
  # ========================================================================= #
  figures <- list()
  if (isTRUE(run_figures)) {
    # ----------------------------------------------------------------------- #
    # Status embedding
    # ----------------------------------------------------------------------- #
    if ("gnrh_status" %in% colnames(md)) {
      figures$status <- GnRHcell::plot_gnrh_embedding(
        object = object,
        group_by = "gnrh_status",
        reduction = reduction,
        plot.ttl = paste0(dataset_label)
      )
      .save_gnrh_plot(
        figures$status,
        file.path(dataset_dir, paste0(dataset_id, "_gnrh_status")),
        width = 6,
        height = 5
      )
    }
    # ----------------------------------------------------------------------- #
    # Stage embedding
    # ----------------------------------------------------------------------- #
    if ("gnrh_stage" %in% colnames(md)) {
      figures$stage <- GnRHcell::plot_gnrh_embedding(
        object = object,
        group_by = "gnrh_stage",
        reduction = reduction,
        plot.ttl = paste0(dataset_label)
      )
      .save_gnrh_plot(
        figures$stage,
        file.path(dataset_dir, paste0(dataset_id, "_gnrh_stage")),
        width = 6,
        height = 5
      )
    }
    # ----------------------------------------------------------------------- #
    # Core GnRH features
    # ----------------------------------------------------------------------- #
    core_features <- .gnrh_features("core")
    core_features <- core_features[core_features %in% colnames(md)]
    if (length(core_features)) {
      figures$core <- GnRHcell::plot_gnrh_feature(
        object = object,
        features = core_features,
        reduction = reduction,
        ncol = 2
      )
      .save_gnrh_plot(
        figures$core,
        file.path(dataset_dir, paste0(dataset_id, "_gnrh_core_features")),
        width = 9,
        height = 8
      )
    }
    # ----------------------------------------------------------------------- #
    # Developmental staging features
    # ----------------------------------------------------------------------- #
    staging_features <- .gnrh_features("staging")
    staging_features <- staging_features[staging_features %in% colnames(md)]
    if (length(staging_features)) {
      figures$staging <- GnRHcell::plot_gnrh_feature(
        object = object,
        features = staging_features,
        reduction = reduction,
        ncol = min(3L, length(staging_features))
      )
      .save_gnrh_plot(
        figures$staging,
        file.path(dataset_dir, paste0(dataset_id, "_gnrh_staging_features")),
        width = 12,
        height = 4.5
      )
    }
    # ----------------------------------------------------------------------- #
    # Status distribution
    # ----------------------------------------------------------------------- #
    if (!is.null(split_by) && "gnrh_status" %in% colnames(md)) {
      figures$status_distribution <- GnRHcell::plot_gnrh_distribution(
        object,
        group.by = "gnrh_status",
        split.by = split_by,
        proportion = TRUE,
        label = FALSE,
        cols = GnRHcell::gnrh_colors("status")
      )
      .save_gnrh_plot(
        figures$status_distribution,
        file.path(dataset_dir, paste0(dataset_id, "_status_distribution")),
        width = 6,
        height = 4
      )
    }
    # ----------------------------------------------------------------------- #
    # Stage distribution
    # ----------------------------------------------------------------------- #
    if (!is.null(split_by) && "gnrh_stage" %in% colnames(md)) {
      figures$stage_distribution <- GnRHcell::plot_gnrh_distribution(
        object,
        group.by = "gnrh_stage",
        split.by = split_by,
        proportion = TRUE,
        label = FALSE,
        cols = GnRHcell::gnrh_colors("stage")
      )
      .save_gnrh_plot(
        figures$stage_distribution,
        file.path(dataset_dir, paste0(dataset_id, "_stage_distribution")),
        width = 6,
        height = 4
      )
    }
  }
  # ========================================================================= #
  # Diagnostic report
  # ========================================================================= #
  report_plot <- NULL
  if (isTRUE(run_report)) {
    report_plot <- tryCatch(
      GnRHcell::gnrh_report(
        object,
        style = "bw"
      ),
      error = function(e) {
        log(
          "Skipping diagnostic report:",
          conditionMessage(e),
          type = "warn"
        )
        NULL
      }
    )
    if (!is.null(report_plot)) {
      .save_gnrh_plot(
        report_plot,
        file.path(dataset_dir, paste0(dataset_id, "_gnrh_report")),
        width = 10,
        height = 5
      )
    }
  }
  # ========================================================================= #
  # Markers
  # ========================================================================= #
  markers <- NULL
  if (isTRUE(run_markers)) {
    default_marker_args <- list(
      object=object,
      group.by="gnrh_status",
      ident.1="pos",
      ident.2=NULL,
      assay=NULL,
      layer="data",
      methods="wilcox",
      coexpr.method="spearman",
      program_col=NULL,
      association_mode="combined",
      min_pct=0.005,
      min_fc=0.10,
      max_padj=0.05,
      coexpr_min=NULL,
      codetect_min_or=NULL,
      codetect_max_fdr=NULL,
      program_min=NULL,
      min_detect=2L,
      exclude_gnrh=FALSE,
      verbose=verbose
    )
  }
    marker_args <- utils::modifyList(default_marker_args,marker_args)
    markers <- do.call(GnRHcell::gnrh_markers,marker_args)
    .save_gnrh_table(markers,file.path(marker_dir,paste0("gnrh_",dataset_id,"_markers")))
    # ========================================================================= #
    # Marker association plots
    # ========================================================================= #
    if (isTRUE(run_figures) && is.data.frame(markers) && nrow(markers)) {
      figures$coexpression <- tryCatch(
        GnRHcell::plot_gnrh_coexpr(
          markers,
          coexp_cutoff=0.10,
          top_n=12L,
          exclude_gnrh=TRUE,
          txtsize=10,
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH co-expression plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$codetection <- tryCatch(
        GnRHcell::plot_gnrh_codetect(
          markers,
          top_n=10L,
          min_or=2,
          max_fdr=0.05,
          min_specificity=0.05,
          exclude_gnrh=TRUE,
          txtsize=10,
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH co-detection plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$detection <- tryCatch(
        GnRHcell::plot_gnrh_detection(
          markers,
          top_n=10L,
          min_specificity=0.05,
          max_padj=0.05,
          exclude_gnrh=TRUE,
          txtsize=10,
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH phenotype-association plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      if (!is.null(figures$coexpression)) {
        .save_gnrh_plot(
          figures$coexpression,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_coexpression")),
          width=6,
          height=6
        )
      }
      if (!is.null(figures$codetection)) {
        .save_gnrh_plot(
          figures$codetection,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_codetection")),
          width=6,
          height=6
        )
      }
      if (!is.null(figures$detection)) {
        .save_gnrh_plot(
          figures$detection,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_detection")),
          width=6,
          height=6
        )
      }
      association_plots <- Filter(
        Negate(is.null),
        list(
          coexpression=figures$coexpression,
          codetection=figures$codetection,
          detection=figures$detection
        )
      )
      if (length(association_plots)) {
        figures$marker_associations <- patchwork::wrap_plots(
          association_plots,
          nrow=1
        )
        .save_gnrh_plot(
          figures$marker_associations,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_marker_associations")),
          width=5.5*length(association_plots),
          height=5.5
        )
      }
    }
    # ========================================================================= #
    # Marker networks
    # ========================================================================= #
    if (isTRUE(run_figures) && is.data.frame(markers) && nrow(markers)) {
      figures$network_coexpression <- tryCatch(
        GnRHcell::plot_network(
          markers,
          mode="coexpression",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize=10,
          style="void"
        ),
        error=function(e) {
          log("Skipping co-expression network:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$network_codetection <- tryCatch(
        GnRHcell::plot_network(
          markers,
          mode="codetection",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize=10,
          style="void"
        ),
        error=function(e) {
          log("Skipping co-detection network:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$network_phenotype <- tryCatch(
        GnRHcell::plot_network(
          markers,
          mode="phenotype",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize=10,
          style="void"
        ),
        error=function(e) {
          log("Skipping phenotype network:",conditionMessage(e),type="warn")
          NULL
        }
      )
      networks <- Filter(
        Negate(is.null),
        list(
          coexpression=figures$network_coexpression,
          codetection=figures$network_codetection,
          phenotype=figures$network_phenotype
        )
      )
      if (length(networks)) {
        purrr::iwalk(
          networks,
          function(p,nm) {
            .save_gnrh_plot(
              p,
              file.path(dataset_dir,paste0(dataset_id,"_gnrh_network_",nm)),
              width=9,
              height=9
            )
          }
        )
        figures$networks <- patchwork::wrap_plots(networks,nrow=1)
        .save_gnrh_plot(
          figures$networks,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_networks")),
          width=8*length(networks),
          height=8
        )
      }
    }
    # ========================================================================= #
    # Marker network
    # ========================================================================= #
    if (isTRUE(run_figures) && is.data.frame(markers) && nrow(markers)) {
      network_mode <- if (any(c("codetect_log2or","codetect_or") %in% colnames(markers))) {
        "codetection"
      } else if (any(c("phenotype_log2or","phenotype_or","specificity","spec") %in% colnames(markers))) {
        "phenotype"
      } else if (any(c("coexpr_cor","coexpr","coexpression","coexpr_pct") %in% colnames(markers))) {
        "coexpression"
      } else {
        NA_character_
      }
      if (!is.na(network_mode)) {
        figures$network <- tryCatch(
          GnRHcell::plot_network(
            markers,
            top_n=30L,
            threshold=0.10,
            mode=network_mode,
            include_gnrh=TRUE,
            txtsize=10,
            style="void"
          ),
          error=function(e) {
            log("Skipping marker network:",conditionMessage(e),type="warn")
            NULL
          }
        )
        if (!is.null(figures$network)) {
          .save_gnrh_plot(
            figures$network,
            file.path(dataset_dir,paste0(dataset_id,"_gnrh_network_",network_mode)),
            width=9,
            height=9
          )
        }
      } else {
        log("Skipping marker network: no usable association column.",type="warn")
      }
    }
    # ========================================================================= #
    # Cleanup
    # ========================================================================= #
    if (isTRUE(clean_object)) {
      invisible(gc())
    }
    # ========================================================================= #
    # Return
    # ========================================================================= #
    list(
      object = object,
      run_info = run_info,
      markers = markers,
      figures = figures,
      report = report_plot,
      reduction = reduction,
      split_by = split_by
    )
  }




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
  #' @param run_gallery Build cross-dataset gallery.
  #' @param objects Optional processed Seurat objects.
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
    run_gallery = FALSE,
    objects = NULL,
    comparison_dir = NULL
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
      .save_gnrh_plot(
        runtime_plot,
        file.path(comparison_dir, "gnrh_runtime_across_datasets"),
        width = 7,
        height = 5
      )
      .save_gnrh_plot(
        detected_plot,
        file.path(comparison_dir, "gnrh_detected_across_datasets"),
        width = 7,
        height = 5
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
      overlap <- gnrh_gene_upset(
        gene_sets = gene_sets,
        venn_title = "Overlap of GnRH-associated markers",
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
      programs <- gnrh_marker_programs(
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
          plot_gnrh_marker_programs(
            programs=programs,
            table="summary",
            type="bar",
            min_genes=1L,
            txtsize=10,
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
        if (!is.null(program_plot)) {
          .save_gnrh_plot(
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
          plot_gnrh_marker_programs(
            programs=programs,
            table="high_confidence",
            type="dot",
            require_coexpr=FALSE,
            top_n=NULL,
            txtsize=9,
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
        if (!is.null(program_high_plot)) {
          .save_gnrh_plot(
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
      conserved_markers <- plot_gnrh_conserved_markers(
        files=marker_files,
        score_col=marker_score_col,
        min_datasets=2L,
        top_n=50L,
        exclude_genes=NULL,
        filename=conserved_pdf
      )
      plot_gnrh_conserved_markers(
        files=marker_files,
        score_col=marker_score_col,
        min_datasets=2L,
        top_n=50L,
        exclude_genes=NULL,
        filename=conserved_png
      )
      message("[INFO] Conserved markers: ",conserved_pdf)
      message("[INFO] Generating dataset-specific marker heatmap...")
      dataset_specific_markers <- plot_gnrh_dataset_specific_markers(
        files=marker_files,
        score_col=marker_score_col,
        top_n_per_dataset=20L,
        max_datasets=2L,
        exclude_genes=NULL,
        filename=dataset_pdf
      )
      plot_gnrh_dataset_specific_markers(
        files=marker_files,
        score_col=marker_score_col,
        top_n_per_dataset=20L,
        max_datasets=2L,
        exclude_genes=NULL,
        filename=dataset_png
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
          output_dir = file.path(comparison_dir, "gallery")
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


  #' Validate a GnRHcell multi-dataset collection
  #'
  #' Performs descriptive cross-dataset validation of GnRH detection,
  #' developmental staging, transcriptomic rescue, migration refinement,
  #' secretory phenotype, and biological evidence.
  #'
  #' @param collection A \code{"gnrh_collection"}.
  #' @param output_dir Validation output directory.
  #' @param positive_classes GnRH-positive evidence classes. Default is
  #'   \code{c("direct", "transcriptomic")}.
  #' @param validation_markers Optional biological validation genes.
  #' @param assay Assay used for expression summaries.
  #' @param layer Expression layer.
  #' @param allow_saved_objects Load saved processed objects when they are not
  #'   retained in memory.
  #' @param write_output Write validation tables.
  #' @param verbose Print progress.
  #'
  #' @return Object of class \code{"gnrh_validation"}.
  #'
  #' @export
  validate_gnrh_collection <- function(
    collection,
    output_dir = file.path(collection$output_dir, "validation"),
    positive_classes = c("direct", "transcriptomic"),
    validation_markers = NULL,
    assay = "RNA",
    layer = "data",
    allow_saved_objects = TRUE,
    write_output = TRUE,
    verbose = TRUE
  ) {
    # ========================================================================= #
    # Validation
    # ========================================================================= #
    if (!inherits(collection, "gnrh_collection")) {
      stop("`collection` must inherit from `gnrh_collection`.", call. = FALSE)
    }
    if (is.null(collection$results) || !length(collection$results)) {
      stop("`collection$results` is empty.", call. = FALSE)
    }
    positive_classes <- unique(as.character(positive_classes))
    if (!length(positive_classes) || anyNA(positive_classes) || any(!nzchar(positive_classes))) {
      stop("`positive_classes` must contain valid class names.", call. = FALSE)
    }
    log <- .msg(verbose)
    log("GNRH COLLECTION VALIDATION", type = "header")
    # ========================================================================= #
    # Helpers
    # ========================================================================= #
    join_metadata <- function(x, metadata) {
      if (!is.data.frame(x) || !nrow(x) || !"id" %in% colnames(x)) {
        return(x)
      }
      dplyr::left_join(x, metadata, by = "id")
    }
    count_true <- function(md, column) {
      if (!column %in% colnames(md)) {
        return(NA_integer_)
      }
      x <- md[[column]]
      if (is.logical(x)) {
        return(sum(x %in% TRUE, na.rm = TRUE))
      }
      if (is.numeric(x)) {
        return(sum(is.finite(x) & x > 0, na.rm = TRUE))
      }
      x <- tolower(trimws(as.character(x)))
      sum(x %in% c("true", "t", "1", "yes", "y", "pos", "positive"), na.rm = TRUE)
    }
    safe_median <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::median(x)
    }
    safe_quantile <- function(x, p) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::quantile(x, p, names = FALSE)
    }
    # ========================================================================= #
    # Dataset metadata
    # ========================================================================= #
    datasets <- tibble::as_tibble(collection$datasets)
    required <- c("id", "label", "species")
    missing <- setdiff(required, colnames(datasets))
    if (length(missing)) {
      stop("Missing dataset metadata: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    datasets$id <- as.character(datasets$id)
    datasets$label <- as.character(datasets$label)
    datasets$species <- as.character(datasets$species)
    if (anyNA(datasets$id) || any(!nzchar(datasets$id)) || anyDuplicated(datasets$id)) {
      stop("`collection$datasets$id` must contain unique non-empty identifiers.", call. = FALSE)
    }
    dataset_metadata <- datasets |>
      dplyr::select(.data$id, .data$label, .data$species)
    # ========================================================================= #
    # Recover processed objects
    # ========================================================================= #
    gnrh_list <- .get_collection_objects(
      collection,
      allow_saved = allow_saved_objects
    )
    missing_objects <- names(gnrh_list)[
      vapply(gnrh_list, is.null, logical(1))
    ]
    if (length(missing_objects)) {
      stop(
        "Processed objects are unavailable for: ",
        paste(missing_objects, collapse = ", "),
        ". Retain them with `clean_objects = FALSE` or use `save_objects = TRUE`.",
        call. = FALSE
      )
    }
    invalid <- !vapply(
      gnrh_list,
      inherits,
      logical(1),
      what = "Seurat"
    )
    if (any(invalid)) {
      stop(
        "Invalid processed object(s): ",
        paste(names(gnrh_list)[invalid], collapse = ", "),
        call. = FALSE
      )
    }
    # ========================================================================= #
    # Required GnRH metadata
    # ========================================================================= #
    required_gnrh <- c(
      "gnrh_status",
      "gnrh_class",
      "gnrh_stage"
    )
    for (id in names(gnrh_list)) {
      md <- gnrh_list[[id]][[]]
      missing <- setdiff(required_gnrh, colnames(md))
      if (length(missing)) {
        stop(
          "Dataset `", id, "` is missing: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
    }
    observed_classes <- unique(unlist(
      lapply(
        gnrh_list,
        function(object) unique(as.character(object[[]]$gnrh_class))
      ),
      use.names = FALSE
    ))
    observed_classes <- observed_classes[
      !is.na(observed_classes) & nzchar(observed_classes)
    ]
    missing_positive_classes <- setdiff(positive_classes, observed_classes)
    if (length(missing_positive_classes)) {
      warning(
        "Positive class(es) not observed in this collection: ",
        paste(missing_positive_classes, collapse = ", "),
        ". Observed classes: ",
        paste(sort(observed_classes), collapse = ", "),
        call. = FALSE
      )
    }
    # ========================================================================= #
    # Validation marker panel
    # ========================================================================= #
    if (is.null(validation_markers)) {
      validation_markers <- unique(c(
        "GNRH1",
        "FEZF1","ISL1","OTX2","SIX3","SIX6","ECEL1",
        "ANOS1","PROKR2","NSMF","ROBO3","SEMA3C","SEMA3F","CXCR4",
        "KISS1R","GNRHR","DOC2B","PTPRN","HCN1",
        "PCSK1","PCSK2","CHGA","CHGB","CPE","SCG2","SCG5","VGF"
      ))
    }
    validation_markers <- unique(as.character(validation_markers))
    # ========================================================================= #
    # Input summary
    # ========================================================================= #
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
      dplyr::left_join(dataset_metadata, by = "id") |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$n_cells,
        .data$n_features
      )
    # ========================================================================= #
    # Detection
    # ========================================================================= #
    detection <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        class <- as.character(md$gnrh_class)
        status <- as.character(md$gnrh_status)
        total <- nrow(md)
        direct <- sum(class == "direct", na.rm = TRUE)
        transcriptomic <- sum(class == "transcriptomic", na.rm = TRUE)
        supported <- sum(class == "supported", na.rm = TRUE)
        positive_class <- sum(class %in% positive_classes, na.rm = TRUE)
        positive_status <- sum(status == "pos", na.rm = TRUE)
        direct_signal <- count_true(md, "gnrh_direct_signal")
        direct_supported <- count_true(md, "gnrh_direct_supported")
        direct_isolated <- count_true(md, "gnrh_direct_isolated")
        reference_positive <- count_true(md, "gnrh_reference_positive")
        confident <- count_true(md, "gnrh_confident")
        candidate_column <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          "gnrh_transcriptomic_candidate"
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          "gnrh_dropout_candidate"
        } else {
          NULL
        }
        transcriptomic_candidates <- if (!is.null(candidate_column)) {
          count_true(md, candidate_column)
        } else {
          NA_integer_
        }
        tibble::tibble(
          id = id,
          n_cells = total,
          direct = direct,
          transcriptomic = transcriptomic,
          supported = supported,
          gnrh_pos = positive_status,
          positive_by_class = positive_class,
          direct_signal = direct_signal,
          direct_supported = direct_supported,
          direct_isolated = direct_isolated,
          reference_positive = reference_positive,
          transcriptomic_candidates = transcriptomic_candidates,
          gnrh_confident = confident,
          pct_gnrh = if (total > 0L) 100 * positive_status / total else NA_real_,
          pct_direct = if (positive_status > 0L) 100 * direct / positive_status else NA_real_,
          pct_transcriptomic = if (positive_status > 0L) 100 * transcriptomic / positive_status else NA_real_,
          pct_supported = if (positive_status > 0L) 100 * supported / positive_status else NA_real_,
          pct_confident = if (positive_status > 0L && !is.na(confident)) 100 * confident / positive_status else NA_real_,
          pct_direct_isolated = if (!is.na(direct_signal) && direct_signal > 0L && !is.na(direct_isolated)) {
            100 * direct_isolated / direct_signal
          } else {
            NA_real_
          },
          pct_transcriptomic_candidates = if (!is.na(transcriptomic_candidates) && total > 0L) {
            100 * transcriptomic_candidates / total
          } else {
            NA_real_
          }
        )
      }
    ) |>
      dplyr::left_join(dataset_metadata, by = "id") |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        dplyr::everything()
      )
    # ========================================================================= #
    # Transcriptomic candidates
    # ========================================================================= #
    transcriptomic_candidates <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        column <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          "gnrh_transcriptomic_candidate"
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          "gnrh_dropout_candidate"
        } else {
          return(tibble::tibble())
        }
        x <- md[
          md[[column]] %in% TRUE,
          ,
          drop = FALSE
        ]
        median_col <- function(column) {
          if (!column %in% colnames(x) || !nrow(x)) return(NA_real_)
          safe_median(x[[column]])
        }
        tibble::tibble(
          id = id,
          n_cells = nrow(x),
          pct_dataset = if (nrow(md)) 100 * nrow(x) / nrow(md) else NA_real_,
          pct_GNRH1_detected = if (nrow(x) && "gnrh_raw" %in% colnames(x)) {
            100 * mean(x$gnrh_raw > 0, na.rm = TRUE)
          } else {
            0
          },
          median_GNRH1 = median_col("gnrh_raw"),
          median_support = median_col("gnrh_support_score_raw"),
          median_identity_primary = median_col("gnrh_identity_primary_hits"),
          median_core = median_col("gnrh_core_hits"),
          median_neuro = median_col("gnrh_neuro_hits"),
          median_knn = median_col("gnrh_knn"),
          median_alternative = median_col("gnrh_alternative_score")
        )
      }
    )
    transcriptomic_candidates <- join_metadata(
      transcriptomic_candidates,
      dataset_metadata
    )
    # ========================================================================= #
    # Status / class consistency
    # ========================================================================= #
    status_class <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        object[[]] |>
          tibble::as_tibble() |>
          dplyr::transmute(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_class = as.character(.data$gnrh_class)
          ) |>
          dplyr::count(
            .data$gnrh_status,
            .data$gnrh_class,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_cells = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    status_class <- join_metadata(status_class, dataset_metadata)
    # ========================================================================= #
    # Logical consistency
    # ========================================================================= #
    classification_consistency <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        status <- as.character(md$gnrh_status)
        class <- as.character(md$gnrh_class)
        raw <- if ("gnrh_raw" %in% colnames(md)) {
          suppressWarnings(as.numeric(md$gnrh_raw))
        } else {
          rep(NA_real_, nrow(md))
        }
        candidate <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          md$gnrh_transcriptomic_candidate %in% TRUE
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          md$gnrh_dropout_candidate %in% TRUE
        } else {
          rep(FALSE, nrow(md))
        }
        expected_positive_class <- class %in% positive_classes
        tibble::tibble(
          id = id,
          n_status_class_discordant_positive = sum(
            status == "pos" & !expected_positive_class,
            na.rm = TRUE
          ),
          n_status_class_discordant_negative = sum(
            status == "neg" & expected_positive_class,
            na.rm = TRUE
          ),
          n_unexplained_positive_without_GNRH1 = sum(
            status == "pos" &
              is.finite(raw) &
              raw <= 0 &
              !candidate,
            na.rm = TRUE
          ),
          n_expected_transcriptomic_positive = sum(
            status == "pos" & candidate,
            na.rm = TRUE
          ),
          n_direct_without_GNRH1 = sum(
            class == "direct" &
              is.finite(raw) &
              raw <= 0,
            na.rm = TRUE
          )
        )
      }
    )
    classification_consistency <- join_metadata(
      classification_consistency,
      dataset_metadata
    )
    # ========================================================================= #
    # Evidence scores
    # ========================================================================= #
    candidate_score_columns <- unique(c(
      "gnrh_score",
      "gnrh_score_raw",
      "gnrh_support_score",
      "gnrh_support_score_raw",
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
      "gnrh_stage_early_score",
      "gnrh_stage_migrating_score",
      "gnrh_stage_mature_score",
      "gnrh_knn"
    ))
    available_columns <- Reduce(
      intersect,
      lapply(
        gnrh_list,
        function(object) colnames(object[[]])
      )
    )
    score_columns <- intersect(
      candidate_score_columns,
      available_columns
    )
    scores <- if (length(score_columns)) {
      purrr::imap_dfr(
        gnrh_list,
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::mutate(
              gnrh_class = as.character(.data$gnrh_class)
            ) |>
            dplyr::filter(
              .data$gnrh_class %in% positive_classes
            ) |>
            dplyr::group_by(
              .data$gnrh_class
            ) |>
            dplyr::summarise(
              n_cells = dplyr::n(),
              dplyr::across(
                dplyr::all_of(score_columns),
                list(
                  median = ~ safe_median(.x),
                  q25 = ~ safe_quantile(.x, 0.25),
                  q75 = ~ safe_quantile(.x, 0.75)
                )
              ),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
    } else {
      tibble::tibble()
    }
    scores <- join_metadata(scores, dataset_metadata)
    # ========================================================================= #
    # Developmental stages
    # ========================================================================= #
    stages <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_stage = as.character(.data$gnrh_stage)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_stage),
            .data$gnrh_stage != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
        x |>
          dplyr::count(
            .data$gnrh_stage,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_positive = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    stages <- join_metadata(stages, dataset_metadata)
    stage_class <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_class = as.character(.data$gnrh_class),
            gnrh_stage = as.character(.data$gnrh_stage)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_stage),
            .data$gnrh_stage != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
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
            pct_class = 100 * .data$n_cells / sum(.data$n_cells)
          ) |>
          dplyr::ungroup() |>
          dplyr::mutate(id = id)
      }
    )
    stage_class <- join_metadata(stage_class, dataset_metadata)
    # ========================================================================= #
    # Secretory phenotype
    # ========================================================================= #
    secretory <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        if (!"gnrh_secretory" %in% colnames(md)) {
          return(tibble::tibble())
        }
        x <- md |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_secretory = as.character(.data$gnrh_secretory)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_secretory),
            .data$gnrh_secretory != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
        x |>
          dplyr::count(
            .data$gnrh_secretory,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_positive = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    secretory <- join_metadata(secretory, dataset_metadata)
    # ========================================================================= #
    # Stage refinement
    # ========================================================================= #
    refinement_required <- c(
      "gnrh_stage_raw",
      "gnrh_stage_reassigned",
      "gnrh_stage_reason"
    )
    has_refinement <- vapply(
      gnrh_list,
      function(object) {
        all(refinement_required %in% colnames(object[[]]))
      },
      logical(1)
    )
    stage_refinement <- tibble::tibble()
    stage_reassignment <- tibble::tibble()
    if (any(has_refinement)) {
      stage_refinement <- purrr::imap_dfr(
        gnrh_list[has_refinement],
        function(object, id) {
          md <- object[[]]
          positive <- as.character(md$gnrh_status) == "pos"
          n_positive <- sum(positive, na.rm = TRUE)
          n_reassigned <- sum(
            positive & md$gnrh_stage_reassigned %in% TRUE,
            na.rm = TRUE
          )
          n_filtered <- sum(
            positive &
              as.character(md$gnrh_stage_reason) == "migration_core_filter",
            na.rm = TRUE
          )
          tibble::tibble(
            id = id,
            n_positive = n_positive,
            n_reassigned = n_reassigned,
            pct_reassigned = if (n_positive > 0L) {
              100 * n_reassigned / n_positive
            } else {
              NA_real_
            },
            n_migration_filtered = n_filtered
          )
        }
      )
      stage_refinement <- join_metadata(
        stage_refinement,
        dataset_metadata
      )
      stage_reassignment <- purrr::imap_dfr(
        gnrh_list[has_refinement],
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::filter(
              as.character(.data$gnrh_status) == "pos"
            ) |>
            dplyr::transmute(
              id = id,
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw),
              gnrh_stage = as.character(.data$gnrh_stage),
              gnrh_stage_reason = as.character(.data$gnrh_stage_reason)
            ) |>
            dplyr::count(
              .data$id,
              .data$gnrh_stage_raw,
              .data$gnrh_stage,
              .data$gnrh_stage_reason,
              name = "n_cells"
            )
        }
      )
      stage_reassignment <- join_metadata(
        stage_reassignment,
        dataset_metadata
      )
    }
    # ========================================================================= #
    # Migration refinement
    # ========================================================================= #
    migration_required <- c(
      "gnrh_stage_raw",
      "gnrh_stage",
      "gnrh_migration_core_hits"
    )
    has_migration <- vapply(
      gnrh_list,
      function(object) {
        all(migration_required %in% colnames(object[[]]))
      },
      logical(1)
    )
    migration_core <- tibble::tibble()
    migration_refinement <- tibble::tibble()
    migration_refinement_summary <- tibble::tibble()
    if (any(has_migration)) {
      migration_core <- purrr::imap_dfr(
        gnrh_list[has_migration],
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::filter(
              as.character(.data$gnrh_status) == "pos"
            ) |>
            dplyr::mutate(
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw)
            ) |>
            dplyr::group_by(
              .data$gnrh_stage_raw
            ) |>
            dplyr::summarise(
              n = dplyr::n(),
              median_hits = safe_median(.data$gnrh_migration_core_hits),
              q25_hits = safe_quantile(.data$gnrh_migration_core_hits, 0.25),
              q75_hits = safe_quantile(.data$gnrh_migration_core_hits, 0.75),
              pct_ge1 = 100 * mean(.data$gnrh_migration_core_hits >= 1, na.rm = TRUE),
              pct_ge2 = 100 * mean(.data$gnrh_migration_core_hits >= 2, na.rm = TRUE),
              pct_ge3 = 100 * mean(.data$gnrh_migration_core_hits >= 3, na.rm = TRUE),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
      migration_core <- join_metadata(
        migration_core,
        dataset_metadata
      )
      migration_refinement <- purrr::imap_dfr(
        gnrh_list[has_migration],
        function(object, id) {
          x <- object[[]] |>
            tibble::as_tibble() |>
            dplyr::mutate(
              gnrh_status = as.character(.data$gnrh_status),
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw),
              gnrh_stage = as.character(.data$gnrh_stage)
            ) |>
            dplyr::filter(
              .data$gnrh_status == "pos",
              .data$gnrh_stage_raw == "migrating"
            ) |>
            dplyr::mutate(
              migration_outcome = ifelse(
                .data$gnrh_stage == "migrating",
                "retained_migrating",
                "reassigned"
              )
            )
          if (!nrow(x)) return(tibble::tibble())
          x |>
            dplyr::group_by(
              .data$migration_outcome,
              .data$gnrh_stage
            ) |>
            dplyr::summarise(
              n_cells = dplyr::n(),
              median_hits = safe_median(.data$gnrh_migration_core_hits),
              mean_hits = mean(.data$gnrh_migration_core_hits, na.rm = TRUE),
              pct_zero = 100 * mean(.data$gnrh_migration_core_hits == 0, na.rm = TRUE),
              pct_ge1 = 100 * mean(.data$gnrh_migration_core_hits >= 1, na.rm = TRUE),
              pct_ge2 = 100 * mean(.data$gnrh_migration_core_hits >= 2, na.rm = TRUE),
              pct_ge3 = 100 * mean(.data$gnrh_migration_core_hits >= 3, na.rm = TRUE),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
      migration_refinement <- join_metadata(
        migration_refinement,
        dataset_metadata
      )
      if (nrow(migration_refinement)) {
        migration_refinement_summary <- migration_refinement |>
          dplyr::group_by(
            .data$id,
            .data$label,
            .data$species
          ) |>
          dplyr::summarise(
            n_raw_migrating = sum(.data$n_cells, na.rm = TRUE),
            n_retained = sum(
              .data$n_cells[.data$migration_outcome == "retained_migrating"],
              na.rm = TRUE
            ),
            n_reassigned = sum(
              .data$n_cells[.data$migration_outcome == "reassigned"],
              na.rm = TRUE
            ),
            n_to_early = sum(
              .data$n_cells[
                .data$migration_outcome == "reassigned" &
                  .data$gnrh_stage == "early"
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
            .groups = "drop"
          ) |>
          dplyr::mutate(
            pct_retained = dplyr::if_else(
              .data$n_raw_migrating > 0,
              100 * .data$n_retained / .data$n_raw_migrating,
              NA_real_
            ),
            pct_reassigned = dplyr::if_else(
              .data$n_raw_migrating > 0,
              100 * .data$n_reassigned / .data$n_raw_migrating,
              NA_real_
            )
          )
      }
    }
    # ========================================================================= #
    # Biological marker validation
    # ========================================================================= #
    biological_markers <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        cells <- rownames(md)[
          as.character(md$gnrh_status) == "pos"
        ]
        if (!length(cells)) return(tibble::tibble())
        expr <- tryCatch(
          .get_expr(
            object,
            assay = assay,
            layer = layer
          ),
          error = function(e) NULL
        )
        if (is.null(expr)) return(tibble::tibble())
        matched <- .match_genes(
          validation_markers,
          rownames(expr)
        )
        if (!length(matched)) return(tibble::tibble())
        cells <- intersect(cells, colnames(expr))
        if (!length(cells)) return(tibble::tibble())
        classes <- as.character(
          md[cells, "gnrh_class", drop = TRUE]
        )
        observed_positive_classes <- intersect(
          positive_classes,
          unique(classes)
        )
        purrr::map_dfr(
          observed_positive_classes,
          function(cl) {
            idx <- classes == cl
            if (!any(idx)) return(tibble::tibble())
            x <- expr[
              matched,
              cells[idx],
              drop = FALSE
            ]
            tibble::tibble(
              id = id,
            gnrh_class = cl,
            gene = rownames(x),
            avg_expression = Matrix::rowMeans(x),
            pct_expressing = 100 * Matrix::rowMeans(x > 0)
          )
        }
      )
    }
  )
  biological_markers <- join_metadata(
    biological_markers,
    dataset_metadata
  )
  # ========================================================================= #
  # Validation outputs
  # ========================================================================= #
  validation_tables <- list(
    input_summary = input_summary,
    detection = detection,
    transcriptomic_candidates = transcriptomic_candidates,
    classification_consistency = classification_consistency,
    status_class = status_class,
    scores = scores,
    stages = stages,
    stage_class = stage_class,
    secretory = secretory,
    stage_refinement = stage_refinement,
    stage_reassignment = stage_reassignment,
    migration_core = migration_core,
    migration_refinement = migration_refinement,
    migration_refinement_summary = migration_refinement_summary,
    biological_markers = biological_markers
  )
  # ========================================================================= #
  # Write outputs
  # ========================================================================= #
  if (isTRUE(write_output)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    purrr::iwalk(
      validation_tables,
      function(x, name) {
        if (is.data.frame(x) && ncol(x)) {
          readr::write_csv(
            x,
            file.path(output_dir, paste0(name, ".csv"))
          )
        }
      }
    )
  }
  # ========================================================================= #
  # Console summary
  # ========================================================================= #
  total_cells <- sum(detection$n_cells, na.rm = TRUE)
  total_positive <- sum(detection$gnrh_pos, na.rm = TRUE)
  has_candidate_counts <- any(
    is.finite(detection$transcriptomic_candidates)
  )
  total_candidates <- if (has_candidate_counts) {
    sum(detection$transcriptomic_candidates, na.rm = TRUE)
  } else {
    NA_integer_
  }
  log(
    sprintf(
      "Validated %d datasets | %s cells | %s GnRH-positive",
      length(gnrh_list),
      format(total_cells, big.mark = ",", trim = TRUE),
      format(total_positive, big.mark = ",", trim = TRUE)
    ),
    type = "info"
  )
  if (is.finite(total_candidates)) {
    log(
      sprintf(
        "Transcriptomic candidates: %s",
        format(total_candidates, big.mark = ",", trim = TRUE)
      ),
      type = "info"
    )
  }
  # ========================================================================= #
  # Consistency summary
  # ========================================================================= #
  violation_columns <- intersect(
    c(
      "n_status_class_discordant_positive",
      "n_status_class_discordant_negative",
      "n_unexplained_positive_without_GNRH1",
      "n_direct_without_GNRH1"
    ),
    colnames(classification_consistency)
  )
  bad_consistency <- if (length(violation_columns)) {
    rowSums(
      as.data.frame(
        classification_consistency[
          ,
          violation_columns,
          drop = FALSE
        ]
      ),
      na.rm = TRUE
    )
  } else {
    rep(0, nrow(classification_consistency))
  }
  if (any(bad_consistency > 0)) {
    log(
      "Classification consistency violations detected.",
      type = "warn"
    )
  } else {
    log(
      "Classification consistency checks passed.",
      type = "done",
      duration = 0
    )
  }
  # ========================================================================= #
  # Return
  # ========================================================================= #
  result <- c(
    list(
      datasets = dataset_metadata
    ),
    validation_tables,
    list(
      parameters = list(
        positive_classes = positive_classes,
        assay = assay,
        layer = layer,
        validation_markers = validation_markers
      ),
      output_dir = if (isTRUE(write_output)) {
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

