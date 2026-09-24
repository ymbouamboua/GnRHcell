
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
  reduction <- gnrh_reduction(object, reduction)
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
  run_info <- .extract_gnrh_run_info(
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
      figures$status <- GnRHcell::gnrh_cellmap(
        object = object,
        group_by = "gnrh_status",
        reduction = reduction,
        plot.ttl = paste0(dataset_label)
      )
      .gnrh_save_plot(
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
      figures$stage <- GnRHcell::gnrh_cellmap(
        object = object,
        group_by = "gnrh_stage",
        reduction = reduction,
        plot.ttl = paste0(dataset_label)
      )
      .gnrh_save_plot(
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
      figures$core <- GnRHcell::gnrh_cellfeat(
        object = object,
        features = core_features,
        reduction = reduction,
        ncol = 2
      )
      .gnrh_save_plot(
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
      figures$staging <- GnRHcell::gnrh_cellfeat(
        object = object,
        features = staging_features,
        reduction = reduction,
        ncol = min(3L, length(staging_features))
      )
      .gnrh_save_plot(
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
      figures$status_distribution <- GnRHcell::gnrh_celldistribution(
        object,
        group.by = "gnrh_status",
        split.by = split_by,
        proportion = TRUE,
        label = FALSE,
        cols = GnRHcell::gnrh_palette("status")
      )
      .gnrh_save_plot(
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
      figures$stage_distribution <- GnRHcell::gnrh_celldistribution(
        object,
        group.by = "gnrh_stage",
        split.by = split_by,
        proportion = TRUE,
        label = FALSE,
        cols = GnRHcell::gnrh_palette("stage")
      )
      .gnrh_save_plot(
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
      .gnrh_save_plot(
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
        GnRHcell::gnrh_coexpr(
          markers,
          coexp_cutoff=0.10,
          top_n=12L,
          exclude_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH co-expression plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$codetection <- tryCatch(
        GnRHcell::gnrh_codetect(
          markers,
          top_n=10L,
          min_or=2,
          max_fdr=0.05,
          min_specificity=0.05,
          exclude_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH co-detection plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$detection <- tryCatch(
        GnRHcell::gnrh_pheno(
          markers,
          top_n=10L,
          min_specificity=0.05,
          max_padj=0.05,
          exclude_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
          style="bw"
        ),
        error=function(e) {
          log("Skipping GnRH phenotype-association plot:",conditionMessage(e),type="warn")
          NULL
        }
      )
      if (!is.null(figures$coexpression)) {
        .gnrh_save_plot(
          figures$coexpression,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_coexpression")),
          width=6,
          height=6
        )
      }
      if (!is.null(figures$codetection)) {
        .gnrh_save_plot(
          figures$codetection,
          file.path(dataset_dir,paste0(dataset_id,"_gnrh_codetection")),
          width=6,
          height=6
        )
      }
      if (!is.null(figures$detection)) {
        .gnrh_save_plot(
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
        .gnrh_save_plot(
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
        GnRHcell::gnrh_network(
          markers,
          mode="coexpression",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
          style="void"
        ),
        error=function(e) {
          log("Skipping co-expression network:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$network_codetection <- tryCatch(
        GnRHcell::gnrh_network(
          markers,
          mode="codetection",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
          style="void"
        ),
        error=function(e) {
          log("Skipping co-detection network:",conditionMessage(e),type="warn")
          NULL
        }
      )
      figures$network_phenotype <- tryCatch(
        GnRHcell::gnrh_network(
          markers,
          mode="phenotype",
          top_n=30L,
          threshold=0.10,
          include_gnrh=TRUE,
          txtsize = getOption("gnrhcell.base_size", 14),
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
            .gnrh_save_plot(
              p,
              file.path(dataset_dir,paste0(dataset_id,"_gnrh_network_",nm)),
              width=9,
              height=9
            )
          }
        )
        figures$networks <- patchwork::wrap_plots(networks,nrow=1)
        .gnrh_save_plot(
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
          GnRHcell::gnrh_network(
            markers,
            top_n=30L,
            threshold=0.10,
            mode=network_mode,
            include_gnrh=TRUE,
            txtsize = getOption("gnrhcell.base_size", 14),
            style="void"
          ),
          error=function(e) {
            log("Skipping marker network:",conditionMessage(e),type="warn")
            NULL
          }
        )
        if (!is.null(figures$network)) {
          .gnrh_save_plot(
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

