# ========================================================================= #
# Selected-dataset figure workflow
# ========================================================================= #

#' Generate a complete GnRHcell figure set for selected datasets
#'
#' Runs per-dataset visualizations and cross-dataset comparisons for a selected
#' subset of a GnRHcell collection.
#'
#' @param collection Result returned by [run_gnrh_collection()].
#' @param selected_ids At least two dataset IDs or unique dataset labels.
#' @param output_dir Directory receiving figures and the result RDS.
#' @param run_individual Generate figures for each selected dataset.
#' @param run_associations Generate marker association plots.
#' @param run_networks Generate marker networks.
#' @param run_gallery Generate the comparison gallery.
#' @param embedding_fields Metadata fields used for embedding plots.
#' @param feature_presets Presets passed to [gnrh_cellfeat()].
#' @param top_conserved Maximum conserved markers in the dot plot.
#' @param top_specific Maximum specific markers selected per dataset.
#' @param top_program Maximum high-confidence genes per program.
#' @param validation_genes Optional external gene signature used for cross-model
#'   validation with [gnrh_signature()].
#' @param validation_name Title describing `validation_genes`.
#' @param show_titles Display main figure titles. Set to `FALSE` when titles
#'   are supplied by a poster or multi-panel layout. Dataset labels, axes,
#'   facet labels, and legends are retained.
#' @param dpi PNG resolution.
#' @param save_objects Logical; include the selected and merged Seurat objects
#'   in `gnrh_selected_dataset_figures.rds`. Default is `FALSE` to
#'   avoid duplicating objects already stored in the collection result. Objects
#'   remain available in the value returned to the current R session.
#' @param save_result Logical; write `gnrh_selected_dataset_figures.rds`.
#'   Default is `FALSE`. This is independent of figure export, which is
#'   always performed by the workflow.
#' @param verbose Display progress messages.
#' @return A list containing datasets, objects, merged GnRH cells, comparisons,
#'   plots, and the output directory.
#' @export
gnrh_compare <- function(
    collection,
    selected_ids,
    output_dir = file.path(collection$output_dir, "selected_dataset_comparison"),
    run_individual = TRUE,
    run_associations = TRUE,
    run_networks = TRUE,
    run_gallery = TRUE,
    embedding_fields = c("gnrh_status", "gnrh_class", "gnrh_stage", "gnrh_confident", "gnrh_secretory"),
    feature_presets = c("core", "modules", "staging", "migration", "secretory", "hits"),
    top_conserved = 25L,
    top_specific = 12L,
    top_program = 8L,
    validation_genes = NULL,
    validation_name = "Published human GnRH signature",
    show_titles = TRUE,
    dpi = 600,
    save_objects = FALSE,
    save_result = FALSE,
    verbose = TRUE
) {
  if (!is.list(collection) || is.null(collection$results) || is.null(collection$datasets)) {
    stop("`collection` must be returned by `run_gnrh_collection()`.", call. = FALSE)
  }
  if (!is.character(selected_ids) || length(selected_ids) < 2L) {
    stop("`selected_ids` must contain at least two dataset IDs or labels.", call. = FALSE)
  }
  if (!is.logical(show_titles) || length(show_titles) != 1L || is.na(show_titles)) {
    stop("`show_titles` must be TRUE or FALSE.", call. = FALSE)
  }
  datasets <- as.data.frame(collection$datasets, stringsAsFactors = FALSE)
  resolve_id <- function(value) {
    hit <- unique(c(which(datasets$id == value), which(datasets$label == value)))
    if (length(hit) != 1L) NA_integer_ else hit[[1L]]
  }
  rows <- vapply(selected_ids, resolve_id, integer(1))
  if (anyNA(rows)) {
    stop("Dataset selection not found or ambiguous: ", paste(selected_ids[is.na(rows)], collapse = ", "), ".", call. = FALSE)
  }
  datasets_sel <- datasets[rows, , drop = FALSE]
  if (anyDuplicated(datasets_sel$id)) stop("Duplicated datasets in `selected_ids`.", call. = FALSE)

  log <- .msg(verbose)
  total_start <- Sys.time()
  step_times <- list()
  log("==== STARTING GnRHcell FIGURE WORKFLOW ====", type = "header")
  log(sprintf("Selected datasets: %d", nrow(datasets_sel)), type = "info")
  for (i in seq_len(nrow(datasets_sel))) {
    log(sprintf("  %d. %s [%s]", i, datasets_sel$label[[i]], datasets_sel$id[[i]]), type = "info")
  }
  log("[1/5] Recovering processed Seurat objects", type = "step")
  step_start <- Sys.time()

  all_objects <- .get_collection_objects(collection, allow_saved = TRUE)
  objects <- all_objects[datasets_sel$id]
  unavailable <- datasets_sel$id[vapply(objects, is.null, logical(1))]
  if (length(unavailable)) stop("Processed object(s) unavailable: ", paste(unavailable, collapse = ", "), ".", call. = FALSE)
  step_times$objects_sec <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  log("Objects recovered.", type = "done", duration = step_times$objects_sec)

  dataset_dir <- file.path(output_dir, "datasets")
  comparison_dir <- file.path(output_dir, "comparisons")
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  individual_plots <- stats::setNames(vector("list", length(objects)), names(objects))
  safe_plot <- function(expr, context) tryCatch(expr, error = function(e) {
    warning(context, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
  without_titles <- function(plot) {
    if (!isTRUE(show_titles) && !is.null(plot)) {
      plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL)
    }
    plot
  }
  save_one <- function(plot, path, width, height) {
    if (!is.null(plot)) gnrh_save(plot, path, width = width, height = height, dpi = dpi)
    invisible(plot)
  }

  log("[2/5] Generating dataset-level figures", type = "step")
  step_start <- Sys.time()
  if (isTRUE(run_individual)) {
    for (i in seq_along(objects)) {
      id <- datasets_sel$id[[i]]
      label <- datasets_sel$label[[i]]
      object <- objects[[i]]
      log(sprintf("  Dataset %d/%d: %s", i, length(objects), label), type = "info")
      reduction <- gnrh_reduction(object, if ("reduction" %in% names(datasets_sel)) datasets_sel$reduction[[i]] else NULL)
      figdir <- file.path(dataset_dir, id)
      dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
      plots <- list(embeddings = list(), features = list(), associations = list(), networks = list())
      for (field in intersect(embedding_fields, names(object[[]]))) {
        plot <- safe_plot(gnrh_cellmap(
          object,
          group_by = field,
          reduction = reduction,
          plot.ttl = if (isTRUE(show_titles)) label else NULL
        ), paste(id, field, sep = " | "))
        save_one(plot, file.path(figdir, paste0(id, "_", field)), 6, 5)
        plots$embeddings[[field]] <- plot
      }
      for (preset in feature_presets) {
        preset_features <- .gnrh_features(preset)
        plot <- safe_plot(gnrh_cellfeat(
          object, preset = preset, reduction = reduction, ncol = 3,
          txtsize = getOption("gnrhcell.base_size", 14), merge.leg = preset %in% c("secretory", "hits")
        ), paste(id, preset, sep = " | "))
        panel_n <- max(1L, length(intersect(preset_features, c(rownames(object), names(object[[]])))))
        save_one(
          plot,
          file.path(figdir, paste0(id, "_gnrh_", preset)),
          max(7, 4.6 * min(3L, panel_n)),
          max(5.8, 5.2 * ceiling(panel_n / 3))
        )
        plots$features[[preset]] <- plot
      }
      markers <- collection$results[[id]]$markers %||% NULL
      marker_file <- file.path(collection$output_dir, "markers", paste0("gnrh_", id, "_markers.tsv"))
      if (is.null(markers) && file.exists(marker_file)) markers <- utils::read.delim(marker_file, check.names = FALSE, stringsAsFactors = FALSE)
      if (is.data.frame(markers) && nrow(markers) && isTRUE(run_associations)) {
        plots$associations <- lapply(Filter(Negate(is.null), list(
          coexpression = safe_plot(gnrh_coexpr(markers, coexp_cutoff = 0.10, top_n = 12L), paste(id, "coexpression")),
          codetection = safe_plot(gnrh_codetect(markers, top_n = 10L, min_or = 2, max_fdr = 0.05, min_specificity = 0.05), paste(id, "codetection")),
          detection = safe_plot(gnrh_pheno(markers, top_n = 10L, min_specificity = 0.05, max_padj = 0.05), paste(id, "detection"))
        )), without_titles)
        for (name in names(plots$associations)) save_one(plots$associations[[name]], file.path(figdir, paste0(id, "_gnrh_", name)), 6, 6)
        if (length(plots$associations)) save_one(patchwork::wrap_plots(plots$associations, nrow = 1), file.path(figdir, paste0(id, "_gnrh_marker_associations")), 5.5 * length(plots$associations), 5.5)
      }
      if (is.data.frame(markers) && nrow(markers) && isTRUE(run_networks)) {
        modes <- c("coexpression", "codetection", "phenotype")
        plots$networks <- stats::setNames(lapply(modes, function(mode) safe_plot(
          gnrh_network(markers, mode = mode, top_n = 30L, threshold = 0.10, txtsize = getOption("gnrhcell.base_size", 14), style = "void"), paste(id, "network", mode)
        )), modes)
        plots$networks <- lapply(Filter(Negate(is.null), plots$networks), without_titles)
        for (name in names(plots$networks)) save_one(plots$networks[[name]], file.path(figdir, paste0(id, "_gnrh_network_", name)), 9, 9)
        if (length(plots$networks)) save_one(
          patchwork::wrap_plots(plots$networks, nrow = 1),
          file.path(figdir, paste0(id, "_gnrh_networks")),
          8 * length(plots$networks), 8
        )
      }
      individual_plots[[id]] <- plots
      invisible(gc())
    }
  } else {
    log("  Dataset-level figures disabled (`run_individual = FALSE`).", type = "info")
  }
  step_times$individual_sec <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  log("Dataset-level figures complete.", type = "done", duration = step_times$individual_sec)

  log("[3/5] Running cross-dataset comparisons", type = "step")
  step_start <- Sys.time()
  comparisons <- compare_gnrh_datasets(
    datasets = datasets_sel, output_dir = collection$output_dir,
    comparison_dir = comparison_dir, run_programs = TRUE,
    run_marker_heatmaps = TRUE, run_gallery = run_gallery, objects = objects,
    show_titles = show_titles
  )
  step_times$comparison_sec <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  log("Cross-dataset comparisons complete.", type = "done", duration = step_times$comparison_sec)

  log("[4/5] Generating comparative GnRH dot plots", type = "step")
  step_start <- Sys.time()
  combined_gnrh <- combine_gnrh_datasets(
    collection, selected_ids = datasets_sel$id, gnrh_only = TRUE,
    group.by = "dataset", label_by = "label"
  )
  common_genes <- Reduce(intersect, lapply(objects, rownames))
  comparison_plots <- list()
  save_dot <- function(features, title, filename, height = 8) {
    if (is.list(features)) {
      features <- lapply(features, intersect, y = common_genes)
      features <- features[lengths(features) > 0L]
    } else features <- intersect(features, common_genes)
    if (!length(features) || (is.list(features) && !length(unlist(features)))) return(NULL)
    n_dot_genes <- length(unique(unlist(features, use.names = FALSE)))
    n_dot_facets <- if (is.list(features)) length(features) else 1L
    dot_dims <- .gnrh_dot_dims(
      labels = datasets_sel$label,
      n_genes = n_dot_genes,
      n_facets = n_dot_facets
    )
    dot_height <- max(height, unname(dot_dims["height"]))
    plot <- gnrh_celldot(
      combined_gnrh, features = features, group.by = "dataset", scale = FALSE,
      dot.scale = 6, flip = TRUE, txtsize = 14, leg.size = 12,
      leg.ttl.size = 12, x.ang = unname(dot_dims["x_angle"]),
      hjust.x = 1, vjust.x = 1,
      title = if (isTRUE(show_titles)) title else NULL
    )
    save_one(
      plot,
      file.path(comparison_dir, filename),
      unname(dot_dims["width"]),
      dot_height
    )
    plot
  }

  stage_modules <- gnrh_stage_modules()
  developmental <- list(
    `GnRH identity` = c("GNRH1", "ISL1", "FEZF1", "SIX3", "SIX6"),
    `Early specification` = stage_modules$early,
    Migration = stage_modules$migrating,
    Maturation = stage_modules$mature,
    Secretion = gnrh_secretory_marker_module()
  )
  comparison_plots$developmental <- save_dot(developmental, "GnRH developmental programs across selected datasets", "gnrh_developmental_programs_dotplot")
  conserved <- comparisons$conserved_markers$conserved %||% NULL
  if (is.data.frame(conserved) && nrow(conserved)) comparison_plots$conserved <- save_dot(
    utils::head(conserved$gene, as.integer(top_conserved)), "Conserved GnRH markers", "conserved_gnrh_markers_dotplot"
  )
  specific <- comparisons$dataset_specific_markers$specific_table %||% NULL
  if (is.data.frame(specific) && nrow(specific)) {
    # Match the heatmap: selected dataset order, then decreasing specificity
    # within each dataset. coord_flip() displays the feature vector bottom-up,
    # so reverse the final vector to preserve that top-to-bottom order.
    specific <- specific[
      order(match(specific$best_dataset, datasets_sel$id), -specific$specificity_score),
      ,
      drop = FALSE
    ]
    genes_by_dataset <- split(specific$gene, specific$best_dataset)
    ordered_groups <- c(
      intersect(datasets_sel$id, names(genes_by_dataset)),
      setdiff(names(genes_by_dataset), datasets_sel$id)
    )
    genes_by_dataset <- genes_by_dataset[ordered_groups]
    genes_by_dataset <- lapply(
      genes_by_dataset,
      utils::head,
      n = as.integer(top_specific)
    )
    heatmap_gene_order <- unique(unlist(genes_by_dataset, use.names = FALSE))
    comparison_plots$specific <- save_dot(
      rev(heatmap_gene_order),
      "Dataset-specific GnRH markers: direct comparison",
      "dataset_specific_gnrh_markers_dotplot",
      9
    )
  }
  high_label_dims <- .gnrh_dot_dims(datasets_sel$label, n_genes = 1L)
  comparison_plots$high_confidence <- without_titles(safe_plot(gnrh_program(
    combined_gnrh, programs = comparisons$programs, table = "high_confidence",
    group.by = "dataset", top_n = top_program, scale = FALSE,
    dot.scale = 6, txtsize = 14, leg.size = 12, leg.ttl.size = 12,
    x.ang = unname(high_label_dims["x_angle"]), hjust.x = 1, vjust.x = 1
  ), "High-confidence marker-program dot plot"))
  high_features <- attr(comparison_plots$high_confidence, "gnrh_features")
  high_n <- if (is.null(high_features)) 12L else length(unique(unlist(high_features, use.names = FALSE)))
  high_dims <- .gnrh_dot_dims(
    labels = datasets_sel$label,
    n_genes = high_n,
    n_facets = max(1L, length(high_features))
  )
  save_one(
    comparison_plots$high_confidence,
    file.path(comparison_dir, "gnrh_high_confidence_marker_programs_dotplot"),
    unname(high_dims["width"]),
    max(8, unname(high_dims["height"]))
  )
  if (!is.null(validation_genes)) {
    comparison_plots$validation_signature <- safe_plot(
      gnrh_signature(
        collection = collection,
        genes = validation_genes,
        selected_ids = datasets_sel$id,
        signature_name = validation_name,
        scale_expression = TRUE,
        show_title = show_titles,
        txtsize = getOption("gnrhcell.base_size", 14)
      ),
      "External signature cross-model validation"
    )
    validation_dims <- attr(comparison_plots$validation_signature, "gnrh_dimensions")
    if (!is.null(comparison_plots$validation_signature)) {
      save_one(
        comparison_plots$validation_signature,
        file.path(comparison_dir, "gnrh_cross_model_signature_dotplot"),
        unname(validation_dims["width"]),
        unname(validation_dims["height"])
      )
      signature_table <- attr(comparison_plots$validation_signature, "gnrh_signature_data")
      if (!is.null(signature_table)) {
        readr::write_tsv(
          signature_table,
          file.path(comparison_dir, "gnrh_cross_model_signature.tsv")
        )
      }
    }
  }
  step_times$dotplots_sec <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  log("Comparative dot plots complete.", type = "done", duration = step_times$dotplots_sec)

  log("[5/5] Saving figure-workflow results", type = "step")
  step_start <- Sys.time()
  result <- list(
    datasets = datasets_sel, objects = objects, combined_gnrh = combined_gnrh,
    comparisons = comparisons,
    plots = list(individual = individual_plots, comparison = comparison_plots),
    output_dir = normalizePath(output_dir, mustWork = FALSE)
  )
  if (isTRUE(save_result)) {
    saved_result <- result
    if (!isTRUE(save_objects)) {
      saved_result$objects <- NULL
      saved_result$combined_gnrh <- NULL
    }
    saveRDS(saved_result, file.path(output_dir, "gnrh_selected_dataset_figures.rds"))
  }
  step_times$save_sec <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  step_times$total_sec <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  log(
    if (isTRUE(save_result)) "RDS result saved." else "RDS result export skipped.",
    type = "done",
    duration = step_times$save_sec
  )
  log("FIGURE WORKFLOW SUMMARY", type = "info")
  log(sprintf("  Datasets: %d | GnRH-positive cells: %s", nrow(datasets_sel), format(ncol(combined_gnrh), big.mark = ",")), type = "info")
  log(sprintf("  Figure-workflow RDS written: %s", if (isTRUE(save_result)) "yes" else "no"), type = "info")
  if (isTRUE(save_result)) log(sprintf("  Seurat objects included in RDS: %s", if (isTRUE(save_objects)) "yes" else "no"), type = "info")
  log(paste0("  Output: ", result$output_dir), type = "info")
  log("==== GnRHcell FIGURE WORKFLOW COMPLETE ====", type = "done", duration = step_times$total_sec)
  invisible(result)
}
