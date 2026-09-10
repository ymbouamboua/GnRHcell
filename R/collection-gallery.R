.build_collection_gallery <- function(objects, datasets, output_dir) {
  # ========================================================================= #
  # Validation
  # ========================================================================= #
  required <- c("id", "label", "species", "reduction")
  missing <- setdiff(required, colnames(datasets))
  if (length(missing)) {
    stop("Missing gallery dataset columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ids <- as.character(datasets$id)
  if (is.null(names(objects))) names(objects) <- ids
  missing_objects <- setdiff(ids, names(objects))
  if (length(missing_objects)) {
    stop("Missing processed objects for gallery: ", paste(missing_objects, collapse = ", "), call. = FALSE)
  }
  objects <- objects[ids]
  datasets <- datasets[match(ids, as.character(datasets$id)), , drop = FALSE]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # ========================================================================= #
  # Helpers
  # ========================================================================= #
  save_gallery_plot <- function(plot, base, width, height, dpi = 320) {
    files <- c(png = paste0(base, ".png"), pdf = paste0(base, ".pdf"))
    ggplot2::ggsave(files[["png"]], plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(files[["pdf"]], plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
    files
  }
  gallery_layout <- function(n) {
    ncol <- min(3L, max(1L, n))
    list(
      ncol = ncol,
      nrow = ceiling(n / ncol),
      width = max(5.5, 4.8 * ncol),
      height = max(4.8, 4.6 * ceiling(n / ncol))
    )
  }
  # ========================================================================= #
  # Collection summary
  # ========================================================================= #
  summary_table <- do.call(
    rbind,
    lapply(seq_along(objects), function(i) {
      md <- objects[[i]][[]]
      positive <- as.character(md$gnrh_status) == "pos"
      data.frame(
        dataset = as.character(datasets$label[[i]]),
        species = as.character(datasets$species[[i]]),
        cells = nrow(md),
        gnrh_detected = sum(positive, na.rm = TRUE),
        detected_percent = 100 * mean(positive, na.rm = TRUE),
        early = sum(positive & as.character(md$gnrh_stage) == "early", na.rm = TRUE),
        migrating = sum(positive & as.character(md$gnrh_stage) == "migrating", na.rm = TRUE),
        post_migratory = sum(positive & as.character(md$gnrh_stage) == "post-migratory", na.rm = TRUE),
        mature = sum(positive & as.character(md$gnrh_stage) == "mature", na.rm = TRUE),
        transitional = sum(positive & as.character(md$gnrh_stage) == "transitional", na.rm = TRUE),
        secretory_supported = if ("gnrh_secretory" %in% colnames(md)) sum(positive & as.character(md$gnrh_secretory) == "supported", na.rm = TRUE) else NA_integer_,
        secretory_limited = if ("gnrh_secretory" %in% colnames(md)) sum(positive & as.character(md$gnrh_secretory) == "limited", na.rm = TRUE) else NA_integer_,
        check.names = FALSE
      )
    })
  )
  rownames(summary_table) <- NULL
  summary_file <- file.path(output_dir, "collection-results.csv")
  utils::write.csv(summary_table, summary_file, row.names = FALSE)
  # ========================================================================= #
  # Status collection map
  # ========================================================================= #
  status_plots <- lapply(seq_along(objects), function(i) {
    object <- objects[[i]]
    md <- object[[]]
    reduction_i <- resolve_reduction(object, as.character(datasets$reduction[[i]]))
    n_neg <- sum(as.character(md$gnrh_status) == "neg", na.rm = TRUE)
    n_pos <- sum(as.character(md$gnrh_status) == "pos", na.rm = TRUE)
    status_labels <- c(
      neg = paste0("Other cells (n = ", format(n_neg, big.mark = ",", trim = TRUE), ")"),
      pos = paste0("GnRH detected (n = ", format(n_pos, big.mark = ",", trim = TRUE), ")")
    )
    plot_gnrh_embedding(
      object = object,
      group_by = "gnrh_status",
      reduction = reduction_i,
      n.cells = FALSE,
      percentage = FALSE,
      axes = TRUE,
      pt.size = 1,
      alpha = 1,
      background_alpha = 1,
      plot.ttl = as.character(datasets$label[[i]]),
      cols = gnrh_colors("status"),
      leg.pos = "bottom",
      leg.dir = "horizontal",
      leg.ncol = 2,
      leg.size = 8,
      item.size = 2.5,
      item.border = FALSE,
      txtsize = 10,
      style = "test"
    ) +
      ggplot2::scale_colour_manual(
        values = gnrh_colors("status"),
        breaks = c("neg", "pos"),
        labels = status_labels,
        drop = FALSE
      ) +
      ggplot2::labs(colour = NULL) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          title = NULL,
          nrow = 1,
          byrow = TRUE,
          override.aes = list(size = 3, alpha = 1)
        )
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 5)),
        aspect.ratio = 1,
        axis.title = ggplot2::element_text(size = 8),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        legend.spacing.x = grid::unit(4, "pt"),
        legend.key.width = grid::unit(10, "pt"),
        legend.key.height = grid::unit(10, "pt"),
        legend.margin = ggplot2::margin(t = 5, r = 0, b = 2, l = 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = ggplot2::margin(4, 4, 4, 4)
      )
  })
  status_plots <- Filter(Negate(is.null), status_plots)
  status_layout <- gallery_layout(length(status_plots))
  umap_plot <- patchwork::wrap_plots(status_plots, ncol = status_layout$ncol) +
    patchwork::plot_annotation(
      title = paste0(
        "GnRHcell detection across ", length(status_plots), " dataset",
        if (length(status_plots) == 1L) "" else "s"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5, margin = ggplot2::margin(b = 6))
      )
    )
  umap_files <- save_gallery_plot(
    umap_plot,
    file.path(output_dir, "collection-umap-status"),
    width = status_layout$width,
    height = status_layout$height + 0.35
  )
  # ========================================================================= #
  # Developmental-stage collection map
  # ========================================================================= #
  stage_col <- "gnrh_stage"
  stage_levels <- c(
    "early", "migrating", "post-migratory", "mature", "transitional"
  )
  stage_palette <- gnrh_colors("stage")[stage_levels]
  stage_display <- c(
    early = "Early",
    migrating = "Migrating",
    "post-migratory" = "Post-migratory",
    mature = "Mature",
    transitional = "Transitional"
  )

  stage_plots <- lapply(seq_along(objects), function(i) {
    object <- objects[[i]]
    md <- object[[]]
    if (!stage_col %in% colnames(md)) return(NULL)
    reduction_i <- resolve_reduction(object, as.character(datasets$reduction[[i]]))
    stage_counts <- table(
      factor(
        as.character(md[[stage_col]]),
        levels = c("non-gnrh", stage_levels)
      )
    )
    stage_labels <- stats::setNames(
      paste0(
        unname(stage_display[stage_levels]),
        " (n = ",
        format(as.integer(stage_counts[stage_levels]), big.mark = ",", trim = TRUE),
        ")"
      ),
      stage_levels
    )
    plot_gnrh_embedding(
      object = object,
      group_by = stage_col,
      reduction = reduction_i,
      n.cells = FALSE,
      percentage = FALSE,
      axes = TRUE,
      pt.size = 1,
      alpha = 1,
      background_alpha = 1,
      plot.ttl = as.character(datasets$label[[i]]),
      cols = stage_palette,
      leg.pos = "bottom",
      leg.dir = "horizontal",
      leg.ncol = min(5L, length(stage_levels)),
      leg.size = 8,
      item.size = 2.5,
      item.border = FALSE,
      txtsize = 10,
      style = "test"
    ) +
      ggplot2::scale_colour_manual(
        values = stage_palette,
        breaks = stage_levels,
        labels = stage_labels,
        drop = FALSE
      ) +
      ggplot2::labs(colour = NULL) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          title = NULL,
          nrow = 1,
          byrow = TRUE,
          override.aes = list(size = 3, alpha = 1)
        )
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 5)),
        aspect.ratio = 1,
        axis.title = ggplot2::element_text(size = 8),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        legend.spacing.x = grid::unit(4, "pt"),
        legend.key.width = grid::unit(10, "pt"),
        legend.key.height = grid::unit(10, "pt"),
        legend.margin = ggplot2::margin(t = 5, r = 0, b = 2, l = 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = ggplot2::margin(4, 4, 4, 4)
      )
  })
  stage_plots <- Filter(Negate(is.null), stage_plots)
  stage_map_plot <- NULL
  stage_map_files <- NULL
  if (length(stage_plots)) {
    stage_layout <- gallery_layout(length(stage_plots))
    stage_map_plot <- patchwork::wrap_plots(stage_plots, ncol = stage_layout$ncol) +
      patchwork::plot_annotation(
        title = paste0(
          "GnRH developmental stages across ", length(stage_plots), " dataset",
          if (length(stage_plots) == 1L) "" else "s"
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5, margin = ggplot2::margin(b = 6))
        )
      )
    stage_map_files <- save_gallery_plot(
      stage_map_plot,
      file.path(output_dir, "collection-map-stage"),
      width = stage_layout$width,
      height = stage_layout$height + 0.35
    )
  }
  # ========================================================================= #
  # Developmental-stage composition
  # ========================================================================= #
  stage_data <- do.call(
    rbind,
    lapply(seq_along(objects), function(i) {
      md <- objects[[i]][[]]
      md <- md[as.character(md$gnrh_status) == "pos", , drop = FALSE]
      counts <- table(factor(as.character(md[[stage_col]]), levels = stage_levels))
      total <- sum(counts)
      data.frame(
        dataset = factor(as.character(datasets$label[[i]]), levels = rev(as.character(datasets$label))),
        stage = factor(names(counts), levels = stage_levels),
        cells = as.integer(counts),
        fraction = if (total > 0L) as.integer(counts) / total else 0
      )
    })
  )
  n_datasets <- length(objects)
  stage_bar_width <- if (n_datasets <= 2L) {
    0.42
  } else if (n_datasets <= 4L) {
    0.52
  } else if (n_datasets <= 8L) {
    0.62
  } else {
    0.72
  }
  stage_plot <- ggplot2::ggplot(
    stage_data,
    ggplot2::aes(.data$dataset, .data$fraction, fill = .data$stage)
  ) +
    ggplot2::geom_col(
      width = stage_bar_width,
      colour = "white",
      linewidth = 0.25,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = stage_palette,
      labels = stage_display,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.005, 0.025))
    ) +
    ggplot2::labs(
      title = "Developmental-stage composition of detected cells",
      x = NULL,
      y = "Detected cells",
      fill = "Stage"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.just = "center",
      plot.margin = ggplot2::margin(8, 12, 6, 8)
    )
  # ========================================================================= #
  # Save stage composition
  # ========================================================================= #
  stage_width <- if (n_datasets <= 3L) 7 else 8.5
  stage_height <- max(
    3.2,
    min(
      8,
      2.25 + 0.58 * n_datasets
    )
  )
  stage_files <- save_gallery_plot(
    stage_plot,
    file.path(output_dir, "collection-stage-composition"),
    width = stage_width,
    height = stage_height
  )
  # ========================================================================= #
  # Return
  # ========================================================================= #
  list(
    summary = summary_table,
    summary_file = summary_file,
    umap_plot = umap_plot,
    umap_file = umap_files,
    stage_map_plot = stage_map_plot,
    stage_map_file = stage_map_files,
    stage_plot = stage_plot,
    stage_file = stage_files
  )
}
# ========================================================================= #
# Load saved gallery objects
# ========================================================================= #
.load_collection_gallery_objects <- function(
    datasets,
    output_dir
) {
  paths <- file.path(
    output_dir,
    "objects",
    paste0(datasets$id, "_gnrh.rds")
  )
  if (!all(file.exists(paths))) {
    return(NULL)
  }
  stats::setNames(
    lapply(paths, readRDS),
    datasets$id
  )
}
