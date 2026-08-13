.build_collection_gallery <- function(objects, datasets, output_dir) {
  required <- c("id", "label", "species", "reduction")
  missing <- setdiff(required, colnames(datasets))
  if (length(missing)) stop("Missing gallery dataset columns: ", paste(missing, collapse = ", "), call. = FALSE)
  ids <- as.character(datasets$id)
  if (is.null(names(objects))) names(objects) <- ids
  missing_objects <- setdiff(ids, names(objects))
  if (length(missing_objects)) stop("Missing processed objects for gallery: ", paste(missing_objects, collapse = ", "), call. = FALSE)
  objects <- objects[ids]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_table <- do.call(rbind, lapply(seq_along(objects), function(i) {
    md <- objects[[i]][[]]; positive <- sum(md$gnrh_status == "pos", na.rm = TRUE)
    data.frame(dataset = as.character(datasets$label[[i]]), species = as.character(datasets$species[[i]]), cells = nrow(md), gnrh_detected = positive, detected_percent = 100 * positive / nrow(md), identity = sum(md$gnrh_stage == "identity", na.rm = TRUE), migrating = sum(md$gnrh_stage == "migrating", na.rm = TRUE), mature = sum(md$gnrh_stage == "mature", na.rm = TRUE), check.names = FALSE)
  }))
  rownames(summary_table) <- NULL
  summary_file <- file.path(output_dir, "collection-results.csv")
  utils::write.csv(summary_table, summary_file, row.names = FALSE)
  embedding_plots <- lapply(seq_along(objects), function(i) {
    plot_gnrh_embedding(
      object = objects[[i]],
      group_by = "gnrh_status",
      reduction = as.character(datasets$reduction[[i]]),
      n.cells = TRUE,
      percentage = FALSE,
      axes = TRUE,
      pt.size = 0.1,
      plot.ttl = as.character(datasets$label[[i]]),
      cols = c(neg = "#D5D8DC", pos = "#D73027"),
      leg.ttl = NULL,
      leg.pos = "bottom",
      leg.dir = "horizontal",
      leg.ncol = 2,
      leg.size = 7,
      item.size = 2.7,
      item.border = FALSE,
      txtsize = 10,
      style = "test"
    ) +
      ggplot2::scale_colour_manual(
        values = c(neg = "#D5D8DC", pos = "#D73027"),
        breaks = c("neg", "pos"),
        labels = c(
          neg = paste0("Other cells\nn = ", format(sum(objects[[i]]$gnrh_status == "neg", na.rm = TRUE), big.mark = ",")),
          pos = paste0("GnRH detected\nn = ", format(sum(objects[[i]]$gnrh_status == "pos", na.rm = TRUE), big.mark = ","))
        ),
        drop = FALSE
      ) +
      ggplot2::labs(
        x = paste0(toupper(resolve_reduction(objects[[i]], as.character(datasets$reduction[[i]]))), "_1"),
        y = paste0(toupper(resolve_reduction(objects[[i]], as.character(datasets$reduction[[i]]))), "_2"),
        colour = NULL
      ) +
      ggplot2::coord_cartesian() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        aspect.ratio = 1,
        axis.title = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.margin = ggplot2::margin(t = -3),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )
  })
  umap_plot <- patchwork::wrap_plots(embedding_plots, ncol = 3) +
    patchwork::plot_layout(widths = rep(1, 3)) +
    patchwork::plot_annotation(
      title = paste0("GnRHcell detection across ", length(objects), " datasets"),
      #subtitle = "Cell numbers are reported separately for every dataset",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11)
      )
    )
  umap_file <- file.path(output_dir, "collection-umap-status.png")
  ggplot2::ggsave(umap_file, umap_plot, width = 10, height = max(6.6, 2.35 * ceiling(length(objects) / 3)), dpi = 320, bg = "white")
  stage_levels <- c("identity", "migrating", "mature")
  stage_data <- do.call(rbind, lapply(seq_along(objects), function(i) {
    md <- objects[[i]][[]]; md <- md[md$gnrh_status == "pos", , drop = FALSE]; counts <- table(factor(md$gnrh_stage, levels = stage_levels))
    data.frame(dataset = factor(as.character(datasets$label[[i]]), levels = rev(as.character(datasets$label))), stage = factor(names(counts), levels = stage_levels), cells = as.integer(counts), fraction = as.integer(counts) / max(1L, sum(counts)))
  }))
  stage_plot <- ggplot2::ggplot(stage_data, ggplot2::aes(.data$dataset, .data$fraction, fill = .data$stage)) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = c(identity = "#3B4CC0", migrating = "#00A6CA", mature = "#F28E2B"), labels = c(identity = "Identity", migrating = "Migrating", mature = "Mature"), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0.005, 0.025))) +
    ggplot2::labs(title = "Developmental-stage composition of detected cells",
                  #subtitle = "Fractions are calculated within identity, migrating, and mature cells",
                  x = NULL, y = "Detected cells", fill = "Stage") +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), axis.text.y = ggplot2::element_text(face = "bold"), legend.position = "bottom", plot.margin = ggplot2::margin(8, 18, 8, 8))
  stage_file <- file.path(output_dir, "collection-stage-composition.png")
  ggplot2::ggsave(stage_file, stage_plot, width = 9, height = max(5.2, 0.62 * length(objects) + 1.8), dpi = 320, bg = "white")
  list(summary = summary_table, summary_file = summary_file, umap_plot = umap_plot, umap_file = umap_file, stage_plot = stage_plot, stage_file = stage_file)
}

.load_collection_gallery_objects <- function(datasets, output_dir) {
  paths <- file.path(output_dir, "objects", paste0(datasets$id, "_gnrh.rds"))
  if (!all(file.exists(paths))) return(NULL)
  stats::setNames(lapply(paths, readRDS), datasets$id)
}
