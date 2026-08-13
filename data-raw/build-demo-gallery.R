# Build the static result assets used by README.Rmd and pkgdown.
#
# Default bundled demo (run from the package root):
#   Rscript data-raw/build-demo-gallery.R
#
# A processed collection can instead call `build_demo_gallery()` directly;
# see the collection example in README.Rmd.

build_demo_gallery <- function(
    objects,
    datasets,
    output_dir = "man/figures",
    file_prefix = "demo"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.", call. = FALSE)
  }

  required <- c("id", "label", "species", "reduction")
  missing <- setdiff(required, colnames(datasets))
  if (length(missing)) {
    stop(
      "Missing dataset columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  ids <- as.character(datasets$id)
  if (is.null(names(objects))) names(objects) <- ids

  missing_objects <- setdiff(ids, names(objects))
  if (length(missing_objects)) {
    stop(
      "Missing processed objects: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  objects <- objects[ids]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summarise_object <- function(object, i) {
    md <- object[[]]
    positive <- sum(md$gnrh_status == "pos", na.rm = TRUE)

    data.frame(
      dataset = as.character(datasets$label[[i]]),
      species = as.character(datasets$species[[i]]),
      cells = nrow(md),
      gnrh_detected = positive,
      detected_in_demo = 100 * positive / nrow(md),
      identity = sum(md$gnrh_stage == "identity", na.rm = TRUE),
      migrating = sum(md$gnrh_stage == "migrating", na.rm = TRUE),
      mature = sum(md$gnrh_stage == "mature", na.rm = TRUE),
      check.names = FALSE
    )
  }

  summary_table <- do.call(
    rbind,
    Map(summarise_object, objects, seq_along(objects))
  )
  rownames(summary_table) <- NULL

  utils::write.csv(
    summary_table,
    file.path(output_dir, paste0(file_prefix, "-results.csv")),
    row.names = FALSE
  )

  extract_umap <- function(object, i) {
    reduction <- as.character(datasets$reduction[[i]])
    available <- names(object@reductions)

    if (!reduction %in% available) {
      fallback <- intersect(c("umap", "umap_scvi"), available)
      if (!length(fallback)) {
        stop(
          "No configured or fallback UMAP reduction for dataset `",
          ids[[i]], "`.",
          call. = FALSE
        )
      }
      warning(
        "Using reduction `", fallback[[1L]], "` for dataset `", ids[[i]],
        "` instead of missing `", reduction, "`.",
        call. = FALSE
      )
      reduction <- fallback[[1L]]
    }

    coords <- as.data.frame(
      object@reductions[[reduction]]@cell.embeddings[, 1:2, drop = FALSE]
    )
    names(coords) <- c("UMAP_1", "UMAP_2")
    md <- object[[]]
    coords$status <- factor(md$gnrh_status, levels = c("neg", "pos"))
    coords$dataset <- factor(
      as.character(datasets$label[[i]]),
      levels = as.character(datasets$label)
    )
    coords
  }

  umap <- do.call(rbind, Map(extract_umap, objects, seq_along(objects)))
  umap <- umap[order(umap$status), ]
  p_umap <- ggplot2::ggplot(
    umap,
    ggplot2::aes(UMAP_1, UMAP_2, colour = status)
  ) +
    ggplot2::geom_point(size = 0.32, alpha = 0.78) +
    ggplot2::facet_wrap(~dataset, ncol = 3) +
    ggplot2::scale_colour_manual(
      values = c(neg = "#D5D8DC", pos = "#D73027"),
      labels = c(neg = "Other cells", pos = "GnRH detected")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("GnRHcell detection across ", length(objects), " datasets"),
      subtitle = "Red cells satisfy the GnRHcell detection model",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(output_dir, paste0(file_prefix, "-umap-status.png")),
    p_umap,
    width = 10,
    height = max(6.6, 2.35 * ceiling(length(objects) / 3)),
    dpi = 320,
    bg = "white"
  )

  stage_levels <- c("identity", "migrating", "mature")
  stage_labels <- c(
    identity = "Identity", migrating = "Migrating",
    mature = "Mature"
  )

  stage_data <- do.call(rbind, lapply(seq_along(objects), function(i) {
    md <- objects[[i]][[]]
    md <- md[md$gnrh_status == "pos", , drop = FALSE]
    counts <- table(factor(md$gnrh_stage, levels = stage_levels))
    data.frame(
      dataset = factor(
        as.character(datasets$label[[i]]),
        levels = rev(as.character(datasets$label))
      ),
      stage = factor(names(counts), levels = stage_levels),
      cells = as.integer(counts),
      fraction = as.integer(counts) / max(1L, nrow(md))
    )
  }))

  p_stage <- ggplot2::ggplot(
    stage_data,
    ggplot2::aes(dataset, fraction, fill = stage)
  ) +
    ggplot2::geom_col(
      width = 0.72,
      colour = "white",
      linewidth = 0.25,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        identity = "#3B4CC0", migrating = "#00A6CA",
        mature = "#F28E2B"
      ),
      labels = stage_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.005, 0.025))
    ) +
    ggplot2::labs(
      title = "Developmental-stage composition of detected cells",
      subtitle = "Fractions are calculated within GnRH-detected cells in each dataset",
      x = NULL, y = "Detected cells", fill = "Stage"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 18, 8, 8)
    )

  ggplot2::ggsave(
    file.path(output_dir, paste0(file_prefix, "-stage-composition.png")),
    p_stage,
    width = 9,
    height = max(5.2, 0.62 * length(objects) + 1.8),
    dpi = 320,
    bg = "white"
  )

  invisible(list(
    summary = summary_table,
    umap_plot = p_umap,
    stage_plot = p_stage
  ))
}


build_bundled_demo_gallery <- function() {
  datasets <- data.frame(
    id = c("hpsc", "hudeca_nose", "human_me", "mouse_hypomap", "human_hypomap"),
    label = c(
      "hPSC-derived GnRH", "Human fetal nose", "Human median eminence",
      "Mouse HypoMap", "Human HypoMap"
    ),
    species = c("Human", "Human", "Human", "Mouse", "Human"),
    reduction = rep("umap", 5L),
    stringsAsFactors = FALSE
  )
  objects <- setNames(
    lapply(datasets$id, function(id) {
      readRDS(file.path("inst", "extdata", paste0(id, ".rds")))
    }),
    datasets$id
  )
  build_demo_gallery(objects, datasets)
}


if (sys.nframe() == 0L) {
  build_bundled_demo_gallery()
}
