# Build the static demo-result assets used by README.Rmd and pkgdown.
# Run from the package root with:
# Rscript data-raw/build-demo-gallery.R

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Packages `ggplot2` and `patchwork` are required.", call. = FALSE)
}

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

dataset_labels <- c(
  hpsc = "hPSC-derived GnRH",
  hudeca_nose = "Human fetal nose",
  human_me = "Human median eminence",
  mouse_hypomap = "Mouse HypoMap",
  human_hypomap = "Human HypoMap"
)

objects <- lapply(
  names(dataset_labels),
  function(id) readRDS(file.path("inst", "extdata", paste0(id, ".rds")))
)
names(objects) <- names(dataset_labels)

summarise_object <- function(object, id) {
  md <- object[[]]
  status <- table(md$gnrh_status)
  positive <- if ("pos" %in% names(status)) as.integer(status[["pos"]]) else 0L

  data.frame(
    dataset = unname(dataset_labels[[id]]),
    species = if (grepl("Mouse", dataset_labels[[id]])) "Mouse" else "Human",
    cells = nrow(md),
    gnrh_detected = positive,
    detected_in_demo = 100 * positive / nrow(md),
    identity = sum(md$gnrh_stage == "identity", na.rm = TRUE),
    migrating = sum(md$gnrh_stage == "migrating", na.rm = TRUE),
    mature = sum(md$gnrh_stage == "mature", na.rm = TRUE),
    secreting = sum(md$gnrh_stage == "secreting", na.rm = TRUE),
    check.names = FALSE
  )
}

summary_table <- do.call(
  rbind,
  Map(summarise_object, objects, names(objects))
)
rownames(summary_table) <- NULL

utils::write.csv(
  summary_table,
  "man/figures/demo-results.csv",
  row.names = FALSE
)

extract_umap <- function(object, id) {
  coords <- as.data.frame(object@reductions$umap@cell.embeddings[, 1:2, drop = FALSE])
  names(coords) <- c("UMAP_1", "UMAP_2")
  md <- object[[]]
  coords$status <- factor(md$gnrh_status, levels = c("neg", "pos"))
  coords$dataset <- factor(
    unname(dataset_labels[[id]]),
    levels = unname(dataset_labels)
  )
  coords
}

umap <- do.call(rbind, Map(extract_umap, objects, names(objects)))
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
    title = "GnRHcell detection across five demo datasets",
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
  "man/figures/demo-umap-status.png",
  p_umap,
  width = 10,
  height = 6.6,
  dpi = 320,
  bg = "white"
)

stage_levels <- c("identity", "migrating", "mature", "secreting")
stage_labels <- c(
  identity = "Identity",
  migrating = "Migrating",
  mature = "Mature",
  secreting = "Secreting"
)

stage_data <- do.call(rbind, lapply(names(objects), function(id) {
  md <- objects[[id]][[]]
  md <- md[md$gnrh_status == "pos", , drop = FALSE]
  counts <- table(factor(md$gnrh_stage, levels = stage_levels))
  data.frame(
    dataset = factor(unname(dataset_labels[[id]]), levels = rev(unname(dataset_labels))),
    stage = factor(names(counts), levels = stage_levels),
    cells = as.integer(counts),
    fraction = as.integer(counts) / max(1L, nrow(md))
  )
}))

p_stage <- ggplot2::ggplot(
  stage_data,
  ggplot2::aes(dataset, fraction, fill = stage)
) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c(
      identity = "#3B4CC0",
      migrating = "#00A6CA",
      mature = "#F28E2B",
      secreting = "#C51B7D"
    ),
    labels = stage_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Developmental-stage composition of detected cells",
    subtitle = "Fractions are calculated within GnRH-detected cells in each demo object",
    x = NULL,
    y = "Detected cells",
    fill = "Stage"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  "man/figures/demo-stage-composition.png",
  p_stage,
  width = 9,
  height = 5.2,
  dpi = 320,
  bg = "white"
)
