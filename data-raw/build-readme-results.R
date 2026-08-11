# Build the static, reproducible results displayed after README workflow steps.
# Run from the package root with:
# Rscript data-raw/build-readme-results.R

required <- c("devtools", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

hpsc <- readRDS("inst/extdata/hpsc.rds")

run_log <- utils::capture.output(
  processed <- run_gnrh(hpsc, verbose = TRUE),
  type = "output"
)
writeLines(run_log, "man/figures/hpsc-run-log.txt")

count_table <- function(x, column) {
  out <- as.data.frame(table(x[[]][[column]], useNA = "ifany"))
  names(out) <- c(column, "cells")
  out
}

utils::write.csv(
  count_table(processed, "gnrh_status"),
  "man/figures/hpsc-status-summary.csv",
  row.names = FALSE
)
utils::write.csv(
  count_table(processed, "gnrh_stage"),
  "man/figures/hpsc-stage-summary.csv",
  row.names = FALSE
)

p_status <- plot_gnrh_embedding(
  processed,
  group_by = "gnrh_status",
  reduction = "umap",
  total.cells = TRUE,
  style = "classic"
) + ggplot2::labs(title = "Detection status")

p_stage <- plot_gnrh_embedding(
  processed,
  group_by = "gnrh_stage",
  reduction = "umap",
  percentage = TRUE,
  cols = gnrh_colors("stage"),
  style = "classic"
) + ggplot2::labs(title = "Developmental stage")

ggplot2::ggsave(
  "man/figures/hpsc-embedding-results.png",
  p_status + p_stage + patchwork::plot_annotation(
    title = "GnRHcell output in the hPSC demo object"
  ),
  width = 11,
  height = 5.4,
  dpi = 320,
  bg = "white"
)

p_feature <- plot_gnrh_feature(
  processed,
  preset = "core",
  reduction = "umap",
  ncol = 2,
  style = "classic"
)

ggplot2::ggsave(
  "man/figures/hpsc-feature-core.png",
  p_feature,
  width = 10,
  height = 7.2,
  dpi = 320,
  bg = "white"
)

p_dot <- plot_gnrh_dot(
  processed,
  features = c("GNRH1", "FEZF1", "ISL1", "ANOS1", "PROKR2", "DCX", "CHGA", "PCSK2"),
  group.by = "gnrh_stage",
  dot.outline = TRUE,
  th.cols = "RdYlBu",
  title = "GnRH developmental programs"
)

ggplot2::ggsave(
  "man/figures/hpsc-program-dotplot.png",
  p_dot,
  width = 9.5,
  height = 5.8,
  dpi = 320,
  bg = "white"
)

p_distribution <- plot_gnrh_distribution(
  processed,
  group.by = "gnrh_stage",
  split.by = "orig.ident",
  proportion = TRUE,
  label = FALSE,
  cols = gnrh_colors("stage")
) + ggplot2::labs(title = "Stage composition by sample")

ggplot2::ggsave(
  "man/figures/hpsc-stage-distribution.png",
  p_distribution,
  width = 9.5,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

md <- processed[[]]
p_score <- ggplot2::ggplot(
  md,
  ggplot2::aes(gnrh_status, gnrh_score, fill = gnrh_status)
) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.8) +
  ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white") +
  ggplot2::scale_fill_manual(values = gnrh_colors("status"), guide = "none") +
  ggplot2::labs(title = "Detection-score separation", x = NULL, y = "GnRH score") +
  ggplot2::theme_classic(base_size = 11)

class_counts <- as.data.frame(table(md$gnrh_class))
names(class_counts) <- c("class", "cells")
p_class <- ggplot2::ggplot(class_counts, ggplot2::aes(class, cells, fill = class)) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(ggplot2::aes(label = cells), vjust = -0.35, size = 3.4) +
  ggplot2::scale_fill_manual(values = gnrh_colors("class"), guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(title = "Detection classes", x = NULL, y = "Cells") +
  ggplot2::theme_classic(base_size = 11)

ggplot2::ggsave(
  "man/figures/hpsc-diagnostic-summary.png",
  p_score + p_class + patchwork::plot_annotation(title = "Internal diagnostic summary"),
  width = 9.5,
  height = 4.8,
  dpi = 320,
  bg = "white"
)

lineage_markers <- find_gnrh_genes(
  processed,
  annotation_col = "ann2",
  control_ident = "GLU",
  donor_col = "orig.ident",
  assay = "RNA",
  layer = "data",
  verbose = FALSE
)

utils::write.csv(
  lineage_markers$candidates,
  "man/figures/hpsc-lineage-markers.csv",
  row.names = FALSE
)

p_coexpr <- plot_gnrh_coexpr(
  lineage_markers$markers,
  coexp_cutoff = 0.30
) + ggplot2::labs(title = "GNRH1 co-expression")

p_network <- plot_network(
  lineage_markers$markers,
  top_n = 25,
  threshold = 0.25
)

ggplot2::ggsave(
  "man/figures/hpsc-lineage-marker-results.png",
  p_coexpr + p_network + patchwork::plot_annotation(title = "GnRH-lineage marker results"),
  width = 11,
  height = 5.6,
  dpi = 320,
  bg = "white"
)

stage_markers <- find_gnrh_stage_markers(
  processed,
  min_cells = 10L,
  verbose = FALSE
)

utils::write.csv(
  stage_markers$markers,
  "man/figures/hpsc-stage-markers.csv",
  row.names = FALSE
)

top_markers <- do.call(rbind, lapply(split(stage_markers$markers, stage_markers$markers$stage), function(x) {
  utils::head(x[order(x$p_val_adj, -x$avg_log2FC), ], 4L)
}))
top_genes <- unique(top_markers$gene)

p_markers <- plot_gnrh_dot(
  processed,
  features = top_genes,
  group.by = "gnrh_stage",
  dot.outline = TRUE,
  th.cols = "RdYlBu",
  title = "Top stage-associated markers"
)

ggplot2::ggsave(
  "man/figures/hpsc-stage-marker-dotplot.png",
  p_markers,
  width = 10,
  height = 6,
  dpi = 320,
  bg = "white"
)
