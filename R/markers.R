# =========================================================
# markers-core.R
# Clean + factorized GnRH marker framework
# =========================================================

# ---------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------

.validate_marker_columns <- function(df, required) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
}


.validate_groups <- function(cells_1, cells_2, min_cells = 5) {

  if (length(cells_1) < min_cells) {
    stop("Too few cells in ident.1", call. = FALSE)
  }

  if (length(cells_2) < min_cells) {
    stop("Too few cells in ident.2", call. = FALSE)
  }

  invisible(TRUE)
}




# ---------------------------------------------------------
# Differential expression
# ---------------------------------------------------------

#' Run marker differential expression
#'
#' @keywords internal
#' @noRd
run_marker_de <- function(
    obj,
    ident.1,
    ident.2 = NULL,
    assay = "RNA",
    layer = "data",
    test.use = "wilcox"
) {

  mk <- Seurat::FindMarkers(
    object = obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    assay = assay,
    slot = layer,
    test.use = test.use,
    logfc.threshold = 0,
    min.pct = 0.01
  )

  mk$gene <- rownames(mk)
  mk$method <- test.use

  mk
}


run_multi_de <- function(
    obj,
    methods,
    ident.1,
    ident.2 = NULL,
    assay = "RNA",
    layer = "data",
    verbose = TRUE
) {

  log <- log_msg(verbose)

  res <- lapply(methods, function(meth) {

    log("Running DE:", meth)

    run_marker_de(
      obj = obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      assay = assay,
      layer = layer,
      test.use = meth
    )
  })

  do.call(rbind, res)
}


# ---------------------------------------------------------
# Coexpression
# ---------------------------------------------------------

.compute_safe_cor <- function(x, y, method = "pearson") {

  out <- suppressWarnings(
    stats::cor(
      x,
      y,
      method = method
    )
  )

  ifelse(is.finite(out), out, 0)
}


compute_coexpression_pearson <- function(X, g) {

  apply(
    as.matrix(X),
    1,
    .compute_safe_cor,
    y = g,
    method = "pearson"
  )
}


compute_coexpression_spearman <- function(X, g) {

  apply(
    as.matrix(X),
    1,
    .compute_safe_cor,
    y = g,
    method = "spearman"
  )
}


compute_coexpression_cosine <- function(X, g) {

  num <- as.vector(X %*% g)

  denom <- sqrt(Matrix::rowSums(X^2)) *
    sqrt(sum(g^2))

  out <- num / denom

  out[!is.finite(out)] <- 0

  out
}


compute_coexpression <- function(
    X,
    g,
    method = c("spearman", "pearson", "cosine")
) {

  method <- match.arg(method)

  switch(
    method,

    spearman =
      compute_coexpression_spearman(X, g),

    pearson =
      compute_coexpression_pearson(X, g),

    cosine =
      compute_coexpression_cosine(X, g)
  )
}


# ---------------------------------------------------------
# Specificity
# ---------------------------------------------------------

#' Compute gene specificity between two cell populations
#'
#' Calculates detection frequency differences between two cell groups.
#'
#' @param expr Sparse or dense expression matrix (genes × cells)
#' @param genes Character vector of genes to evaluate
#' @param cells_1 Cell IDs for group 1
#' @param cells_2 Cell IDs for group 2
#'
#' @return A data.frame with gene-level specificity metrics:
#' pct1, pct2, and spec (difference)
#'
#' @export
compute_specificity <- function(
    expr,
    genes,
    cells_1,
    cells_2
) {

  genes <- intersect(genes, rownames(expr))

  m1 <- expr[genes, cells_1, drop = FALSE]
  m2 <- expr[genes, cells_2, drop = FALSE]

  pct1 <- Matrix::rowSums(m1 > 0) / ncol(m1)
  pct2 <- Matrix::rowSums(m2 > 0) / ncol(m2)

  data.frame(
    gene = genes,
    pct1 = as.numeric(pct1),
    pct2 = as.numeric(pct2),
    spec = as.numeric(pct1 - pct2),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------
# Scoring
# ---------------------------------------------------------

compute_marker_score <- function(df) {

  df$score <- with(
    df,
    avg_log2FC *
      pmax(spec, 0) *
      sqrt(pct1) *
      coexpr *
      log1p(method_support)
  )

  df
}


filter_marker_table <- function(
    df,
    min_pct = 0.01,
    min_logfc = 0.1,
    max_padj = 0.05,
    coexpr_thresh = 0.3,
    coexpr_filter = c("none", "only", "annotate")
) {

  coexpr_filter <- match.arg(coexpr_filter)

  keep <- (
    df$pct1 >= min_pct &
      df$avg_log2FC >= min_logfc &
      df$p_val_adj <= max_padj
  )

  res <- df[keep, , drop = FALSE]

  res$coexpressed_flag <-
    res$coexpr >= coexpr_thresh

  if (coexpr_filter == "only") {

    res <- res[
      res$coexpressed_flag,
      ,
      drop = FALSE
    ]
  }

  res
}


rank_markers <- function(df) {

  df <- df[!duplicated(df$gene), , drop = FALSE]

  df <- df[order(-df$score), , drop = FALSE]

  rownames(df) <- df$gene

  df
}


# ---------------------------------------------------------
# Utilities
# ---------------------------------------------------------

detect_gnrh_gene <- function(expr) {

  candidates <- c(
    "GNRH1",
    "Gnrh1",
    "gnrh1",
    "LHRH"
  )

  hit <- rownames(expr)[
    toupper(rownames(expr)) %in%
      toupper(candidates)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


prepare_marker_groups <- function(
    obj,
    group.by,
    ident.1,
    ident.2 = NULL
) {

  meta <- obj@meta.data

  group_vec <- meta[[group.by]]

  names(group_vec) <- colnames(obj)

  cells_1 <- names(group_vec)[group_vec == ident.1]

  cells_2 <- if (is.null(ident.2)) {
    names(group_vec)[group_vec != ident.1]
  } else {
    names(group_vec)[group_vec == ident.2]
  }

  list(
    group_vec = group_vec,
    cells_1 = cells_1,
    cells_2 = cells_2
  )
}


# ---------------------------------------------------------
# Aggregation
# ---------------------------------------------------------

#' Aggregate differential expression results across methods
#'
#' Combines DE results from multiple statistical methods into
#' a single consensus table per gene.
#'
#' @param markers_all Data.frame containing DE results with columns:
#' gene, avg_log2FC, p_val_adj, method
#'
#' @return A data.frame with aggregated gene-level statistics.
#'
#' @export
aggregate_marker_de <- function(markers_all) {

  required <- c("gene", "avg_log2FC", "p_val_adj", "method")

  .validate_marker_columns(markers_all, required)

  markers_all <- markers_all[!is.na(markers_all$gene), , drop = FALSE]

  split_markers <- split(markers_all, markers_all$gene)

  out <- lapply(split_markers, function(df) {
    data.frame(
      gene = df$gene[1],
      avg_log2FC = mean(df$avg_log2FC, na.rm = TRUE),
      p_val_adj = min(df$p_val_adj, na.rm = TRUE),
      method_support = length(unique(df$method)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}



# ---------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------

#' GnRH marker detection pipeline
#'
#' Runs a full differential expression + specificity + coexpression
#' scoring pipeline for identifying GnRH-related marker genes in a
#' Seurat object.
#'
#' The pipeline performs:
#' \itemize{
#'   \item Differential expression across groups
#'   \item Gene filtering based on detection thresholds
#'   \item Specificity scoring between cell populations
#'   \item Coexpression with reference GnRH gene
#'   \item Composite marker scoring and ranking
#' }
#'
#' @param obj A Seurat object.
#' @param group.by Metadata column defining groups.
#' @param ident.1 Identity/class of target group.
#' @param ident.2 Optional second group; defaults to all others.
#' @param assay Assay to use (default "RNA").
#' @param layer Data layer in Seurat assay (default "data").
#' @param min_cells Minimum number of cells required per group.
#' @param gnrh_gene Optional reference GnRH gene. If NULL, auto-detected.
#' @param coexpr_thresh Coexpression threshold for filtering.
#' @param min_detect Minimum detection count per gene.
#' @param coexpr.method Method for coexpression ("spearman", "cosine", "pearson").
#' @param methods Character vector of DE methods (e.g. "wilcox").
#' @param coexpr_filter How to handle coexpression filtering:
#'   "none", "only", or "annotate".
#' @param verbose Logical; print progress messages.
#'
#' @return A data.frame of ranked marker genes.
#'
#' @export
gnrh_markers <- function(
    obj,
    group.by = "gnrh_status",
    ident.1 = TRUE,
    ident.2 = NULL,
    assay = "RNA",
    layer = "data",
    min_cells = 5,
    gnrh_gene = NULL,
    coexpr_thresh = 0.3,
    min_detect = 10,
    coexpr.method = c("spearman", "cosine", "pearson"),
    methods = c("wilcox"),
    coexpr_filter = c("none", "only", "annotate"),
    verbose = TRUE
) {

  stopifnot(inherits(obj, "Seurat"))

  coexpr.method <- match.arg(coexpr.method)

  coexpr_filter <- match.arg(coexpr_filter)

  log <- log_msg(verbose)

  log("=== GNRH MARKER PIPELINE START ===")

  # -------------------------------------------------------
  # Expression matrix
  # -------------------------------------------------------

  expr <- SeuratObject::GetAssayData(
    obj,
    assay = assay,
    layer = layer
  )

  if (!inherits(expr, "dgCMatrix")) {
    expr <- methods::as(expr, "dgCMatrix")
  }

  expr <- expr[, colnames(obj), drop = FALSE]

  rownames(expr) <- make.unique(rownames(expr))

  # -------------------------------------------------------
  # Groups
  # -------------------------------------------------------

  grp <- prepare_marker_groups(
    obj = obj,
    group.by = group.by,
    ident.1 = ident.1,
    ident.2 = ident.2
  )

  .validate_groups(
    grp$cells_1,
    grp$cells_2,
    min_cells = min_cells
  )

  SeuratObject::Idents(obj) <- grp$group_vec

  # -------------------------------------------------------
  # GNRH gene
  # -------------------------------------------------------

  if (is.null(gnrh_gene)) {
    gnrh_gene <- detect_gnrh_gene(expr)
  }

  if (is.na(gnrh_gene)) {
    stop("No GNRH1-like gene found", call. = FALSE)
  }

  log("Reference gene:", gnrh_gene)

  # -------------------------------------------------------
  # DE
  # -------------------------------------------------------

  markers_all <- run_multi_de(
    obj = obj,
    methods = methods,
    ident.1 = ident.1,
    ident.2 = ident.2,
    assay = assay,
    layer = layer,
    verbose = verbose
  )

  de <- aggregate_marker_de(markers_all)

  # -------------------------------------------------------
  # Detection filter
  # -------------------------------------------------------

  keep <- Matrix::rowSums(expr != 0) >= min_detect

  expr <- expr[keep, , drop = FALSE]

  de <- de[
    de$gene %in% rownames(expr),
    ,
    drop = FALSE
  ]

  # -------------------------------------------------------
  # Specificity
  # -------------------------------------------------------

  spec <- compute_specificity(
    expr = expr,
    genes = de$gene,
    cells_1 = grp$cells_1,
    cells_2 = grp$cells_2
  )

  idx <- match(de$gene, spec$gene)

  de$pct1 <- spec$pct1[idx]
  de$pct2 <- spec$pct2[idx]
  de$spec <- spec$spec[idx]

  # -------------------------------------------------------
  # Coexpression
  # -------------------------------------------------------

  g <- as.numeric(
    expr[gnrh_gene, grp$cells_1]
  )

  X <- expr[, grp$cells_1, drop = FALSE]

  coexpr <- compute_coexpression(
    X = X,
    g = g,
    method = coexpr.method
  )

  de$coexpr <- coexpr[de$gene]

  # -------------------------------------------------------
  # Score + filter + rank
  # -------------------------------------------------------

  de <- compute_marker_score(de)

  de <- filter_marker_table(
    df = de,
    coexpr_thresh = coexpr_thresh,
    coexpr_filter = coexpr_filter
  )

  de <- rank_markers(de)

  # -------------------------------------------------------
  # Cleanup
  # -------------------------------------------------------

  de <- de[
    !is.na(de$gene) &
      nzchar(de$gene) &
      is.finite(de$score),
    ,
    drop = FALSE
  ]

  rownames(de) <- make.unique(as.character(de$gene))

  # -------------------------------------------------------
  # Summary
  # -------------------------------------------------------

  summary_list <- list(
    n_genes = nrow(de),
    top_gene = if (nrow(de) > 0) de$gene[1] else NA_character_,
    mean_coexpr = mean(de$coexpr, na.rm = TRUE),
    mean_score = mean(de$score, na.rm = TRUE)
  )

  params_list <- list(
    assay = assay,
    layer = layer,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min_cells = min_cells,
    min_detect = min_detect,
    coexpr_thresh = coexpr_thresh,
    coexpr_filter = coexpr_filter,
    coexpr_method = coexpr.method
  )

  res <- new_gnrh_markers(
    markers = de,
    params = params_list,
    summary = summary_list,
    gene = gnrh_gene,
    method = coexpr.method
  )

  obj@misc$gnrh_markers <- res

  log("Final markers:", nrow(res$markers))

  log("=== DONE ===")

  res$markers
}


# ---------------------------------------------------------
# S3
# ---------------------------------------------------------

new_gnrh_markers <- function(
    markers,
    params,
    summary,
    gene,
    method
) {

  structure(
    list(
      markers = markers,
      params = params,
      summary = summary,
      gene = gene,
      method = method
    ),
    class = "gnrh_markers"
  )
}


#' Summarize GnRH marker results
#'
#' @param object A gnrh_markers object
#' @param n Number of top genes to display
#' @param ... Additional arguments (unused)
#'
#' @return A data.frame of top-ranked markers
#'
#' @export
summary.gnrh_markers <- function(object, n = 10, ...) {

  head(
    object$markers[
      order(-object$markers$score),
      c(
        "gene",
        "avg_log2FC",
        "coexpr",
        "score"
      )
    ],
    n
  )
}


#' Print GnRH marker results
#'
#' @param x A gnrh_markers object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns input object
#'
#' @export
print.gnrh_markers <- function(x, ...) {

  cat("\n")
  cat("GnRH Marker Analysis\n")
  cat("--------------------\n")

  cat("Reference gene :", x$gene, "\n")
  cat("Coexpression   :", x$method, "\n")
  cat("Markers        :", nrow(x$markers), "\n")

  if (!is.null(x$summary$top_gene)) {
    cat("Top marker     :", x$summary$top_gene, "\n")
  }

  cat("\n")

  invisible(x)
}
