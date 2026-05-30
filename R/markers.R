# =========================================================
# GnRHcell marker discovery
# =========================================================

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

.safe_cor <- function(x, y, method = "spearman") {

  out <- suppressWarnings(
    stats::cor(x, y, method = method)
  )

  ifelse(is.finite(out), out, 0)
}


.find_gnrh_gene <- function(genes) {

  candidates <- c(
    "GNRH1",
    "LHRH"
  )

  hit <- genes[
    toupper(genes) %in%
      toupper(candidates)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


.get_cells <- function(
    object,
    group.by,
    ident.1,
    ident.2 = NULL
) {

  grp <- object[[group.by]][, 1]

  names(grp) <- colnames(object)

  cells_1 <- names(grp)[grp == ident.1]

  cells_2 <- if (is.null(ident.2)) {
    names(grp)[grp != ident.1]
  } else {
    names(grp)[grp == ident.2]
  }

  list(
    cells_1 = cells_1,
    cells_2 = cells_2,
    groups = grp
  )
}


# ---------------------------------------------------------
# Core metrics
# ---------------------------------------------------------

.compute_specificity <- function(
    expr,
    genes,
    cells_1,
    cells_2
) {

  genes <- intersect(genes, rownames(expr))

  m1 <- expr[genes, cells_1, drop = FALSE]
  m2 <- expr[genes, cells_2, drop = FALSE]

  pct1 <- Matrix::rowMeans(m1 > 0)
  pct2 <- Matrix::rowMeans(m2 > 0)

  data.frame(
    gene = genes,
    pct1 = as.numeric(pct1),
    pct2 = as.numeric(pct2),
    spec = as.numeric(pct1 - pct2),
    stringsAsFactors = FALSE
  )
}


.compute_coexpr <- function(
    expr,
    gnrh_gene,
    cells,
    method = "spearman"
) {

  g <- as.numeric(
    expr[gnrh_gene, cells]
  )

  X <- expr[, cells, drop = FALSE]

  apply(
    as.matrix(X),
    1,
    .safe_cor,
    y = g,
    method = method
  )
}


.compute_score <- function(df) {

  df$score <- with(
    df,
    avg_log2FC *
      pmax(spec, 0) *
      sqrt(pmax(pct1, 0)) *
      pmax(coexpr, 0) *
      log1p(method_support)
  )

  df
}


# ---------------------------------------------------------
# DE aggregation
# ---------------------------------------------------------

.aggregate_de <- function(res) {

  split_res <- split(res, res$gene)

  out <- lapply(split_res, function(x) {

    data.frame(
      gene = x$gene[1],
      avg_log2FC = mean(
        x$avg_log2FC,
        na.rm = TRUE
      ),
      p_val_adj = min(
        x$p_val_adj,
        na.rm = TRUE
      ),
      method_support = length(
        unique(x$method)
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}


# ---------------------------------------------------------
# Main function
# ---------------------------------------------------------

#' Detect positive GnRH lineage markers
#'
#' @param object Seurat object
#' @param group.by Metadata column
#' @param ident.1 Positive group
#' @param ident.2 Optional comparison group
#' @param assay Assay
#' @param layer Layer
#' @param methods DE methods
#' @param coexpr.method Correlation method
#' @param min_pct Minimum detection fraction
#' @param min_fc Minimum log2FC
#' @param max_padj Maximum adjusted p-value
#' @param coexpr_min Minimum coexpression
#' @param min_detect Minimum detected cells
#' @param verbose Print progress
#'
#' @return Ranked marker table
#' @export
gnrh_markers <- function(
    object,
    group.by = "gnrh_status",
    ident.1 = "pos",
    ident.2 = NULL,
    assay = NULL,
    layer = "data",
    methods = "wilcox",
    coexpr.method = "spearman",
    min_pct = 0.01,
    min_fc = 0.25,
    max_padj = 0.05,
    coexpr_min = 0.15,
    min_detect = 3,
    verbose = TRUE,
    ...
) {

  stopifnot(inherits(object, "Seurat"))

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  log <- .msg(verbose)

  log("GnRH marker discovery", type = "header")

  expr <- SeuratObject::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )

  if (!inherits(expr, "dgCMatrix")) {
    expr <- methods::as(expr, "dgCMatrix")
  }

  rownames(expr) <- make.unique(rownames(expr))

  gnrh_gene <- .find_gnrh_gene(rownames(expr))

  if (is.na(gnrh_gene)) {
    stop("GNRH1 not found.", call. = FALSE)
  }

  grp <- .get_cells(
    object = object,
    group.by = group.by,
    ident.1 = ident.1,
    ident.2 = ident.2
  )

  if (length(grp$cells_1) < 5) {
    stop("Too few positive cells.", call. = FALSE)
  }

  if (length(grp$cells_2) < 5) {
    stop("Too few comparison cells.", call. = FALSE)
  }

  SeuratObject::Idents(object) <- grp$groups

  all_res <- lapply(methods, function(m) {

    log("Running DE method:", m, type = "step")

    mk <- Seurat::FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      assay = assay,
      test.use = m,
      ...
    )

    if (!nrow(mk)) {
      return(NULL)
    }

    mk$gene <- rownames(mk)
    mk$method <- m
    mk
  })

  all_res <- all_res[!vapply(all_res, is.null, logical(1))]

  if (!length(all_res)) {
    stop("No markers detected by Seurat::FindMarkers().", call. = FALSE)
  }

  de <- .aggregate_de(do.call(rbind, all_res))

  keep <- Matrix::rowSums(expr > 0) >= min_detect

  expr <- expr[keep, , drop = FALSE]

  de <- de[
    de$gene %in% rownames(expr),
    ,
    drop = FALSE
  ]

  if (!nrow(de)) {
    stop("No markers passed min_detect filtering.", call. = FALSE)
  }

  spec <- .compute_specificity(
    expr = expr,
    genes = de$gene,
    cells_1 = grp$cells_1,
    cells_2 = grp$cells_2
  )

  idx <- match(de$gene, spec$gene)

  de$pct1 <- spec$pct1[idx]
  de$pct2 <- spec$pct2[idx]
  de$spec <- spec$spec[idx]

  coexpr <- .compute_coexpr(
    expr = expr,
    gnrh_gene = gnrh_gene,
    cells = grp$cells_1,
    method = coexpr.method
  )

  de$coexpr <- coexpr[de$gene]
  de$coexpr[is.na(de$coexpr)] <- 0

  de$coexpr_flag <- de$coexpr >= coexpr_min

  de <- .compute_score(de)

  de <- de[
    de$pct1 >= min_pct &
      de$avg_log2FC >= min_fc &
      de$p_val_adj <= max_padj &
      de$coexpr_flag,
    ,
    drop = FALSE
  ]

  de <- de[order(-de$score), , drop = FALSE]
  de <- de[!duplicated(de$gene), , drop = FALSE]

  rownames(de) <- de$gene

  object@misc$gnrh_markers <- de

  log("Markers detected:", nrow(de), type = "done")

  if (nrow(de) > 0) {
    log("Top marker:", de$gene[1], type = "info")
  }

  de
}
