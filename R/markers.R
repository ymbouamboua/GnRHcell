# ========================================================= #
# GnRHcell marker discovery
# ========================================================= #

# --------------------------------------------------------- #
# Helpers
# --------------------------------------------------------- #

.safe_cor <- function(x, y, method = "spearman") {

  out <- suppressWarnings(
    stats::cor(x, y, method = method)
  )

  ifelse(is.finite(out), out, 0)
}


.get_gnrh_gene <- function(genes) {

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


# --------------------------------------------------------- #
# Core metrics
# --------------------------------------------------------- #

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


# --------------------------------------------------------- #
# DE aggregation
# --------------------------------------------------------- #

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


# --------------------------------------------------------- #
# Main function
# --------------------------------------------------------- #

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
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers()}.
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

  gnrh_gene <- .get_gnrh_gene(rownames(expr))

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

  # Canonical Seurat-compatible column names
  de$`pct.1` <- spec$pct1[idx]
  de$`pct.2` <- spec$pct2[idx]
  de$specificity <- spec$spec[idx]

  # Backward-compatible aliases
  de$pct1 <- de$`pct.1`
  de$pct2 <- de$`pct.2`
  de$spec <- de$specificity

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


#' Identify specific and reproducible GnRH-expressed genes
#'
#' Identifies genes enriched in GnRH-positive cells relative to matched
#' control cells, measures their same-cell co-detection with \code{GNRH1},
#' and optionally evaluates their reproducibility across donors.
#'
#' This function provides a high-level interface to
#' \code{\link{gnrh_markers}}. The latter performs differential-expression
#' testing and expression-correlation analysis, while
#' \code{find_gnrh_genes} additionally selects matched controls, calculates
#' the percentage of \code{GNRH1}-positive cells expressing each marker,
#' excludes generic genes, and ranks candidates by specificity and donor
#' support.
#'
#' @param object A Seurat object containing GnRH classification metadata and
#'   normalized expression data.
#' @param status_col Character. Metadata column containing GnRH-positive and
#'   GnRH-negative classifications. Default is \code{"gnrh_status"}.
#' @param positive Character. Value in \code{status_col} identifying
#'   GnRH-positive cells. Default is \code{"pos"}.
#' @param negative Character. Value in \code{status_col} identifying
#'   GnRH-negative cells. Default is \code{"neg"}.
#' @param confidence_col Optional character. Metadata column used to restrict
#'   the positive population to confident GnRH calls. If \code{NULL}, no
#'   confidence filter is applied.
#' @param confidence_value Value or values in \code{confidence_col} defining
#'   confident cells. Default is \code{"pos"}.
#' @param annotation_col Character. Metadata column containing broad cell-type
#'   annotations used to select matched controls. Default is \code{"ann1"}.
#' @param control_ident Character vector. Values in \code{annotation_col}
#'   defining the control population. Default is \code{"Neuronal"}.
#' @param donor_col Optional character. Metadata column identifying biological
#'   donors or samples. If supplied, candidate detection is evaluated
#'   independently in each donor.
#' @param assay Character. Seurat assay used for marker discovery and
#'   co-expression analysis. Default is \code{"RNA"}.
#' @param layer Character. Assay layer containing normalized expression values.
#'   Default is \code{"data"}.
#' @param coexpr_gene Character. Reference gene used for same-cell
#'   co-expression analysis. Default is \code{"GNRH1"}.
#' @param methods Character vector of differential-expression methods passed
#'   to \code{\link{gnrh_markers}}. Default is \code{"wilcox"}.
#' @param min_fc Numeric. Minimum average log2 fold change required for final
#'   candidates. Default is \code{0.5}.
#' @param min_pct Numeric. Minimum detection fraction among GnRH-positive cells.
#'   Default is \code{0.20}.
#' @param max_control_pct Numeric. Maximum detection fraction allowed among
#'   control cells. Default is \code{0.10}.
#' @param min_coexpr_pct Numeric. Minimum percentage of
#'   \code{coexpr_gene}-positive GnRH cells in which a candidate must be
#'   detected. Expressed on a 0--100 scale. Default is \code{20}.
#' @param max_padj Numeric. Maximum adjusted p-value allowed for candidate
#'   markers. Default is \code{0.05}.
#' @param min_donor_pct Numeric. Minimum within-donor detection fraction used
#'   to count a candidate as supported by a donor. Default is \code{0.10}.
#' @param min_cells_donor Integer. Minimum number of selected GnRH cells
#'   required for a donor to be evaluated. Default is \code{3}.
#' @param generic_genes Character vector of generic neuronal, housekeeping, or
#'   broadly expressed genes excluded from the final candidate table.
#' @param seed Integer. Random seed used for differential-expression analysis.
#'   Default is \code{1234}.
#' @param verbose Logical. Display progress messages. Default is \code{TRUE}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{candidates}}{
#'     Filtered and ranked candidate GnRH genes. This table includes
#'     differential-expression statistics, GnRH-cell and control-cell
#'     detection fractions, correlation-based co-expression,
#'     same-cell co-detection counts and percentages, specificity score,
#'     and donor support when available.
#'   }
#'   \item{\code{markers}}{
#'     The complete marker table returned by \code{\link{gnrh_markers}},
#'     augmented with same-cell co-detection metrics.
#'   }
#'   \item{\code{donor_detection}}{
#'     A gene-by-donor matrix containing within-donor detection fractions.
#'     Returns \code{NULL} when \code{donor_col} is not supplied or no donor
#'     is evaluable.
#'   }
#'   \item{\code{donor_cells}}{
#'     Number of selected GnRH cells per donor, or \code{NULL}.
#'   }
#'   \item{\code{gnrh_cells}}{
#'     Cell barcodes assigned to the GnRH-positive group.
#'   }
#'   \item{\code{control_cells}}{
#'     Cell barcodes assigned to the matched control group.
#'   }
#'   \item{\code{parameters}}{
#'     Main selection thresholds used for the analysis.
#'   }
#' }
#'
#' @details
#' Candidate genes are required to be positively enriched in GnRH cells,
#' detected in at least \code{min_pct} of GnRH cells, detected in no more than
#' \code{max_control_pct} of matched controls, co-detected with
#' \code{coexpr_gene} in at least \code{min_coexpr_pct} percent of reference
#' gene-positive cells, and significant at \code{max_padj}.
#'
#' The candidate specificity score is:
#' \deqn{
#' logFC \times (pct.1 - pct.2) \times log(1 + coexpression\ percentage)
#' }
#'
#' Correlation-based co-expression and same-cell co-detection measure distinct
#' properties. Correlation evaluates coordinated expression variation, whereas
#' co-detection records whether both transcripts are detected in the same
#' cells. Neither metric alone demonstrates a direct regulatory or molecular
#' interaction.
#'
#' Donor support should be interpreted cautiously when few GnRH cells are
#' available. Donors represented by fewer than \code{min_cells_donor} selected
#' cells are excluded from donor-level evaluation.
#'
#' For Seurat v5 objects containing multiple layers for the selected assay,
#' layers should be joined with \code{Seurat::JoinLayers()} before calling this
#' function.
#'
#' @examples
#' \dontrun{
#' obj <- Seurat::JoinLayers(obj, assay = "RNA")
#'
#' result <- find_gnrh_genes(
#'   object = obj,
#'   status_col = "gnrh_status",
#'   positive = "pos",
#'   negative = "neg",
#'   confidence_col = "gnrh_confident",
#'   confidence_value = "pos",
#'   annotation_col = "ann1",
#'   control_ident = "Neuronal",
#'   donor_col = "status"
#' )
#'
#' head(result$candidates, 20)
#' result$donor_cells
#' result$donor_detection
#'
#' # Run without confidence or donor filtering
#' result <- find_gnrh_genes(
#'   object = obj,
#'   annotation_col = "ann1",
#'   control_ident = "Neuronal"
#' )
#' }
#'
#' @seealso
#' \code{\link{gnrh_markers}},
#' \code{\link{gnrh_marker_programs}},
#' \code{\link{gnrh_stage_modules}}
#'
#' @family marker discovery
#' @export
find_gnrh_genes <- function(
    object,
    status_col = "gnrh_status",
    positive = "pos",
    negative = "neg",
    confidence_col = NULL,
    confidence_value = "pos",
    annotation_col = "ann1",
    control_ident = "Neuronal",
    donor_col = NULL,
    assay = "RNA",
    layer = "data",
    coexpr_gene = "GNRH1",
    methods = "wilcox",
    min_fc = 0.5,
    min_pct = 0.20,
    max_control_pct = 0.10,
    min_coexpr_pct = 20,
    max_padj = 0.05,
    min_donor_pct = 0.10,
    min_cells_donor = 3,
    generic_genes = c(
      "TUBB3", "RBFOX3", "GAPDH", "ACTB", "MALAT1",
      "GAD1", "GAD2", "SLC17A7", "SNAP25", "SYT1"
    ),
    seed = 1234,
    verbose = TRUE
) {
  stopifnot(inherits(object, "Seurat"))
  set.seed(seed)

  md <- object[[]]

  required <- c(status_col, annotation_col)

  if (!is.null(confidence_col))
    required <- c(required, confidence_col)

  if (!is.null(donor_col))
    required <- c(required, donor_col)

  missing <- setdiff(required, colnames(md))

  if (length(missing)) {
    stop(
      "Missing metadata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!coexpr_gene %in% rownames(object)) {
    stop("Gene not found: ", coexpr_gene, call. = FALSE)
  }

  # Select GnRH and matched neuronal controls

  positive_index <- !is.na(md[[status_col]]) &
    md[[status_col]] == positive

  if (!is.null(confidence_col)) {
    positive_index <- positive_index &
      !is.na(md[[confidence_col]]) &
      md[[confidence_col]] %in% confidence_value
  }

  control_index <- !is.na(md[[status_col]]) &
    md[[status_col]] == negative &
    !is.na(md[[annotation_col]]) &
    md[[annotation_col]] %in% control_ident

  gnrh_cells <- rownames(md)[positive_index]
  control_cells <- rownames(md)[control_index]

  if (length(gnrh_cells) < 5L) {
    stop("Fewer than five GnRH cells selected.", call. = FALSE)
  }

  if (length(control_cells) < 5L) {
    stop("Fewer than five control cells selected.", call. = FALSE)
  }

  object$gnrh_comparison <- "Other"
  object$gnrh_comparison[gnrh_cells] <- "GnRH"
  object$gnrh_comparison[control_cells] <- "Control"

  Seurat::DefaultAssay(object) <- assay

  # Reuse existing GnRHcell marker engine

  markers <- gnrh_markers(
    object = object,
    group.by = "gnrh_comparison",
    ident.1 = "GnRH",
    ident.2 = "Control",
    assay = assay,
    layer = layer,
    methods = methods,
    min_pct = 0.05,
    min_fc = 0.25,
    max_padj = max_padj,
    coexpr_min = -Inf,
    verbose = verbose,
    only.pos = TRUE,
    logfc.threshold = 0.25
  )

  if (!nrow(markers)) {
    stop("No positive GnRH markers detected.", call. = FALSE)
  }

  # Ensure canonical names for older gnrh_markers versions
  if (!"pct.1" %in% colnames(markers))
    markers$`pct.1` <- markers$pct1

  if (!"pct.2" %in% colnames(markers))
    markers$`pct.2` <- markers$pct2

  # Same-cell co-detection with GNRH1

  expression <- SeuratObject::LayerData(
    object,
    assay = assay,
    layer = layer
  )[unique(c(coexpr_gene, markers$gene)), gnrh_cells, drop = FALSE]

  gnrh1_positive <- expression[coexpr_gene, ] > 0

  if (!any(gnrh1_positive)) {
    stop(
      "No selected cells express ",
      coexpr_gene,
      ".",
      call. = FALSE
    )
  }

  detection <- expression[
    markers$gene,
    gnrh1_positive,
    drop = FALSE
  ] > 0

  # Preserve correlation produced by gnrh_markers()
  markers$coexpr_cor <- markers$coexpr
  markers$coexpr_n <- Matrix::rowSums(detection)
  markers$coexpr_pct <- 100 * Matrix::rowMeans(detection)
  markers$coexpr_flag <- markers$coexpr_pct >= min_coexpr_pct

  # Select specific candidates

  candidates <- markers[
    markers$avg_log2FC >= min_fc &
      markers$`pct.1` >= min_pct &
      markers$`pct.2` <= max_control_pct &
      markers$coexpr_flag &
      markers$p_val_adj <= max_padj,
    ,
    drop = FALSE
  ]

  candidates <- candidates[
    !toupper(candidates$gene) %in% toupper(generic_genes),
    ,
    drop = FALSE
  ]

  candidates$specificity_score <- with(
    candidates,
    avg_log2FC *
      (`pct.1` - `pct.2`) *
      log1p(coexpr_pct)
  )

  donor_detection <- NULL
  donor_cells <- NULL

  # Donor reproducibility

  if (!is.null(donor_col) && nrow(candidates)) {
    donor <- factor(md[gnrh_cells, donor_col, drop = TRUE])
    donor_cells <- table(donor)

    evaluable <- names(donor_cells)[
      donor_cells >= min_cells_donor
    ]

    if (length(evaluable)) {
      candidate_expression <- expression[
        candidates$gene,
        ,
        drop = FALSE
      ]

      donor_detection <- vapply(
        evaluable,
        \(d) {
          Matrix::rowMeans(
            candidate_expression[
              , donor == d,
              drop = FALSE
            ] > 0
          )
        },
        numeric(nrow(candidate_expression))
      )

      if (is.null(dim(donor_detection))) {
        donor_detection <- matrix(
          donor_detection,
          ncol = 1L,
          dimnames = list(
            rownames(candidate_expression),
            evaluable
          )
        )
      } else {
        rownames(donor_detection) <-
          rownames(candidate_expression)
      }

      candidates$n_donors <- rowSums(
        donor_detection[
          candidates$gene,
          ,
          drop = FALSE
        ] >= min_donor_pct
      )

      candidates$n_donors_evaluable <- length(evaluable)
    }
  }

  if (!"n_donors" %in% colnames(candidates))
    candidates$n_donors <- NA_integer_

  if (!"n_donors_evaluable" %in% colnames(candidates))
    candidates$n_donors_evaluable <- NA_integer_

  candidates <- candidates[
    order(
      candidates$n_donors,
      candidates$specificity_score,
      decreasing = TRUE,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  list(
    candidates = candidates,
    markers = markers,
    donor_detection = donor_detection,
    donor_cells = donor_cells,
    gnrh_cells = gnrh_cells,
    control_cells = control_cells,
    parameters = list(
      positive = positive,
      negative = negative,
      control_ident = control_ident,
      min_fc = min_fc,
      min_pct = min_pct,
      max_control_pct = max_control_pct,
      min_coexpr_pct = min_coexpr_pct,
      max_padj = max_padj
    )
  )
}
