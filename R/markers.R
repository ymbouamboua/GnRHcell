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




#' Identify stage-specific markers in GnRH-lineage cells
#'
#' Finds genes enriched in each GnRH developmental stage relative to the other
#' GnRH-positive stages. Unlike `gnrh_markers()`, this function does not require
#' correlation or co-detection with `GNRH1`, because those filters favor general
#' GnRH-lineage genes rather than stage-specific programs.
#'
#' @param object A Seurat object processed with `run_gnrh()`.
#' @param stage_col Metadata column containing developmental stages.
#' @param stages Stages to test. By default, `identity`, `migrating`, and
#'   `mature` are tested when present.
#' @param class_col Metadata column defining GnRH detection classes. Set to
#'   `NULL` to select cells using `stage_col` alone.
#' @param positive_classes Detection classes considered GnRH-lineage positive.
#' @param exclude_stages Stage values excluded from all comparisons.
#' @param assay Assay used for differential expression.
#' @param layer Normalized-expression layer used for donor concordance.
#' @param test_use Differential-expression test passed to
#'   `Seurat::FindMarkers()`.
#' @param min_cells Minimum number of cells required in both the target stage
#'   and its reference group.
#' @param min_pct Minimum expression fraction passed to `FindMarkers()`.
#' @param min_log2fc Minimum positive average log2 fold change.
#' @param max_padj Maximum adjusted p-value.
#' @param min_specificity Minimum `pct.1 - pct.2` detection difference.
#' @param donor_col Optional biological replicate column. When supplied, the
#'   direction of the stage effect is evaluated independently within donors.
#' @param min_cells_donor Minimum cells required in both comparison groups for
#'   a donor to be evaluable for a stage.
#' @param min_donor_support Minimum fraction of evaluable donors in which the
#'   mean normalized expression effect is positive.
#' @param min_donors Minimum number of evaluable donors required before donor
#'   support is used as a candidate filter.
#' @param max_cells_per_ident Optional maximum cells sampled from each identity
#'   by `FindMarkers()`. Useful for very large atlases.
#' @param exclude_pattern Optional regular expression for genes to exclude from
#'   the candidate table, for example mitochondrial or ribosomal genes.
#' @param seed Random seed used by differential-expression subsampling.
#' @param verbose Display progress messages.
#' @param ... Additional arguments passed to `Seurat::FindMarkers()`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{candidates}{Filtered and ranked stage-specific markers.}
#'   \item{markers}{Complete positive marker results for all tested stages.}
#'   \item{donor_effects}{Long table of within-donor mean-expression effects.}
#'   \item{stage_counts}{Numbers of selected cells per stage.}
#'   \item{cells}{Cell barcodes used in the analysis.}
#'   \item{parameters}{Principal analysis parameters.}
#' }
#'
#' @details
#' Each stage is tested against the union of the other selected GnRH stages
#' (one-versus-rest). The reported `specificity` is `pct.1 - pct.2`. The ranking
#' score combines positive fold change, detection specificity, statistical
#' significance, and—when available—donor concordance.
#'
#' Cell-level differential expression is useful for marker discovery but does
#' not replace a replicate-aware pseudobulk analysis for formal inference.
#' `donor_support` should therefore be used to prioritize reproducible markers,
#' while final publication claims should be confirmed with pseudobulk counts.
#'
#' @examples
#' \dontrun{
#' stage_markers <- find_gnrh_stage_markers(
#'   object = wang,
#'   stage_col = "gnrh_stage",
#'   donor_col = "orig.ident",
#'   stages = c("identity", "migrating", "mature")
#' )
#'
#' head(stage_markers$candidates)
#' subset(stage_markers$candidates, stage == "migrating")
#'
#' # Faster discovery in a very large atlas
#' stage_markers <- find_gnrh_stage_markers(
#'   object = human_hypomap,
#'   donor_col = "orig.ident",
#'   max_cells_per_ident = 5000
#' )
#' }
#'
#' @seealso `gnrh_markers()`, `find_gnrh_genes()`
#' @family marker discovery
#' @export
find_gnrh_stage_markers <- function(
    object,
    stage_col = "gnrh_stage",
    stages = c("identity", "migrating", "mature"),
    class_col = "gnrh_class",
    positive_classes = c("direct", "supported"),
    exclude_stages = c("non-gnrh", "secreting"),
    assay = NULL,
    layer = "data",
    test_use = "wilcox",
    min_cells = 10L,
    min_pct = 0.10,
    min_log2fc = 0.25,
    max_padj = 0.05,
    min_specificity = 0.05,
    donor_col = NULL,
    min_cells_donor = 3L,
    min_donor_support = 0.60,
    min_donors = 2L,
    max_cells_per_ident = Inf,
    exclude_pattern = NULL,
    seed = 1234L,
    verbose = TRUE,
    ...) {

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  metadata <- object[[]]
  required <- c(stage_col, class_col, donor_col)
  required <- required[!vapply(required, is.null, logical(1))]
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop("Missing metadata column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Assay `", assay, "` was not found.", call. = FALSE)
  }

  stage_values <- as.character(metadata[[stage_col]])
  keep <- !is.na(stage_values) & nzchar(stage_values) &
    !stage_values %in% exclude_stages
  if (!is.null(class_col)) {
    keep <- keep & !is.na(metadata[[class_col]]) &
      metadata[[class_col]] %in% positive_classes
  }

  cells <- rownames(metadata)[keep]
  observed <- unique(stage_values[keep])
  stages <- intersect(stages, observed)
  if (length(stages) < 2L) {
    stop(
      "At least two requested stages must be present among selected GnRH cells. Present: ",
      paste(observed, collapse = ", "), call. = FALSE
    )
  }

  # Restrict the reference pool to the explicitly selected stages so that each
  # comparison is target stage versus the other biologically ordered stages.
  cells <- cells[stage_values[match(cells, rownames(metadata))] %in% stages]
  stage_factor <- factor(as.character(metadata[cells, stage_col, drop = TRUE]), levels = stages)
  stage_counts <- table(stage_factor)
  too_small <- names(stage_counts)[stage_counts < min_cells]
  if (length(too_small)) {
    stop(
      "Stage(s) below `min_cells = ", min_cells, "`: ",
      paste(paste0(too_small, " (n=", stage_counts[too_small], ")"), collapse = ", "),
      call. = FALSE
    )
  }

  analysis_object <- subset(object, cells = cells)
  analysis_object[[".gnrh_stage_test"]] <- stage_factor
  SeuratObject::Idents(analysis_object) <- ".gnrh_stage_test"
  set.seed(seed)

  run_stage <- function(stage) {
    if (verbose) message("Testing stage `", stage, "` versus other GnRH stages...")

    reference <- setdiff(stages, stage)
    result <- Seurat::FindMarkers(
      object = analysis_object,
      ident.1 = stage,
      ident.2 = reference,
      assay = assay,
      test.use = test_use,
      only.pos = TRUE,
      min.pct = min_pct,
      logfc.threshold = min_log2fc,
      max.cells.per.ident = max_cells_per_ident,
      random.seed = seed,
      verbose = FALSE,
      ...
    )

    if (!nrow(result)) return(NULL)
    result$gene <- rownames(result)
    result$stage <- stage
    result$reference <- paste(reference, collapse = "+")

    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(result))
    if (!length(fc_col)) stop("FindMarkers result has no average log-fold-change column.", call. = FALSE)
    result$avg_log2FC <- result[[fc_col[[1L]]]]
    result
  }

  marker_list <- lapply(stages, run_stage)
  marker_list <- marker_list[!vapply(marker_list, is.null, logical(1))]
  if (!length(marker_list)) stop("No stage markers were detected.", call. = FALSE)
  markers <- dplyr::bind_rows(marker_list)
  markers$specificity <- markers$`pct.1` - markers$`pct.2`

  donor_effects <- data.frame()
  donor_summary <- data.frame()

  if (!is.null(donor_col)) {
    expression <- SeuratObject::LayerData(
      analysis_object,
      assay = assay,
      layer = layer
    )
    genes <- intersect(unique(markers$gene), rownames(expression))
    expression <- expression[genes, cells, drop = FALSE]
    donor <- as.character(metadata[cells, donor_col, drop = TRUE])
    stage_by_cell <- as.character(metadata[cells, stage_col, drop = TRUE])

    donor_rows <- lapply(stages, function(stage) {
      evaluable <- unique(donor[!is.na(donor)])
      evaluable <- evaluable[vapply(evaluable, function(id) {
        idx <- donor == id
        sum(idx & stage_by_cell == stage) >= min_cells_donor &&
          sum(idx & stage_by_cell %in% setdiff(stages, stage)) >= min_cells_donor
      }, logical(1))]

      if (!length(evaluable)) return(NULL)
      stage_genes <- intersect(markers$gene[markers$stage == stage], genes)

      dplyr::bind_rows(lapply(evaluable, function(id) {
        target <- donor == id & stage_by_cell == stage
        reference <- donor == id & stage_by_cell %in% setdiff(stages, stage)
        effect <- Matrix::rowMeans(expression[stage_genes, target, drop = FALSE]) -
          Matrix::rowMeans(expression[stage_genes, reference, drop = FALSE])
        data.frame(
          stage = stage, gene = stage_genes, donor = id,
          mean_difference = unname(effect),
          n_target = sum(target), n_reference = sum(reference)
        )
      }))
    })

    donor_effects <- dplyr::bind_rows(donor_rows)
    if (nrow(donor_effects)) {
      donor_summary <- donor_effects |>
        dplyr::group_by(.data$stage, .data$gene) |>
        dplyr::summarise(
          n_donors = dplyr::n(),
          donors_positive = sum(.data$mean_difference > 0),
          donor_support = mean(.data$mean_difference > 0),
          median_donor_effect = stats::median(.data$mean_difference),
          .groups = "drop"
        )
      markers <- dplyr::left_join(markers, donor_summary, by = c("stage", "gene"))
    }
  }

  for (column in c("n_donors", "donors_positive", "donor_support", "median_donor_effect")) {
    if (!column %in% names(markers)) markers[[column]] <- NA_real_
  }

  candidates <- markers |>
    dplyr::filter(
      .data$avg_log2FC >= min_log2fc,
      .data$p_val_adj <= max_padj,
      .data$specificity >= min_specificity
    )

  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    candidates <- candidates[!grepl(exclude_pattern, candidates$gene, ignore.case = TRUE), , drop = FALSE]
  }

  # Apply donor filtering only when enough donors were evaluable. Markers with
  # insufficient donor information are retained but clearly flagged.
  candidates$donor_evaluable <- !is.na(candidates$n_donors) & candidates$n_donors >= min_donors
  candidates$donor_reproducible <- candidates$donor_evaluable &
    candidates$donor_support >= min_donor_support
  if (!is.null(donor_col)) {
    candidates <- candidates[!candidates$donor_evaluable | candidates$donor_reproducible, , drop = FALSE]
  }

  significance <- pmin(-log10(pmax(candidates$p_val_adj, .Machine$double.xmin)), 50)
  donor_weight <- ifelse(
    is.na(candidates$donor_support), 1,
    0.5 + candidates$donor_support
  )
  candidates$stage_score <- candidates$avg_log2FC *
    pmax(candidates$specificity, 0) * significance * donor_weight
  candidates <- candidates |>
    dplyr::arrange(.data$stage, dplyr::desc(.data$stage_score), .data$p_val_adj)

  list(
    candidates = candidates,
    markers = markers,
    donor_effects = donor_effects,
    stage_counts = data.frame(stage = names(stage_counts), n_cells = as.integer(stage_counts)),
    cells = cells,
    parameters = list(
      stages = stages, positive_classes = positive_classes,
      assay = assay, layer = layer, test_use = test_use,
      min_cells = min_cells, min_pct = min_pct,
      min_log2fc = min_log2fc, max_padj = max_padj,
      min_specificity = min_specificity,
      donor_col = donor_col, min_donor_support = min_donor_support
    )
  )
}

