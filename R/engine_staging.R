# =========================================================
# staging_modules.R
# =========================================================

# ---------------------------------------------------------
# Build developmental stage gene sets
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
build_stage_gene_sets <- function(genes) {

  map_genes <- function(features) {
    na.omit(
      match_gene_symbols(
        target = features,
        pool = genes
      )
    )
  }

  list(

    progenitor = map_genes(c(
      "FEZF1", "OTX2", "SIX3",
      "SIX6", "HES1", "SOX2",
      "NES", "VIM"
    )),

    migrating = map_genes(c(
      "ANOS1", "PROKR2", "DCX",
      "L1CAM", "SEMA3A", "ROBO1",
      "NRP1"
    )),

    transitional = map_genes(c(
      "GNRH1", "ISL1", "DLX5",
      "DLX6", "SCG2", "SYP"
    )),

    mature = map_genes(c(
      "GNRH1", "PCSK1", "PCSK2",
      "CHGA", "CHGB", "CPE",
      "VGF", "RAB3A"
    )),

    secreting = map_genes(c(
      "GNRH1", "KISS1R", "TAC3",
      "TACR3", "PCSK1", "SCG2",
      "VGF", "SYP"
    ))
  )
}


# =========================================================
# staging_scores.R
# =========================================================

# ---------------------------------------------------------
# Compute stage module scores
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
compute_stage_module_scores <- function(
    expr,
    modules
) {

  lib <- Matrix::colSums(expr)

  score_module <- function(gset) {

    if (length(gset) == 0) {
      return(rep(0, ncol(expr)))
    }

    mat <- expr[gset, , drop = FALSE]

    total_expr <- Matrix::colSums(mat)

    detected_n <- Matrix::colSums(mat > 0)
    detected_n[detected_n == 0] <- 1

    avg_expr <- total_expr / detected_n

    library_normalize(
      x = avg_expr,
      lib = lib
    )
  }

  scores <- lapply(modules, score_module)

  scores <- as.data.frame(
    scores,
    check.names = FALSE
  )

  rownames(scores) <- colnames(expr)

  scores
}


# ---------------------------------------------------------
# Assign dominant developmental stage
# ---------------------------------------------------------

#' @keywords internal
#' @noRd
assign_cell_stage <- function(scores) {

  apply(
    scores,
    1,
    function(x) {
      names(which.max(x))
    }
  )
}


# =========================================================
# staging_pipeline.R
# =========================================================

# ---------------------------------------------------------
# Stage GnRH cells
# ---------------------------------------------------------

#' Stage GnRH-expressing cells
#'
#' Assigns developmental stages using predefined
#' GnRH lineage gene programs.
#'
#' @param object A Seurat object.
#' @param assay Assay name.
#' @param layer Expression layer.
#' @param verbose Print progress messages.
#'
#' @return Updated Seurat object.
#'
#' @export
stage_gnrh_cells <- function(
    object,
    assay = "RNA",
    layer = "counts",
    verbose = TRUE
) {

  stopifnot(inherits(object, "Seurat"))

  log <- log_msg(verbose)

  log("==== GNRH STAGING START ====")

  # -------------------------------------------------------
  # Expression matrix
  # -------------------------------------------------------

  expr <- extract_expression_matrix(
    object = object,
    assay = assay,
    layer = layer
  )

  # -------------------------------------------------------
  # Stage gene programs
  # -------------------------------------------------------

  modules <- build_stage_gene_sets(
    genes = rownames(expr)
  )

  # -------------------------------------------------------
  # Compute module scores
  # -------------------------------------------------------

  scores <- compute_stage_module_scores(
    expr = expr,
    modules = modules
  )

  # -------------------------------------------------------
  # Assign stages
  # -------------------------------------------------------

  stage <- assign_cell_stage(scores)

  # -------------------------------------------------------
  # Non-GnRH cells
  # -------------------------------------------------------

  if ("gnrh_class" %in% colnames(object[[]])) {

    stage[
      object$gnrh_class == "gnrh_neg"
    ] <- "other"
  }

  # -------------------------------------------------------
  # Store metadata
  # -------------------------------------------------------

  object$gnrh_stage <- factor(
    stage,
    levels = c(
      "progenitor",
      "migrating",
      "transitional",
      "mature",
      "secreting",
      "other"
    )
  )

  # -------------------------------------------------------
  # Store scores
  # -------------------------------------------------------

  object$gnrh_prog_score <- scores$progenitor
  object$gnrh_mig_score  <- scores$migrating
  object$gnrh_trans_score <- scores$transitional
  object$gnrh_mat_score  <- scores$mature
  object$gnrh_sec_score  <- scores$secreting

  # -------------------------------------------------------
  # Store modules
  # -------------------------------------------------------

  object@misc$gnrh_stage_modules <- modules

  log(
    "Stage distribution:",
    paste(
      names(table(object$gnrh_stage)),
      table(object$gnrh_stage),
      collapse = " | "
    )
  )

  log("==== GNRH STAGING DONE ====")

  object
}
