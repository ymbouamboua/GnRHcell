# =========================================================
# pipeline_main.R
# =========================================================

# ---------------------------------------------------------
# Run complete GnRH analysis pipeline
# ---------------------------------------------------------

#' Run complete GnRH analysis workflow
#'
#' Executes:
#' \itemize{
#'   \item GnRH cell detection
#'   \item developmental staging
#'   \item diagnostic analysis
#' }
#'
#' @param object A Seurat object.
#' @param detect Run GnRH cell detection.
#' @param stage Run developmental stage assignment.
#' @param diagnostics Run diagnostic metrics and plots.
#' @param verbose Print progress messages.
#' @param ... Additional arguments passed to
#'   \code{detect_gnrh_cells()}.
#'
#' @return Updated Seurat object.
#'
#' @export
run_gnrh <- function(
    object,
    detect = TRUE,
    stage = TRUE,
    diagnostics = TRUE,
    verbose = TRUE,
    ...
) {

  stopifnot(inherits(object, "Seurat"))

  log <- log_msg(verbose)

  log("==== RUNNING GNRH PIPELINE ====")

  # -------------------------------------------------------
  # GnRH detection
  # -------------------------------------------------------

  if (isTRUE(detect)) {

    log("[1] Detecting GnRH cells...")

    object <- detect_gnrh_cells(
      object = object,
      verbose = verbose,
      ...
    )
  }

  # -------------------------------------------------------
  # Developmental staging
  # -------------------------------------------------------

  if (isTRUE(stage)) {

    log("[2] Assigning developmental stages...")

    object <- stage_gnrh_cells(
      object = object,
      verbose = verbose
    )
  }

  # -------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------

  if (isTRUE(diagnostics)) {

    log("[3] Running diagnostics...")

    object <- run_detection_diagnostics(
      object = object,
      verbose = verbose
    )
  }

  # -------------------------------------------------------
  # Summary
  # -------------------------------------------------------

  if ("gnrh_class" %in% colnames(object[[]])) {

    class_tab <- table(object$gnrh_class)

    log(
      "Detected classes:",
      paste(
        names(class_tab),
        class_tab,
        collapse = " | "
      )
    )
  }

  if ("gnrh_stage" %in% colnames(object[[]])) {

    stage_tab <- table(object$gnrh_stage)

    log(
      "Detected stages:",
      paste(
        names(stage_tab),
        stage_tab,
        collapse = " | "
      )
    )
  }

  log("==== GNRH PIPELINE DONE ====")

  object
}
