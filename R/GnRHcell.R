#' Run complete GnRHcell analysis pipeline
#'
#' Executes the full GnRHcell workflow for identification,
#' developmental staging, diagnostics, and runtime reporting.
#'
#' The pipeline can perform:
#' \itemize{
#'   \item GnRH neuron detection via \code{\link{detect_gnrh}}
#'   \item developmental stage assignment via \code{\link{stage_gnrh}}
#'   \item diagnostic analysis via \code{\link{gnrh_diagnostics}}
#' }
#'
#' Runtime metrics, dataset summaries, analysis parameters,
#' and session information are stored in
#' \code{object@misc$gnrh$run_info}.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param detect Logical; run GnRH detection.
#' Default is \code{TRUE}.
#' @param stage Logical; run developmental staging.
#' Default is \code{TRUE}.
#' @param diagnostics Logical; run diagnostic analysis.
#' Default is \code{TRUE}.
#' @param verbose Logical; print progress messages.
#' Default is \code{TRUE}.
#' @param ... Additional arguments passed to
#' \code{\link{detect_gnrh}}.
#'
#' @return A Seurat object updated with GnRH detection,
#' developmental staging, diagnostics, and runtime metadata.
#'
#' @details
#' Workflow steps:
#' \enumerate{
#'   \item GnRH detection using transcriptomic and marker-based scoring
#'   \item developmental state assignment
#'   \item diagnostic performance evaluation
#'   \item runtime and summary reporting
#' }
#'
#' Metadata added may include:
#' \describe{
#'   \item{\code{gnrh_status}}{Binary GnRH classification.}
#'   \item{\code{gnrh_class}}{Internal classification labels.}
#'   \item{\code{gnrh_truth}}{High-confidence truth labels.}
#'   \item{\code{gnrh_stage}}{Developmental stage assignments.}
#' }
#'
#' Stored runtime metadata:
#' \describe{
#'   \item{\code{detect_sec}}{Detection runtime in seconds.}
#'   \item{\code{stage_sec}}{Staging runtime in seconds.}
#'   \item{\code{diagnostics_sec}}{Diagnostics runtime in seconds.}
#'   \item{\code{total_sec}}{Total pipeline runtime.}
#' }
#'
#' If specific steps are disabled, only selected workflow components
#' are executed.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{stage_gnrh}},
#' \code{\link{gnrh_diagnostics}},
#' \code{\link{extract_gnrh_run_info}}
#'
#' @examples
#' \dontrun{
#' obj <- run_gnrh(seurat_obj)
#'
#' obj <- run_gnrh(
#'   seurat_obj,
#'   detect = TRUE,
#'   stage = TRUE,
#'   diagnostics = TRUE
#' )
#' }
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

  log <- .msg(verbose)

  object <- .init_gnrh_misc(object)

  # IMPORTANT
  step_times <- list()
  total_start <- Sys.time()

  log <- .msg()

  log("==== STARTING GnRHcell PIPELINE ====", type = "header")

  # detect
  if (detect) {
    log("[1/3] Detecting GnRH cells", type = "step")


    t0 <- Sys.time()

    object <- detect_gnrh(
      object,
      verbose = verbose,
      ...
    )

    step_times$detect_sec <- as.numeric(
      difftime(Sys.time(), t0, units = "secs")
    )
    log("Detection complete.", type = "done",
        duration = step_times$detect_sec)
  }

  # stage
  if (stage) {
    log("[2/3] Assigning developmental stages", type = "step")

    t0 <- Sys.time()

    object <- stage_gnrh(
      object,
      verbose = verbose
    )

    step_times$stage_sec <- as.numeric(
      difftime(Sys.time(), t0, units = "secs")
    )
  }

  # diagnostics
  if (diagnostics) {
    log("[3/3] Running diagnostics", type = "info")


    t0 <- Sys.time()

    object <- gnrh_diagnostics(
      object,
      verbose = verbose
    )

    step_times$diagnostics_sec <- as.numeric(
      difftime(Sys.time(), t0, units = "secs")
    )
    log("Assigning stages complete.", type = "done",
        duration = step_times$detect_sec)
  }


  # total
  step_times$total_sec <- as.numeric(
    difftime(Sys.time(), total_start, units = "secs")
  )

  object <- .add_run_info(
    object,
    step_times = step_times,
    params = list(...)
  )

  # summary
  md <- object[[]]


  log("PIPELINE SUMMARY", type = "info")


  if ("gnrh_status" %in% colnames(md)) {

    log("Status:")

    status_tab <- table(md$gnrh_status)

    for (nm in names(status_tab)) {
      log(sprintf("  %s: %d", nm, status_tab[[nm]]))
    }
  }

  if ("gnrh_truth" %in% colnames(md)) {

    log("Truth:")

    truth_tab <- table(md$gnrh_truth)

    for (nm in names(truth_tab)) {
      log(sprintf("  %s: %d", nm, truth_tab[[nm]]))
    }
  }

  if ("gnrh_stage" %in% colnames(md)) {

    log("Stage:")

    stage_tab <- table(md$gnrh_stage)

    for (nm in names(stage_tab)) {
      log(sprintf("  %s: %d", nm, stage_tab[[nm]]))
    }
  }

  log("==== GnRHcell PIPELINE COMPLETE ====", type = "done")

  object
}
