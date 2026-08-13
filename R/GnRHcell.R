#' Run complete GnRHcell analysis pipeline
#'
#' Executes the full GnRHcell workflow for identification,
#' developmental staging, diagnostics, and runtime reporting.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param detect Logical; run GnRH detection. Default is \code{TRUE}.
#' @param stage Logical; run developmental staging. Default is \code{TRUE}.
#' @param diagnostics Logical; run diagnostic analysis. Default is \code{TRUE}.
#' @param verbose Logical; print progress messages. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{detect_gnrh}}.
#'
#' @return A Seurat object updated with GnRH detection,
#' developmental staging, diagnostics, and runtime metadata.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{stage_gnrh}},
#' \code{\link{gnrh_diagnostics}},
#' \code{\link{extract_gnrh_run_info}}
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

  stopifnot(
    inherits(
      object,
      "Seurat"
    )
  )

  object <- .init_gnrh_misc(
    object
  )

  log <- .msg(
    verbose
  )

  step_times <- list()
  total_start <- Sys.time()

  log(
    "==== STARTING GnRHcell PIPELINE ====",
    type = "header"
  )

  # --------------------------------------------------------------------------- #
  # Step 1: detection
  # --------------------------------------------------------------------------- #

  if (isTRUE(detect)) {

    log(
      "[1/3] Detecting GnRH cells",
      type = "step"
    )

    t0 <- Sys.time()

    object <- detect_gnrh(
      object,
      verbose = verbose,
      ...
    )

    step_times$detect_sec <- as.numeric(
      difftime(
        Sys.time(),
        t0,
        units = "secs"
      )
    )

    log(
      "Detection complete.",
      type = "done",
      duration = step_times$detect_sec
    )
  }

  # --------------------------------------------------------------------------- #
  # Step 2: staging
  # --------------------------------------------------------------------------- #

  if (isTRUE(stage)) {

    log(
      "[2/3] Assigning developmental stages",
      type = "step"
    )

    t0 <- Sys.time()

    object <- stage_gnrh(
      object,
      verbose = FALSE
    )

    step_times$stage_sec <- as.numeric(
      difftime(
        Sys.time(),
        t0,
        units = "secs"
      )
    )

    log(
      "Staging complete.",
      type = "done",
      duration = step_times$stage_sec
    )
  }

  # --------------------------------------------------------------------------- #
  # Step 3: diagnostics
  # --------------------------------------------------------------------------- #

  if (isTRUE(diagnostics)) {

    log(
      "[3/3] Running diagnostics",
      type = "step"
    )

    t0 <- Sys.time()

    object <- gnrh_diagnostics(
      object,
      verbose = FALSE
    )

    step_times$diagnostics_sec <- as.numeric(
      difftime(
        Sys.time(),
        t0,
        units = "secs"
      )
    )

    log(
      "Diagnostics complete.",
      type = "done",
      duration = step_times$diagnostics_sec
    )
  }

  # --------------------------------------------------------------------------- #
  # Total runtime
  # --------------------------------------------------------------------------- #

  step_times$total_sec <- as.numeric(
    difftime(
      Sys.time(),
      total_start,
      units = "secs"
    )
  )

  object <- .add_run_info(
    object,
    step_times = step_times,
    params = list(...)
  )

  # --------------------------------------------------------------------------- #
  # Pipeline summary
  # --------------------------------------------------------------------------- #

  md <- object[[]]

  log(
    "PIPELINE SUMMARY",
    type = "info"
  )

  if (
    "gnrh_status" %in%
    colnames(md)
  ) {

    log(
      "Status:"
    )

    status_tab <- table(
      md$gnrh_status,
      useNA = "no"
    )

    for (nm in names(status_tab)) {

      log(
        sprintf(
          "  %s: %d",
          nm,
          status_tab[[nm]]
        )
      )
    }
  }

  if (all(c("gnrh_direct_signal", "gnrh_direct_supported",
            "gnrh_direct_isolated") %in% colnames(md))) {
    log("Detection evidence:")
    log(sprintf(
      "  GNRH1 direct signal: %d",
      sum(md$gnrh_direct_signal, na.rm = TRUE)
    ))
    log(sprintf(
      "  identity-supported direct: %d",
      sum(md$gnrh_direct_supported, na.rm = TRUE)
    ))
    log(sprintf(
      "  isolated GNRH1 signal: %d",
      sum(md$gnrh_direct_isolated, na.rm = TRUE)
    ))
  }

  if (
    "gnrh_stage" %in%
    colnames(md)
  ) {

    log(
      "Stage:"
    )

    stage_tab <- table(
      md$gnrh_stage,
      useNA = "no"
    )

    for (nm in names(stage_tab)) {

      log(
        sprintf(
          "  %s: %d",
          nm,
          stage_tab[[nm]]
        )
      )
    }
  }

  if (
    "gnrh_secretory" %in%
    colnames(md)
  ) {

    log(
      "Secretory:"
    )

    secretory_tab <- table(
      md$gnrh_secretory,
      useNA = "no"
    )

    for (nm in names(secretory_tab)) {

      log(
        sprintf(
          "  %s: %d",
          nm,
          secretory_tab[[nm]]
        )
      )
    }
  }

  # --------------------------------------------------------------------------- #
  # Complete
  # --------------------------------------------------------------------------- #

  log(
    "==== GnRHcell PIPELINE COMPLETE ====",
    type = "done",
    duration = step_times$total_sec
  )

  object
}
