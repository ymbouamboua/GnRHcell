#' Run complete GnRHcell analysis pipeline
#'
#' Executes the complete GnRHcell workflow for GnRH detection,
#' developmental staging, diagnostics, and runtime reporting.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param detect Logical. Run GnRH detection. Default is \code{TRUE}.
#' @param stage Logical. Run developmental staging. Default is \code{TRUE}.
#' @param diagnostics Logical. Run diagnostic analysis. Default is \code{TRUE}.
#' @param detect_args Named list of additional arguments passed to
#'   \code{\link{detect_gnrh}}.
#' @param stage_args Named list of additional arguments passed to
#'   \code{\link{stage_gnrh}}.
#' @param diagnostic_args Named list of additional arguments passed to
#'   \code{\link{gnrh_diagnostics}}.
#' @param summary_level Character scalar controlling the final pipeline summary.
#'   \code{"concise"} prints four compact lines, \code{"detailed"} prints the
#'   full diagnostic summary, and \code{"none"} suppresses the summary.
#'   Default is \code{"concise"}.
#' @param verbose Logical. Print progress messages. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{detect_gnrh}} for
#'   backward compatibility. Values in \code{detect_args} take precedence.
#'
#' @return A Seurat object updated with GnRH detection, developmental staging,
#' diagnostics, and runtime metadata.
#'
#' @details
#' Detection, developmental staging, and diagnostics can be configured
#' independently through \code{detect_args}, \code{stage_args}, and
#' \code{diagnostic_args}.
#'
#' When developmental staging is available, the pipeline summary reports both
#' dominant developmental stages and stage-resolution status
#' (\code{resolved} or \code{transitional}).
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
    detect_args = list(),
    stage_args = list(),
    diagnostic_args = list(),
    summary_level = c("concise", "detailed", "none"),
    verbose = TRUE,
    ...
) {

  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  summary_level <- match.arg(summary_level)

  # ------------------------------------------------------------------------- #
  # Validate argument lists
  # ------------------------------------------------------------------------- #

  for (x in c(
    "detect_args",
    "stage_args",
    "diagnostic_args"
  )) {
    value <- get(x)

    if (!is.list(value)) {
      stop(
        "`", x, "` must be a named list.",
        call. = FALSE
      )
    }

    if (
      length(value) &&
      (
        is.null(names(value)) ||
        any(!nzchar(names(value)))
      )
    ) {
      stop(
        "`", x, "` must be a named list.",
        call. = FALSE
      )
    }
  }

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

  # ------------------------------------------------------------------------- #
  # Resolve detection arguments
  # ------------------------------------------------------------------------- #

  legacy_detect_args <- list(...)

  duplicate_detect_args <- intersect(
    names(legacy_detect_args),
    names(detect_args)
  )

  if (length(duplicate_detect_args)) {
    legacy_detect_args[
      duplicate_detect_args
    ] <- NULL
  }

  detect_call_args <- c(
    legacy_detect_args,
    detect_args
  )

  detect_call_args$verbose <- NULL
  stage_args$verbose <- NULL
  diagnostic_args$verbose <- NULL

  # ------------------------------------------------------------------------- #
  # Determine active steps
  # ------------------------------------------------------------------------- #

  active_steps <- c(
    detect = isTRUE(detect),
    stage = isTRUE(stage),
    diagnostics = isTRUE(diagnostics)
  )

  n_steps <- sum(
    active_steps
  )

  current_step <- 0L

  step_label <- function(name) {
    current_step <<-
      current_step + 1L

    sprintf(
      "[%d/%d] %s",
      current_step,
      n_steps,
      name
    )
  }

  # ------------------------------------------------------------------------- #
  # Detection
  # ------------------------------------------------------------------------- #

  if (isTRUE(detect)) {

    log(
      step_label(
        "Detecting GnRH cells"
      ),
      type = "step"
    )

    t0 <- Sys.time()

    object <- do.call(
      detect_gnrh,
      c(
        list(
          object = object,
          verbose = verbose
        ),
        detect_call_args
      )
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
      duration =
        step_times$detect_sec
    )
  }

  # ------------------------------------------------------------------------- #
  # Staging
  # ------------------------------------------------------------------------- #

  if (isTRUE(stage)) {

    if (
      !"gnrh_status" %in%
      colnames(
        object[[]]
      )
    ) {
      stop(
        "Developmental staging requires GnRH detection metadata. ",
        "Run with `detect = TRUE` or provide an object previously processed ",
        "with `detect_gnrh()`.",
        call. = FALSE
      )
    }

    log(
      step_label(
        "Assigning developmental stages"
      ),
      type = "step"
    )

    t0 <- Sys.time()

    object <- do.call(
      stage_gnrh,
      c(
        list(
          object = object,
          verbose = FALSE
        ),
        stage_args
      )
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
      duration =
        step_times$stage_sec
    )
  }

  # ------------------------------------------------------------------------- #
  # Diagnostics
  # ------------------------------------------------------------------------- #

  if (isTRUE(diagnostics)) {

    if (
      !"gnrh_status" %in%
      colnames(
        object[[]]
      )
    ) {
      stop(
        "Diagnostics require GnRH detection metadata. ",
        "Run with `detect = TRUE` or provide an object previously processed ",
        "with `detect_gnrh()`.",
        call. = FALSE
      )
    }

    log(
      step_label(
        "Running diagnostics"
      ),
      type = "step"
    )

    t0 <- Sys.time()

    object <- do.call(
      gnrh_diagnostics,
      c(
        list(
          object = object,
          verbose = FALSE
        ),
        diagnostic_args
      )
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
      duration =
        step_times$diagnostics_sec
    )
  }

  # ------------------------------------------------------------------------- #
  # Runtime
  # ------------------------------------------------------------------------- #

  step_times$total_sec <- as.numeric(
    difftime(
      Sys.time(),
      total_start,
      units = "secs"
    )
  )

  # ------------------------------------------------------------------------- #
  # Store pipeline parameters
  # ------------------------------------------------------------------------- #

  params <- list(
    detect =
      isTRUE(detect),

    stage =
      isTRUE(stage),

    diagnostics =
      isTRUE(diagnostics),

    detect_args =
      detect_call_args,

    stage_args =
      stage_args,

    diagnostic_args =
      diagnostic_args
  )

  object <- .add_run_info(
    object,
    step_times = step_times,
    params = params
  )

  # ------------------------------------------------------------------------- #
  # Pipeline summary
  # ------------------------------------------------------------------------- #

  md <- object[[]]

  count_true <- function(column) {
    if (
      !column %in%
      colnames(md)
    ) {
      return(
        NA_integer_
      )
    }

    sum(
      md[[column]] %in% TRUE,
      na.rm = TRUE
    )
  }

  print_table <- function(
    x,
    title,
    exclude = NULL,
    order = NULL
  ) {
    x <- as.character(x)

    x <- x[
      !is.na(x)
    ]

    if (!is.null(exclude)) {
      x <- x[
        !x %in% exclude
      ]
    }

    if (!length(x)) {
      return(
        invisible(NULL)
      )
    }

    tab <- table(
      x,
      useNA = "no"
    )

    if (!is.null(order)) {
      ordered_names <- c(
        intersect(order, names(tab)),
        setdiff(names(tab), order)
      )
      tab <- tab[ordered_names]
    }

    log(
      paste0(
        title,
        ":"
      )
    )

    for (nm in names(tab)) {
      log(
        sprintf(
          "  %s: %d",
          nm,
          tab[[nm]]
        )
      )
    }

    invisible(tab)
  }

  if (!identical(summary_level, "none")) {
    log(
      "PIPELINE SUMMARY",
      type = "info"
    )
  }

  # ------------------------------------------------------------------------- #
  # Concise summary
  # ------------------------------------------------------------------------- #

  if (identical(summary_level, "concise")) {
    total <- nrow(md)
    positive <- if ("gnrh_status" %in% colnames(md)) {
      sum(as.character(md$gnrh_status) == "pos", na.rm = TRUE)
    } else {
      0L
    }
    confident <- count_true("gnrh_confident")
    positive_pct <- if (total > 0L) 100 * positive / total else 0

    log(sprintf(
      "  Cells: %s | GnRH+: %s (%.1f%%) | high confidence: %s",
      format(total, big.mark = ","),
      format(positive, big.mark = ","),
      positive_pct,
      format(confident, big.mark = ",")
    ))

    direct <- if ("gnrh_class" %in% colnames(md)) {
      sum(as.character(md$gnrh_class) == "direct", na.rm = TRUE)
    } else {
      count_true("gnrh_direct_supported")
    }
    transcriptomic <- if ("gnrh_class" %in% colnames(md)) {
      sum(as.character(md$gnrh_class) == "supported", na.rm = TRUE)
    } else {
      0L
    }
    isolated <- count_true("gnrh_direct_isolated")

    log(sprintf(
      "  Evidence: direct %s | transcriptomic %s | isolated GNRH1 signal %s",
      format(direct, big.mark = ","),
      format(transcriptomic, big.mark = ","),
      format(isolated, big.mark = ",")
    ))

    if ("gnrh_stage" %in% colnames(md) && positive > 0L) {
      stages <- table(factor(
        as.character(md$gnrh_stage),
        levels = c("identity", "migrating", "mature")
      ))
      log(sprintf(
        "  Stage: identity %d | migrating %d | mature %d",
        stages[["identity"]],
        stages[["migrating"]],
        stages[["mature"]]
      ))
    }

    if ("gnrh_stage_resolution" %in% colnames(md) && positive > 0L) {
      resolution <- as.character(md$gnrh_stage_resolution)
      positive_cells <- as.character(md$gnrh_status) == "pos"
      resolved <- sum(positive_cells & resolution == "resolved", na.rm = TRUE)
      transitional <- sum(
        positive_cells & resolution == "transitional",
        na.rm = TRUE
      )
      log(sprintf(
        "  Stage resolution: resolved %d (%.1f%%) | transitional %d (%.1f%%)",
        resolved,
        100 * resolved / positive,
        transitional,
        100 * transitional / positive
      ))
    }
  }

  if (identical(summary_level, "detailed")) {

  # ------------------------------------------------------------------------- #
  # Status
  # ------------------------------------------------------------------------- #

  if (
    "gnrh_status" %in%
    colnames(md)
  ) {
    print_table(
      md$gnrh_status,
      "Status"
    )
  }

  # ------------------------------------------------------------------------- #
  # Classification
  # ------------------------------------------------------------------------- #

  if (
    "gnrh_class" %in%
    colnames(md)
  ) {
    print_table(
      md$gnrh_class,
      "Class"
    )
  }

  # ------------------------------------------------------------------------- #
  # Detection evidence
  # ------------------------------------------------------------------------- #

  evidence_labels <- c(
    gnrh_direct_signal =
      "GNRH1 direct signal",

    gnrh_direct_supported =
      "identity-supported direct",

    gnrh_direct_isolated =
      "isolated GNRH1 signal",

    gnrh_reference_positive =
      "high-specificity reference",

    gnrh_transcriptomic_candidate =
      "transcriptomic candidates",

    gnrh_confident =
      "high-confidence GnRH"
  )

  available_evidence <- intersect(
    names(evidence_labels),
    colnames(md)
  )

  if (length(available_evidence)) {

    log(
      "Detection evidence:"
    )

    for (
      column in
      available_evidence
    ) {
      log(
        sprintf(
          "  %s: %d",
          evidence_labels[[column]],
          count_true(column)
        )
      )
    }
  }

  # ------------------------------------------------------------------------- #
  # Developmental stage
  # ------------------------------------------------------------------------- #

  if (
    "gnrh_stage" %in%
    colnames(md)
  ) {
    print_table(
      md$gnrh_stage,
      "Stage",
      exclude = "non-gnrh",
      order = c(
        "identity",
        "migrating",
        "mature"
      )
    )
  }

  # ------------------------------------------------------------------------- #
  # Stage resolution
  # ------------------------------------------------------------------------- #

  if (
    "gnrh_stage_resolution" %in%
    colnames(md)
  ) {

    resolution <- as.character(
      md$gnrh_stage_resolution
    )

    if (
      "gnrh_status" %in%
      colnames(md)
    ) {
      resolution <- resolution[
        !is.na(md$gnrh_status) &
          as.character(
            md$gnrh_status
          ) == "pos"
      ]
    } else {
      resolution <- resolution[
        resolution != "non-gnrh"
      ]
    }

    resolution <- resolution[
      !is.na(resolution) &
        resolution != "non-gnrh"
    ]

    if (length(resolution)) {

      resolution_tab <- table(
        resolution,
        useNA = "no"
      )

      total_resolution <- sum(
        resolution_tab
      )

      log(
        "Stage resolution:"
      )

      for (
        nm in
        names(resolution_tab)
      ) {
        n <- as.integer(
          resolution_tab[[nm]]
        )

        pct <- if (
          total_resolution > 0
        ) {
          100 *
            n /
            total_resolution
        } else {
          NA_real_
        }

        log(
          sprintf(
            "  %s: %d (%.1f%%)",
            nm,
            n,
            pct
          )
        )
      }
    }
  }

  # ------------------------------------------------------------------------- #
  # Resolution by developmental stage
  # ------------------------------------------------------------------------- #

  if (
    all(
      c(
        "gnrh_stage",
        "gnrh_stage_resolution"
      ) %in%
      colnames(md)
    )
  ) {

    stage_value <- as.character(
      md$gnrh_stage
    )

    resolution_value <- as.character(
      md$gnrh_stage_resolution
    )

    valid <-
      !is.na(stage_value) &
      !is.na(resolution_value) &
      stage_value != "non-gnrh" &
      resolution_value != "non-gnrh"

    if (any(valid)) {

      stage_resolution_tab <- table(
        stage_value[valid],
        resolution_value[valid]
      )

      resolved_column <- if (
        "resolved" %in%
        colnames(
          stage_resolution_tab
        )
      ) {
        stage_resolution_tab[
          ,
          "resolved"
        ]
      } else {
        rep(
          0,
          nrow(
            stage_resolution_tab
          )
        )
      }

      stage_totals <- rowSums(
        stage_resolution_tab
      )

      log(
        "Stage confidence:"
      )

      stage_order <- c(
        "identity",
        "migrating",
        "mature"
      )

      stage_order <- c(
        intersect(
          stage_order,
          rownames(stage_resolution_tab)
        ),
        setdiff(
          rownames(stage_resolution_tab),
          stage_order
        )
      )

      for (nm in stage_order) {

        n_resolved <- as.integer(
          resolved_column[[nm]]
        )

        total_stage <- as.integer(
          stage_totals[[nm]]
        )

        pct_resolved <- if (
          total_stage > 0
        ) {
          100 *
            n_resolved /
            total_stage
        } else {
          NA_real_
        }

        log(
          sprintf(
            "  %s: %d/%d resolved (%.1f%%)",
            nm,
            n_resolved,
            total_stage,
            pct_resolved
          )
        )
      }
    }
  }

  # ------------------------------------------------------------------------- #
  # Secretory state
  # ------------------------------------------------------------------------- #

  if (
    "gnrh_secretory" %in%
    colnames(md)
  ) {
    print_table(
      md$gnrh_secretory,
      "Secretory",
      exclude = "non-gnrh"
    )
  }

  # ------------------------------------------------------------------------- #
  # External validation
  # ------------------------------------------------------------------------- #

  validation <-
    object@misc$gnrh$validation

  if (
    !is.null(validation) &&
    is.list(validation)
  ) {

    metric <- function(name) {
      x <- validation[[name]]

      if (
        is.null(x) ||
        length(x) != 1L ||
        !is.finite(x)
      ) {
        return(
          NA_real_
        )
      }

      as.numeric(x)
    }

    auc <- metric("auc")
    auprc <- metric("auprc")
    precision <- metric("precision")
    recall <- metric("recall")
    F1 <- metric("F1")

    false_positives_per_100k <-
      metric(
        "false_positives_per_100k"
      )

    available_validation <- any(
      is.finite(
        c(
          auc,
          auprc,
          precision,
          recall,
          F1,
          false_positives_per_100k
        )
      )
    )

    if (available_validation) {

      log(
        "Validation:"
      )

      if (is.finite(auc)) {
        log(
          sprintf(
            "  ROC AUC: %.3f",
            auc
          )
        )
      }

      if (is.finite(auprc)) {
        log(
          sprintf(
            "  PR AUC: %.3f",
            auprc
          )
        )
      }

      if (is.finite(precision)) {
        log(
          sprintf(
            "  Precision: %.3f",
            precision
          )
        )
      }

      if (is.finite(recall)) {
        log(
          sprintf(
            "  Recall: %.3f",
            recall
          )
        )
      }

      if (is.finite(F1)) {
        log(
          sprintf(
            "  F1: %.3f",
            F1
          )
        )
      }

      if (
        is.finite(
          false_positives_per_100k
        )
      ) {
        log(
          sprintf(
            "  False positives / 100k: %.1f",
            false_positives_per_100k
          )
        )
      }
    }
  }
  }

  # ------------------------------------------------------------------------- #
  # Complete
  # ------------------------------------------------------------------------- #

  log(
    "==== GnRHcell PIPELINE COMPLETE ====",
    type = "done",
    duration =
      step_times$total_sec
  )

  object
}
