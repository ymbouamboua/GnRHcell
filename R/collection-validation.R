  #' Validate a GnRHcell multi-dataset collection
  #'
  #' Performs descriptive cross-dataset validation of GnRH detection,
  #' developmental staging, transcriptomic rescue, migration refinement,
  #' secretory phenotype, and biological evidence.
  #'
  #' @param collection A \code{"gnrh_collection"}.
  #' @param output_dir Validation output directory.
  #' @param positive_classes GnRH-positive evidence classes. Default is
  #'   \code{c("direct", "transcriptomic")}.
  #' @param validation_markers Optional biological validation genes.
  #' @param assay Assay used for expression summaries.
  #' @param layer Expression layer.
  #' @param allow_saved_objects Load saved processed objects when they are not
  #'   retained in memory.
  #' @param write_output Write validation tables.
  #' @param verbose Print progress.
  #'
  #' @return Object of class \code{"gnrh_validation"}.
  #'
  #' @export
  validate_gnrh_collection <- function(
    collection,
    output_dir = file.path(collection$output_dir, "validation"),
    positive_classes = c("direct", "transcriptomic"),
    validation_markers = NULL,
    assay = "RNA",
    layer = "data",
    allow_saved_objects = TRUE,
    write_output = TRUE,
    verbose = TRUE
  ) {
    # ========================================================================= #
    # Validation
    # ========================================================================= #
    if (!inherits(collection, "gnrh_collection")) {
      stop("`collection` must inherit from `gnrh_collection`.", call. = FALSE)
    }
    if (is.null(collection$results) || !length(collection$results)) {
      stop("`collection$results` is empty.", call. = FALSE)
    }
    positive_classes <- unique(as.character(positive_classes))
    if (!length(positive_classes) || anyNA(positive_classes) || any(!nzchar(positive_classes))) {
      stop("`positive_classes` must contain valid class names.", call. = FALSE)
    }
    log <- .msg(verbose)
    log("GNRH COLLECTION VALIDATION", type = "header")
    # ========================================================================= #
    # Helpers
    # ========================================================================= #
    join_metadata <- function(x, metadata) {
      if (!is.data.frame(x) || !nrow(x) || !"id" %in% colnames(x)) {
        return(x)
      }
      dplyr::left_join(x, metadata, by = "id")
    }
    count_true <- function(md, column) {
      if (!column %in% colnames(md)) {
        return(NA_integer_)
      }
      x <- md[[column]]
      if (is.logical(x)) {
        return(sum(x %in% TRUE, na.rm = TRUE))
      }
      if (is.numeric(x)) {
        return(sum(is.finite(x) & x > 0, na.rm = TRUE))
      }
      x <- tolower(trimws(as.character(x)))
      sum(x %in% c("true", "t", "1", "yes", "y", "pos", "positive"), na.rm = TRUE)
    }
    safe_median <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::median(x)
    }
    safe_quantile <- function(x, p) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::quantile(x, p, names = FALSE)
    }
    # ========================================================================= #
    # Dataset metadata
    # ========================================================================= #
    datasets <- tibble::as_tibble(collection$datasets)
    required <- c("id", "label", "species")
    missing <- setdiff(required, colnames(datasets))
    if (length(missing)) {
      stop("Missing dataset metadata: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    datasets$id <- as.character(datasets$id)
    datasets$label <- as.character(datasets$label)
    datasets$species <- as.character(datasets$species)
    if (anyNA(datasets$id) || any(!nzchar(datasets$id)) || anyDuplicated(datasets$id)) {
      stop("`collection$datasets$id` must contain unique non-empty identifiers.", call. = FALSE)
    }
    dataset_metadata <- datasets |>
      dplyr::select(.data$id, .data$label, .data$species)
    # ========================================================================= #
    # Recover processed objects
    # ========================================================================= #
    gnrh_list <- .get_collection_objects(
      collection,
      allow_saved = allow_saved_objects
    )
    missing_objects <- names(gnrh_list)[
      vapply(gnrh_list, is.null, logical(1))
    ]
    if (length(missing_objects)) {
      stop(
        "Processed objects are unavailable for: ",
        paste(missing_objects, collapse = ", "),
        ". Retain them with `clean_objects = FALSE` or use `save_objects = TRUE`.",
        call. = FALSE
      )
    }
    invalid <- !vapply(
      gnrh_list,
      inherits,
      logical(1),
      what = "Seurat"
    )
    if (any(invalid)) {
      stop(
        "Invalid processed object(s): ",
        paste(names(gnrh_list)[invalid], collapse = ", "),
        call. = FALSE
      )
    }
    # ========================================================================= #
    # Required GnRH metadata
    # ========================================================================= #
    required_gnrh <- c(
      "gnrh_status",
      "gnrh_class",
      "gnrh_stage"
    )
    for (id in names(gnrh_list)) {
      md <- gnrh_list[[id]][[]]
      missing <- setdiff(required_gnrh, colnames(md))
      if (length(missing)) {
        stop(
          "Dataset `", id, "` is missing: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
    }
    observed_classes <- unique(unlist(
      lapply(
        gnrh_list,
        function(object) unique(as.character(object[[]]$gnrh_class))
      ),
      use.names = FALSE
    ))
    observed_classes <- observed_classes[
      !is.na(observed_classes) & nzchar(observed_classes)
    ]
    missing_positive_classes <- setdiff(positive_classes, observed_classes)
    if (length(missing_positive_classes)) {
      warning(
        "Positive class(es) not observed in this collection: ",
        paste(missing_positive_classes, collapse = ", "),
        ". Observed classes: ",
        paste(sort(observed_classes), collapse = ", "),
        call. = FALSE
      )
    }
    # ========================================================================= #
    # Validation marker panel
    # ========================================================================= #
    if (is.null(validation_markers)) {
      validation_markers <- unique(c(
        "GNRH1",
        "FEZF1","ISL1","OTX2","SIX3","SIX6","ECEL1",
        "ANOS1","PROKR2","NSMF","ROBO3","SEMA3C","SEMA3F","CXCR4",
        "KISS1R","GNRHR","DOC2B","PTPRN","HCN1",
        "PCSK1","PCSK2","CHGA","CHGB","CPE","SCG2","SCG5","VGF"
      ))
    }
    validation_markers <- unique(as.character(validation_markers))
    # ========================================================================= #
    # Input summary
    # ========================================================================= #
    input_summary <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        tibble::tibble(
          id = id,
          n_cells = ncol(object),
          n_features = nrow(object)
        )
      }
    ) |>
      dplyr::left_join(dataset_metadata, by = "id") |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        .data$n_cells,
        .data$n_features
      )
    # ========================================================================= #
    # Detection
    # ========================================================================= #
    detection <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        class <- as.character(md$gnrh_class)
        status <- as.character(md$gnrh_status)
        total <- nrow(md)
        direct <- sum(class == "direct", na.rm = TRUE)
        transcriptomic <- sum(class == "transcriptomic", na.rm = TRUE)
        supported <- sum(class == "supported", na.rm = TRUE)
        positive_class <- sum(class %in% positive_classes, na.rm = TRUE)
        positive_status <- sum(status == "pos", na.rm = TRUE)
        direct_signal <- count_true(md, "gnrh_direct_signal")
        direct_supported <- count_true(md, "gnrh_direct_supported")
        direct_isolated <- count_true(md, "gnrh_direct_isolated")
        reference_positive <- count_true(md, "gnrh_reference_positive")
        confident <- count_true(md, "gnrh_confident")
        candidate_column <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          "gnrh_transcriptomic_candidate"
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          "gnrh_dropout_candidate"
        } else {
          NULL
        }
        transcriptomic_candidates <- if (!is.null(candidate_column)) {
          count_true(md, candidate_column)
        } else {
          NA_integer_
        }
        tibble::tibble(
          id = id,
          n_cells = total,
          direct = direct,
          transcriptomic = transcriptomic,
          supported = supported,
          gnrh_pos = positive_status,
          positive_by_class = positive_class,
          direct_signal = direct_signal,
          direct_supported = direct_supported,
          direct_isolated = direct_isolated,
          reference_positive = reference_positive,
          transcriptomic_candidates = transcriptomic_candidates,
          gnrh_confident = confident,
          pct_gnrh = if (total > 0L) 100 * positive_status / total else NA_real_,
          pct_direct = if (positive_status > 0L) 100 * direct / positive_status else NA_real_,
          pct_transcriptomic = if (positive_status > 0L) 100 * transcriptomic / positive_status else NA_real_,
          pct_supported = if (positive_status > 0L) 100 * supported / positive_status else NA_real_,
          pct_confident = if (positive_status > 0L && !is.na(confident)) 100 * confident / positive_status else NA_real_,
          pct_direct_isolated = if (!is.na(direct_signal) && direct_signal > 0L && !is.na(direct_isolated)) {
            100 * direct_isolated / direct_signal
          } else {
            NA_real_
          },
          pct_transcriptomic_candidates = if (!is.na(transcriptomic_candidates) && total > 0L) {
            100 * transcriptomic_candidates / total
          } else {
            NA_real_
          }
        )
      }
    ) |>
      dplyr::left_join(dataset_metadata, by = "id") |>
      dplyr::select(
        .data$id,
        .data$label,
        .data$species,
        dplyr::everything()
      )
    # ========================================================================= #
    # Transcriptomic candidates
    # ========================================================================= #
    transcriptomic_candidates <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        column <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          "gnrh_transcriptomic_candidate"
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          "gnrh_dropout_candidate"
        } else {
          return(tibble::tibble())
        }
        x <- md[
          md[[column]] %in% TRUE,
          ,
          drop = FALSE
        ]
        median_col <- function(column) {
          if (!column %in% colnames(x) || !nrow(x)) return(NA_real_)
          safe_median(x[[column]])
        }
        tibble::tibble(
          id = id,
          n_cells = nrow(x),
          pct_dataset = if (nrow(md)) 100 * nrow(x) / nrow(md) else NA_real_,
          pct_GNRH1_detected = if (nrow(x) && "gnrh_raw" %in% colnames(x)) {
            100 * mean(x$gnrh_raw > 0, na.rm = TRUE)
          } else {
            0
          },
          median_GNRH1 = median_col("gnrh_raw"),
          median_support = median_col("gnrh_support_score_raw"),
          median_identity_primary = median_col("gnrh_identity_primary_hits"),
          median_core = median_col("gnrh_core_hits"),
          median_neuro = median_col("gnrh_neuro_hits"),
          median_knn = median_col("gnrh_knn"),
          median_alternative = median_col("gnrh_alternative_score")
        )
      }
    )
    transcriptomic_candidates <- join_metadata(
      transcriptomic_candidates,
      dataset_metadata
    )
    # ========================================================================= #
    # Status / class consistency
    # ========================================================================= #
    status_class <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        object[[]] |>
          tibble::as_tibble() |>
          dplyr::transmute(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_class = as.character(.data$gnrh_class)
          ) |>
          dplyr::count(
            .data$gnrh_status,
            .data$gnrh_class,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_cells = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    status_class <- join_metadata(status_class, dataset_metadata)
    # ========================================================================= #
    # Logical consistency
    # ========================================================================= #
    classification_consistency <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        status <- as.character(md$gnrh_status)
        class <- as.character(md$gnrh_class)
        raw <- if ("gnrh_raw" %in% colnames(md)) {
          suppressWarnings(as.numeric(md$gnrh_raw))
        } else {
          rep(NA_real_, nrow(md))
        }
        candidate <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) {
          md$gnrh_transcriptomic_candidate %in% TRUE
        } else if ("gnrh_dropout_candidate" %in% colnames(md)) {
          md$gnrh_dropout_candidate %in% TRUE
        } else {
          rep(FALSE, nrow(md))
        }
        expected_positive_class <- class %in% positive_classes
        tibble::tibble(
          id = id,
          n_status_class_discordant_positive = sum(
            status == "pos" & !expected_positive_class,
            na.rm = TRUE
          ),
          n_status_class_discordant_negative = sum(
            status == "neg" & expected_positive_class,
            na.rm = TRUE
          ),
          n_unexplained_positive_without_GNRH1 = sum(
            status == "pos" &
              is.finite(raw) &
              raw <= 0 &
              !candidate,
            na.rm = TRUE
          ),
          n_expected_transcriptomic_positive = sum(
            status == "pos" & candidate,
            na.rm = TRUE
          ),
          n_direct_without_GNRH1 = sum(
            class == "direct" &
              is.finite(raw) &
              raw <= 0,
            na.rm = TRUE
          )
        )
      }
    )
    classification_consistency <- join_metadata(
      classification_consistency,
      dataset_metadata
    )
    # ========================================================================= #
    # Evidence scores
    # ========================================================================= #
    candidate_score_columns <- unique(c(
      "gnrh_score",
      "gnrh_score_raw",
      "gnrh_support_score",
      "gnrh_support_score_raw",
      "gnrh_identity_score",
      "gnrh_migration_score",
      "gnrh_neuro_score",
      "gnrh_hormone_score",
      "gnrh_guidance_score",
      "gnrh_alternative_score",
      "gnrh_identity_primary_hits",
      "gnrh_identity_supportive_hits",
      "gnrh_core_hits",
      "gnrh_migration_primary_hits",
      "gnrh_migration_supportive_hits",
      "gnrh_mig_hits",
      "gnrh_neuro_primary_hits",
      "gnrh_neuro_supportive_hits",
      "gnrh_neuro_hits",
      "gnrh_migration_core_hits",
      "gnrh_stage_early_score",
      "gnrh_stage_migrating_score",
      "gnrh_stage_mature_score",
      "gnrh_knn"
    ))
    available_columns <- Reduce(
      intersect,
      lapply(
        gnrh_list,
        function(object) colnames(object[[]])
      )
    )
    score_columns <- intersect(
      candidate_score_columns,
      available_columns
    )
    scores <- if (length(score_columns)) {
      purrr::imap_dfr(
        gnrh_list,
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::mutate(
              gnrh_class = as.character(.data$gnrh_class)
            ) |>
            dplyr::filter(
              .data$gnrh_class %in% positive_classes
            ) |>
            dplyr::group_by(
              .data$gnrh_class
            ) |>
            dplyr::summarise(
              n_cells = dplyr::n(),
              dplyr::across(
                dplyr::all_of(score_columns),
                list(
                  median = ~ safe_median(.x),
                  q25 = ~ safe_quantile(.x, 0.25),
                  q75 = ~ safe_quantile(.x, 0.75)
                )
              ),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
    } else {
      tibble::tibble()
    }
    scores <- join_metadata(scores, dataset_metadata)
    # ========================================================================= #
    # Developmental stages
    # ========================================================================= #
    stages <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_stage = as.character(.data$gnrh_stage)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_stage),
            .data$gnrh_stage != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
        x |>
          dplyr::count(
            .data$gnrh_stage,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_positive = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    stages <- join_metadata(stages, dataset_metadata)
    stage_class <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        x <- object[[]] |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_class = as.character(.data$gnrh_class),
            gnrh_stage = as.character(.data$gnrh_stage)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_stage),
            .data$gnrh_stage != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
        x |>
          dplyr::count(
            .data$gnrh_class,
            .data$gnrh_stage,
            name = "n_cells"
          ) |>
          dplyr::group_by(
            .data$gnrh_class
          ) |>
          dplyr::mutate(
            pct_class = 100 * .data$n_cells / sum(.data$n_cells)
          ) |>
          dplyr::ungroup() |>
          dplyr::mutate(id = id)
      }
    )
    stage_class <- join_metadata(stage_class, dataset_metadata)
    # ========================================================================= #
    # Secretory phenotype
    # ========================================================================= #
    secretory <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        if (!"gnrh_secretory" %in% colnames(md)) {
          return(tibble::tibble())
        }
        x <- md |>
          tibble::as_tibble() |>
          dplyr::mutate(
            gnrh_status = as.character(.data$gnrh_status),
            gnrh_secretory = as.character(.data$gnrh_secretory)
          ) |>
          dplyr::filter(
            .data$gnrh_status == "pos",
            !is.na(.data$gnrh_secretory),
            .data$gnrh_secretory != "non-gnrh"
          )
        if (!nrow(x)) return(tibble::tibble())
        x |>
          dplyr::count(
            .data$gnrh_secretory,
            name = "n_cells"
          ) |>
          dplyr::mutate(
            id = id,
            pct_positive = 100 * .data$n_cells / sum(.data$n_cells)
          )
      }
    )
    secretory <- join_metadata(secretory, dataset_metadata)
    # ========================================================================= #
    # Stage refinement
    # ========================================================================= #
    refinement_required <- c(
      "gnrh_stage_raw",
      "gnrh_stage_reassigned",
      "gnrh_stage_reason"
    )
    has_refinement <- vapply(
      gnrh_list,
      function(object) {
        all(refinement_required %in% colnames(object[[]]))
      },
      logical(1)
    )
    stage_refinement <- tibble::tibble()
    stage_reassignment <- tibble::tibble()
    if (any(has_refinement)) {
      stage_refinement <- purrr::imap_dfr(
        gnrh_list[has_refinement],
        function(object, id) {
          md <- object[[]]
          positive <- as.character(md$gnrh_status) == "pos"
          n_positive <- sum(positive, na.rm = TRUE)
          n_reassigned <- sum(
            positive & md$gnrh_stage_reassigned %in% TRUE,
            na.rm = TRUE
          )
          n_filtered <- sum(
            positive &
              as.character(md$gnrh_stage_reason) == "migration_core_filter",
            na.rm = TRUE
          )
          tibble::tibble(
            id = id,
            n_positive = n_positive,
            n_reassigned = n_reassigned,
            pct_reassigned = if (n_positive > 0L) {
              100 * n_reassigned / n_positive
            } else {
              NA_real_
            },
            n_migration_filtered = n_filtered
          )
        }
      )
      stage_refinement <- join_metadata(
        stage_refinement,
        dataset_metadata
      )
      stage_reassignment <- purrr::imap_dfr(
        gnrh_list[has_refinement],
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::filter(
              as.character(.data$gnrh_status) == "pos"
            ) |>
            dplyr::transmute(
              id = id,
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw),
              gnrh_stage = as.character(.data$gnrh_stage),
              gnrh_stage_reason = as.character(.data$gnrh_stage_reason)
            ) |>
            dplyr::count(
              .data$id,
              .data$gnrh_stage_raw,
              .data$gnrh_stage,
              .data$gnrh_stage_reason,
              name = "n_cells"
            )
        }
      )
      stage_reassignment <- join_metadata(
        stage_reassignment,
        dataset_metadata
      )
    }
    # ========================================================================= #
    # Migration refinement
    # ========================================================================= #
    migration_required <- c(
      "gnrh_stage_raw",
      "gnrh_stage",
      "gnrh_migration_core_hits"
    )
    has_migration <- vapply(
      gnrh_list,
      function(object) {
        all(migration_required %in% colnames(object[[]]))
      },
      logical(1)
    )
    migration_core <- tibble::tibble()
    migration_refinement <- tibble::tibble()
    migration_refinement_summary <- tibble::tibble()
    if (any(has_migration)) {
      migration_core <- purrr::imap_dfr(
        gnrh_list[has_migration],
        function(object, id) {
          object[[]] |>
            tibble::as_tibble() |>
            dplyr::filter(
              as.character(.data$gnrh_status) == "pos"
            ) |>
            dplyr::mutate(
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw)
            ) |>
            dplyr::group_by(
              .data$gnrh_stage_raw
            ) |>
            dplyr::summarise(
              n = dplyr::n(),
              median_hits = safe_median(.data$gnrh_migration_core_hits),
              q25_hits = safe_quantile(.data$gnrh_migration_core_hits, 0.25),
              q75_hits = safe_quantile(.data$gnrh_migration_core_hits, 0.75),
              pct_ge1 = 100 * mean(.data$gnrh_migration_core_hits >= 1, na.rm = TRUE),
              pct_ge2 = 100 * mean(.data$gnrh_migration_core_hits >= 2, na.rm = TRUE),
              pct_ge3 = 100 * mean(.data$gnrh_migration_core_hits >= 3, na.rm = TRUE),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
      migration_core <- join_metadata(
        migration_core,
        dataset_metadata
      )
      migration_refinement <- purrr::imap_dfr(
        gnrh_list[has_migration],
        function(object, id) {
          x <- object[[]] |>
            tibble::as_tibble() |>
            dplyr::mutate(
              gnrh_status = as.character(.data$gnrh_status),
              gnrh_stage_raw = as.character(.data$gnrh_stage_raw),
              gnrh_stage = as.character(.data$gnrh_stage)
            ) |>
            dplyr::filter(
              .data$gnrh_status == "pos",
              .data$gnrh_stage_raw == "migrating"
            ) |>
            dplyr::mutate(
              migration_outcome = ifelse(
                .data$gnrh_stage == "migrating",
                "retained_migrating",
                "reassigned"
              )
            )
          if (!nrow(x)) return(tibble::tibble())
          x |>
            dplyr::group_by(
              .data$migration_outcome,
              .data$gnrh_stage
            ) |>
            dplyr::summarise(
              n_cells = dplyr::n(),
              median_hits = safe_median(.data$gnrh_migration_core_hits),
              mean_hits = mean(.data$gnrh_migration_core_hits, na.rm = TRUE),
              pct_zero = 100 * mean(.data$gnrh_migration_core_hits == 0, na.rm = TRUE),
              pct_ge1 = 100 * mean(.data$gnrh_migration_core_hits >= 1, na.rm = TRUE),
              pct_ge2 = 100 * mean(.data$gnrh_migration_core_hits >= 2, na.rm = TRUE),
              pct_ge3 = 100 * mean(.data$gnrh_migration_core_hits >= 3, na.rm = TRUE),
              .groups = "drop"
            ) |>
            dplyr::mutate(id = id)
        }
      )
      migration_refinement <- join_metadata(
        migration_refinement,
        dataset_metadata
      )
      if (nrow(migration_refinement)) {
        migration_refinement_summary <- migration_refinement |>
          dplyr::group_by(
            .data$id,
            .data$label,
            .data$species
          ) |>
          dplyr::summarise(
            n_raw_migrating = sum(.data$n_cells, na.rm = TRUE),
            n_retained = sum(
              .data$n_cells[.data$migration_outcome == "retained_migrating"],
              na.rm = TRUE
            ),
            n_reassigned = sum(
              .data$n_cells[.data$migration_outcome == "reassigned"],
              na.rm = TRUE
            ),
            n_to_early = sum(
              .data$n_cells[
                .data$migration_outcome == "reassigned" &
                  .data$gnrh_stage == "early"
              ],
              na.rm = TRUE
            ),
            n_to_mature = sum(
              .data$n_cells[
                .data$migration_outcome == "reassigned" &
                  .data$gnrh_stage == "mature"
              ],
              na.rm = TRUE
            ),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            pct_retained = dplyr::if_else(
              .data$n_raw_migrating > 0,
              100 * .data$n_retained / .data$n_raw_migrating,
              NA_real_
            ),
            pct_reassigned = dplyr::if_else(
              .data$n_raw_migrating > 0,
              100 * .data$n_reassigned / .data$n_raw_migrating,
              NA_real_
            )
          )
      }
    }
    # ========================================================================= #
    # Biological marker validation
    # ========================================================================= #
    biological_markers <- purrr::imap_dfr(
      gnrh_list,
      function(object, id) {
        md <- object[[]]
        cells <- rownames(md)[
          as.character(md$gnrh_status) == "pos"
        ]
        if (!length(cells)) return(tibble::tibble())
        expr <- tryCatch(
          .get_expr(
            object,
            assay = assay,
            layer = layer
          ),
          error = function(e) NULL
        )
        if (is.null(expr)) return(tibble::tibble())
        matched <- .match_genes(
          validation_markers,
          rownames(expr)
        )
        if (!length(matched)) return(tibble::tibble())
        cells <- intersect(cells, colnames(expr))
        if (!length(cells)) return(tibble::tibble())
        classes <- as.character(
          md[cells, "gnrh_class", drop = TRUE]
        )
        observed_positive_classes <- intersect(
          positive_classes,
          unique(classes)
        )
        purrr::map_dfr(
          observed_positive_classes,
          function(cl) {
            idx <- classes == cl
            if (!any(idx)) return(tibble::tibble())
            x <- expr[
              matched,
              cells[idx],
              drop = FALSE
            ]
            tibble::tibble(
              id = id,
            gnrh_class = cl,
            gene = rownames(x),
            avg_expression = Matrix::rowMeans(x),
            pct_expressing = 100 * Matrix::rowMeans(x > 0)
          )
        }
      )
    }
  )
  biological_markers <- join_metadata(
    biological_markers,
    dataset_metadata
  )
  # ========================================================================= #
  # Validation outputs
  # ========================================================================= #
  validation_tables <- list(
    input_summary = input_summary,
    detection = detection,
    transcriptomic_candidates = transcriptomic_candidates,
    classification_consistency = classification_consistency,
    status_class = status_class,
    scores = scores,
    stages = stages,
    stage_class = stage_class,
    secretory = secretory,
    stage_refinement = stage_refinement,
    stage_reassignment = stage_reassignment,
    migration_core = migration_core,
    migration_refinement = migration_refinement,
    migration_refinement_summary = migration_refinement_summary,
    biological_markers = biological_markers
  )
  # ========================================================================= #
  # Write outputs
  # ========================================================================= #
  if (isTRUE(write_output)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    purrr::iwalk(
      validation_tables,
      function(x, name) {
        if (is.data.frame(x) && ncol(x)) {
          readr::write_csv(
            x,
            file.path(output_dir, paste0(name, ".csv"))
          )
        }
      }
    )
  }
  # ========================================================================= #
  # Console summary
  # ========================================================================= #
  total_cells <- sum(detection$n_cells, na.rm = TRUE)
  total_positive <- sum(detection$gnrh_pos, na.rm = TRUE)
  has_candidate_counts <- any(
    is.finite(detection$transcriptomic_candidates)
  )
  total_candidates <- if (has_candidate_counts) {
    sum(detection$transcriptomic_candidates, na.rm = TRUE)
  } else {
    NA_integer_
  }
  log(
    sprintf(
      "Validated %d datasets | %s cells | %s GnRH-positive",
      length(gnrh_list),
      format(total_cells, big.mark = ",", trim = TRUE),
      format(total_positive, big.mark = ",", trim = TRUE)
    ),
    type = "info"
  )
  if (is.finite(total_candidates)) {
    log(
      sprintf(
        "Transcriptomic candidates: %s",
        format(total_candidates, big.mark = ",", trim = TRUE)
      ),
      type = "info"
    )
  }
  # ========================================================================= #
  # Consistency summary
  # ========================================================================= #
  violation_columns <- intersect(
    c(
      "n_status_class_discordant_positive",
      "n_status_class_discordant_negative",
      "n_unexplained_positive_without_GNRH1",
      "n_direct_without_GNRH1"
    ),
    colnames(classification_consistency)
  )
  bad_consistency <- if (length(violation_columns)) {
    rowSums(
      as.data.frame(
        classification_consistency[
          ,
          violation_columns,
          drop = FALSE
        ]
      ),
      na.rm = TRUE
    )
  } else {
    rep(0, nrow(classification_consistency))
  }
  if (any(bad_consistency > 0)) {
    log(
      "Classification consistency violations detected.",
      type = "warn"
    )
  } else {
    log(
      "Classification consistency checks passed.",
      type = "done",
      duration = 0
    )
  }
  # ========================================================================= #
  # Return
  # ========================================================================= #
  result <- c(
    list(
      datasets = dataset_metadata
    ),
    validation_tables,
    list(
      parameters = list(
        positive_classes = positive_classes,
        assay = assay,
        layer = layer,
        validation_markers = validation_markers
      ),
      output_dir = if (isTRUE(write_output)) {
        normalizePath(
          output_dir,
          mustWork = FALSE
        )
      } else {
        NULL
      }
    )
  )
  structure(
    result,
    class = "gnrh_validation"
  )
}

