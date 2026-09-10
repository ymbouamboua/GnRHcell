test_that("detection exposes coherent signal and identity-supported statuses", {

  # ------------------------------------------------------------------------- #
  # Resolve test dataset
  # ------------------------------------------------------------------------- #

  installed_path <- system.file(
    "extdata",
    "hpsc.rds",
    package = "GnRHcell"
  )

  source_candidates <- c(
    file.path("inst", "extdata", "hpsc.rds"),
    testthat::test_path("..", "..", "inst", "extdata", "hpsc.rds")
  )

  source_path <- source_candidates[
    file.exists(source_candidates)
  ]

  path <- if (
    nzchar(installed_path) &&
    file.exists(installed_path)
  ) {
    installed_path
  } else if (length(source_path)) {
    source_path[[1L]]
  } else {
    ""
  }

  if (!nzchar(path) || !file.exists(path)) {
    testthat::fail(
      paste0(
        "Test dataset `hpsc.rds` was not found in ",
        "`inst/extdata` or in the installed GnRHcell package."
      )
    )
  }

  # ------------------------------------------------------------------------- #
  # Load test dataset
  # ------------------------------------------------------------------------- #

  object <- readRDS(path)

  testthat::expect_s4_class(
    object,
    "Seurat"
  )

  # ------------------------------------------------------------------------- #
  # Run GnRH detection
  # ------------------------------------------------------------------------- #

  object <- detect_gnrh(
    object,
    verbose = FALSE
  )

  metadata <- object[[]]

  # ------------------------------------------------------------------------- #
  # Detection parameters
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    "gnrh_params" %in%
      names(object@misc)
  )

  min_umi <-
    object@misc$gnrh_params$min_umi

  testthat::expect_true(
    length(min_umi) == 1L &&
      is.finite(min_umi) &&
      min_umi > 0
  )

  # ------------------------------------------------------------------------- #
  # Required detection metadata
  # ------------------------------------------------------------------------- #

  required_columns <- c(
    "gnrh_status",
    "gnrh_class",
    "gnrh_raw",

    "gnrh_direct_signal",
    "gnrh_direct_supported",
    "gnrh_direct_isolated",

    "gnrh_identity_moderate"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )

  if (length(missing_columns)) {
    testthat::fail(
      paste0(
        "Missing detection metadata: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }

  # ------------------------------------------------------------------------- #
  # Normalize logical metadata
  # ------------------------------------------------------------------------- #

  direct_signal <-
    metadata$gnrh_direct_signal %in% TRUE

  direct_supported <-
    metadata$gnrh_direct_supported %in% TRUE

  direct_isolated <-
    metadata$gnrh_direct_isolated %in% TRUE

  identity_moderate <-
    metadata$gnrh_identity_moderate %in% TRUE

  status <-
    as.character(
      metadata$gnrh_status
    )

  class <-
    as.character(
      metadata$gnrh_class
    )

  raw <-
    as.numeric(
      metadata$gnrh_raw
    )

  observed_positive <-
    status == "pos"

  expected_positive <-
    class %in%
    c(
      "direct",
      "supported"
    )

  # ------------------------------------------------------------------------- #
  # Positive calls require detectable GNRH1
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      !observed_positive |
        raw > 0,
      na.rm = TRUE
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct calls require direct GNRH1 signal
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      class != "direct" |
        direct_signal,
      na.rm = TRUE
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct calls require identity-supported direct evidence
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      class != "direct" |
        direct_supported,
      na.rm = TRUE
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct-supported cells are direct-signal cells
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      !direct_supported |
        direct_signal
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct-isolated cells are direct-signal cells
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      !direct_isolated |
        direct_signal
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct-supported and direct-isolated are mutually exclusive
  # ------------------------------------------------------------------------- #

  testthat::expect_false(
    any(
      direct_supported &
        direct_isolated
    )
  )

  # ------------------------------------------------------------------------- #
  # Direct signal is partitioned into supported or isolated evidence
  # ------------------------------------------------------------------------- #

  testthat::expect_identical(
    direct_signal,
    direct_supported |
      direct_isolated
  )

  # ------------------------------------------------------------------------- #
  # Direct-isolated signals remain GnRH-negative
  # ------------------------------------------------------------------------- #

  testthat::expect_false(
    any(
      direct_isolated &
        observed_positive
    )
  )

  # ------------------------------------------------------------------------- #
  # Supported calls require detectable but sub-threshold GNRH1
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      class != "supported" |
        (
          raw > 0 &
            raw < min_umi
        ),
      na.rm = TRUE
    )
  )

  # ------------------------------------------------------------------------- #
  # Supported calls must not overlap direct signal
  # ------------------------------------------------------------------------- #

  testthat::expect_false(
    any(
      class == "supported" &
        direct_signal,
      na.rm = TRUE
    )
  )

  # ------------------------------------------------------------------------- #
  # Positive status corresponds exactly to positive detection classes
  # ------------------------------------------------------------------------- #

  testthat::expect_identical(
    observed_positive,
    expected_positive
  )

  # ------------------------------------------------------------------------- #
  # Positive cells require moderate GnRH identity evidence
  # ------------------------------------------------------------------------- #

  testthat::expect_true(
    all(
      !observed_positive |
        identity_moderate
    )
  )

  # ------------------------------------------------------------------------- #
  # Transcriptomic candidates remain diagnostic only
  # ------------------------------------------------------------------------- #

  candidate_col <- if (
    "gnrh_transcriptomic_candidate" %in%
    colnames(metadata)
  ) {
    "gnrh_transcriptomic_candidate"

  } else if (
    "gnrh_dropout_candidate" %in%
    colnames(metadata)
  ) {
    "gnrh_dropout_candidate"

  } else {
    NULL
  }

  if (!is.null(candidate_col)) {

    candidate <-
      metadata[[candidate_col]] %in%
      TRUE

    # Transcriptomic candidates cannot be GnRH-positive.
    testthat::expect_false(
      any(
        candidate &
          observed_positive
      )
    )

    # Transcriptomic candidates must have no detected GNRH1.
    testthat::expect_true(
      all(
        !candidate |
          raw == 0,
        na.rm = TRUE
      )
    )
  }

  # ------------------------------------------------------------------------- #
  # Classification levels remain stable
  # ------------------------------------------------------------------------- #

  testthat::expect_setequal(
    levels(
      metadata$gnrh_status
    ),
    c(
      "neg",
      "pos"
    )
  )

  testthat::expect_setequal(
    levels(
      metadata$gnrh_class
    ),
    c(
      "neg",
      "supported",
      "direct"
    )
  )
})
