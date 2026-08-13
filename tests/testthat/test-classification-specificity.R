test_that("isolated GNRH1 signal is not classified as a GnRH neuron", {
  result <- GnRHcell:::.classify(
    raw = c(3, 3), norm = c(2, 2), score = c(2, 2),
    support_score = c(0, 2),
    hits = list(core = c(0L, 2L), mig = c(0L, 0L), neuro = c(0L, 0L)),
    lib = c(1000, 1000), min_umi = 2, min_counts = 500,
    expr_thr = 1, identity_strong = c(FALSE, FALSE),
    identity_moderate = c(FALSE, TRUE),
    independent_support = c(FALSE, TRUE)
  )

  expect_equal(result$direct_signal, c(TRUE, TRUE))
  expect_equal(result$direct_isolated, c(TRUE, FALSE))
  expect_equal(result$direct_supported, c(FALSE, TRUE))
  expect_equal(result$status, c("neg", "pos"))
  expect_equal(result$class, c("neg", "direct"))
})

test_that("detection exposes signal and identity-supported statuses", {
  object_file <- system.file(
    "extdata",
    "hpsc.rds",
    package = "GnRHcell",
    mustWork = TRUE
  )
  object <- readRDS(object_file)
  object <- run_gnrh(object, verbose = FALSE)
  metadata <- object[[]]

  expect_true("gnrh_signal_status" %in% names(metadata))
  expect_equal(
    as.character(metadata$gnrh_signal_status),
    ifelse(metadata$gnrh_direct_signal, "signal", "no_signal")
  )
  expect_true(all((metadata$gnrh_status == "pos") <= metadata$gnrh_direct_signal))
  expect_true(all(metadata$gnrh_direct_isolated <= (metadata$gnrh_status == "neg")))
})
