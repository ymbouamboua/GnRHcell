test_that("normalized module scores are not library-normalized twice", {
  expr <- Matrix::Matrix(
    matrix(
      c(
        1, 2,
        3, 4,
        0, 5
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("A", "B", "C"), c("cell1", "cell2"))
    ),
    sparse = TRUE
  )

  result <- .score_modules(
    expr,
    modules = list(program = c("A", "B")),
    expression_weight = 1,
    detection_weight = 0,
    input_type = "normalized"
  )

  expect_equal(
    as.numeric(result$score$program),
    c(2, 3)
  )
  expect_identical(result$input_type, "normalized")
})


test_that("count module scores retain library-size normalization", {
  expr <- Matrix::Matrix(
    matrix(
      c(10, 20, 30, 40),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("A", "B"), c("cell1", "cell2"))
    ),
    sparse = TRUE
  )

  result <- .score_modules(
    expr,
    modules = list(program = c("A", "B")),
    expression_weight = 1,
    detection_weight = 0,
    input_type = "counts"
  )

  expected <- log1p(c(20, 30) / c(40, 60) * 10000)
  expect_equal(as.numeric(result$score$program), expected)
})


test_that("orthogonal developmental states preserve uncertainty", {
  state <- .derive_developmental_state(
    raw_stage = c("early", "migrating", "mature", "mature", "migrating", NA),
    margin = c(0.2, 0.2, 0.2, 0.2, 0.05, NA),
    positive = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    migration_filtered = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    secretory_supported = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    stage_margin = 0.1
  )

  expect_identical(
    state,
    c(
      "early",
      "migrating",
      "post-migratory",
      "mature",
      "transitional",
      "non-gnrh"
    )
  )
})


test_that("migration filtering cannot yield a confident raw stage", {
  scores <- data.frame(
    early = 0.3,
    migrating = 0.8,
    mature = 0.2
  )

  assignment <- assign_stage(
    scores,
    migration_core_hits = 0,
    min_migration_hits = 2
  )

  expect_identical(assignment$raw_stage, "migrating")
  expect_true(assignment$migration_filtered)
  expect_identical(assignment$stage, "early")
})
