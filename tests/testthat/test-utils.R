# Core utility tests

test_that("auto_raster works", {
  expect_true(.auto_raster(1e6))
  expect_false(.auto_raster(1e4))
})


test_that("plot_defaults returns list", {
  obj <- matrix(1, ncol = 10)
  out <- .plot_defaults(obj)

  expect_type(out, "list")
  expect_true("raster" %in% names(out))
})

test_that("gnrh_colors returns expected keys", {
  x <- gnrh_colors("status")
  expect_true(all(c("neg", "pos") %in% names(x)))
})
