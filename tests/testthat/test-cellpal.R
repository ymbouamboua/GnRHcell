test_that("cellpal returns correct number of colors", {
  cols <- cellpal(n = 5, preset = "base")

  expect_type(cols, "character")
  expect_equal(length(cols), 5)
})
