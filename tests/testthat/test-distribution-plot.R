test_that("distribution plots adapt to the number of samples", {
  counts <- Matrix::Matrix(1, nrow = 2, ncol = 12, sparse = TRUE,
    dimnames = list(c("GNRH1", "ACTB"), paste0("cell", 1:12)))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- rep(c("D22", "D24", "D26"), each = 4)
  object$gnrh_stage <- rep(c("non-gnrh", "identity", "migrating", "mature"), 3)

  plot <- plot_gnrh_distribution(
    object, group.by = "gnrh_stage", split.by = "sample",
    proportion = TRUE, label = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(attr(plot, "n_samples"), 3L)
  expect_equal(unname(attr(plot, "recommended_size")["height"]), 5.2)
  expect_equal(plot$theme$legend.position, "bottom")
  y_scale <- plot$scales$get_scales("y")
  expect_equal(y_scale$breaks, seq(0, 1, 0.25))
  expect_equal(y_scale$labels(c(0, 0.5, 1)), c("0.00", "0.50", "100"))
})

test_that("bar gap is independent of bar width", {
  counts <- Matrix::Matrix(1, nrow = 2, ncol = 12, sparse = TRUE,
    dimnames = list(c("GNRH1", "ACTB"), paste0("cell", 1:12)))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- rep(c("D22", "D24", "D26"), each = 4)
  object$gnrh_stage <- rep(c("non-gnrh", "identity", "migrating", "mature"), 3)

  plot <- plot_gnrh_distribution(
    object, group.by = "gnrh_stage", split.by = "sample",
    bar.width = 0.3, bar.gap = 0.2, adaptive = FALSE, label = FALSE
  )
  positions <- sort(unique(plot$data$.split_position))
  expect_equal(diff(positions), c(0.5, 0.5))
  expect_equal(plot$layers[[1]]$aes_params$width, 0.3)
  expect_equal(plot$scales$get_scales("x")$limits, c(0.75, 2.25))
})
