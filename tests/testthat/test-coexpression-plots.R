test_that("coexpression plot titles are bold and centered", {
  markers <- data.frame(
    gene = c("GNRH1", "ECEL1", "DLX5"),
    coexpr = c(0.9, 0.5, 0.4),
    score = c(1, 0.5, 0.4),
    p_val_adj = c(1e-20, 1e-10, 1e-5)
  )

  coexpr <- plot_gnrh_coexpr(markers, coexp_cutoff = 0.3)
  network <- plot_network(markers, top_n = 3, threshold = 0.3)

  expect_equal(coexpr$theme$plot.title$face, "bold")
  expect_equal(coexpr$theme$plot.title$hjust, 0.5)
  expect_equal(network$theme$plot.title$face, "bold")
  expect_equal(network$theme$plot.title$hjust, 0.5)
})
