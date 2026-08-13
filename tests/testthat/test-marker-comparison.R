test_that("GnRH marker scores merge with automatic score detection", {
  a <- data.frame(gene = c("GNRH1", "ECEL1"), avg_log2FC = c(9, 2))
  b <- data.frame(gene = c("gnrh1", "PAX6"), avg_log2FC = c(8, 3))
  matrix <- merge_gnrh_marker_scores(
    list(A = a, B = b), dataset_names = c("A", "B")
  )
  expect_equal(matrix["GNRH1", ], c(A = 9, B = 8))
  expect_true(is.na(matrix["ECEL1", "B"]))
  expect_equal(attr(matrix, "score_columns"), c(A = "avg_log2FC", B = "avg_log2FC"))
})

test_that("GnRH marker score validation is informative", {
  bad <- data.frame(gene = "GNRH1", p_val = 0.01)
  expect_error(merge_gnrh_marker_scores(list(A = bad)), "No supported score column")
})

test_that("named marker file vectors are accepted", {
  directory <- tempfile("gnrh-markers-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  table <- data.frame(gene = "GNRH1", avg_log2FC = 9)
  utils::write.table(table, file.path(directory, "a.tsv"), sep = "\t",
    row.names = FALSE, quote = FALSE)
  matrix <- merge_gnrh_marker_scores(c(A = "a.tsv"), dir = directory)
  expect_equal(unname(matrix["GNRH1", "A"]), 9)
})

test_that("row scaling preserves missing markers and displays singletons", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  matrix <- rbind(
    singleton = c(A = 3, B = NA, C = NA),
    shared = c(A = 1, B = 3, C = NA)
  )
  result <- GnRHcell:::.gnrh_marker_heatmap(matrix, scale_rows = TRUE)
  expect_equal(unname(result$plot_matrix["singleton", "A"]), 2)
  expect_true(all(is.na(result$plot_matrix["singleton", c("B", "C")])))
  expect_true(is.na(result$plot_matrix["shared", "C"]))
  expect_true(all(is.finite(result$plot_matrix["shared", c("A", "B")])))
})

test_that("marker heatmap clustering tolerates missing markers", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  matrix <- rbind(
    GNRH1 = c(A = 9, B = 8, C = 7),
    ECEL1 = c(A = 3, B = NA, C = 2),
    OTX2 = c(A = NA, B = 4, C = 1)
  )
  result <- GnRHcell:::.gnrh_marker_heatmap(
    matrix, scale_rows = TRUE, cluster_rows = TRUE,
    cluster_columns = TRUE
  )
  expect_true(all(is.finite(result$clustering_matrix)))
  expect_s4_class(result$heatmap, "Heatmap")
  expect_silent(ComplexHeatmap::draw(result$heatmap))
})
