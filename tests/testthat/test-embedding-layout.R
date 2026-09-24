test_that("embedding uses the standard free panel aspect", {
  counts <- Matrix::Matrix(
    1,
    nrow = 2,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(c("GNRH1", "ACTB"), paste0("cell", 1:4))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  embedding <- matrix(
    c(0, 0, 1, 3, 2, 6, 3, 9),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  object$gnrh_stage <- c("non-gnrh", "early", "migrating", "mature")

  plot <- gnrh_cellmap(object, group_by = "gnrh_stage")

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$coordinates, "CoordCartesian")
  expect_false(inherits(plot$coordinates, "CoordFixed"))
})

test_that("void embedding style removes every axis component", {
  counts <- Matrix::Matrix(
    1,
    nrow = 2,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(c("GNRH1", "ACTB"), paste0("cell", 1:4))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  embedding <- matrix(
    c(0, 0, 1, 3, 2, 6, 3, 9),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  object$gnrh_status <- c("neg", "neg", "pos", "pos")

  plot <- gnrh_cellmap(
    object,
    group_by = "gnrh_status",
    style = "void"
  )
  built_theme <- ggplot2::ggplot_build(plot)@plot@theme

  expect_s3_class(built_theme$axis.title.x, "element_blank")
  expect_s3_class(built_theme$axis.title.y, "element_blank")
  expect_s3_class(built_theme$axis.text.x, "element_blank")
  expect_s3_class(built_theme$axis.text.y, "element_blank")
  expect_s3_class(built_theme$axis.ticks.x, "element_blank")
  expect_s3_class(built_theme$axis.ticks.y, "element_blank")
})

test_that("gnrh_theme centrally enforces a complete void style", {
  plot <- ggplot2::ggplot(
    data.frame(x = 1:2, y = 1:2),
    ggplot2::aes(.data$x, .data$y)
  ) +
    ggplot2::geom_point() +
    gnrh_theme(
      style = "void",
      ticks = TRUE,
      border = TRUE,
      grid.major = TRUE,
      x.ttl = TRUE,
      y.ttl = TRUE
    )
  built_theme <- ggplot2::ggplot_build(plot)@plot@theme

  expect_s3_class(built_theme$axis.title.x, "element_blank")
  expect_s3_class(built_theme$axis.title.y, "element_blank")
  expect_s3_class(built_theme$axis.text.x, "element_blank")
  expect_s3_class(built_theme$axis.text.y, "element_blank")
  expect_s3_class(built_theme$axis.ticks.x, "element_blank")
  expect_s3_class(built_theme$axis.ticks.y, "element_blank")
  expect_s3_class(built_theme$panel.grid.major, "element_blank")
  expect_s3_class(built_theme$panel.border, "element_blank")
})
