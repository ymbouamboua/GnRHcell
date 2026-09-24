test_that("external signature validation accepts dataset labels and returns plot data", {
  make_object <- function(prefix, shift = 0) {
    counts <- Matrix::Matrix(
      c(
        4 + shift, 2, 0,
        1, 3 + shift, 0,
        2, 1, 1
      ),
      nrow = 3,
      byrow = TRUE,
      sparse = TRUE,
      dimnames = list(
        c("GNRH1", "ISL1", "ACTB"),
        paste0(prefix, "_", 1:3)
      )
    )
    object <- SeuratObject::CreateSeuratObject(counts = counts)
    object <- suppressMessages(Seurat::NormalizeData(object, verbose = FALSE))
    object$gnrh_status <- c("pos", "pos", "neg")
    object
  }

  collection <- list(
    datasets = data.frame(
      id = c("nose", "hypo"),
      label = c("Human fetal nose", "Human fetal hypothalamus"),
      species = c("Human", "Human"),
      stringsAsFactors = FALSE
    ),
    results = list(
      nose = list(object = make_object("nose")),
      hypo = list(object = make_object("hypo", 2))
    ),
    output_dir = tempdir()
  )

  plot <- gnrh_signature(
    collection = collection,
    genes = c("GNRH1", "ISL1", "MISSING"),
    selected_ids = c("Human fetal nose", "Human fetal hypothalamus"),
    show_title = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$title)
  data <- attr(plot, "gnrh_signature_data")
  expect_equal(nrow(data), 4L)
  expect_setequal(as.character(data$gene), c("GNRH1", "ISL1"))
  expect_true(all(is.finite(data$pct_detected)))
  expect_true(all(is.finite(data$scaled_mean_expr)))
  expect_equal(
    attr(plot, "gnrh_missing_genes")$nose,
    "MISSING"
  )
  expect_equal(
    names(attr(plot, "gnrh_dimensions")),
    c("width", "height")
  )
  expect_identical(attr(plot, "gnrh_renderer"), "gnrh_celldot")
  expect_true(all(c("avg.exp", "pct.exp", "features.plot", "id") %in% names(plot$data)))
  expect_equal(plot$theme$legend.position, "right")
  expect_equal(plot$theme$axis.text.x$angle, 60)
})

test_that("external signature validation selects human datasets by default", {
  counts <- Matrix::Matrix(
    c(4, 2, 1, 1),
    nrow = 2,
    sparse = TRUE,
    dimnames = list(c("GNRH1", "ACTB"), c("cell1", "cell2"))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object <- suppressMessages(Seurat::NormalizeData(object, verbose = FALSE))
  object$gnrh_status <- c("pos", "neg")
  collection <- list(
    datasets = data.frame(
      id = "human",
      label = "Human model",
      species = "Human",
      stringsAsFactors = FALSE
    ),
    results = list(human = list(object = object)),
    output_dir = tempdir()
  )

  plot <- gnrh_signature(collection, genes = "GNRH1")
  expect_equal(attr(plot, "gnrh_datasets")$id, "human")
})
