test_that("marker-program dot plot selects and groups ranked genes", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(
    c(5, 4, 0, 0, 3, 2, 1, 0, 0, 1, 4, 5),
    nrow = 3,
    dimnames = list(c("GNRH1", "FEZF1", "ROBO3"), paste0("cell", 1:4))
  ), sparse = TRUE)
  object <- Seurat::CreateSeuratObject(counts)
  object$dataset <- rep(c("nose", "hypothalamus"), each = 2)
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  marker_programs <- list(
    high_confidence = data.frame(
      gene = c("FEZF1", "FEZF1", "ROBO3", "ROBO3", "MISSING"),
      program = c("early", "early", "migrating", "migrating", "early"),
      integrated_score = c(2, 1, 3, 2, 10),
      n_datasets = c(2, 2, 2, 2, 2),
      stringsAsFactors = FALSE
    )
  )

  plot <- gnrh_program(
    object,
    programs = marker_programs,
    group.by = "dataset",
    top_n = 1
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(
    attr(plot, "gnrh_features"),
    list(Early = "FEZF1", Migrating = "ROBO3")
  )

  reversed <- gnrh_program(
    object,
    programs = marker_programs,
    group.by = "dataset",
    top_n = 1,
    program_order = c("migrating", "early")
  )
  expect_equal(names(attr(reversed, "gnrh_features")), c("Migrating", "Early"))
})

test_that("named dot-plot groups tolerate overlapping genes", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(5, 4, 3, 2, 1, 1),
      nrow = 2,
      dimnames = list(c("FEZF1", "SIX3"), paste0("cell", 1:3))
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts)
  object$dataset <- "test"
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  expect_no_error(
    gnrh_celldot(
      object,
      features = list(
        `GnRH identity` = c("FEZF1", "SIX3"),
        Early = c("FEZF1", "SIX3")
      ),
      group.by = "dataset",
      scale = FALSE
    )
  )
})

test_that("marker-program dot plot validates its inputs", {
  skip_if_not_installed("Seurat")

  object <- suppressWarnings(Seurat::CreateSeuratObject(Matrix::Matrix(
      matrix(
        c(1, 1, 0, 1),
        nrow = 2,
        dimnames = list(c("GNRH1", "FEZF1"), c("a", "b"))
      ),
      sparse = TRUE
    )))

  expect_error(
    gnrh_program(object, list(), group.by = "orig.ident"),
    "find_gnrh_programs"
  )
  expect_error(
    gnrh_program(
      object,
      list(high_confidence = data.frame(gene = "MISSING", program = "early")),
      group.by = "orig.ident"
    ),
    "No marker-program genes"
  )
})

test_that("collection input combines user-selected datasets", {
  skip_if_not_installed("Seurat")

  make_object <- function(prefix) {
    counts <- Matrix::Matrix(
      matrix(
        c(4, 3, 1, 0, 2, 1),
        nrow = 2,
        dimnames = list(c("FEZF1", "ROBO3"), paste0(prefix, 1:3))
      ),
      sparse = TRUE
    )
    object <- Seurat::CreateSeuratObject(counts)
    object$gnrh_status <- c("pos", "pos", "neg")
    Seurat::NormalizeData(object, verbose = FALSE)
  }
  marker_table <- data.frame(
    gene = c("FEZF1", "ROBO3"),
    program = c("early", "migrating"),
    integrated_score = c(2, 1),
    n_datasets = 2,
    stringsAsFactors = FALSE
  )
  collection <- list(
    datasets = data.frame(
      id = c("nose", "hypo"),
      label = c("Fetal nose", "Fetal hypothalamus")
    ),
    results = list(
      nose = list(object = make_object("n")),
      hypo = list(object = make_object("h"))
    ),
    comparisons = list(programs = list(high_confidence = marker_table)),
    output_dir = tempdir()
  )

  combined <- combine_gnrh_datasets(collection, c("nose", "Fetal hypothalamus"))
  expect_s4_class(combined, "Seurat")
  expect_equal(sort(unique(combined$dataset)), c("Fetal hypothalamus", "Fetal nose"))
  expect_equal(sort(as.integer(table(combined$dataset))), c(2L, 2L))

  plot <- gnrh_program(
    collection,
    selected_ids = c("nose", "hypo"),
    group.by = "dataset",
    top_n = 1
  )
  expect_s3_class(plot, "ggplot")
})
