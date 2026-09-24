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

test_that("cellplot visualization core supports adaptive points and gradients", {
  expect_equal(.auto_raster(3e5), TRUE)
  expect_lt(.auto_pt_size(1e5), .auto_pt_size(1e3))
  expect_length(.cellplot_gradient("bwr", n = 17), 17L)
  expect_equal(
    .cellplot_gradient(colors = c("black", "white"), n = 2),
    c("#000000", "#FFFFFF")
  )
})

test_that("gnrh_theme exposes cellplot typography and layout controls", {
  th <- gnrh_theme(
    axis.size = 9,
    axis.ttl.size = 11,
    ttl.size = 16,
    facet.size = 12,
    leg.key.size = 0.6,
    aspect.ratio = 1.2
  )
  expect_equal(th$axis.text.x$size, 9)
  expect_equal(th$axis.title$size, 11)
  expect_equal(th$plot.title$size, 16)
  expect_equal(th$aspect.ratio, 1.2)
  expect_equal(th$plot.title.position, "panel")

  whole_plot <- gnrh_theme(ttl.scope = "plot")
  expect_equal(whole_plot$plot.title.position, "plot")
})

test_that("gnrh_palette returns expected keys", {
  x <- gnrh_palette("status")
  expect_true(all(c("neg", "pos") %in% names(x)))
})


test_that("figure workflow helpers are public", {
  exports <- getNamespaceExports("GnRHcell")
  expect_true("gnrh_reduction" %in% exports)
  expect_true("gnrh_save" %in% exports)
})

test_that("dot-plot dimensions adapt to identities, genes, and labels", {
  two <- .gnrh_dot_dims(c("Fetal nose", "Fetal hypothalamus"), n_genes = 20)
  four <- .gnrh_dot_dims(
    c("Fetal nose", "Fetal hypothalamus", "Adult ME", "Adult hypothalamus"),
    n_genes = 20
  )
  more_genes <- .gnrh_dot_dims(c("Fetal nose", "Fetal hypothalamus"), n_genes = 40)

  expect_lt(unname(two["width"]), unname(four["width"]))
  expect_lt(unname(two["height"]), unname(more_genes["height"]))
  expect_true(unname(two["width"]) < 10)
})

test_that("gnrh_save saves vector and raster formats", {
  plot <- ggplot2::ggplot(
    data.frame(x = 1:2, y = 1:2),
    ggplot2::aes(.data$x, .data$y)
  ) + ggplot2::geom_point()
  pdf_file <- tempfile(fileext = ".pdf")
  png_file <- tempfile(fileext = ".png")

  expect_no_error(
    gnrh_save(plot, pdf_file, width = 4, height = 3)
  )
  expect_no_error(
    gnrh_save(plot, png_file, width = 4, height = 3, dpi = 100)
  )
  expect_true(file.exists(pdf_file))
  expect_true(file.exists(png_file))
})

test_that("gnrh_save saves tables and serialized R objects", {
  table <- data.frame(gene = c("GNRH1", "ISL1"), score = c(2.1, 1.4))
  csv_file <- tempfile(fileext = ".csv")
  tsv_base <- tempfile()
  rds_file <- tempfile(fileext = ".rds")
  rda_file <- tempfile(fileext = ".RData")

  expect_equal(gnrh_save(table, csv_file), csv_file)
  expect_equal(gnrh_save(table, tsv_base), paste0(tsv_base, ".tsv"))
  expect_equal(gnrh_save(table, rds_file), rds_file)
  expect_equal(gnrh_save(table, rda_file, object_name = "markers"), rda_file)
  expect_true(all(file.exists(c(csv_file, paste0(tsv_base, ".tsv"), rds_file, rda_file))))
  expect_equal(readRDS(rds_file), table)

  loaded <- new.env(parent = emptyenv())
  expect_equal(load(rda_file, envir = loaded), "markers")
  expect_equal(loaded$markers, table)
})

test_that("gnrh_save infers defaults and validates mismatched formats", {
  object <- list(stage = c("early", "mature"))
  base <- tempfile()
  expect_equal(gnrh_save(object, base), paste0(base, ".rds"))
  expect_equal(readRDS(paste0(base, ".rds")), object)
  expect_error(
    gnrh_save(object, tempfile(fileext = ".rds"), format = "csv"),
    "does not match"
  )
})
