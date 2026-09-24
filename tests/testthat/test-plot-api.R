test_that("concise plotting API is exported", {
  expected <- c(
    "gnrh_cellmap", "gnrh_cellfeat", "gnrh_celldot",
    "gnrh_pheno", "gnrh_cellmark", "gnrh_program_plot",
    "gnrh_program", "gnrh_cellcount", "gnrh_specificity",
    "gnrh_celldistribution", "gnrh_hits", "gnrh_network",
    "gnrh_coexpr", "gnrh_codetect", "gnrh_conserved",
    "gnrh_specific", "gnrh_compare", "gnrh_signature",
    "gnrh_runtime", "gnrh_upset", "gnrh_save",
    "gnrh_theme", "gnrh_palette", "gnrh_reduction"
  )
  exports <- getNamespaceExports("GnRHcell")
  expect_true(all(expected %in% exports))
  expect_true("find_gnrh_programs" %in% exports)
  expect_false(any(grepl("^plot_gnrh_|^plot_class_counts$|^plot_network$", exports)))
  expect_false(any(c(
    "gnrh_celltheme", "gnrh_colors", "gnrh_gene_upset",
    "save_gnrh_plot", "resolve_reduction",
    "gnrh_cellconserved", "gnrh_cellspecific",
    "gnrh_cellupset", "gnrh_cellspecificity",
    "gnrh_cellprogramdot", "gnrh_cellnetwork",
    "gnrh_cellcodetect", "gnrh_cellcoexpr",
    "gnrh_geneconserved", "gnrh_genespecific",
    "gnrh_cellcompare", "gnrh_celldetect",
    "gnrh_cellruntime", "gnrh_cellprogram",
    "gnrh_detect", "gnrh_cellpalette", "gnrh_cellsave",
    "gnrh_cellreduction", "gnrh_cellhits", "gnrh_cellsignature",
    "gnrh_geneupset", "gnrh_marker_program", "gnrh_marker_programs",
    "validate_input", "extract_gnrh_run_info", "merge_gnrh_marker_scores"
  ) %in% exports))
})

test_that("concise core plotting names forward their arguments", {
  expect_error(gnrh_cellmap(object = matrix(1)), "Seurat")
  expect_error(gnrh_cellfeat(object = matrix(1), features = "GNRH1"), "Seurat")
  expect_error(gnrh_celldot(object = matrix(1), features = "ISL1"), "Seurat")
})
