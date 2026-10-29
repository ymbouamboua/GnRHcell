plotfix_fixture <- function() {
  set.seed(110)
  genes <- unique(c("GNRH1","ISL1","GAPDH","ACTB",unlist(gnrh_stage_modules(include_secretory=TRUE))))
  m <- matrix(rpois(length(genes)*24,2),length(genes),dimnames=list(genes,paste0("c",1:24)))
  m["GNRH1",] <- c(10,5,2,1,0,0,6,rep(0,17))
  obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m,sparse=TRUE))
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  coords <- cbind(UMAP_1=seq_len(24),UMAP_2=seq_len(24)*3)
  rownames(coords) <- colnames(obj)
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(coords,key="UMAP_",assay="RNA")
  run_gnrh(obj,min_reference_cells=2,min_counts=1,verbose=FALSE)
}

collect_report_labels <- function(p) {
  labels <- unlist(lapply(p$scales$scales,function(sc) if (is.character(sc$labels)) unname(sc$labels) else NULL))
  if (inherits(p,"patchwork")) labels <- c(labels,unlist(lapply(p$patches$plots,collect_report_labels)))
  labels
}

test_that("report legends retain detection and stage names with their counts", {
  obj <- plotfix_fixture()
  p <- gnrh_report(obj,verbose=FALSE)
  labels <- collect_report_labels(p)
  expect_true(any(grepl("^neg\\n",labels)))
  expect_true(any(grepl("^pos\\n",labels)))
  present <- unique(as.character(obj$gnrh_stage))
  for (stage in setdiff(present,"non-gnrh")) expect_true(any(startsWith(labels,paste0(stage,"\n"))))
  expect_true(any(grepl("^non-GnRH\\n",labels)))
})

test_that("default zoom uses positive cells and preserves square embedding distances", {
  obj <- plotfix_fixture()
  q <- gnrh_qc(obj,qc.columns=character(),zoom.quantile=1)
  positive <- colnames(obj)[obj$gnrh_status=="pos"]
  expect_setequal(q$data$zoom$cells,positive)
  expect_identical(q$parameters$zoom.by,"positive")
  expect_equal(diff(q$data$zoom$xlim),diff(q$data$zoom$ylim))
  expect_equal(diff(q$data$zoom$context_xlim),diff(q$data$zoom$context_ylim))
  xy <- SeuratObject::Embeddings(obj,"umap")[positive,,drop=FALSE]
  expect_true(all(xy[,1]>=q$data$zoom$xlim[1] & xy[,1]<=q$data$zoom$xlim[2]))
  expect_true(all(xy[,2]>=q$data$zoom$ylim[1] & xy[,2]<=q$data$zoom$ylim[2]))
  q <- gnrh_qc(obj,zoom.by="confident",qc.columns=character())
  expect_setequal(q$data$zoom$cells,colnames(obj)[obj$gnrh_confident])
})

test_that("explicit zoom can select background cells and candidate mode is retained", {
  obj <- plotfix_fixture()
  negative <- colnames(obj)[obj$gnrh_status=="neg"]
  q <- gnrh_qc(obj,zoom.cells=negative[1],max_background=1,qc.columns=character())
  expect_identical(q$data$zoom$cells,negative[1])
  expect_true(negative[1] %in% q$data$plotted$cell)
  expect_lt(diff(q$data$zoom$xlim),diff(q$data$zoom$context_xlim))
  q <- gnrh_qc(obj,zoom.by="candidates",qc.columns=character())
  expect_setequal(q$data$zoom$cells,q$data$cells$cell[q$data$cells$category!="negative"])
  expect_error(gnrh_qc(obj,zoom.cells="unknown"),"zoom.cells")
  expect_error(gnrh_qc(obj,zoom.padding=-1),"zoom.padding")
})


test_that("central zoom reports cropped cells without changing the global evidence", {
  obj <- plotfix_fixture()
  obj$gnrh_status <- factor(rep("pos",24),levels=c("neg","pos"))
  obj$gnrh_class <- factor(rep("direct",24),levels=c("neg","supported","direct"))
  coords <- SeuratObject::Embeddings(obj,"umap")
  coords[24,] <- c(1000,1000)
  obj[["umap"]] <- NULL
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(coords,key="UMAP_",assay="RNA")
  q <- gnrh_qc(obj,qc.columns=character(),zoom.quantile=.8,zoom.padding=0)
  expect_true("c24" %in% q$data$zoom$cells_outside)
  expect_true("c24" %in% q$data$plotted$cell)
  expect_equal(length(q$data$zoom$cells_inside)+length(q$data$zoom$cells_outside),24)
  expect_equal(nrow(q$data$cells),24)
  expect_error(gnrh_qc(obj,zoom.quantile=0),"zoom.quantile")
})
