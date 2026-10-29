qc_fixture <- function() {
  set.seed(110)
  genes <- unique(c("GNRH1","ISL1","GAPDH","ACTB",unlist(gnrh_stage_modules(include_secretory=TRUE))))
  m <- matrix(rpois(length(genes)*24,2),length(genes),dimnames=list(genes,paste0("c",1:24)))
  m["GNRH1",] <- c(10,5,2,1,0,0,6,rep(0,17))
  obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m,sparse=TRUE))
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  obj$sample <- rep(c("A","B",NA),8)
  obj$percent_mito <- seq(0,10,length.out=24)
  obj$scDblFinder_score <- seq(0,1,length.out=24)
  coords <- matrix(rnorm(48),24,dimnames=list(colnames(obj),c("UMAP_1","UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(coords,key="UMAP_",assay="RNA")
  run_gnrh(obj,min_reference_cells=2,min_counts=1,verbose=FALSE)
}

test_that("QC reads archived fields, retains all candidates and counts every sample", {
  obj <- qc_fixture()
  original <- obj
  q <- gnrh_qc(obj,sample.by="sample",max_background=2)
  expect_s3_class(q$plot,"patchwork")
  expect_equal(sum(q$data$samples$n_cells),ncol(obj))
  expect_equal(sum(q$data$samples$n_positive),sum(obj$gnrh_status=="pos"))
  expect_true("<missing sample>" %in% q$data$samples$sample)
  expect_equal(q$data$samples$percent_positive,100*q$data$samples$n_positive/q$data$samples$n_cells)
  expect_true(all(c("nCount_RNA","nFeature_RNA","percent_mito","scDblFinder_score") %in% q$parameters$qc.columns))
  expect_true(all(q$data$cells$cell[q$data$cells$category!="negative"] %in% q$data$plotted$cell))
  expect_lte(sum(q$data$plotted$category=="negative"),2)
  expect_identical(obj,original)
  f <- tempfile(fileext=".png")
  expect_no_error(ggplot2::ggsave(f,q$plot,width=16,height=16,dpi=90))
  unlink(f)
})

test_that("QC uses current cell order, preserves RNG and handles optional panels", {
  obj <- qc_fixture()
  set.seed(912);before <- .Random.seed
  q <- gnrh_qc(obj,max_background=1)
  expect_identical(.Random.seed,before)
  small <- obj[,c("c7","c2","c10")]
  q <- gnrh_qc(small,reduction=NULL,qc.columns=character())
  expect_identical(q$data$cells$cell,colnames(small))
  expect_equal(sum(q$data$samples$n_cells),3)
  expect_null(q$plots$technical)
  expect_null(q$plots$embedding)
  q <- gnrh_qc(obj,reduction="missing",qc.columns=character())
  expect_true(any(grepl("unavailable",q$notes)))
  expect_error(gnrh_qc(obj,sample.by="missing"),"sample.by")
  expect_error(gnrh_qc(obj,qc.columns="sample"),"numeric")
  expect_error(gnrh_qc(obj,dims=c(1,1)),"dims")
  expect_error(gnrh_qc(obj,max_background=0),"max_background")
  expect_error(gnrh_qc(SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(c(1,2,3,4),2,2),sparse=TRUE))),"Missing GnRH")
})

test_that("stability reruns detection and returns exact positive frequencies", {
  obj <- qc_fixture();original <- obj
  grid <- data.frame(min_umi=c(1,2,4),supported_q=c(.5,.6,.8))
  s <- gnrh_stability(obj,grid=grid,sample.by="sample",verbose=FALSE)
  expect_s3_class(s$plot,"patchwork")
  expect_equal(nrow(s$runs),3)
  expect_identical(s$cells$cell,colnames(obj))
  manual <- matrix(FALSE,ncol(obj),3)
  for (i in 1:3) {
    fit <- detect_gnrh(obj,min_umi=grid$min_umi[i],supported_q=grid$supported_q[i],
      min_reference_cells=2,min_counts=1,verbose=FALSE)
    manual[,i] <- fit$gnrh_status=="pos"
    expect_equal(s$runs$n_positive[i],sum(manual[,i]))
  }
  expect_equal(s$cells$positive_frequency,rowMeans(manual))
  expect_equal(sum(s$cells$positive_frequency)*3,sum(s$runs$n_positive))
  expect_true(all(s$runs$jaccard>=0 & s$runs$jaccard<=1))
  expect_equal(s$parameters$baseline_disagreements,0)
  expect_equal(nrow(s$samples),9)
  expect_true(anyNA(s$samples$sample))
  expect_identical(obj,original)
  f <- tempfile(fileext=".png")
  expect_no_error(ggplot2::ggsave(f,s$plot,width=11,height=12,dpi=90))
  unlink(f)
})

test_that("stability handles subsets, duplicate grid rows and zero positives", {
  obj <- qc_fixture()
  small <- obj[,c("c7","c1","c4","c8")]
  s <- gnrh_stability(small,grid=data.frame(k=c(2,2)),verbose=FALSE)
  expect_identical(s$cells$cell,colnames(small))
  expect_equal(nrow(s$runs),1)
  s <- gnrh_stability(obj,grid=data.frame(min_umi=2),detect_args=list(min_counts=1e9),verbose=FALSE)
  expect_equal(s$runs$n_positive,0)
  expect_equal(s$runs$jaccard,1)
  expect_true(all(s$cells$positive_frequency==0))
  expect_null(s$plots$frequency)
  expect_error(gnrh_stability(obj,grid=data.frame(supported_q=1),verbose=FALSE),"candidate_q")
  expect_error(gnrh_stability(obj,grid=data.frame(k=0),verbose=FALSE),"positive integers")
  expect_error(gnrh_stability(obj,grid=data.frame(unknown=1),verbose=FALSE),"grid")
  expect_error(gnrh_stability(obj,detect_args=list(unknown=1),verbose=FALSE),"detect_args")
  expect_error(gnrh_stability(obj,grid=data.frame(min_umi=2),detect_args=list(assay="missing"),verbose=FALSE),"Baseline failed")
})
