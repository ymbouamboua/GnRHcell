review_object <- function(n=3L) {
  counts <- matrix(rep(c(0,600),n),nrow=2,
                   dimnames=list(c("GNRH1","ISL1"),paste0("cell",seq_len(n))))
  SeuratObject::CreateSeuratObject(Matrix::Matrix(counts,sparse=TRUE))
}

test_that("detailed summaries handle exclusions for multiple cells", {
  object <- detect_gnrh(review_object(),verbose=FALSE)
  object$gnrh_stage <- factor(c("early","non-gnrh","mature"))
  object$gnrh_secretory <- factor(c("limited","non-gnrh","supported"))
  expect_no_error(result <- run_gnrh(object,detect=FALSE,stage=FALSE,
                                    diagnostics=FALSE,summary_level="detailed",
                                    verbose=FALSE))
  expect_identical(result[[]],object[[]])
})

test_that("detection rejects parameters that invalidate count gates", {
  object <- review_object()
  for (x in list(0,-1,NA_real_,Inf,1.5,c(1,2),"2")) {
    expect_error(detect_gnrh(object,min_umi=x,verbose=FALSE),"min_umi")
  }
  invalid <- list(min_counts=-1,k=0,min_reference_cells=1.5,mad_factor=-1,
                  scale_factor=0,max_alternative=2,supported_q=-0.1,
                  candidate_q=1.1,dims=c(1,NA_real_))
  for (nm in names(invalid)) {
    args <- c(list(object=object,verbose=FALSE),setNames(list(invalid[[nm]]),nm))
    expect_error(do.call(detect_gnrh,args),nm)
  }
  expect_error(detect_gnrh(object,supported_q=0.8,candidate_q=0.7,
                           verbose=FALSE),"candidate_q")
  for (umi in c(1L,2L)) {
    result <- detect_gnrh(object,min_umi=umi,verbose=FALSE)
    expect_true(all(as.character(result$gnrh_status)=="neg"))
  }
})

review_diagnostics <- function(score,truth,n_thresholds=200L) {
  object <- detect_gnrh(review_object(length(score)),verbose=FALSE)
  object$test_score <- score
  gnrh_diagnostics(object,truth=truth,predictor="test_score",
                   n_thresholds=n_thresholds,verbose=FALSE)@misc$gnrh
}

test_that("average precision covers perfect rankings and endpoints", {
  one <- review_diagnostics(c(3,2,1),c(TRUE,FALSE,FALSE))
  expect_equal(one$auprc,1)
  expect_equal(one$pr_curve[1,],data.frame(recall=0,precision=1))
  multiple <- review_diagnostics(c(4,3,2,1),c(TRUE,TRUE,FALSE,FALSE))
  expect_equal(multiple$auprc,1)
})

test_that("average precision groups ties and is independent of the grid", {
  # Positive ranks 1 and 3: AP = (1 + 2/3)/2.
  x <- review_diagnostics(c(4,3,2,1),c(TRUE,FALSE,TRUE,FALSE),2L)
  y <- review_diagnostics(c(4,3,2,1),c(TRUE,FALSE,TRUE,FALSE),200L)
  expect_equal(x$auprc,5/6)
  expect_equal(x$auprc,y$auprc)
  tied <- review_diagnostics(c(3,3,1),c(TRUE,FALSE,FALSE))
  permuted <- review_diagnostics(c(3,3,1),c(FALSE,TRUE,FALSE))
  expect_equal(tied$auprc,0.5)
  expect_equal(tied$auprc,permuted$auprc)
  constant <- review_diagnostics(rep(1,4),c(TRUE,FALSE,TRUE,FALSE))
  expect_equal(constant$auprc,0.5)
})
