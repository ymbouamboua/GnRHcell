# ----------------------------------------------------------------------- #
# Module scoring
# ----------------------------------------------------------------------- #
#' Compute module scores and detection metrics
#' @keywords internal
#' @noRd
.score_modules <- function(expr,modules,scale_factor=10000,expression_weight=0.5,detection_weight=0.5,input_type=c("counts","normalized")) {
  input_type <- match.arg(input_type)
  if (!length(modules)) stop("`modules` must contain at least one gene set.",call.=FALSE)
  if (length(expression_weight)!=1L || !is.finite(expression_weight) || expression_weight<0) stop("`expression_weight` must be a non-negative number.",call.=FALSE)
  if (length(detection_weight)!=1L || !is.finite(detection_weight) || detection_weight<0) stop("`detection_weight` must be a non-negative number.",call.=FALSE)
  if (expression_weight+detection_weight<=0) stop("At least one module-score weight must be positive.",call.=FALSE)
  w <- expression_weight+detection_weight
  expression_weight <- expression_weight/w
  detection_weight <- detection_weight/w
  lib <- if (input_type=="counts") pmax(Matrix::colSums(expr),1) else NULL
  modules <- lapply(modules,intersect,y=rownames(expr))
  score_fun <- function(g) {
    if (!length(g)) return(rep(0,ncol(expr)))
    x <- Matrix::colMeans(expr[g,,drop=FALSE])
    if (input_type=="counts") log1p(x/lib*scale_factor) else as.numeric(x)
  }
  hit_fun <- function(g) if (!length(g)) integer(ncol(expr)) else as.integer(Matrix::colSums(expr[g,,drop=FALSE]>0))
  frac_fun <- function(g) if (!length(g)) rep(0,ncol(expr)) else Matrix::colSums(expr[g,,drop=FALSE]>0)/length(g)
  to_df <- function(x) {
    x <- as.data.frame(x,check.names=FALSE)
    rownames(x) <- colnames(expr)
    x
  }
  score <- to_df(lapply(modules,score_fun))
  hits <- to_df(lapply(modules,hit_fun))
  fraction <- to_df(lapply(modules,frac_fun))
  integrated <- as.data.frame(expression_weight*as.matrix(score)+detection_weight*as.matrix(fraction),check.names=FALSE)
  rownames(integrated) <- colnames(expr)
  list(score=score,hits=hits,fraction=fraction,integrated=integrated,genes=modules,weights=c(expression=expression_weight,detection=detection_weight),input_type=input_type)
}
# ----------------------------------------------------------------------- #
# Ambient RNA
# ----------------------------------------------------------------------- #
#' Estimate ambient RNA background signal
#' @keywords internal
#' @noRd
.ambient <- function(expr,gene,lib,max_library=100) {
  low <- is.finite(lib) & lib>0 & lib<max_library
  if (!any(low)) return(0)
  x <- as.numeric(Matrix::rowMeans(expr[,low,drop=FALSE])[gene])
  if (!length(x) || !is.finite(x)) 0 else x
}
# ----------------------------------------------------------------------- #
# KNN support
# ----------------------------------------------------------------------- #
#' Compute neighborhood support signal
#' @keywords internal
#' @noRd
.knn_signal <- function(object,signal,reduction="pca",dims=1:20,k=20) {
  n <- length(signal)
  if (!reduction %in% names(object@reductions)) return(rep(0,n))
  emb <- Seurat::Embeddings(object,reduction=reduction)
  dims <- dims[dims>=1L & dims<=ncol(emb)]
  if (!length(dims)) return(rep(0,n))
  emb <- emb[,dims,drop=FALSE]
  if (nrow(emb)!=n) stop("Signal length does not match reduction.",call.=FALSE)
  if (n<=1L) return(rep(0,n))
  k <- min(as.integer(k),n-1L)
  if (k<1L) return(rep(0,n))
  signal[!is.finite(signal)] <- 0
  nn <- FNN::get.knn(emb,k=k)
  rowMeans(matrix(signal[nn$nn.index],nrow=n,ncol=k))
}
