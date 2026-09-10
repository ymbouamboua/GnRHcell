# ========================================================================= #
# Input validation
# ========================================================================= #
#' Validate input Seurat object for GnRH analysis
#' @param object A Seurat object.
#' @param assay Assay to validate; defaults to `DefaultAssay(object)`.
#' @param required_layers Required layers: counts, data, or scale.data.
#' @param auto_normalize Generate missing `data` with `NormalizeData()`.
#' @param min_cells Minimum number of cells.
#' @param min_features Minimum number of features.
#' @param verbose Print progress messages.
#' @return A validated Seurat object, optionally normalized.
#' @export
validate_input <- function(object,assay=NULL,required_layers=c("counts","data"),auto_normalize=FALSE,min_cells=1L,min_features=1L,verbose=TRUE) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  log <- .msg(verbose)
  assay <- assay %||% Seurat::DefaultAssay(object)
  if (length(assay)!=1L || is.na(assay) || !nzchar(assay)) stop("`assay` must be a single valid assay name.",call.=FALSE)
  assays <- SeuratObject::Assays(object)
  if (!assay %in% assays) stop("Assay `",assay,"` not found. Available assays: ",paste(assays,collapse=", "),call.=FALSE)
  required_layers <- unique(as.character(required_layers))
  required_layers <- required_layers[!is.na(required_layers) & nzchar(required_layers)]
  supported <- c("counts","data","scale.data")
  invalid <- setdiff(required_layers,supported)
  if (length(invalid)) stop("Unsupported required layer(s): ",paste(invalid,collapse=", "),". Supported values: ",paste(supported,collapse=", "),".",call.=FALSE)
  .int <- function(x,nm) {
    if (length(x)!=1L || is.na(x) || !is.finite(x) || x<1 || x!=floor(x)) stop("`",nm,"` must be an integer >= 1.",call.=FALSE)
    as.integer(x)
  }
  min_cells <- .int(min_cells,"min_cells")
  min_features <- .int(min_features,"min_features")
  if (ncol(object)<min_cells) stop("Object contains ",ncol(object)," cells; at least ",min_cells," required.",call.=FALSE)
  if (nrow(object[[assay]])<min_features) stop("Assay `",assay,"` contains ",nrow(object[[assay]])," features; at least ",min_features," required.",call.=FALSE)
  log("Using assay:",assay)
  get_layers <- function(x) {
    if (inherits(x,"Assay5")) return(SeuratObject::Layers(x))
    slots <- methods::slotNames(x)
    out <- character()
    for (nm in intersect(c("counts","data","scale.data"),slots)) {
      z <- methods::slot(x,nm)
      if (nrow(z)>0L && ncol(z)>0L) out <- c(out,nm)
    }
    out
  }
  obj_assay <- object[[assay]]
  available <- get_layers(obj_assay)
  if (inherits(obj_assay,"Assay5")) {
    counts_layers <- grep("^counts(\\.|$)",available,value=TRUE)
    data_layers <- grep("^data(\\.|$)",available,value=TRUE)
    if (length(counts_layers)>1L || length(data_layers)>1L) warning("Assay `",assay,"` contains multiple expression layers. Consider `Seurat::JoinLayers()` before GnRHcell analysis.",call.=FALSE)
  }
  if ("data" %in% required_layers && !"data" %in% available) {
    if (!isTRUE(auto_normalize)) stop("Normalized `data` layer is required but missing from assay `",assay,"`. Run `Seurat::NormalizeData()` or use `auto_normalize = TRUE`.",call.=FALSE)
    if (!"counts" %in% available) stop("Cannot normalize assay `",assay,"` because the `counts` layer is unavailable.",call.=FALSE)
    log("Normalized data missing; running NormalizeData().",type="step")
    object <- Seurat::NormalizeData(object,assay=assay,verbose=FALSE)
    available <- get_layers(object[[assay]])
  }
  missing <- setdiff(required_layers,available)
  if (length(missing)) stop("Assay `",assay,"` is missing required layer(s): ",paste(missing,collapse=", "),". Available layers: ",paste(available,collapse=", "),".",call.=FALSE)
  for (layer in required_layers) {
    mat <- tryCatch(Seurat::GetAssayData(object,assay=assay,layer=layer),error=function(e) NULL)
    if (is.null(mat)) stop("Unable to retrieve layer `",layer,"` from assay `",assay,"`.",call.=FALSE)
    if (!nrow(mat) || !ncol(mat)) stop("Layer `",layer,"` in assay `",assay,"` is empty.",call.=FALSE)
    if (ncol(mat)!=ncol(object)) stop("Layer `",layer,"` contains ",ncol(mat)," cells but the Seurat object contains ",ncol(object),".",call.=FALSE)
  }
  invisible(object)
}
