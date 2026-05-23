#' Validate input Seurat object for GnRH analysis
#'
#' Validates that the input Seurat object contains the required
#' expression layers for GnRH analysis.
#'
#' This function checks compatibility with both Seurat v5
#' (\code{Assay5}) and earlier assay formats, ensuring that
#' raw counts and normalized expression data are available.
#'
#' If normalized data are missing and \code{auto_normalize = TRUE},
#' normalization is performed automatically using
#' \code{Seurat::NormalizeData()}.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param assay Assay name to validate. If \code{NULL}, the default
#' assay of the Seurat object is used.
#' @param auto_normalize Logical; automatically normalize the assay
#' if normalized data are missing. Default is \code{FALSE}.
#' @param verbose Logical; print progress messages.
#' Default is \code{TRUE}.
#'
#' @return A validated Seurat object, optionally normalized if
#' requested.
#'
#' @details
#' Validation checks:
#' \itemize{
#'   \item presence of raw counts (\code{counts})
#'   \item presence of normalized expression data (\code{data})
#'   \item compatibility with Seurat v5 \code{Assay5} objects
#'   \item compatibility with legacy Seurat assay objects
#' }
#'
#' If counts are missing, the function stops with an error.
#'
#' If normalized data are missing:
#' \itemize{
#'   \item with \code{auto_normalize = FALSE}: error
#'   \item with \code{auto_normalize = TRUE}: normalization is performed
#' }
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link[Seurat:NormalizeData]{Seurat::NormalizeData}}
#'
#' @export
validate_input <- function(
    object,
    assay = NULL,
    auto_normalize = FALSE,
    verbose = TRUE
) {

  stopifnot(inherits(object, "Seurat"))

  log <- .msg(verbose)

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  log("Using assay:", assay)

  obj_assay <- object[[assay]]

  # Assay5 (Seurat v5)
  if (inherits(obj_assay, "Assay5")) {

    layer_names <- SeuratObject::Layers(obj_assay)

    has_counts <- "counts" %in% layer_names
    has_data <- "data" %in% layer_names

  } else {

    # Old Assay
    slots <- slotNames(obj_assay)

    has_counts <- "counts" %in% slots &&
      nrow(obj_assay@counts) > 0

    has_data <- "data" %in% slots &&
      nrow(obj_assay@data) > 0
  }

  if (!has_counts) {
    stop("counts layer/slot missing")
  }

  if (!has_data) {

    if (!auto_normalize) {
      stop("normalized data missing")
    }

    log("Running NormalizeData()")

    object <- Seurat::NormalizeData(
      object,
      assay = assay,
      verbose = FALSE
    )
  }

  return(object)
}
