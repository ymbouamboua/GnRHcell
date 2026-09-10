#' Validate input Seurat object for GnRH analysis
#'
#' Validates that a Seurat object contains the assay and expression layers
#' required by GnRHcell.
#'
#' The function supports Seurat v5 \code{Assay5} objects and legacy Seurat
#' assays. Required expression layers can be specified explicitly, allowing
#' detection workflows to require raw counts only while staging or marker
#' workflows can additionally require normalized data.
#'
#' When normalized expression is required but unavailable,
#' \code{auto_normalize = TRUE} runs \code{Seurat::NormalizeData()}.
#'
#' @param object A Seurat object.
#' @param assay Assay to validate. If \code{NULL}, the default assay is used.
#' @param required_layers Character vector containing required expression
#'   layers. Default is \code{c("counts", "data")}.
#' @param auto_normalize Logical. Automatically generate the \code{"data"}
#'   layer with \code{Seurat::NormalizeData()} when required and missing.
#'   Default is \code{FALSE}.
#' @param min_cells Minimum number of cells required. Default is \code{1}.
#' @param min_features Minimum number of features required. Default is
#'   \code{1}.
#' @param verbose Logical. Print progress messages.
#'
#' @return A validated Seurat object, optionally normalized.
#'
#' @details
#' For \code{\link{detect_gnrh}}, validation should normally use
#' \code{required_layers = "counts"}.
#'
#' Functions based on normalized expression, such as developmental staging or
#' marker discovery, should require \code{c("counts", "data")} or
#' \code{"data"} as appropriate.
#'
#' For Seurat v5 assays containing multiple count or data layers,
#' \code{Seurat::JoinLayers()} may be required before analysis.
#'
#' @seealso
#' \code{\link{detect_gnrh}},
#' \code{\link{stage_gnrh}},
#' \code{\link[Seurat:NormalizeData]{Seurat::NormalizeData}}
#'
#' @export
validate_input <- function(
    object,
    assay = NULL,
    required_layers = c("counts", "data"),
    auto_normalize = FALSE,
    min_cells = 1L,
    min_features = 1L,
    verbose = TRUE
) {
  if (!inherits(object, "Seurat"))
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )

  log <- .msg(verbose)

  assay <- assay %||%
    Seurat::DefaultAssay(object)

  if (
    length(assay) != 1L ||
    is.na(assay) ||
    !nzchar(assay)
  ) {
    stop(
      "`assay` must be a single valid assay name.",
      call. = FALSE
    )
  }

  if (!assay %in% SeuratObject::Assays(object))
    stop(
      "Assay `",
      assay,
      "` not found. Available assays: ",
      paste(
        SeuratObject::Assays(object),
        collapse = ", "
      ),
      call. = FALSE
    )

  required_layers <- unique(
    as.character(
      required_layers
    )
  )

  required_layers <- required_layers[
    !is.na(required_layers) &
      nzchar(required_layers)
  ]

  supported_layers <- c(
    "counts",
    "data",
    "scale.data"
  )

  invalid_layers <- setdiff(
    required_layers,
    supported_layers
  )

  if (length(invalid_layers))
    stop(
      "Unsupported required layer(s): ",
      paste(
        invalid_layers,
        collapse = ", "
      ),
      ". Supported values are: ",
      paste(
        supported_layers,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )

  min_cells <- as.integer(min_cells)
  min_features <- as.integer(min_features)

  if (
    length(min_cells) != 1L ||
    is.na(min_cells) ||
    min_cells < 1L
  ) {
    stop(
      "`min_cells` must be >= 1.",
      call. = FALSE
    )
  }

  if (
    length(min_features) != 1L ||
    is.na(min_features) ||
    min_features < 1L
  ) {
    stop(
      "`min_features` must be >= 1.",
      call. = FALSE
    )
  }

  if (ncol(object) < min_cells)
    stop(
      "Object contains ",
      ncol(object),
      " cells; at least ",
      min_cells,
      " required.",
      call. = FALSE
    )

  if (nrow(object[[assay]]) < min_features)
    stop(
      "Assay `",
      assay,
      "` contains ",
      nrow(object[[assay]]),
      " features; at least ",
      min_features,
      " required.",
      call. = FALSE
    )

  log(
    "Using assay:",
    assay
  )

  obj_assay <- object[[assay]]

  # ------------------------------------------------------------------------- #
  # Determine available layers
  # ------------------------------------------------------------------------- #

  if (inherits(obj_assay, "Assay5")) {
    available_layers <-
      SeuratObject::Layers(
        obj_assay
      )

  } else {
    available_layers <- character(0)

    if (
      "counts" %in% methods::slotNames(obj_assay) &&
      nrow(obj_assay@counts) > 0L &&
      ncol(obj_assay@counts) > 0L
    ) {
      available_layers <- c(
        available_layers,
        "counts"
      )
    }

    if (
      "data" %in% methods::slotNames(obj_assay) &&
      nrow(obj_assay@data) > 0L &&
      ncol(obj_assay@data) > 0L
    ) {
      available_layers <- c(
        available_layers,
        "data"
      )
    }

    if (
      "scale.data" %in% methods::slotNames(obj_assay) &&
      nrow(obj_assay@scale.data) > 0L &&
      ncol(obj_assay@scale.data) > 0L
    ) {
      available_layers <- c(
        available_layers,
        "scale.data"
      )
    }
  }

  # ------------------------------------------------------------------------- #
  # Seurat v5 multi-layer warning
  # ------------------------------------------------------------------------- #

  if (inherits(obj_assay, "Assay5")) {
    count_layers <- grep(
      "^counts(\\.|$)",
      available_layers,
      value = TRUE
    )

    data_layers <- grep(
      "^data(\\.|$)",
      available_layers,
      value = TRUE
    )

    if (
      length(count_layers) > 1L ||
      length(data_layers) > 1L
    ) {
      warning(
        "Assay `",
        assay,
        "` contains multiple expression layers. ",
        "Consider `Seurat::JoinLayers()` before GnRHcell analysis.",
        call. = FALSE
      )
    }
  }

  # ------------------------------------------------------------------------- #
  # Missing normalized data
  # ------------------------------------------------------------------------- #

  data_required <-
    "data" %in%
    required_layers

  data_available <-
    "data" %in%
    available_layers

  if (
    data_required &&
    !data_available
  ) {
    if (!isTRUE(auto_normalize))
      stop(
        "Normalized `data` layer is required but missing from assay `",
        assay,
        "`. Run `Seurat::NormalizeData()` or use ",
        "`auto_normalize = TRUE`.",
        call. = FALSE
      )

    counts_available <-
      "counts" %in%
      available_layers

    if (!counts_available)
      stop(
        "Cannot normalize assay `",
        assay,
        "` because the `counts` layer is unavailable.",
        call. = FALSE
      )

    log(
      "Normalized data missing; running NormalizeData().",
      type = "step"
    )

    object <- Seurat::NormalizeData(
      object,
      assay = assay,
      verbose = FALSE
    )

    obj_assay <- object[[assay]]

    available_layers <- if (
      inherits(obj_assay, "Assay5")
    ) {
      SeuratObject::Layers(
        obj_assay
      )
    } else {
      unique(
        c(
          available_layers,
          "data"
        )
      )
    }
  }

  # ------------------------------------------------------------------------- #
  # Validate all required layers
  # ------------------------------------------------------------------------- #

  missing_layers <- setdiff(
    required_layers,
    available_layers
  )

  if (length(missing_layers))
    stop(
      "Assay `",
      assay,
      "` is missing required layer(s): ",
      paste(
        missing_layers,
        collapse = ", "
      ),
      ". Available layers: ",
      paste(
        available_layers,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )

  # ------------------------------------------------------------------------- #
  # Validate dimensions of required matrices
  # ------------------------------------------------------------------------- #

  for (layer in required_layers) {
    mat <- tryCatch(
      Seurat::GetAssayData(
        object,
        assay = assay,
        layer = layer
      ),
      error = function(e)
        NULL
    )

    if (is.null(mat))
      stop(
        "Unable to retrieve layer `",
        layer,
        "` from assay `",
        assay,
        "`.",
        call. = FALSE
      )

    if (
      nrow(mat) == 0L ||
      ncol(mat) == 0L
    ) {
      stop(
        "Layer `",
        layer,
        "` in assay `",
        assay,
        "` is empty.",
        call. = FALSE
      )
    }

    if (ncol(mat) != ncol(object))
      stop(
        "Layer `",
        layer,
        "` contains ",
        ncol(mat),
        " cells but the Seurat object contains ",
        ncol(object),
        ".",
        call. = FALSE
      )
  }

  invisible(object)
}
