#' GnRH module weights
#'
#' Internal weights used to combine curated GnRH marker modules.
#'
#' @return Named numeric vector of module weights.
#'
#' @keywords internal
#' @noRd
.gnrh_module_weights <- function() {
  
  c(
    identity_primary = 3.0,
    identity_supportive = 1.5,
    
    migration_primary = 1.0,
    migration_supportive = 0.5,
    
    neuroendocrine_primary = 1.5,
    neuroendocrine_supportive = 0.75,
    
    hormone_supportive = 0.25,
    
    guidance_environment = 0
  )
}



#' Compute a weighted GnRH module score
#'
#' @param scores Named list or data frame containing module scores.
#' @param weights Named numeric vector of module weights.
#'
#' @return Numeric vector containing the weighted GnRH score.
#'
#' @keywords internal
#' @noRd
.weight_gnrh_scores <- function(
    scores,
    weights = .gnrh_module_weights()
) {
  
  common <- intersect(
    names(weights),
    names(scores)
  )
  
  if (length(common) == 0L) {
    stop(
      "No GnRH module scores matched the supplied weights.",
      call. = FALSE
    )
  }
  
  w <- weights[common]
  
  keep <- is.finite(w) & w != 0
  
  common <- common[keep]
  w <- w[keep]
  
  if (length(common) == 0L) {
    stop(
      "All matched GnRH module weights are zero or non-finite.",
      call. = FALSE
    )
  }
  
  score_matrix <- do.call(
    cbind,
    scores[common]
  )
  
  score_matrix <- as.matrix(score_matrix)
  
  weighted <- sweep(
    score_matrix,
    MARGIN = 2,
    STATS = w,
    FUN = "*"
  )
  
  rowSums(
    weighted,
    na.rm = TRUE
  ) / sum(abs(w))
}

