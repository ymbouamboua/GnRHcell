
#' Curated GnRH gene modules
#'
#' Defines literature-informed gene programs associated with GnRH neuronal
#' identity, migration, neuroendocrine maturation, hormonal responsiveness,
#' and extracellular guidance.
#'
#' @param genes Character vector containing genes available in the expression
#'   matrix.
#'
#' @return A nested named list containing matched gene symbols.
#'
#' @keywords internal
#' @noRd
.gnrh_modules <- function(genes) {
  modules <- list(
    identity = list(
      primary = c(
        #"GNRH1",
        "FEZF1",
        "ISL1"
      ),
      supportive = c(
        "OTX2",
        "SIX3",
        "SIX6",
        "DLX1",
        "DLX2",
        "DLX5",
        "DLX6"
      )
    ),
    migration = list(
      primary = c(
        "PROKR2",
        "NRP1",
        "NRP2",
        "PLXNA1",
        "ROBO1",
        "ROBO2",
        "ROBO3",
        "CXCR4",
        "NSMF"
      ),
      supportive = c(
        "L1CAM",
        "DCX"
      )
    ),
    neuroendocrine = list(
      primary = c(
        "KISS1R",
        "GNRHR",
        "PCSK1",
        "PCSK2"
      ),
      supportive = c(
        "CPE",
        "SCG2",
        "CHGA",
        "CHGB",
        "VGF",
        "SYP",
        "RAB3A"
      )
    ),
    hormone = list(
      supportive = c(
        "ESR1",
        "ESR2",
        "PGR",
        "AR"
      )
    ),
    guidance_environment = list(
      supportive = c(
        "ANOS1",
        "PROK2",
        "SEMA3A",
        "SEMA3C",
        "SEMA3E",
        "SEMA3F",
        "SEMA7A",
        "SLIT1",
        "SLIT2",
        "SLIT3",
        "CXCL12"
      )
    )
  )
  lapply(
    modules,
    function(module) {
      lapply(
        module,
        function(x) .match_genes(x, genes)
      )
    }
  )
}


#' Alternative neuronal identity modules
#'
#' Defines neuronal and neuroendocrine programs that may partially overlap
#' with GnRH-associated transcriptional programs.
#'
#' These modules are used only to guard against false-positive dropout rescue
#' and do not penalize cells with direct \code{GNRH1} evidence.
#'
#' @param genes Character vector containing genes available in the expression
#'   matrix.
#'
#' @return A named list containing matched gene symbols for each alternative
#'   neuronal program.
#'
#' @keywords internal
#' @noRd
.gnrh_alternative_modules <- function(genes) {
  modules <- list(
    kndy = c(
      "KISS1",
      "TAC3",
      "PDYN"
    ),
    pomc = c(
      "POMC"
    ),
    agrp_npy = c(
      "AGRP",
      "NPY"
    ),
    avp = c(
      "AVP"
    ),
    oxt = c(
      "OXT"
    ),
    crh = c(
      "CRH"
    ),
    trh = c(
      "TRH"
    ),
    sst = c(
      "SST"
    )
  )
  lapply(
    modules,
    .match_genes,
    genes = genes
  )
}



#' Score alternative neuronal identity programs
#'
#' Computes per-cell scores for alternative neuronal or neuroendocrine
#' programs and returns the strongest alternative identity signal.
#'
#' @param expr Gene-by-cell expression matrix.
#' @param modules Named list containing alternative marker gene sets.
#'
#' @return A list containing the maximum alternative score per cell,
#'   maximum marker hit counts, individual program scores, and program-specific
#'   hit counts.
#'
#' @keywords internal
#' @noRd
.score_gnrh_alternatives <- function(
    expr,
    modules
) {
  n <- ncol(expr)
  if (!length(modules)) {
    return(
      list(
        score = rep(0, n),
        hits = integer(n),
        scores = matrix(
          numeric(0),
          nrow = n,
          ncol = 0
        ),
        hit_matrix = matrix(
          integer(0),
          nrow = n,
          ncol = 0
        )
      )
    )
  }
  score_one <- function(x) {
    if (!length(x)) {
      return(rep(0, n))
    }
    Matrix::colMeans(
      expr[
        x,
        ,
        drop = FALSE
      ] > 0
    )
  }
  hits_one <- function(x) {
    if (!length(x)) {
      return(integer(n))
    }
    as.integer(
      Matrix::colSums(
        expr[
          x,
          ,
          drop = FALSE
        ] > 0
      )
    )
  }
  scores <- do.call(
    cbind,
    lapply(
      modules,
      score_one
    )
  )
  hit_matrix <- do.call(
    cbind,
    lapply(
      modules,
      hits_one
    )
  )
  colnames(scores) <- names(modules)
  colnames(hit_matrix) <- names(modules)
  score <- apply(
    scores,
    1,
    max,
    na.rm = TRUE
  )
  hits <- apply(
    hit_matrix,
    1,
    max,
    na.rm = TRUE
  )
  list(
    score = as.numeric(score),
    hits = as.integer(hits),
    scores = scores,
    hit_matrix = hit_matrix
  )
}



#' Build GnRH developmental stage modules
#'
#' Constructs predefined marker modules for developmental staging
#' of GnRH neurons from single-cell RNA-seq data.
#'
#' The modules represent major biological states of GnRH neuron
#' development, including identity specification, migration,
#' maturation, and secretory activity.
#'
#' These modules are used internally by staging functions such as
#' \code{\link{stage_gnrh}} to score GnRH cells across developmental
#' programs.
#'
#' @param genes Character vector of gene names present in the
#' expression matrix.
#'
#' @return A named list of character vectors containing matched genes:
#' \describe{
#'   \item{\code{identity}}{GnRH lineage identity and specification markers.}
#'   \item{\code{migrating}}{Migration, axon-guidance, adhesion, and cytoskeletal markers.}
#'   \item{\code{mature}}{Neuroendocrine maturation and synaptic function markers.}
#'   \item{\code{secreting}}{Peptide processing, dense-core vesicle, and secretion markers.}
#' }
#'
#' @keywords internal
#' @noRd
build_stage_modules <- function(genes) {

  list(
    identity = .match_genes(c(
      "GNRH1",
      "ISL1",
      "SIX3", "SIX6",
      "DLX1","DLX2", "DLX5", "DLX6",
      "OTX2",
      "KISS1R",
      "PBX3",
      "RASD1",
      "RMST",
      "MIAT"
    ), genes),

    migrating = .match_genes(c(
      "GNRH1",
      "SEMA3C", "SEMA3A", "SEMA3F",
      "NRP1", "NRP2",
      "ROBO1", "ROBO2", "ROBO3",
      "L1CAM", "NCAM1", "CNTN2",
      "STMN2", "STMN3",
      "GAP43", "MAP1B",
      "RIPOR2", "SPOCK1", "UNC5D",
      "PLXNA3", "SLIT1"
    ), genes),

    mature = .match_genes(c(
      "GNRH1",
      "KISS1R",
      "ISL1",
      "DOC2B",
      "PTPRN",
      "BAIAP3",
      "ECEL1",
      "SCG2", "SCG5",
      "SYP", "VAMP2", "SNAP25",
      "RAB3A", "RAB3B", "RAB3C",
      "HCN1", "NALCN", "SCN3A",
      "GAD1"
    ), genes),

    secreting = .match_genes(c(
      "GNRH1",
      "PCSK1", "PCSK2", "CPE",
      "SCG2", "SCG5",
      "CHGA", "CHGB", "VGF",
      "PTPRN", "BAIAP3", "DOC2B",
      "SYP", "VAMP2", "SNAP25", "STX1A",
      "SYT1", "RAB3A",
      "FOS", "JUN", "EGR1",
      "ESR1","PGR","AR"
    ), genes)
  )
}




