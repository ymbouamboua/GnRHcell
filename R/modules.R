
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
        "FEZF1",
        "ISL1",
        "SIX6",
        "ECEL1"
      ),
      supportive = c(
        "OTX2",
        "SIX3",
        "DLX1",
        "DLX2",
        "DLX5",
        "DLX6",
        "PBX3",
        "ARX",
        "FOXG1"
      ),
      contextual = c(
        "GAD2",
        "HESX1",
        "FGFR1",
        "RBFOX1",
        "MYT1L",
        "BCL11B"
      )
    ),

    migration = list(
      primary = c(
        "PROKR2",
        "NSMF",
        "ROBO3"
      ),
      supportive = c(
        "NRP1",
        "NRP2",
        "ROBO1",
        "ROBO2",
        "CXCR4",
        "L1CAM",
        "DCX"
      ),
      contextual = c(
        "DCC",
        "CNTN1",
        "ITGAV",
        "ACKR3",
        "PLXNA1",
        "PLXNA2",
        "PLXNA3",
        "PLXNA4",
        "RELN",
        "DSCAM"
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
      ),
      contextual = c(
        "SCN2A",
        "TAC1",
        "PTPRN2"
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
        "NTN1",
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
        .match_genes,
        genes = genes
      )
    }
  )
}




#' Alternative neuronal identity modules
#'
#' Defines neuronal and neuroendocrine programs that may partially overlap
#' with GnRH-associated transcriptional programs.
#'
#' These modules do not reject cells with detectable \code{GNRH1}. They are
#' used only when identifying exploratory GNRH1-negative transcriptomic
#' candidates.
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
    kndy = c("KISS1", "TAC3", "PDYN"),
    pomc = "POMC",
    agrp_npy = c("AGRP", "NPY"),
    avp = "AVP",
    oxt = "OXT",
    crh = "CRH",
    trh = "TRH",
    sst = "SST"
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
#'
#' @keywords internal
#' @noRd
.score_gnrh_alternatives <- function(expr, modules) {
  n <- ncol(expr)

  if (!length(modules)) {
    return(list(
      score = rep(0, n),
      hits = integer(n),
      strong = rep(FALSE, n),
      scores = matrix(numeric(0), nrow = n),
      hit_matrix = matrix(integer(0), nrow = n)
    ))
  }

  hits_one <- function(g) {
    if (!length(g)) return(integer(n))

    as.integer(
      Matrix::colSums(
        expr[g, , drop = FALSE] > 0
      )
    )
  }

  score_one <- function(g) {
    if (!length(g)) return(rep(0, n))

    hits <- hits_one(g)

    # Avoid score=1 from a singleton merely because one gene was detected.
    hits / max(length(g), 2L)
  }

  strong_one <- function(g) {
    if (!length(g)) return(rep(FALSE, n))

    hits <- hits_one(g)

    if (length(g) == 1L) {
      hits >= 1L
    } else {
      hits >= min(2L, length(g))
    }
  }

  scores <- do.call(cbind, lapply(modules, score_one))
  hit_matrix <- do.call(cbind, lapply(modules, hits_one))
  strong_matrix <- do.call(cbind, lapply(modules, strong_one))

  colnames(scores) <- names(modules)
  colnames(hit_matrix) <- names(modules)
  colnames(strong_matrix) <- names(modules)

  list(
    score = apply(scores, 1, max, na.rm = TRUE),
    hits = apply(hit_matrix, 1, max, na.rm = TRUE),
    strong = apply(strong_matrix, 1, any),
    scores = scores,
    hit_matrix = hit_matrix,
    strong_matrix = strong_matrix
  )
}


#' Build GnRH developmental stage modules
#'
#' Constructs predefined marker modules for developmental staging
#' of GnRH neurons from single-cell RNA-seq data.
#'
#' The modules represent major biological states of GnRH neuron
#' development, including identity specification, migration and
#' maturation activity.
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
#' }
#'
#' @keywords internal
#' @noRd
.build_stage_modules <- function(genes) {
  list(
    identity = .match_genes(
      c(
        "FEZF1",
        "ISL1",
        "SIX6",
        "SIX3",
        "OTX2",
        "DLX1",
        "DLX2",
        "DLX5",
        "DLX6",
        "PBX3",
        "RASD1",
        "RMST",
        "ECEL1"
      ),
      genes
    ),

    migrating = .match_genes(
      c(
        "PROKR2",
        "NSMF",
        "SEMA3C",
        "SEMA3F",
        "ROBO2",
        "ROBO3",
        "RIPOR2",
        "PLXNA3",
        "SLIT1",
        "CXCR4"
      ),
      genes
    ),

    mature = .match_genes(
      c(
        "KISS1R",
        "DOC2B",
        "PTPRN",
        "BAIAP3",
        "SCG2",
        "SCG5",
        "HCN1",
        "NALCN"
      ),
      genes
    )
  )
}



#' @keywords internal
#' @noRd
.build_migration_core <- function(genes) {
  .match_genes(
    c(
      "PROKR2",
      "NSMF",
      "SEMA3C",
      "SEMA3F",
      "ROBO2",
      "ROBO3",
      "RIPOR2",
      "PLXNA3",
      "SLIT1"
    ),
    genes
  )
}


#' Build GnRH secretory module
#'
#' Constructs core and supportive neuroendocrine secretory programs used
#' independently from developmental stage classification.
#'
#' @param genes Character vector containing genes available in the expression
#'   matrix.
#'
#' @return A named list containing matched core and supportive secretory genes.
#'
#' @keywords internal
#' @noRd
.build_secretory_module <- function(genes) {
  list(
    core = .match_genes(
      c(
        "PCSK1",
        "PCSK2",
        "CHGA",
        "CHGB"
      ),
      genes
    ),

    supportive = .match_genes(
      c(
        "CPE",
        "VGF",
        "SCG2",
        "SCG5",
        "PTPRN",
        "BAIAP3",
        "DOC2B"
      ),
      genes
    )
  )
}
