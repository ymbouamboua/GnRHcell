# ==============================================================================
# MODULE V1
# ==============================================================================

# .gnrh_modules <- function(genes) {
#
#   list(
#
#     core = .match_genes(c(
#       "GNRH1","FEZF1","ISL1",
#       "OTX2","SIX3","SIX6",
#       "DLX1","DLX2","DLX5","DLX6"
#     ), genes),
#
#     mig = .match_genes(c(
#       "ANOS1","PROKR2","PROK2",
#       "NRP1","NRP2","SEMA3A",
#       "SEMA3C","SEMA3F",
#       "ROBO1","ROBO2",
#       "L1CAM","DCX"
#     ), genes),
#
#     neuro = .match_genes(c(
#       "KISS1R","TAC3","TACR3",
#       "GNRHR","PCSK1","PCSK2",
#       "SCG2","CHGA","CHGB",
#       "CPE","VGF","SYP","RAB3A"
#       #"ESR1","PGR","AR" # hormone
#     ), genes)
#   )
# }


# build_stage_modules <- function(genes) {
#
#   list(
#     # # Neurogenesis / placode progenitors
#     # # Key studies: Wray et al., 1989; Forni & Wray 2012; Wang et al., 2022
#     # neurogenesis = .match_genes(c(
#     #   "SOX2",   # stem cell identity, olfactory placode progenitors
#     #   "PAX6",   # olfactory placode patterning http://dx.doi.org/10.1016/j.mce.2017.02.030
#     #   "TUJ1",   # http://dx.doi.org/10.1016/j.mce.2017.02.030
#     #   "FOXG1",  # forebrain/placode neural progenitors
#     #   "SIX3",   # cranial placode development
#     #   "EYA1",   # placode development
#     #   "HES1",   # progenitor maintenance
#     #   "ASCL1",  # neuronal differentiation
#     #   "NEUROG1",# proneural TF
#     #   "NEUROD1",# neuronal differentiation
#     #   "NOTCH1", # progenitor maintenance
#     #   "DLL1"    # Notch ligand, maintains progenitors
#     # ), genes),
#
#     # Specification / early GnRH lineage TFs
#     # References: Taroc et al., 2017; Wray 2010; Wang et al., 2022
#     identity = .match_genes(c(
#       "GNRH1",  # hallmark of GnRH neurons
#       "ISL1",   # transcription factor, early GnRH specification
#       "DLX5",   # TF, anterior forebrain and GnRH lineage
#       "DLX6",   # TF, forebrain/GnRH neuron specification
#       "LHX2",   # TF involved in forebrain and olfactory placode
#       "LHX3",   # TF implicated in hypothalamic neuron differentiation
#       "OTX2",   # TF, GnRH neuron specification
#       "FEZF1",  # TF, olfactory placode/GnRH specification
#       "POU3F2", # TF, neuronal differentiation
#       "NKX2-1"  # hypothalamic fate, early GnRH progenitors
#     ),genes),
#
#     # Migrating / axon guidance and migration machinery
#     # References: Forni & Wray 2015; Taroc et al., 2017; Wang et al., 2022
#     migrating = .match_genes(c(
#       "GNRH1",
#       "ISL1","DLX5","DLX6","FEZF1","OTX2",
#       "ANOS1","PROK2","PROKR2", # Kallmann syndrome genes
#       "SEMA3A","SEMA3F",        # axon guidance cues
#       "ROBO1","ROBO3",          # axon guidance receptors
#       "CNTN2","NCAM1",          # cell adhesion molecules
#       "NRP1","NRP2",            # guidance receptors
#       "HS6ST1",                 # heparan sulfate modification, migration
#       "FGFR1","FGF8","IL17RD",  # FGFR pathway, migration
#       "TUBB3","STMN2","STMN3","GAP43","DCX" # cytoskeleton, migration
#     ),genes),
#
#     # Mature / neuroendocrine function
#     # References: Wray 2010; Wang et al., 2022
#     mature = .match_genes(c(
#       "GNRH1",
#       "PCSK1",
#       "SCG2",
#       "VGF",
#       "TAC3",
#       "GNRHR", # receptor for GnRH
#       "KISS1R",# receptor for kisspeptin
#       "AR",     # "ESR1", steroid hormone receptors
#       "SYP","SYT1",    # synaptic vesicle machinery
#       "CHGA","CHGB"    # neurosecretory vesicles
#     ), genes),
#
#     secreting = .match_genes(c(
#       "GNRH1",
#       # peptide processing / maturation
#       "PCSK1", "PCSK2", "CPE",
#       # dense-core vesicle / regulated secretory pathway
#       "CHGA", "CHGB", "SCG2", "VGF",
#       # synaptic vesicle docking, trafficking, release
#       "SYP", "SYT1", "SNAP25", "STX1A", "VAMP2",
#       "RAB3A", "RIMS1", "UNC13A",
#       # activity-dependent / secretory response
#       "FOS", "JUN", "EGR1",
#       # reproductive neuropeptide responsiveness
#       "KISS1R", "TAC3", "TACR3"
#     ), genes)
#
#   )
#
# }


# ==============================================================================
# MODULE V2
# ==============================================================================


#' Build GnRH detection marker modules
#'
#' Internal helper that maps predefined GnRH-associated marker
#' genes to the genes available in the input expression matrix.
#'
#' The returned modules are used by \code{\link{detect_gnrh}}
#' to support GnRH neuron detection through core identity,
#' migration, and neuroendocrine marker signals.
#'
#' @param genes Character vector of gene names present in the
#' expression matrix.
#'
#' @return A named list of character vectors containing matched genes:
#' \describe{
#'   \item{\code{core}}{Core GnRH lineage and identity markers.}
#'   \item{\code{mig}}{Migration and axon-guidance markers.}
#'   \item{\code{neuro}}{Neuroendocrine and neuronal function markers.}
#' }
#'
#' @keywords internal
#' @noRd
.gnrh_modules <- function(genes) {

  list(
    core = .match_genes(c(
      "GNRH1",
      "ISL1",
      "SIX3", "SIX6",
      "DLX1","DLX2","DLX5", "DLX6",
      "OTX2",
      "KISS1R"
    ), genes),

    mig = .match_genes(c(
      "ANOS1", "PROK2", "PROKR2",
      "FGFR1", "FGF8", "IL17RD", "HS6ST1",
      "SEMA3A", "SEMA3C", "SEMA3F",
      "NRP1", "NRP2",
      "ROBO1", "ROBO2", "ROBO3",
      "L1CAM", "NCAM1", "CNTN2",
      "GAP43", "STMN2", "STMN3", "PLXNA3" ,"SLIT1",
      "MAP1B", "RIPOR2", "SPOCK1", "UNC5D"
    ), genes),

    neuro = .match_genes(c(
      "GNRH1",
      "KISS1R", "TAC3", "TACR3", "TAC1",
      "GNRHR",
      "PCSK1", "PCSK2", "CPE",
      "SCG2", "SCG5", "CHGA", "CHGB", "VGF",
      "SYP", "VAMP2", "SNAP25",
      "RAB3A", "RAB3B", "RAB3C",
      "PTPRN", "BAIAP3", "DOC2B",
      "ECEL1", "RASD1",
      "HCN1", "NALCN", "SCN3A"
    ), genes)
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




