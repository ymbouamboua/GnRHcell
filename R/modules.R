# ----------------------------------------------------------------------- #
# Curated GnRH modules
# ----------------------------------------------------------------------- #
#' Curated GnRH gene modules
#' @keywords internal
#' @noRd
.gnrh_modules <- function(genes) {
  modules <- list(
    identity=list(
      primary=c("FEZF1","ISL1","SIX6","ECEL1"),
      supportive=c("OTX2","SIX3","DLX1","DLX2","DLX5","DLX6","PBX3","ARX","FOXG1"),
      contextual=c("GAD2","HESX1","FGFR1","RBFOX1","MYT1L","BCL11B")
    ),
    migration=list(
      primary=c("PROKR2","NSMF","ROBO3"),
      supportive=c("NRP1","NRP2","ROBO1","ROBO2","CXCR4","L1CAM","DCX"),
      contextual=c("DCC","CNTN1","ITGAV","ACKR3","PLXNA1","PLXNA2","PLXNA3","PLXNA4","RELN","DSCAM")
    ),
    neuroendocrine=list(
      primary=c("KISS1R","GNRHR","PCSK1","PCSK2"),
      supportive=c("CPE","SCG2","CHGA","CHGB","VGF","SYP","RAB3A"),
      contextual=c("SCN2A","TAC1","PTPRN2")
    ),
    hormone=list(supportive=c("ESR1","ESR2","PGR","AR")),
    guidance_environment=list(supportive=c("ANOS1","PROK2","NTN1","SEMA3A","SEMA3C","SEMA3E","SEMA3F","SEMA7A","SLIT1","SLIT2","SLIT3","CXCL12"))
  )
  lapply(modules,function(x) lapply(x,.match_genes,genes=genes))
}
# ----------------------------------------------------------------------- #
# Alternative identities
# ----------------------------------------------------------------------- #
#' Alternative neuronal identity modules
#' @keywords internal
#' @noRd
.gnrh_alternative_modules <- function(genes) {
  modules <- list(
    kndy=c("KISS1","TAC3","PDYN"),
    pomc="POMC",
    agrp_npy=c("AGRP","NPY"),
    avp="AVP",
    oxt="OXT",
    crh="CRH",
    trh="TRH",
    sst="SST"
  )
  lapply(modules,.match_genes,genes=genes)
}
# ----------------------------------------------------------------------- #
# Alternative scoring
# ----------------------------------------------------------------------- #
#' Score alternative neuronal identity programs
#' @keywords internal
#' @noRd
.score_gnrh_alternatives <- function(expr,modules) {
  n <- ncol(expr)
  if (!length(modules)) return(list(score=rep(0,n),hits=integer(n),strong=rep(FALSE,n),scores=matrix(numeric(0),nrow=n),hit_matrix=matrix(integer(0),nrow=n),strong_matrix=matrix(logical(0),nrow=n)))
  hits_one <- function(g) if (!length(g)) integer(n) else as.integer(Matrix::colSums(expr[g,,drop=FALSE]>0))
  score_one <- function(g) if (!length(g)) rep(0,n) else hits_one(g)/max(length(g),2L)
  strong_one <- function(g) {
    if (!length(g)) return(rep(FALSE,n))
    h <- hits_one(g)
    if (length(g)==1L) h>=1L else h>=min(2L,length(g))
  }
  scores <- do.call(cbind,lapply(modules,score_one))
  hits <- do.call(cbind,lapply(modules,hits_one))
  strong <- do.call(cbind,lapply(modules,strong_one))
  colnames(scores) <- colnames(hits) <- colnames(strong) <- names(modules)
  list(
    score=apply(scores,1,max,na.rm=TRUE),
    hits=apply(hits,1,max,na.rm=TRUE),
    strong=apply(strong,1,any),
    scores=scores,
    hit_matrix=hits,
    strong_matrix=strong
  )
}
# ----------------------------------------------------------------------- #
# Developmental stage modules
# ----------------------------------------------------------------------- #
#' Build GnRH developmental stage modules
#' @keywords internal
#' @noRd
.build_stage_modules <- function(genes) {
  list(
    early=.match_genes(c("FEZF1","SIX6","SIX3","OTX2","DLX1","DLX2","DLX5","DLX6"),genes),
    migrating=.match_genes(c("PROKR2","NSMF","SEMA3C","SEMA3F","ROBO2","ROBO3","RIPOR2","PLXNA3","SLIT1","CXCR4"),genes),
    mature=.match_genes(c("KISS1R","DOC2B","PTPRN","BAIAP3","SCG2","SCG5","HCN1","NALCN"),genes)
  )
}
# ----------------------------------------------------------------------- #
# Migration core
# ----------------------------------------------------------------------- #
#' Build GnRH migration-core module
#' @keywords internal
#' @noRd
.build_migration_core <- function(genes) {
  .match_genes(c("PROKR2","NSMF","SEMA3C","SEMA3F","ROBO2","ROBO3","RIPOR2","PLXNA3","SLIT1"),genes)
}
# ----------------------------------------------------------------------- #
# Secretory module
# ----------------------------------------------------------------------- #
#' Build GnRH secretory module
#' @keywords internal
#' @noRd
.build_secretory_module <- function(genes) {
  list(
    core=.match_genes(c("PCSK1","PCSK2","CHGA","CHGB"),genes),
    supportive=.match_genes(c("CPE","VGF","SCG2","SCG5","PTPRN","BAIAP3","DOC2B"),genes)
  )
}
