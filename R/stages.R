#' Assign dominant developmental stage
#' @keywords internal
#' @noRd
assign_stage <- function(scores,migration_core_hits=NULL,min_migration_hits=2L) {
  scores <- as.matrix(scores)
  req <- c("early","migrating","mature")
  miss <- setdiff(req,colnames(scores))
  if (length(miss)) stop("Missing developmental-stage score(s): ",paste(miss,collapse=", "),call.=FALSE)
  scores <- scores[,req,drop=FALSE]
  if (!nrow(scores)) return(list(stage=character(),raw_stage=character(),top_score=numeric(),second_score=numeric(),margin=numeric(),migration_filtered=logical()))
  finite <- is.finite(scores)
  x <- scores
  x[!finite] <- -Inf
  noinfo <- rowSums(finite)==0L
  top_i <- max.col(x,ties.method="first")
  raw <- colnames(x)[top_i]
  raw[noinfo] <- NA_character_
  top <- x[cbind(seq_len(nrow(x)),top_i)]
  second <- apply(x,1,function(z) { z <- sort(z,decreasing=TRUE); if (length(z)<2L || !is.finite(z[2])) NA_real_ else z[2] })
  top[noinfo] <- NA_real_
  margin <- top-second
  stage <- raw
  filtered <- rep(FALSE,nrow(x))
  if (!is.null(migration_core_hits)) {
    if (length(migration_core_hits)!=nrow(x)) stop("`migration_core_hits` must contain one value per cell.",call.=FALSE)
    hits <- suppressWarnings(as.numeric(migration_core_hits))
    hits[!is.finite(hits)] <- 0
    filtered <- !is.na(stage) & stage=="migrating" & hits<min_migration_hits
    if (any(filtered)) {
      alt <- x[filtered,c("early","mature"),drop=FALSE]
      stage[filtered] <- colnames(alt)[max.col(alt,ties.method="first")]
    }
  }
  list(stage=stage,raw_stage=raw,top_score=top,second_score=second,margin=margin,migration_filtered=filtered)
}
# ----------------------------------------------------------------------- #
# Developmental-state interpretation
# ----------------------------------------------------------------------- #
#' Derive a conservative developmental state
#' @keywords internal
#' @noRd
.derive_developmental_state <- function(raw_stage,margin,positive,migration_filtered,secretory_supported,stage_margin) {
  n <- length(raw_stage)
  if (!all(vapply(list(margin,positive,migration_filtered,secretory_supported),length,integer(1))==n)) stop("Developmental-state inputs must have equal lengths.",call.=FALSE)
  state <- rep("undetermined",n)
  state[!positive] <- "non-gnrh"
  uncertain <- positive & (is.na(raw_stage) | !is.finite(margin) | margin<stage_margin | migration_filtered)
  state[uncertain] <- "transitional"
  ok <- positive & !uncertain
  state[ok & raw_stage=="early"] <- "early"
  state[ok & raw_stage=="migrating"] <- "migrating"
  state[ok & raw_stage=="mature" & !secretory_supported] <- "post-migratory"
  state[ok & raw_stage=="mature" & secretory_supported] <- "mature"
  state
}
# ----------------------------------------------------------------------- #
# Stage GnRH cells
# ----------------------------------------------------------------------- #
#' Stage GnRH lineage cells
#' @param object A Seurat object previously processed with `detect_gnrh()`.
#' @param assay Assay used for developmental and secretory scoring.
#' @param layer Expression layer used for module scoring.
#' @param min_migration_hits Minimum migration-core hits.
#' @param min_secretory_core_hits Minimum core secretory hits.
#' @param min_secretory_supportive_hits Minimum supportive secretory hits.
#' @param expression_weight Weight assigned to module expression.
#' @param detection_weight Weight assigned to module detection fraction.
#' @param stage_margin Minimum top-vs-second score margin.
#' @param verbose Print progress messages.
#' @return A Seurat object with GnRH developmental-state metadata.
#' @export
stage_gnrh <- function(object,assay="RNA",layer="data",min_migration_hits=2L,min_secretory_core_hits=1L,min_secretory_supportive_hits=2L,expression_weight=0.5,detection_weight=0.5,stage_margin=0.10,verbose=TRUE) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  .int <- function(x,nm) {
    if (length(x)!=1L || is.na(x) || !is.finite(x) || x<0 || x!=floor(x)) stop("`",nm,"` must be a single non-negative integer.",call.=FALSE)
    as.integer(x)
  }
  min_migration_hits <- .int(min_migration_hits,"min_migration_hits")
  min_secretory_core_hits <- .int(min_secretory_core_hits,"min_secretory_core_hits")
  min_secretory_supportive_hits <- .int(min_secretory_supportive_hits,"min_secretory_supportive_hits")
  if (length(stage_margin)!=1L || !is.finite(stage_margin) || stage_margin<0) stop("`stage_margin` must be a single non-negative number.",call.=FALSE)
  log <- .msg(verbose)
  object <- validate_input(object,assay=assay,required_layers=layer,auto_normalize=FALSE,verbose=FALSE)
  md <- object[[]]
  if (!"gnrh_status" %in% colnames(md)) stop("`gnrh_status` is missing. Run `detect_gnrh()` first.",call.=FALSE)
  positive <- !is.na(md$gnrh_status) & as.character(md$gnrh_status)=="pos"
  negative <- !positive
  log(sprintf("Staging %s GnRH-positive cells",format(sum(positive),big.mark=",")))
  expr <- .get_expr(object,assay=assay,layer=layer)
  genes <- rownames(expr)
  req <- c("early","migrating","mature")
  modules <- .build_stage_modules(genes)
  miss <- setdiff(req,names(modules))
  if (length(miss)) stop("Developmental-stage modules are incomplete. Missing: ",paste(miss,collapse=", "),call.=FALSE)
  dev_modules <- modules[req]
  mod <- .score_modules(expr=expr,modules=dev_modules,expression_weight=expression_weight,detection_weight=detection_weight,input_type=if (grepl("^counts($|\\.)",layer)) "counts" else "normalized")
  expr_score <- as.data.frame(mod$score,check.names=FALSE)
  frac_score <- as.data.frame(mod$fraction,check.names=FALSE)
  stage_score <- as.data.frame(mod$integrated,check.names=FALSE)
  mig_core <- .build_migration_core(genes)
  mig_hits <- if (length(mig_core)) as.integer(Matrix::colSums(expr[mig_core,,drop=FALSE]>0)) else integer(ncol(expr))
  a <- assign_stage(stage_score,mig_hits,min_migration_hits)
  raw <- a$raw_stage
  stage <- a$stage
  reassigned <- !is.na(raw) & !is.na(stage) & raw!=stage
  reason <- rep("max_score",length(stage))
  reason[is.na(raw)] <- "insufficient_stage_signal"
  reason[a$migration_filtered] <- "migration_core_filter"
  confident <- positive & !is.na(stage) & !a$migration_filtered & is.finite(a$margin) & a$margin>=stage_margin
  resolution <- ifelse(negative,"non-gnrh",ifelse(confident,"resolved","transitional"))
  stage[negative] <- "non-gnrh"
  reassigned[negative] <- FALSE
  confident[negative] <- FALSE
  reason[negative] <- "non_gnrh"
  sec <- .build_secretory_module(genes)
  if (!all(c("core","supportive") %in% names(sec))) stop("Secretory module must contain `core` and `supportive` gene sets.",call.=FALSE)
  core_hits <- if (length(sec$core)) as.integer(Matrix::colSums(expr[sec$core,,drop=FALSE]>0)) else integer(ncol(expr))
  sup_hits <- if (length(sec$supportive)) as.integer(Matrix::colSums(expr[sec$supportive,,drop=FALSE]>0)) else integer(ncol(expr))
  sec_hits <- core_hits+sup_hits
  sec_supported <- positive & core_hits>=min_secretory_core_hits & (core_hits>=max(2L,min_secretory_core_hits) | sup_hits>=min_secretory_supportive_hits)
  sec_state <- ifelse(negative,"non-gnrh",ifelse(sec_supported,"supported","limited"))
  dev_state <- .derive_developmental_state(raw_stage=raw,margin=a$margin,positive=positive,migration_filtered=a$migration_filtered,secretory_supported=sec_supported,stage_margin=stage_margin)
  dev_index <- stage_score$mature-stage_score$migrating
  # Keep the dominant-program diagnostic restricted to the detected lineage.
  # Otherwise every background cell is misleadingly labelled as a GnRH stage.
  raw_display <- raw
  raw_display[negative] <- "non-gnrh"
  object$gnrh_stage_raw <- factor(
    raw_display,
    levels = c(req, "non-gnrh")
  )
  object$gnrh_stage_reassigned <- reassigned
  object$gnrh_stage_reason <- factor(reason,levels=c("max_score","migration_core_filter","insufficient_stage_signal","non_gnrh"))
  object$gnrh_stage_score <- a$top_score
  object$gnrh_stage_second_score <- a$second_score
  object$gnrh_stage_margin <- a$margin
  object$gnrh_stage_confident <- confident
  object$gnrh_stage_resolution <- factor(resolution,levels=c("resolved","transitional","non-gnrh"))
  object$gnrh_stage <- factor(dev_state,levels=c("early","migrating","post-migratory","mature","transitional","undetermined","non-gnrh"))
  object$gnrh_developmental_index <- dev_index
  object$gnrh_early_specification_score <- stage_score$early
  object$gnrh_lineage_identity_score <- if ("gnrh_identity_score" %in% colnames(md)) as.numeric(md$gnrh_identity_score) else stage_score$early
  obsolete <- c("gnrh_stage_legacy","gnrh_developmental_state","gnrh_stage_identity_score","gnrh_stage_identity_expression","gnrh_stage_identity_fraction")
  for (nm in intersect(obsolete,colnames(object[[]]))) object[[nm]] <- NULL
  for (nm in req) {
    object[[paste0("gnrh_stage_",nm,"_score")]] <- stage_score[[nm]]
    object[[paste0("gnrh_stage_",nm,"_expression")]] <- expr_score[[nm]]
    object[[paste0("gnrh_stage_",nm,"_fraction")]] <- frac_score[[nm]]
  }
  object$gnrh_migration_core_hits <- mig_hits
  object$gnrh_secretory_core_hits <- core_hits
  object$gnrh_secretory_supportive_hits <- sup_hits
  object$gnrh_secretory_hits <- sec_hits
  object$gnrh_secretory_supported <- sec_supported
  object$gnrh_secretory <- factor(sec_state,levels=c("limited","supported","non-gnrh"))
  pars <- list(assay=assay,layer=layer,score_input_type=mod$input_type,expression_weight=expression_weight,detection_weight=detection_weight,stage_margin=stage_margin,min_migration_hits=min_migration_hits,min_secretory_core_hits=min_secretory_core_hits,min_secretory_supportive_hits=min_secretory_supportive_hits)
  object@misc$gnrh_stage_modules <- dev_modules
  object@misc$gnrh_migration_core <- mig_core
  object@misc$gnrh_secretory_module <- sec
  object@misc$gnrh_stage_parameters <- pars
  if (is.null(object@misc$gnrh)) object@misc$gnrh <- list()
  object@misc$gnrh$stage <- list(modules=dev_modules,migration_core=mig_core,secretory_module=sec,parameters=pars)
  object
}
