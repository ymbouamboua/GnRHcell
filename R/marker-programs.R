# ============================================================================= #
# Marker-program helpers
# ============================================================================= #
#' Safely extract table column
#' @keywords internal
#' @noRd
.marker_safe_col <- function(x,column,default=NA) if (column %in% colnames(x)) x[[column]] else rep(default,nrow(x))
#' Convert arbitrary values to logical
#' @keywords internal
#' @noRd
.marker_as_logical <- function(x) {
  out <- if (is.logical(x)) x else if (is.numeric(x)) is.finite(x) & x>0 else tolower(trimws(as.character(x))) %in% c("true","t","yes","y","1","positive","pos")
  out[is.na(out)] <- FALSE
  out
}


#' Compute marker evidence score
#' @keywords internal
#' @noRd
.marker_evidence_score <- function(df) {
  for (nm in c("final_score","marker_score","specificity_score")) {
    if (nm %in% colnames(df)) {
      x <- suppressWarnings(as.numeric(df[[nm]]))
      x[!is.finite(x)] <- 0
      if (any(x>0)) return(x)
    }
  }
  fc <- suppressWarnings(as.numeric(.marker_safe_col(df,"avg_log2FC",0)))
  pct1 <- suppressWarnings(as.numeric(.marker_safe_col(df,"pct.1",0)))
  pct2 <- suppressWarnings(as.numeric(.marker_safe_col(df,"pct.2",0)))
  padj <- suppressWarnings(as.numeric(.marker_safe_col(df,"p_val_adj",1)))
  fc[!is.finite(fc)] <- 0
  pct1[!is.finite(pct1)] <- 0
  pct2[!is.finite(pct2)] <- 0
  padj[!is.finite(padj)] <- 1
  sig <- pmin(-log10(pmax(padj,.Machine$double.xmin)),50)
  out <- pmax(fc,0)*pmax(pct1-pct2,0)*sqrt(pmax(pct1,0))*sig
  out[!is.finite(out)] <- 0
  out
}


#' Assign genes to curated programs
#' @keywords internal
#' @noRd
.assign_marker_programs <- function(genes,developmental_modules,secretory_module) {
  genes <- unique(toupper(as.character(genes)))
  modules <- c(developmental_modules,list(secretory=secretory_module))
  dimension <- c(stats::setNames(rep("developmental",length(developmental_modules)),names(developmental_modules)),secretory="functional")
  rows <- lapply(names(modules),function(program) {
    hit <- intersect(genes,toupper(modules[[program]]))
    if (!length(hit)) return(NULL)
    data.frame(gene=hit,program=program,program_dimension=dimension[[program]],stringsAsFactors=FALSE)
  })
  rows <- Filter(Negate(is.null),rows)
  if (!length(rows)) return(data.frame(gene=character(),program=character(),program_dimension=character(),stringsAsFactors=FALSE))
  do.call(rbind,rows)
}


# ============================================================================= #
# GnRH marker programs
# ============================================================================= #
#' Classify GnRH marker candidates into biological programs
#' @export
gnrh_marker_programs <- function(files,results,outdir,modules=gnrh_stage_modules(),secretory_module=gnrh_secretory_marker_module(),coexpr_col="coexpr_flag",gene_col="gene",min_medium_score=0.25,min_high_score=1,min_high_datasets=2L,min_medium_datasets=1L,include_gnrh1=TRUE,write_output=TRUE) {
  .validate_named_files(files,"files")
  if (!is.list(results)) stop("`results` must be returned by `gnrh_gene_upset()`.",call.=FALSE)
  if (!is.list(modules) || !length(modules) || is.null(names(modules))) stop("`modules` must be a named list.",call.=FALSE)
  min_high_datasets <- as.integer(min_high_datasets)
  min_medium_datasets <- as.integer(min_medium_datasets)
  if (min_high_datasets<1L) stop("`min_high_datasets` must be >= 1.",call.=FALSE)
  if (min_medium_datasets<1L) stop("`min_medium_datasets` must be >= 1.",call.=FALSE)
  miss <- files[!file.exists(unname(files))]
  if (length(miss)) stop("Missing marker files: ",paste(unname(miss),collapse=", "),call.=FALSE)
  refs <- gnrh_stage_gene_references()
  known <- unique(toupper(c(refs$gene,"GNRH1")))
  generic <- toupper(c("TUBB","TUBB2A","TUBB2B","TUBB3","TUBA1A","TUBA1B","EEF1A2","ACTG1","PKM","ALDOA","GAPDH","HSP90AB1","HSPA5","HSPH1","MALAT1","MIAT","GNAS","SNAP25","SYT1","RBFOX3"))
  overlap_sets <- results$gene_sets %||% NULL
  if (!is.null(overlap_sets)) overlap_sets <- lapply(overlap_sets,function(x) unique(toupper(as.character(x))))
  unique_sets <- lapply(results$unique %||% list(),function(x) unique(toupper(as.character(x))))
  build_one <- function(dataset) {
    f <- files[[dataset]]
    x <- utils::read.delim(f,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
    if (!gene_col %in% colnames(x)) stop("Missing gene column `",gene_col,"` in ",dataset,".",call.=FALSE)
    x[[gene_col]] <- toupper(trimws(as.character(x[[gene_col]])))
    x <- x[!is.na(x[[gene_col]]) & nzchar(x[[gene_col]]),,drop=FALSE]
    if (!nrow(x)) return(NULL)
    names(x)[names(x)==gene_col] <- "gene"
    x$dataset <- dataset
    x$is_unique <- x$gene %in% (unique_sets[[dataset]] %||% character())
    x$coexpr_GNRH1 <- if (!is.null(coexpr_col) && coexpr_col %in% names(x)) .marker_as_logical(x[[coexpr_col]]) else if ("coexpr_pct" %in% names(x)) {
      z <- suppressWarnings(as.numeric(x$coexpr_pct)); is.finite(z) & z>0
    } else NA
    x$marker_evidence_score <- .marker_evidence_score(x)
    x
  }
  marker_list <- Filter(Negate(is.null),lapply(names(files),build_one))
  if (!length(marker_list)) stop("No marker tables could be loaded.",call.=FALSE)
  all_markers <- do.call(rbind,marker_list)
  rownames(all_markers) <- NULL
  genes <- unique(all_markers$gene)
  support <- vapply(genes,function(g) if (!is.null(overlap_sets)) sum(vapply(overlap_sets,function(x) g %in% x,logical(1))) else length(unique(all_markers$dataset[all_markers$gene==g])),integer(1))
  names(support) <- genes
  all_markers$n_datasets <- support[all_markers$gene]
  all_markers$dataset_fraction <- all_markers$n_datasets/length(files)
  all_markers$reproducible <- all_markers$n_datasets>=min_high_datasets
  all_markers$context_specific <- all_markers$n_datasets==1L
  all_markers$known_status <- ifelse(all_markers$gene %in% generic,"generic_neuronal",ifelse(all_markers$gene %in% known,"known_GnRH_or_developmental",ifelse(all_markers$context_specific,"context_specific_candidate","candidate_marker")))
  program_df <- .assign_marker_programs(all_markers$gene,modules,secretory_module)
  if (nrow(program_df)) candidate_table <- merge(all_markers,program_df,by="gene",all.x=TRUE,sort=FALSE) else {
    candidate_table <- all_markers
    candidate_table$program <- NA_character_
    candidate_table$program_dimension <- NA_character_
  }
  candidate_table$program[is.na(candidate_table$program)] <- "unassigned"
  candidate_table$program_dimension[is.na(candidate_table$program_dimension)] <- "unassigned"
  candidate_table$is_reference_gene <- candidate_table$gene=="GNRH1"
  if (!isTRUE(include_gnrh1)) candidate_table <- candidate_table[!candidate_table$is_reference_gene,,drop=FALSE]
  score <- candidate_table$marker_evidence_score
  not_generic <- candidate_table$known_status!="generic_neuronal"
  assigned <- candidate_table$program!="unassigned"
  high <- not_generic & assigned & candidate_table$n_datasets>=min_high_datasets & score>=min_high_score
  medium <- not_generic & candidate_table$n_datasets>=min_medium_datasets & score>=min_medium_score
  candidate_table$confidence_level <- "low"
  candidate_table$confidence_level[medium] <- "medium"
  candidate_table$confidence_level[high] <- "high"
  candidate_table$confidence_level <- factor(candidate_table$confidence_level,levels=c("low","medium","high"),ordered=TRUE)
  candidate_table$integrated_score <- score*(0.5+candidate_table$dataset_fraction)*ifelse(assigned,1.25,1)*ifelse(candidate_table$coexpr_GNRH1 %in% TRUE,1.10,1)*ifelse(not_generic,1,0.25)
  candidate_table <- candidate_table[order(-as.integer(candidate_table$confidence_level),-candidate_table$n_datasets,-candidate_table$integrated_score,candidate_table$gene),,drop=FALSE]
  rownames(candidate_table) <- NULL
  high_confidence <- candidate_table[candidate_table$confidence_level=="high" & !candidate_table$is_reference_gene,,drop=FALSE]
  context_specific <- candidate_table[candidate_table$context_specific & candidate_table$known_status!="generic_neuronal" & !candidate_table$is_reference_gene,,drop=FALSE]
  gene_summary <- stats::aggregate(cbind(n_datasets,marker_evidence_score,integrated_score)~gene,candidate_table,max,na.rm=TRUE)
  programs <- stats::aggregate(program~gene,candidate_table,function(x) paste(sort(unique(x[x!="unassigned"])),collapse=", "))
  gene_summary <- merge(gene_summary,programs,by="gene",all.x=TRUE,sort=FALSE)
  gene_summary <- gene_summary[order(-gene_summary$n_datasets,-gene_summary$integrated_score),,drop=FALSE]
  ss <- candidate_table[candidate_table$confidence_level %in% c("medium","high") & !candidate_table$is_reference_gene,,drop=FALSE]
  if (nrow(ss)) {
    summary_table <- ss |>
      dplyr::group_by(.data$program_dimension,.data$program,.data$confidence_level) |>
      dplyr::summarise(n_genes=dplyr::n_distinct(.data$gene),genes=paste(sort(unique(.data$gene)),collapse=", "),.groups="drop")
    split_table <- split(ss$gene,list(ss$program_dimension,ss$program),drop=TRUE)
    split_table <- lapply(split_table,function(x) sort(unique(x)))
  } else {
    summary_table <- data.frame(program_dimension=character(),program=character(),confidence_level=character(),n_genes=integer(),genes=character(),stringsAsFactors=FALSE)
    split_table <- list()
  }
  dataset_summary <- candidate_table |>
    dplyr::group_by(.data$dataset) |>
    dplyr::summarise(
      n_markers=dplyr::n_distinct(.data$gene),
      n_high=dplyr::n_distinct(.data$gene[.data$confidence_level=="high"]),
      n_medium=dplyr::n_distinct(.data$gene[.data$confidence_level=="medium"]),
      n_context_specific=dplyr::n_distinct(.data$gene[.data$context_specific]),
      .groups="drop"
    )
  if (isTRUE(write_output)) {
    out_tables <- file.path(outdir,"tables")
    dir.create(out_tables,recursive=TRUE,showWarnings=FALSE)
    write_tsv <- function(x,f) utils::write.table(x,file.path(out_tables,f),sep="\t",quote=FALSE,row.names=FALSE)
    write_tsv(candidate_table,"gnrh_candidate_marker_table.tsv")
    write_tsv(high_confidence,"gnrh_high_confidence_candidate_markers.tsv")
    write_tsv(context_specific,"gnrh_context_specific_candidate_markers.tsv")
    write_tsv(gene_summary,"gnrh_marker_cross_dataset_summary.tsv")
    write_tsv(summary_table,"gnrh_candidate_marker_summary_by_program.tsv")
    write_tsv(dataset_summary,"gnrh_candidate_marker_summary_by_dataset.tsv")
  }
  list(
    candidate_table=candidate_table,
    high_confidence=high_confidence,
    context_specific=context_specific,
    gene_summary=gene_summary,
    summary=summary_table,
    dataset_summary=dataset_summary,
    split_by_program=split_table,
    split_by_dataset_program=split_table,
    modules=list(developmental=modules,secretory=secretory_module),
    parameters=list(min_medium_score=min_medium_score,min_high_score=min_high_score,min_medium_datasets=min_medium_datasets,min_high_datasets=min_high_datasets)
  )
}


# ============================================================================= #
# Curated programs
# ============================================================================= #
#' Default GnRH developmental marker programs
#' @param include_secretory Include secretory module for backward compatibility.
#' @return Named list of developmental gene programs.
#' @export
gnrh_stage_modules <- function(include_secretory=FALSE) {
  modules <- list(
    early=c("FEZF1","OTX2","SIX3","SIX6","DLX1","DLX2","DLX5","DLX6"),
    migrating=c("ANOS1","PROKR2","PROK2","NSMF","DCX","L1CAM","SEMA3A",
                "SEMA3C","SEMA3F","NRP1","NRP2","ROBO1","ROBO2","ROBO3",
                "SLIT1","UNC5D","RIPOR2","PTPRO","PLXNA3","NFASC","ADGRV1",
                "CDH22","CTNNA2","FREM1","CXCR4","GDNF","TMEM131L","DLX6OS1"),
    mature=c("KISS1R","GNRHR","DOC2B","PTPRN","BAIAP3","SCG2","SCG5","VGF",
             "HCN1","NALCN","SCN3A","SCN9A","KCNMB2","CACNA1B","CHRNB4","CHRNA3")
  )
  modules <- lapply(modules,function(x) unique(toupper(x)))
  if (isTRUE(include_secretory)) modules$secretory <- gnrh_secretory_marker_module()
  modules
}

#' GnRH secretory marker program
#' @return Character vector of secretory genes.
#' @export
gnrh_secretory_marker_module <- function() {
  unique(toupper(c("PCSK1","PCSK2","CPE","CHGA","CHGB","SCG2",
                   "SCG5","VGF","SYP","VAMP2","RAB3A","RAB3B",
                   "RAB3C","BAIAP3","PTPRN","DOC2B","CADPS","PCSK1N","RAB27B")))
}


#' GnRH developmental and functional marker gene references
#' @return Gene-level reference table.
#' @export
gnrh_stage_gene_references <- function() {
  modules <- c(gnrh_stage_modules(),list(secretory=gnrh_secretory_marker_module()))
  dimension <- c(early="developmental",migrating="developmental",mature="developmental",secretory="functional")
  evidence <- c(early="early_specification",migrating="migration_guidance",mature="neuroendocrine_maturation",secretory="secretory_machinery")
  refs <- c(
    early="Wray 2010; Stevenson et al. 2013; Cho et al. 2019",
    migrating="Schwanzel-Fukuda & Pfaff 1989; Wray 2010; Cho et al. 2019; Taroc et al. 2019",
    mature="Wray 2010; Stevenson et al. 2013",
    secretory="Wray 2010; general neuroendocrine secretory machinery"
  )
  out <- do.call(rbind,lapply(names(modules),function(program) {
    data.frame(gene=modules[[program]],program=program,dimension=dimension[[program]],evidence_class=evidence[[program]],references=refs[[program]],stringsAsFactors=FALSE)
  }))
  rownames(out) <- NULL
  unique(out)
}
