# =========================================================
# marker-programs.R
# GnRHcell marker program discovery
# =========================================================

#' Default GnRH developmental marker programs
#'
#' Returns the curated developmental gene modules used by
#' \code{gnrh_marker_programs()} to classify candidate GnRH markers
#' into major developmental programs.
#'
#' The modules represent four key biological stages of the GnRH
#' neuronal lineage:
#'
#' \itemize{
#'   \item \strong{identity}: specification and developmental identity
#'   genes.
#'   \item \strong{migrating}: genes involved in neuronal migration and
#'   axon guidance.
#'   \item \strong{mature}: markers associated with differentiated GnRH
#'   neurons.
#'   \item \strong{secreting}: genes involved in neuropeptide processing
#'   and secretion.
#' }
#'
#' All genes are returned in uppercase to facilitate cross-species
#' comparisons and marker matching.
#'
#' @return
#' A named list containing four character vectors:
#' \describe{
#'   \item{identity}{Developmental identity and specification genes.}
#'   \item{migrating}{Migration and guidance-associated genes.}
#'   \item{mature}{Markers of differentiated GnRH neurons.}
#'   \item{secreting}{Genes involved in secretory function.}
#' }
#'
#' @details
#' These modules were manually curated from the GnRH developmental
#' literature and are intended for exploratory marker discovery rather
#' than strict cell-state annotation.
#'
#' Users may modify the returned gene sets and provide a customized
#' module list to \code{gnrh_marker_programs()}.
#'
#' @examples
#' modules <- gnrh_stage_modules()
#'
#' names(modules)
#'
#' modules$identity
#'
#' length(modules$migrating)
#'
#' custom_modules <- modules
#' custom_modules$identity <- unique(
#'   c(custom_modules$identity, "FEZF2")
#' )
#'
#' @seealso
#' \code{\link{gnrh_marker_programs}}
#'
#' @references
#' Wray S. Development of gonadotropin-releasing hormone-1 neurons.
#' Front Neuroendocrinol. 2010.
#'
#' Schwanzel-Fukuda M, Pfaff DW. Origin of luteinizing hormone-releasing
#' hormone neurons. Nature. 1989.
#'
#' @export
#'
gnrh_stage_modules <- function() {
  list(
    identity = toupper(c(
      "FEZF1", "SOX2", "HES1", "OTX2", "SIX3", "SIX6",
      "ISL1", "DLX1", "DLX2", "DLX5", "DLX6",
      "ARX", "FOXG1", "RAX", "ISL2", "MYT1",
      "KLF7", "ZIC2", "ZIC4", "ZIC5",
      "ZNF483", "EBF3", "MEIS1"
    )),
    migrating = toupper(c(
      "ANOS1", "PROKR2", "PROK2", "DCX", "L1CAM",
      "SEMA3A", "SEMA3C", "SEMA3F",
      "NRP1", "NRP2", "ROBO1", "ROBO2",
      "SLIT1", "UNC5D", "RIPOR2",
      "PTPRO", "PLXNA3", "NFASC", "ADGRV1",
      "CDH22", "CTNNA2", "NALCN", "FREM1",
      "CXCR4", "GDNF", "TMEM131L", "DLX6OS1"
    )),
    mature = toupper(c(
      "GNRH1", "KISS1R", "GNRHR",
      "PCSK1", "PCSK2", "CPE",
      "CHGA", "CHGB", "SCG2", "SCG5",
      "SYP", "RAB3A", "RAB3B", "RAB3C",
      "VGF", "PTPRN", "DOC2B", "CADPS",
      "HCN1", "SCN3A", "SCN9A",
      "KCNMB2", "CACNA1B", "CHRNB4", "CHRNA3"
    )),
    secreting = toupper(c(
      "PCSK1", "PCSK2", "CPE",
      "SYP", "VAMP2", "RAB3A", "RAB3B", "RAB3C",
      "SCG2", "SCG5", "BAIAP3", "PTPRN",
      "TAC3", "TAC1", "KISS1R",
      "DOC2B", "CADPS", "SLC18A1",
      "PCSK1N", "RAB27B", "PLD5"
    ))
  )
}




#' Classify GnRH marker candidates into developmental programs
#'
#' Combines marker tables, dataset-specific marker overlap results, GNRH1
#' co-expression evidence, specificity scores, and curated developmental
#' programs to identify candidate GnRH markers associated with identity,
#' migration, maturation, or secretion.
#'
#' This function is designed to be used downstream of marker discovery and
#' overlap analysis. It reads marker tables generated for multiple datasets,
#' intersects dataset-specific unique markers with GNRH1-coexpressed genes,
#' assigns candidates to curated GnRH developmental programs, and reports
#' confidence levels based on co-expression, specificity, and program
#' membership.
#'
#' @param files Named character vector. Names correspond to dataset names and
#'   values correspond to marker table file names located in
#'   \code{file.path(outdir, "tables")}.
#' @param results A list returned by \code{gene_upset()}, containing at least
#'   \code{results$unique}, a named list of dataset-specific unique genes.
#' @param outdir Character. Output directory containing a \code{tables/}
#'   subdirectory with marker tables.
#' @param modules Named list of developmental gene modules. By default, uses
#'   \code{gnrh_stage_modules()}.
#' @param coexpr_col Character. Name of the column indicating GNRH1
#'   co-expression. Default is \code{"coexpr_flag"}.
#' @param gene_col Character. Name of the gene column in marker tables.
#'   Default is \code{"gene"}.
#' @param min_medium_score Numeric. Minimum specificity score required for
#'   medium confidence. Default is \code{0.25}.
#' @param min_high_score Numeric. Minimum specificity score required for
#'   high confidence. Default is \code{1}.
#' @param write_output Logical. If \code{TRUE}, writes output tables to
#'   \code{file.path(outdir, "tables")}. Default is \code{TRUE}.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{\code{candidate_table}}{Data frame containing all evaluated marker
#'   candidates with dataset, gene, uniqueness, known marker status,
#'   developmental program, GNRH1 co-expression status, specificity score, and
#'   confidence level.}
#'   \item{\code{high_confidence}}{Subset of \code{candidate_table} classified
#'   as high-confidence candidate markers.}
#'   \item{\code{summary}}{Summary table of candidate genes grouped by dataset,
#'   developmental program, and confidence level.}
#'   \item{\code{split_by_dataset_program}}{List of candidate genes split by
#'   dataset and developmental program.}
#' }
#'
#' @details
#' The specificity score is computed as:
#' \deqn{
#' avg\_log2FC \times (pct.1 - pct.2) \times -log10(p\_val\_adj + \epsilon)
#' }
#' where \eqn{\epsilon = 1e-300}. Missing columns are handled safely and
#' replaced with conservative default values.
#'
#' Confidence levels are assigned as follows:
#' \itemize{
#'   \item \code{high}: unique, GNRH1-coexpressed, assigned to a developmental
#'   program, specificity score greater than or equal to
#'   \code{min_high_score}, and not classified as a generic neuronal marker.
#'   \item \code{medium}: unique, GNRH1-coexpressed, specificity score greater
#'   than or equal to \code{min_medium_score}, and not classified as generic.
#'   \item \code{low}: all remaining candidates.
#' }
#'
#' @examples
#' \dontrun{
#' files <- c(
#'   "HuDeCa Nose" = "gnrh_nose_markers.tsv",
#'   "HuDeCa Hypo" = "gnrh_hudeca_hypo_markers.tsv",
#'   "HPSC Wang 2022" = "gnrh_wang_markers.tsv",
#'   "Human HypoMap" = "gnrh_human_hypomap_markers.tsv",
#'   "Mouse HypoMap" = "gnrh_mouse_hypomap_markers.tsv",
#'   "Mouse Amato 2024" = "gnrh_mouse_amato_markers.tsv",
#'   "Mouse POA" = "gnrh_mouse_poa_markers.tsv",
#'   "Mouse MBH" = "gnrh_mouse_mbh_markers.tsv"
#' )
#'
#' marker_programs <- gnrh_marker_programs(
#'   files = files,
#'   results = overlap_results,
#'   outdir = "results/gnrh",
#'   write_output = TRUE
#' )
#'
#' marker_programs$candidate_table
#' marker_programs$high_confidence
#' marker_programs$summary
#' marker_programs$split_by_dataset_program
#' }
#'
#' @seealso
#' \code{\link{gnrh_stage_modules}}, \code{\link{gene_upset}},
#' \code{\link{gnrh_markers}}
#'
#' @export
gnrh_marker_programs <- function(
    files,
    results,
    outdir,
    modules = gnrh_stage_modules(),
    coexpr_col = "coexpr_flag",
    gene_col = "gene",
    min_medium_score = 0.25,
    min_high_score = 1,
    write_output = TRUE
) {

  if (is.null(names(files)) || any(names(files) == "")) {
    stop("'files' must be a named character vector.", call. = FALSE)
  }

  if (!is.list(results) || is.null(results$unique)) {
    stop("'results' must be a list containing results$unique.", call. = FALSE)
  }

  if (!is.list(modules) || is.null(names(modules))) {
    stop("'modules' must be a named list.", call. = FALSE)
  }

  known_gnrh_markers <- toupper(c(
    "GNRH1", "ISL1", "SIX3", "SIX6",
    "DLX1", "DLX2", "DLX5", "DLX6",
    "OTX2", "KISS1R", "TAC1", "TAC2", "TAC3",
    "SEMA3C", "RIPOR2", "UNC5D",
    "PROKR2", "ANOS1", "L1CAM",
    "GAD1", "GAD2"
  ))

  generic_neuronal <- toupper(c(
    "TUBB", "TUBB2A", "TUBB2B", "TUBB3",
    "TUBA1A", "TUBA1B", "EEF1A2",
    "ACTG1", "PKM", "ALDOA",
    "GAPDH", "HSP90AB1", "HSPA5",
    "HSPH1", "MALAT1", "MIAT", "GNAS"
  ))

  candidate_novel <- toupper(c(
    "CADPS", "PTPRO", "RAB3C", "ADGRV1",
    "PLXNA3", "SLC18A1", "RAX", "ISL2",
    "SCG2", "SCG5", "BAIAP3", "PCSK1N",
    "CFAP206", "PDE11A", "RAB27B", "KCTD8",
    "TMEM131L", "GDNF", "HCN1", "SCN3A",
    "SCN9A", "KCNMB2", "CNGB1", "CHRNB4",
    "CHRNA3", "S100Z", "HAP1", "PLD5",
    "FREM1", "MEIS1", "EBF3", "DLX6OS1"
  ))

  classify_known_status <- function(gene) {
    gene <- toupper(gene)

    if (gene %in% known_gnrh_markers) {
      return("known_GnRH_or_developmental")
    }

    if (gene %in% generic_neuronal) {
      return("generic_neuronal")
    }

    if (gene %in% candidate_novel) {
      return("candidate_novel")
    }

    "unknown_or_context_specific"
  }

  safe_col <- function(df, col, default = NA) {
    if (col %in% colnames(df)) {
      return(df[[col]])
    }
    rep(default, nrow(df))
  }

  compute_specificity_score <- function(df) {
    avg_log2FC <- suppressWarnings(
      as.numeric(safe_col(df, "avg_log2FC", 0))
    )
    pct1 <- suppressWarnings(
      as.numeric(safe_col(df, "pct.1", 0))
    )
    pct2 <- suppressWarnings(
      as.numeric(safe_col(df, "pct.2", 0))
    )
    padj <- suppressWarnings(
      as.numeric(safe_col(df, "p_val_adj", 1))
    )

    score <- avg_log2FC * (pct1 - pct2) * -log10(padj + 1e-300)
    score[is.na(score)] <- 0
    score
  }

  assign_programs <- function(genes, modules) {
    genes <- unique(toupper(genes))

    out <- lapply(names(modules), function(program) {
      hits <- intersect(genes, toupper(modules[[program]]))

      if (length(hits) == 0) {
        return(NULL)
      }

      data.frame(
        gene = hits,
        program = program,
        stringsAsFactors = FALSE
      )
    })

    out <- Filter(Negate(is.null), out)

    if (length(out) == 0) {
      return(data.frame(
        gene = character(0),
        program = character(0),
        stringsAsFactors = FALSE
      ))
    }

    do.call(rbind, out)
  }

  empty_candidate_table <- function() {
    data.frame(
      dataset = character(0),
      gene = character(0),
      is_unique = logical(0),
      known_status = character(0),
      program = character(0),
      coexpr_GNRH1 = logical(0),
      specificity_score = numeric(0),
      confidence_level = character(0),
      stringsAsFactors = FALSE
    )
  }

  build_one <- function(dataset) {
    marker_file <- file.path(outdir, "tables", files[[dataset]])

    if (!file.exists(marker_file)) {
      warning("Missing file: ", marker_file, call. = FALSE)
      return(NULL)
    }

    markers <- utils::read.delim(
      marker_file,
      sep = "\t",
      stringsAsFactors = FALSE
    )

    if (!gene_col %in% colnames(markers)) {
      warning("Missing gene column in: ", dataset, call. = FALSE)
      return(NULL)
    }

    markers[[gene_col]] <- toupper(markers[[gene_col]])

    if (!coexpr_col %in% colnames(markers)) {
      warning("Missing co-expression column in: ", dataset, call. = FALSE)
      markers[[coexpr_col]] <- FALSE
    }

    unique_genes <- toupper(results$unique[[dataset]])

    if (is.null(unique_genes)) {
      unique_genes <- character(0)
    }

    keep_genes <- unique(c(unique_genes, "GNRH1"))

    df <- markers[markers[[gene_col]] %in% keep_genes, , drop = FALSE]

    if (nrow(df) == 0) {
      return(NULL)
    }

    colnames(df)[colnames(df) == gene_col] <- "gene"

    df$dataset <- dataset
    df$is_unique <- df$gene %in% unique_genes
    df$coexpr_GNRH1 <- df[[coexpr_col]] == TRUE
    df$known_status <- vapply(
      df$gene,
      classify_known_status,
      character(1)
    )

    df$specificity_score <- compute_specificity_score(df)

    program_df <- assign_programs(df$gene, modules)

    if (nrow(program_df) == 0) {
      df$program <- "unassigned"
      expanded <- df
    } else {
      expanded <- merge(
        df,
        program_df,
        by = "gene",
        all.x = TRUE,
        sort = FALSE
      )

      expanded$program[is.na(expanded$program)] <- "unassigned"
    }

    expanded$confidence_level <- "low"

    is_medium <- (
      expanded$is_unique == TRUE &
        expanded$coexpr_GNRH1 == TRUE &
        expanded$specificity_score >= min_medium_score &
        expanded$known_status != "generic_neuronal"
    )

    is_high <- (
      expanded$is_unique == TRUE &
        expanded$coexpr_GNRH1 == TRUE &
        expanded$program != "unassigned" &
        expanded$specificity_score >= min_high_score &
        expanded$known_status %in% c(
          "candidate_novel",
          "known_GnRH_or_developmental",
          "unknown_or_context_specific"
        )
    )

    expanded$confidence_level[is_medium] <- "medium"
    expanded$confidence_level[is_high] <- "high"

    expanded[, c(
      "dataset",
      "gene",
      "is_unique",
      "known_status",
      "program",
      "coexpr_GNRH1",
      "specificity_score",
      "confidence_level"
    )]
  }

  candidate_list <- lapply(names(files), build_one)
  candidate_list <- Filter(Negate(is.null), candidate_list)

  if (length(candidate_list) == 0) {
    candidate_table <- empty_candidate_table()
  } else {
    candidate_table <- do.call(rbind, candidate_list)
  }

  if (nrow(candidate_table) > 0) {
    candidate_table <- candidate_table[
      order(
        candidate_table$dataset,
        candidate_table$program,
        candidate_table$confidence_level,
        -candidate_table$specificity_score
      ),
      ,
      drop = FALSE
    ]

    rownames(candidate_table) <- NULL
  }

  high_confidence <- candidate_table[
    candidate_table$confidence_level == "high",
    ,
    drop = FALSE
  ]

  candidate_for_summary <- candidate_table[
    candidate_table$confidence_level %in% c("high", "medium") &
      candidate_table$is_unique == TRUE &
      candidate_table$coexpr_GNRH1 == TRUE,
    ,
    drop = FALSE
  ]

  if (nrow(candidate_for_summary) == 0) {
    summary_table <- data.frame(
      dataset = character(0),
      program = character(0),
      confidence_level = character(0),
      gene = character(0),
      n_genes = integer(0),
      stringsAsFactors = FALSE
    )
    split_table <- list()
  } else {
    summary_table <- stats::aggregate(
      gene ~ dataset + program + confidence_level,
      data = candidate_for_summary,
      FUN = function(x) paste(sort(unique(x)), collapse = ", ")
    )

    summary_table$n_genes <- vapply(
      strsplit(summary_table$gene, ", "),
      length,
      integer(1)
    )

    summary_table <- summary_table[
      order(summary_table$dataset, summary_table$program),
      ,
      drop = FALSE
    ]

    split_table <- split(
      candidate_for_summary$gene,
      list(candidate_for_summary$dataset, candidate_for_summary$program),
      drop = TRUE
    )
  }

  if (write_output) {
    out_tables <- file.path(outdir, "tables")

    if (!dir.exists(out_tables)) {
      dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
    }

    utils::write.table(
      candidate_table,
      file = file.path(out_tables, "gnrh_candidate_marker_table.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      high_confidence,
      file = file.path(out_tables, "gnrh_high_confidence_candidate_markers.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      summary_table,
      file = file.path(out_tables, "gnrh_candidate_marker_summary_by_program.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  list(
    candidate_table = candidate_table,
    high_confidence = high_confidence,
    summary = summary_table,
    split_by_dataset_program = split_table
  )
}
