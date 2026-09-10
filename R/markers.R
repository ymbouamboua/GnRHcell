# ============================================================================= #
# GnRHcell marker discovery utilities
# ============================================================================= #

# ----------------------------------------------------------------------------- #
# Internal helpers
# ----------------------------------------------------------------------------- #

#' Fast row-wise correlation
#' @keywords internal
#' @noRd
.fast_row_cor <- function(x,y,method=c("spearman","pearson")) {
  method <- match.arg(method)
  y <- as.numeric(y)
  ok <- is.finite(y)
  x <- x[,ok,drop=FALSE]
  y <- y[ok]
  if (length(y)<3L || !nrow(x)) return(rep(0,nrow(x)))
  if (method=="pearson") {
    n <- length(y)
    sx <- Matrix::rowSums(x)
    sxx <- Matrix::rowSums(x^2)
    sy <- sum(y)
    syy <- sum(y^2)
    sxy <- as.numeric(x%*%y)
    num <- sxy-sx*sy/n
    den <- sqrt(pmax(sxx-sx^2/n,0)*pmax(syy-sy^2/n,0))
    out <- num/den
    out[!is.finite(out)] <- 0
    return(out)
  }
  yr <- rank(y,ties.method="average")
  xr <- t(vapply(seq_len(nrow(x)),function(i) rank(as.numeric(x[i,]),ties.method="average"),numeric(ncol(x))))
  yr <- yr-mean(yr)
  xr <- xr-rowMeans(xr)
  out <- as.numeric(xr%*%yr)/sqrt(rowSums(xr^2)*sum(yr^2))
  out[!is.finite(out)] <- 0
  out
}

#' Resolve GnRH gene
#' @keywords internal
#' @noRd
.get_gnrh_gene <- function(genes) {
  hit <- .match_genes(c("GNRH1","Gnrh1","LHRH"),genes)
  if (!length(hit)) return(NA_character_)
  hit[[1L]]
}

#' Resolve comparison cells
#' @keywords internal
#' @noRd
.get_cells <- function(object,group.by,ident.1,ident.2=NULL) {
  md <- object[[]]
  if (!group.by %in% colnames(md)) stop("Metadata column `",group.by,"` not found.",call.=FALSE)
  grp <- md[[group.by]]
  names(grp) <- rownames(md)
  cells_1 <- names(grp)[!is.na(grp) & grp %in% ident.1]
  cells_2 <- if (is.null(ident.2)) names(grp)[!is.na(grp) & !grp %in% ident.1] else names(grp)[!is.na(grp) & grp %in% ident.2]
  list(cells_1=cells_1,cells_2=cells_2,groups=grp)
}

#' Compute marker specificity
#' @keywords internal
#' @noRd
.compute_specificity <- function(expr,genes,cells_1,cells_2) {
  genes <- intersect(genes,rownames(expr))
  if (!length(genes)) return(data.frame())
  m1 <- expr[genes,cells_1,drop=FALSE]
  m2 <- expr[genes,cells_2,drop=FALSE]
  pct1 <- Matrix::rowMeans(m1>0)
  pct2 <- Matrix::rowMeans(m2>0)
  data.frame(
    gene=genes,
    pct1=as.numeric(pct1),
    pct2=as.numeric(pct2),
    specificity=as.numeric(pct1-pct2),
    stringsAsFactors=FALSE
  )
}

#' Compute gene-wise coexpression with GnRH
#' @keywords internal
#' @noRd
.compute_coexpr <- function(expr,genes,gnrh_gene,cells,method="spearman") {
  genes <- intersect(genes,rownames(expr))
  cells <- intersect(cells,colnames(expr))
  if (!length(genes) || length(cells)<3L) return(setNames(numeric(0),character(0)))
  x <- expr[genes,cells,drop=FALSE]
  y <- as.numeric(expr[gnrh_gene,cells,drop=TRUE])
  out <- .fast_row_cor(x,y,method=method)
  setNames(out,genes)
}

#' Compute GnRH program association
#' @keywords internal
#' @noRd
.compute_program_association <- function(object,expr,genes,cells,program_col,method="spearman") {
  genes <- intersect(genes,rownames(expr))
  cells <- intersect(cells,colnames(expr))
  if (!length(genes) || length(cells)<3L) return(setNames(numeric(0),character(0)))
  meta <- object[[]]
  cells <- intersect(cells,rownames(meta))
  if (length(cells)<3L) return(setNames(rep(NA_real_,length(genes)),genes))
  y <- suppressWarnings(as.numeric(meta[cells,program_col]))
  x <- expr[genes,cells,drop=FALSE]
  out <- .fast_row_cor(x,y,method=method)
  setNames(out,genes)
}

#' Compute binary GnRH co-detection enrichment
#' @keywords internal
#' @noRd
.compute_codetect <- function(expr,genes,gnrh_gene,cells,pseudocount=0.5) {
  genes <- intersect(genes,rownames(expr))
  cells <- intersect(cells,colnames(expr))
  if (!length(genes) || !length(cells)) {
    return(data.frame(gene=character(),odds_ratio=numeric(),log2_or=numeric(),p_value=numeric(),p_adj=numeric()))
  }
  x <- expr[genes,cells,drop=FALSE]>0
  g <- as.numeric(expr[gnrh_gene,cells,drop=TRUE])>0
  n1 <- sum(g)
  n0 <- length(g)-n1
  if (!n1 || !n0) {
    return(data.frame(gene=genes,odds_ratio=NA_real_,log2_or=NA_real_,p_value=NA_real_,p_adj=NA_real_,row.names=NULL))
  }
  a <- as.numeric(Matrix::rowSums(x[,g,drop=FALSE]))
  b <- as.numeric(Matrix::rowSums(x[,!g,drop=FALSE]))
  c <- n1-a
  d <- n0-b
  aa <- a+pseudocount
  bb <- b+pseudocount
  cc <- c+pseudocount
  dd <- d+pseudocount
  odds_ratio <- (aa*dd)/(bb*cc)
  log_or <- log(odds_ratio)
  se <- sqrt(1/aa+1/bb+1/cc+1/dd)
  z <- log_or/se
  p_value <- stats::pnorm(z,lower.tail=FALSE)
  data.frame(
    gene=genes,
    odds_ratio=odds_ratio,
    log2_or=log_or/log(2),
    p_value=p_value,
    p_adj=stats::p.adjust(p_value,method="BH"),
    row.names=NULL,
    stringsAsFactors=FALSE
  )
}

#' Compute GnRH marker ranking score
#' @keywords internal
#' @noRd
.compute_marker_score <- function(df) {
  significance <- pmin(-log10(pmax(df$p_val_adj,.Machine$double.xmin)),50)
  method_weight <- log1p(pmax(df$method_support,1))
  df$marker_score <- pmax(df$avg_log2FC,0)*pmax(df$specificity,0)*sqrt(pmax(df$`pct.1`,0))*significance*method_weight
  df
}

#' Aggregate differential-expression results
#' @keywords internal
#' @noRd
.aggregate_de <- function(res) {
  if (!nrow(res)) return(data.frame())
  out <- lapply(split(res,res$gene),function(x) {
    fc_col <- intersect(c("avg_log2FC","avg_logFC"),colnames(x))
    if (!length(fc_col)) stop("Differential-expression results contain no logFC column.",call.=FALSE)
    fc <- x[[fc_col[[1L]]]]
    padj <- x$p_val_adj
    data.frame(
      gene=x$gene[[1L]],
      avg_log2FC=mean(fc,na.rm=TRUE),
      p_val_adj=if (any(is.finite(padj))) min(padj,na.rm=TRUE) else NA_real_,
      method_support=length(unique(x$method)),
      stringsAsFactors=FALSE
    )
  })
  out <- do.call(rbind,out)
  rownames(out) <- NULL
  out
}


#' Compute phenotype association
#' @keywords internal
#' @noRd
.compute_phenotype_association <- function(expr,genes,cells_1,cells_2,pseudocount=0.5) {
  genes <- intersect(genes,rownames(expr))
  cells_1 <- intersect(cells_1,colnames(expr))
  cells_2 <- intersect(cells_2,colnames(expr))
  if (!length(genes) || !length(cells_1) || !length(cells_2)) {
    return(data.frame(
      gene=character(),
      phenotype_or=numeric(),
      phenotype_log2or=numeric(),
      phenotype_delta=numeric(),
      phenotype_p=numeric(),
      phenotype_fdr=numeric()
    ))
  }
  x1 <- expr[genes,cells_1,drop=FALSE]>0
  x2 <- expr[genes,cells_2,drop=FALSE]>0
  n1 <- length(cells_1)
  n2 <- length(cells_2)
  a <- as.numeric(Matrix::rowSums(x1))
  c <- as.numeric(Matrix::rowSums(x2))
  b <- n1-a
  d <- n2-c
  aa <- a+pseudocount
  bb <- b+pseudocount
  cc <- c+pseudocount
  dd <- d+pseudocount
  phenotype_or <- (aa*dd)/(bb*cc)
  log_or <- log(phenotype_or)
  se <- sqrt(1/aa+1/bb+1/cc+1/dd)
  z <- log_or/se
  phenotype_p <- 2*stats::pnorm(abs(z),lower.tail=FALSE)
  phenotype_delta <- a/n1-c/n2
  data.frame(
    gene=genes,
    phenotype_or=phenotype_or,
    phenotype_log2or=log_or/log(2),
    phenotype_delta=phenotype_delta,
    phenotype_p=phenotype_p,
    phenotype_fdr=stats::p.adjust(phenotype_p,method="BH"),
    row.names=NULL,
    stringsAsFactors=FALSE
  )
}


# ------------------------------------------------------------------------- #
# GnRH marker discovery
# ------------------------------------------------------------------------- #

#' Discover and rank GnRH-associated marker genes
#'
#' Identifies genes enriched in a target GnRH population relative to a
#' comparison population and ranks them using differential expression,
#' detection specificity, GnRH co-expression, GnRH co-detection, and optional
#' correlation with a GnRH transcriptional program.
#'
#' @param object A Seurat object.
#' @param group.by Character scalar giving the metadata column defining target
#'   and comparison groups. Default is \code{"gnrh_status"}.
#' @param ident.1 Character vector defining the target identity or identities.
#'   Default is \code{"pos"}.
#' @param ident.2 Optional character vector defining the comparison identity or
#'   identities. If \code{NULL}, all cells not belonging to \code{ident.1} are
#'   used as the comparison group.
#' @param assay Optional character scalar giving the assay used for differential
#'   expression and association calculations. If \code{NULL}, the default
#'   Seurat assay is used.
#' @param layer Character scalar giving the expression layer used for
#'   specificity and association calculations. Default is \code{"data"}.
#' @param methods Character vector of differential-expression methods passed to
#'   \code{Seurat::FindMarkers()}. Default is \code{"wilcox"}.
#' @param coexpr.method Character scalar giving the correlation method used for
#'   quantitative GnRH co-expression and program association. One of
#'   \code{"pearson"} or \code{"spearman"}. Default is \code{"pearson"}.
#' @param program_col Optional character scalar giving a metadata column
#'   containing a GnRH transcriptional-program score. If \code{NULL} and
#'   \code{"gnrh_identity_score"} exists, that column is used automatically.
#' @param association_mode Character scalar specifying how GnRH association is
#'   incorporated into marker ranking. One of \code{"combined"},
#'   \code{"auto"}, \code{"coexpression"}, \code{"codetection"}, or
#'   \code{"program"}. Default is \code{"combined"}.
#' @param min_pct Minimum detection fraction in the target population required
#'   for candidate retention. Default is \code{0.005}.
#' @param min_fc Minimum average log2 fold change required for candidate
#'   retention. Default is \code{0.10}.
#' @param max_padj Maximum adjusted P value allowed for candidate markers.
#'   Default is \code{0.05}.
#' @param coexpr_min Optional minimum correlation with GnRH expression required
#'   for marker retention.
#' @param codetect_min_or Optional minimum odds ratio for binary co-detection
#'   with the GnRH transcript.
#' @param codetect_max_fdr Optional maximum FDR for GnRH co-detection
#'   enrichment.
#' @param program_min Optional minimum correlation with the GnRH transcriptional
#'   program.
#' @param min_detect Minimum number of cells in which a gene must be detected
#'   before marker evaluation. Default is \code{2L}.
#' @param exclude_gnrh Logical. Whether to exclude the GnRH reference gene
#'   itself from the returned marker table. Default is \code{TRUE}.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers()}.
#'
#' @return A data frame of ranked GnRH-associated marker genes. The returned
#'   table may contain:
#' \describe{
#'   \item{\code{gene}}{Gene symbol.}
#'   \item{\code{avg_log2FC}}{Average log2 fold change across DE methods.}
#'   \item{\code{p_val_adj}}{Best adjusted P value across DE methods.}
#'   \item{\code{method_support}}{Number of DE methods supporting the gene.}
#'   \item{\code{pct.1}, \code{pct.2}}{Detection fractions in target and
#'     comparison cells.}
#'   \item{\code{specificity}}{Detection difference
#'     \code{pct.1 - pct.2}.}
#'   \item{\code{phenotype_or}}{Detection odds ratio for target versus
#'     comparison cells.}
#'   \item{\code{phenotype_log2or}}{Log2 phenotype odds ratio.}
#'   \item{\code{phenotype_delta}}{Difference in detection frequency between
#'     target and comparison cells.}
#'   \item{\code{phenotype_p}, \code{phenotype_fdr}}{Phenotype-association
#'     significance statistics.}
#'   \item{\code{coexpr_cor}}{Correlation between candidate-gene and GnRH
#'     expression within target cells.}
#'   \item{\code{codetect_or}}{Odds ratio for binary co-detection with the GnRH
#'     transcript.}
#'   \item{\code{codetect_log2or}}{Log2 co-detection odds ratio.}
#'   \item{\code{codetect_p}, \code{codetect_fdr}}{Co-detection significance
#'     statistics.}
#'   \item{\code{program_cor}}{Correlation with the GnRH transcriptional
#'     program when available.}
#'   \item{\code{association_score}}{Scaled GnRH-association score used to
#'     refine marker ranking.}
#'   \item{\code{marker_score}}{Final marker-ranking score.}
#' }
#'
#' @details
#' Differential-expression results from all requested \code{methods} are first
#' aggregated by gene. Candidate genes are then filtered using minimum
#' detection, fold change, target-cell detection frequency, and adjusted
#' P value before the more computationally intensive association metrics are
#' calculated.
#'
#' Quantitative co-expression measures correlation between each candidate gene
#' and the GnRH reference transcript within target cells. Binary co-detection
#' evaluates whether candidate-gene detection is enriched in cells in which the
#' GnRH transcript is detected. When \code{program_col} is available,
#' candidate-gene expression can additionally be correlated with the GnRH
#' transcriptional-program score.
#'
#' With \code{association_mode = "combined"}, the ranking association score is
#' based on scaled co-detection and program association when a program score is
#' available, otherwise on co-detection alone. Quantitative co-expression is
#' still calculated so that it remains available for downstream diagnostics and
#' visualization.
#'
#' With \code{association_mode = "auto"}, the function evaluates the information
#' content of the GnRH transcript. Low-information expression profiles favor
#' co-detection or combined scoring, whereas sufficiently variable GnRH
#' expression favors quantitative co-expression.
#'
#' Pearson correlation is the default because it permits efficient matrix-based
#' calculation on large single-cell datasets. Spearman correlation is supported
#' but requires rank transformation and is therefore more computationally
#' expensive.
#'
#' The final marker score combines positive fold change, target-cell
#' specificity, target detection frequency, statistical significance, method
#' support, and the selected GnRH association score.
#'
#' @examples
#' \dontrun{
#' markers <- gnrh_markers(
#'   object = wang,
#'   group.by = "gnrh_status",
#'   ident.1 = "pos"
#' )
#'
#' head(markers)
#'
#' markers <- gnrh_markers(
#'   object = wang,
#'   group.by = "gnrh_status",
#'   ident.1 = "pos",
#'   association_mode = "combined",
#'   coexpr.method = "pearson"
#' )
#'
#' strict <- gnrh_markers(
#'   object = wang,
#'   ident.1 = "pos",
#'   min_fc = 0.5,
#'   min_pct = 0.1,
#'   codetect_min_or = 2
#' )
#' }
#'
#' @seealso
#' \code{\link{find_gnrh_genes}},
#' \code{\link{find_gnrh_stage_markers}},
#' \code{\link{run_gnrh}}
#'
#' @family marker discovery
#' @export
gnrh_markers <- function(
    object,
    group.by="gnrh_status",
    ident.1="pos",
    ident.2=NULL,
    assay=NULL,
    layer="data",
    methods="wilcox",
    coexpr.method="pearson",
    program_col=NULL,
    association_mode=c("combined","auto","coexpression","codetection","program"),
    min_pct=0.005,
    min_fc=0.10,
    max_padj=0.05,
    coexpr_min=NULL,
    codetect_min_or=NULL,
    codetect_max_fdr=NULL,
    program_min=NULL,
    min_detect=2L,
    exclude_gnrh=TRUE,
    verbose=TRUE,
    ...
) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  association_mode <- match.arg(association_mode)
  coexpr.method <- match.arg(coexpr.method,c("spearman","pearson"))
  assay <- assay %||% Seurat::DefaultAssay(object)
  log <- .msg(verbose)
  log("GnRH marker discovery",type="header")
  expr <- .get_expr(object,assay=assay,layer=layer)
  genes <- rownames(expr)
  gnrh_gene <- .get_gnrh_gene(genes)
  if (is.na(gnrh_gene)) stop("GNRH1/Gnrh1 not found.",call.=FALSE)
  grp <- .get_cells(object=object,group.by=group.by,ident.1=ident.1,ident.2=ident.2)
  if (length(grp$cells_1)<5L) stop("Too few target cells: ",length(grp$cells_1),".",call.=FALSE)
  if (length(grp$cells_2)<5L) stop("Too few comparison cells: ",length(grp$cells_2),".",call.=FALSE)
  SeuratObject::Idents(object) <- grp$groups
  # ------------------------------------------------------------------------- #
  # Differential expression
  # ------------------------------------------------------------------------- #
  all_res <- lapply(methods,function(m) {
    log("Running DE method:",m,type="step")
    mk <- Seurat::FindMarkers(
      object=object,ident.1=ident.1,ident.2=ident.2,assay=assay,
      test.use=m,verbose=FALSE,...
    )
    if (!nrow(mk)) return(NULL)
    mk$gene <- rownames(mk)
    mk$method <- m
    mk
  })
  all_res <- all_res[!vapply(all_res,is.null,logical(1))]
  if (!length(all_res)) stop("No markers detected by Seurat::FindMarkers().",call.=FALSE)
  de <- .aggregate_de(do.call(rbind,all_res))
  # ------------------------------------------------------------------------- #
  # Minimum detection
  # ------------------------------------------------------------------------- #
  detected <- Matrix::rowSums(expr>0)
  keep_genes <- names(detected)[detected>=min_detect]
  de <- de[de$gene %in% keep_genes,,drop=FALSE]
  if (!nrow(de)) stop("No markers passed `min_detect` filtering.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Specificity
  # ------------------------------------------------------------------------- #
  spec <- .compute_specificity(expr=expr,genes=de$gene,cells_1=grp$cells_1,cells_2=grp$cells_2)
  idx <- match(de$gene,spec$gene)
  de$`pct.1` <- spec$pct1[idx]
  de$`pct.2` <- spec$pct2[idx]
  de$specificity <- spec$specificity[idx]
  de$pct1 <- de$`pct.1`
  de$pct2 <- de$`pct.2`
  de$spec <- de$specificity
  # ------------------------------------------------------------------------- #
  # Early candidate filtering
  # ------------------------------------------------------------------------- #
  keep <- is.finite(de$avg_log2FC) &
    de$avg_log2FC>=min_fc &
    is.finite(de$`pct.1`) &
    de$`pct.1`>=min_pct &
    is.finite(de$p_val_adj) &
    de$p_val_adj<=max_padj
  de <- de[keep,,drop=FALSE]
  if (!nrow(de)) stop("No markers passed primary filtering.",call.=FALSE)
  if (isTRUE(exclude_gnrh)) de <- de[toupper(de$gene)!=toupper(gnrh_gene),,drop=FALSE]
  if (!nrow(de)) stop("No markers remained after excluding GNRH1/Gnrh1.",call.=FALSE)
  log("Candidate markers after primary filtering:",nrow(de),type="info")
  # ------------------------------------------------------------------------- #
  # Association setup
  # ------------------------------------------------------------------------- #
  assoc_cells <- unique(c(grp$cells_1,grp$cells_2))
  gnrh_vals <- as.numeric(expr[gnrh_gene,assoc_cells,drop=TRUE])
  gnrh_pos_vals <- gnrh_vals[is.finite(gnrh_vals) & gnrh_vals>0]
  gnrh_var <- if (length(gnrh_pos_vals)>1L) stats::var(gnrh_pos_vals) else NA_real_
  low_information <- length(gnrh_pos_vals)<30L ||
    length(unique(gnrh_pos_vals))<3L ||
    !is.finite(gnrh_var) ||
    gnrh_var<0.25
  if (is.null(program_col) && "gnrh_identity_score" %in% colnames(object[[]])) program_col <- "gnrh_identity_score"
  has_program <- !is.null(program_col) &&
    length(program_col)==1L &&
    !is.na(program_col) &&
    nzchar(program_col) &&
    program_col %in% colnames(object[[]])
  resolved_mode <- association_mode
  if (association_mode=="auto") {
    if (low_information && has_program) resolved_mode <- "combined"
    else if (low_information) resolved_mode <- "codetection"
    else resolved_mode <- "coexpression"
  }
  if (resolved_mode=="program" && !has_program) stop("`association_mode = 'program'` requires a valid `program_col`.",call.=FALSE)
  if (resolved_mode=="combined" && !has_program) warning("No valid GnRH program column found; `combined` will use co-detection only.",call.=FALSE)
  log("Association mode:",resolved_mode,type="info")
  log("Correlation method:",coexpr.method,type="info")
  # ------------------------------------------------------------------------- #
  # Required association metrics
  # ------------------------------------------------------------------------- #
  need_coexpr <- resolved_mode %in% c("coexpression","combined") || !is.null(coexpr_min)
  need_codetect <- resolved_mode %in% c("codetection","combined") || !is.null(codetect_min_or) || !is.null(codetect_max_fdr)
  need_program <- (resolved_mode %in% c("program","combined") || !is.null(program_min)) && has_program
  # ------------------------------------------------------------------------- #
  # GnRH phenotype association
  # ------------------------------------------------------------------------- #
  log("Computing phenotype association",type="step")
  phenotype <- .compute_phenotype_association(expr=expr,genes=de$gene,cells_1=grp$cells_1,cells_2=grp$cells_2)
  idx <- match(de$gene,phenotype$gene)
  de$phenotype_or <- phenotype$phenotype_or[idx]
  de$phenotype_log2or <- phenotype$phenotype_log2or[idx]
  de$phenotype_delta <- phenotype$phenotype_delta[idx]
  de$phenotype_p <- phenotype$phenotype_p[idx]
  de$phenotype_fdr <- phenotype$phenotype_fdr[idx]
  # ------------------------------------------------------------------------- #
  # Initialize association metrics
  # ------------------------------------------------------------------------- #
  de$coexpr_cor <- NA_real_
  de$coexpr <- NA_real_
  de$codetect_or <- NA_real_
  de$codetect_log2or <- NA_real_
  de$codetect_p <- NA_real_
  de$codetect_fdr <- NA_real_
  de$program_cor <- NA_real_
  # ------------------------------------------------------------------------- #
  # GNRH1 quantitative co-expression
  # ------------------------------------------------------------------------- #
  if (need_coexpr) {
    log("Computing GNRH1 co-expression",type="step")
    coexpr <- .compute_coexpr(expr=expr,genes=de$gene,gnrh_gene=gnrh_gene,cells=grp$cells_1,method=coexpr.method)
    de$coexpr_cor <- unname(coexpr[de$gene])
    de$coexpr_cor[!is.finite(de$coexpr_cor)] <- 0
    de$coexpr <- de$coexpr_cor
  }
  # ------------------------------------------------------------------------- #
  # Binary GNRH1 co-detection
  # ------------------------------------------------------------------------- #
  if (need_codetect) {
    log("Computing GNRH1 co-detection",type="step")
    codetect <- .compute_codetect(expr=expr,genes=de$gene,gnrh_gene=gnrh_gene,cells=assoc_cells)
    idx <- match(de$gene,codetect$gene)
    de$codetect_or <- codetect$odds_ratio[idx]
    de$codetect_log2or <- codetect$log2_or[idx]
    de$codetect_p <- codetect$p_value[idx]
    de$codetect_fdr <- codetect$p_adj[idx]
  }
  # ------------------------------------------------------------------------- #
  # GnRH program association
  # ------------------------------------------------------------------------- #
  if (need_program) {
    log("Computing GnRH program association",type="step")
    program_cor <- .compute_program_association(
      object=object,expr=expr,genes=de$gene,cells=assoc_cells,
      program_col=program_col,method=coexpr.method
    )
    de$program_cor <- unname(program_cor[de$gene])
    de$program_cor[!is.finite(de$program_cor)] <- 0
  }
  # ------------------------------------------------------------------------- #
  # Association scaling
  # ------------------------------------------------------------------------- #
  scale01 <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    if (!any(is.finite(x))) return(rep(0,length(x)))
    rng <- range(x,na.rm=TRUE)
    if (!all(is.finite(rng)) || diff(rng)==0) return(rep(0,length(x)))
    out <- (x-rng[1])/diff(rng)
    out[!is.finite(out)] <- 0
    out
  }
  coexpr_s <- scale01(pmax(de$coexpr_cor,0))
  codetect_s <- scale01(pmax(de$codetect_log2or,0))
  program_s <- scale01(pmax(de$program_cor,0))
  # ------------------------------------------------------------------------- #
  # Association score
  # ------------------------------------------------------------------------- #
  de$association_score <- switch(
    resolved_mode,
    coexpression=coexpr_s,
    codetection=codetect_s,
    program=program_s,
    combined=if (has_program && need_program) rowMeans(cbind(codetect_s,program_s),na.rm=TRUE) else codetect_s
  )
  de$association_score[!is.finite(de$association_score)] <- 0
  de$association_mode <- resolved_mode
  # ------------------------------------------------------------------------- #
  # Optional association filtering
  # ------------------------------------------------------------------------- #
  keep <- rep(TRUE,nrow(de))
  if (!is.null(coexpr_min)) keep <- keep & is.finite(de$coexpr_cor) & de$coexpr_cor>=coexpr_min
  if (!is.null(codetect_min_or)) keep <- keep & is.finite(de$codetect_or) & de$codetect_or>=codetect_min_or
  if (!is.null(codetect_max_fdr)) keep <- keep & is.finite(de$codetect_fdr) & de$codetect_fdr<=codetect_max_fdr
  if (!is.null(program_min)) keep <- keep & is.finite(de$program_cor) & de$program_cor>=program_min
  de <- de[keep,,drop=FALSE]
  # ------------------------------------------------------------------------- #
  # Ranking
  # ------------------------------------------------------------------------- #
  if (nrow(de)) {
    de <- .compute_marker_score(de)
    de$association_score[!is.finite(de$association_score)] <- 0
    de$marker_score <- de$marker_score*(1+de$association_score)
    de <- de[order(-de$marker_score,de$p_val_adj),,drop=FALSE]
    de <- de[!duplicated(de$gene),,drop=FALSE]
    rownames(de) <- de$gene
  }
  # ------------------------------------------------------------------------- #
  # Store attributes
  # ------------------------------------------------------------------------- #
  attr(de,"association_mode") <- resolved_mode
  attr(de,"gnrh_low_information") <- low_information
  attr(de,"gnrh_gene") <- gnrh_gene
  attr(de,"coexpr_method") <- coexpr.method
  log("Markers detected:",nrow(de),type="done")
  if (nrow(de)) log("Top marker:",de$gene[[1L]],type="info")
  de
}



# ============================================================================= #
# Identify specific and reproducible GnRH-expressed genes
# ============================================================================= #

#' Identify specific and reproducible GnRH-expressed genes
#'
#' Identifies genes enriched in GnRH-positive cells relative to matched negative
#' neuronal controls, then prioritizes candidates using expression specificity,
#' same-cell co-detection with the GnRH transcript, statistical significance,
#' and optional donor-level reproducibility.
#'
#' @param object A Seurat object containing GnRHcell classifications.
#' @param status_col Character scalar giving the GnRH status metadata column.
#'   Default is \code{"gnrh_status"}.
#' @param positive Character vector defining GnRH-positive status labels.
#'   Default is \code{"pos"}.
#' @param negative Character vector defining GnRH-negative status labels.
#'   Default is \code{"neg"}.
#' @param confidence_col Optional character scalar giving a metadata column used
#'   to restrict the positive population to confident GnRH cells.
#' @param confidence_value Value or values in \code{confidence_col} considered
#'   confident. Default is \code{TRUE}.
#' @param annotation_col Character scalar giving the broad cell-type annotation
#'   used to select matched negative controls. Default is \code{"ann1"}.
#' @param control_ident Character vector giving the identities retained as
#'   matched controls. Default is \code{"Neuronal"}.
#' @param donor_col Optional character scalar identifying biological donors,
#'   samples, or replicates.
#' @param assay Character scalar giving the assay used for differential
#'   expression. Default is \code{"RNA"}.
#' @param layer Character scalar giving the normalized-expression layer used for
#'   co-detection and donor-level calculations. Default is \code{"data"}.
#' @param coexpr_gene Character scalar giving the GnRH reference gene used for
#'   same-cell co-detection. Default is \code{"GNRH1"}.
#' @param methods Character vector of differential-expression methods passed to
#'   \code{\link{gnrh_markers}}. Default is \code{"wilcox"}.
#' @param min_fc Minimum average log2 fold change required for candidate
#'   retention. Default is \code{0.5}.
#' @param min_pct Minimum fraction of GnRH cells expressing a candidate gene.
#'   Default is \code{0.20}.
#' @param max_control_pct Maximum fraction of matched control cells expressing
#'   a candidate gene. Default is \code{0.10}.
#' @param min_coexpr_pct Minimum percentage of GNRH-positive GnRH cells in which
#'   the candidate gene must be co-detected. Default is \code{20}.
#' @param max_padj Maximum adjusted P value allowed for candidate markers.
#'   Default is \code{0.05}.
#' @param min_donor_pct Minimum within-donor detection fraction required for a
#'   donor to support a candidate gene. Default is \code{0.10}.
#' @param min_cells_donor Minimum number of GnRH cells required for a donor to
#'   be evaluable. Default is \code{3L}.
#' @param exclude_transitional_controls Logical. Whether to exclude cells with
#'   direct isolated GnRH signal, transcriptomic candidate status, dropout
#'   candidate status, or direct GnRH signal from the matched control
#'   population. Default is \code{TRUE}.
#' @param generic_genes Character vector of generic neuronal, housekeeping, or
#'   broadly expressed genes excluded from the final candidate table.
#' @param seed Integer random seed. Default is \code{1234L}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{candidates}}{Filtered and ranked candidate GnRH markers.}
#'   \item{\code{markers}}{Complete marker table returned by
#'     \code{\link{gnrh_markers}} before the final candidate filters.}
#'   \item{\code{donor_detection}}{Gene-by-donor matrix of within-donor
#'     detection fractions when \code{donor_col} is supplied.}
#'   \item{\code{donor_cells}}{Number of selected GnRH cells per donor.}
#'   \item{\code{gnrh_cells}}{Cell barcodes assigned to the GnRH population.}
#'   \item{\code{control_cells}}{Cell barcodes assigned to matched controls.}
#'   \item{\code{parameters}}{Principal parameters used for marker discovery.}
#' }
#'
#' @details
#' GnRH-positive cells are contrasted against negative cells restricted to
#' \code{control_ident}. Differential-expression discovery is performed with
#' \code{\link{gnrh_markers}} using permissive internal thresholds so that the
#' final filters can be applied consistently by this function.
#'
#' Same-cell co-detection is calculated among GnRH cells in which
#' \code{coexpr_gene} itself is detected. Candidate genes must satisfy the
#' requested fold-change, target detection, control detection, adjusted
#' P-value, and co-detection thresholds.
#'
#' Candidate specificity is ranked using positive fold change, the difference
#' in detection frequency between GnRH and control cells, target-cell detection
#' frequency, co-detection with the GnRH transcript, and statistical
#' significance.
#'
#' When \code{donor_col} is supplied, donor support is defined as the fraction
#' of evaluable donors in which a candidate is detected in at least
#' \code{min_donor_pct} of selected GnRH cells. Donor support modifies the
#' final ranking score but is not used as a strict exclusion criterion.
#'
#' The GnRH reference gene itself is excluded from the final candidate table
#' because it defines the target lineage rather than representing a newly
#' discovered marker.
#'
#' @examples
#' \dontrun{
#' genes <- find_gnrh_genes(
#'   object = wang,
#'   annotation_col = "ann1",
#'   control_ident = "Neuron",
#'   donor_col = "orig.ident"
#' )
#'
#' head(genes$candidates)
#' }
#'
#' @seealso
#' \code{\link{gnrh_markers}},
#' \code{\link{find_gnrh_stage_markers}}
#'
#' @family marker discovery
#' @export
find_gnrh_genes <- function(
    object,status_col="gnrh_status",positive="pos",negative="neg",
    confidence_col=NULL,confidence_value=TRUE,annotation_col="ann1",
    control_ident="Neuronal",donor_col=NULL,assay="RNA",layer="data",
    coexpr_gene="GNRH1",methods="wilcox",min_fc=0.5,min_pct=0.20,
    max_control_pct=0.10,min_coexpr_pct=20,max_padj=0.05,
    min_donor_pct=0.10,min_cells_donor=3L,
    exclude_transitional_controls=TRUE,
    generic_genes=c("TUBB3","RBFOX3","GAPDH","ACTB","MALAT1","GAD1","GAD2","SLC17A7","SNAP25","SYT1"),
    seed=1234L,verbose=TRUE
) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  set.seed(seed)
  log <- .msg(verbose)
  md <- object[[]]
  required <- c(status_col,annotation_col,if (!is.null(confidence_col)) confidence_col,if (!is.null(donor_col)) donor_col)
  missing <- setdiff(required,colnames(md))
  if (length(missing)) stop("Missing metadata column(s): ",paste(missing,collapse=", "),call.=FALSE)
  expression <- .get_expr(object,assay=assay,layer=layer)
  coexpr_gene_match <- .match_genes(coexpr_gene,rownames(expression))
  if (!length(coexpr_gene_match)) stop("Gene not found: ",coexpr_gene,call.=FALSE)
  coexpr_gene <- coexpr_gene_match[[1L]]
  # ------------------------------------------------------------------------- #
  # Select populations
  # ------------------------------------------------------------------------- #
  positive_index <- !is.na(md[[status_col]]) & as.character(md[[status_col]]) %in% positive
  if (!is.null(confidence_col)) positive_index <- positive_index & !is.na(md[[confidence_col]]) & md[[confidence_col]] %in% confidence_value
  control_index <- !is.na(md[[status_col]]) & as.character(md[[status_col]]) %in% negative & !is.na(md[[annotation_col]]) & md[[annotation_col]] %in% control_ident
  if (isTRUE(exclude_transitional_controls)) {
    if ("gnrh_direct_isolated" %in% colnames(md)) control_index <- control_index & !(md$gnrh_direct_isolated %in% TRUE)
    candidate_col <- if ("gnrh_transcriptomic_candidate" %in% colnames(md)) "gnrh_transcriptomic_candidate" else if ("gnrh_dropout_candidate" %in% colnames(md)) "gnrh_dropout_candidate" else NULL
    if (!is.null(candidate_col)) control_index <- control_index & !(md[[candidate_col]] %in% TRUE)
    if ("gnrh_direct_signal" %in% colnames(md)) control_index <- control_index & !(md$gnrh_direct_signal %in% TRUE)
  }
  gnrh_cells <- rownames(md)[positive_index]
  control_cells <- rownames(md)[control_index]
  if (length(gnrh_cells)<5L) stop("Fewer than five GnRH cells selected.",call.=FALSE)
  if (length(control_cells)<5L) stop("Fewer than five matched control cells selected.",call.=FALSE)
  log(sprintf("Selected %d GnRH and %d control cells",length(gnrh_cells),length(control_cells)))
  # ------------------------------------------------------------------------- #
  # Marker discovery
  # ------------------------------------------------------------------------- #
  object$gnrh_comparison <- "Other"
  object$gnrh_comparison[gnrh_cells] <- "GnRH"
  object$gnrh_comparison[control_cells] <- "Control"
  markers <- gnrh_markers(
    object=object,group.by="gnrh_comparison",ident.1="GnRH",ident.2="Control",
    assay=assay,layer=layer,methods=methods,
    min_pct=min(0.05,min_pct),min_fc=min(0.25,min_fc),max_padj=max_padj,
    coexpr_min=NULL,exclude_gnrh=FALSE,verbose=verbose,
    only.pos=TRUE,logfc.threshold=0
  )
  if (!nrow(markers)) stop("No positive GnRH markers detected.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Same-cell GNRH co-detection
  # ------------------------------------------------------------------------- #
  genes <- intersect(unique(c(coexpr_gene,markers$gene)),rownames(expression))
  expr_pos <- expression[genes,gnrh_cells,drop=FALSE]
  gnrh_detected <- expr_pos[coexpr_gene,,drop=TRUE]>0
  n_gnrh_detected <- sum(gnrh_detected)
  if (n_gnrh_detected>0L) {
    detection <- expr_pos[markers$gene,gnrh_detected,drop=FALSE]>0
    markers$coexpr_n <- Matrix::rowSums(detection)
    markers$coexpr_pct <- 100*Matrix::rowMeans(detection)
  } else {
    markers$coexpr_n <- 0L
    markers$coexpr_pct <- 0
  }
  markers$coexpr_flag <- markers$coexpr_pct>=min_coexpr_pct
  # ------------------------------------------------------------------------- #
  # Candidate filtering
  # ------------------------------------------------------------------------- #
  candidates <- markers[
    is.finite(markers$avg_log2FC) &
      markers$avg_log2FC>=min_fc &
      markers$`pct.1`>=min_pct &
      markers$`pct.2`<=max_control_pct &
      is.finite(markers$p_val_adj) &
      markers$p_val_adj<=max_padj &
      markers$coexpr_flag,
    ,drop=FALSE
  ]
  candidates <- candidates[toupper(candidates$gene)!=toupper(coexpr_gene),,drop=FALSE]
  if (length(generic_genes)) candidates <- candidates[!toupper(candidates$gene) %in% toupper(generic_genes),,drop=FALSE]
  significance <- pmin(-log10(pmax(candidates$p_val_adj,.Machine$double.xmin)),50)
  candidates$specificity_score <- candidates$avg_log2FC*pmax(candidates$`pct.1`-candidates$`pct.2`,0)*sqrt(pmax(candidates$`pct.1`,0))*log1p(candidates$coexpr_pct)*significance
  # ------------------------------------------------------------------------- #
  # Donor reproducibility
  # ------------------------------------------------------------------------- #
  donor_detection <- NULL
  donor_cells <- NULL
  if (!is.null(donor_col) && nrow(candidates)) {
    donor <- as.character(md[gnrh_cells,donor_col,drop=TRUE])
    donor_cells <- table(donor,useNA="no")
    evaluable <- names(donor_cells)[donor_cells>=min_cells_donor]
    if (length(evaluable)) {
      candidate_expression <- expression[candidates$gene,gnrh_cells,drop=FALSE]
      donor_detection <- vapply(evaluable,function(d) Matrix::rowMeans(candidate_expression[,donor==d,drop=FALSE]>0),numeric(nrow(candidate_expression)))
      if (is.null(dim(donor_detection))) donor_detection <- matrix(donor_detection,ncol=1L,dimnames=list(rownames(candidate_expression),evaluable))
      else {
        rownames(donor_detection) <- rownames(candidate_expression)
        colnames(donor_detection) <- evaluable
      }
      candidates$n_donors <- rowSums(donor_detection[candidates$gene,,drop=FALSE]>=min_donor_pct)
      candidates$n_donors_evaluable <- length(evaluable)
      candidates$donor_support <- candidates$n_donors/candidates$n_donors_evaluable
    }
  }
  if (!"n_donors" %in% colnames(candidates)) candidates$n_donors <- NA_integer_
  if (!"n_donors_evaluable" %in% colnames(candidates)) candidates$n_donors_evaluable <- NA_integer_
  if (!"donor_support" %in% colnames(candidates)) candidates$donor_support <- NA_real_
  # ------------------------------------------------------------------------- #
  # Final ranking
  # ------------------------------------------------------------------------- #
  donor_weight <- ifelse(is.na(candidates$donor_support),1,0.5+candidates$donor_support)
  candidates$final_score <- candidates$specificity_score*donor_weight
  candidates <- candidates[order(-candidates$final_score,-candidates$donor_support,candidates$p_val_adj,na.last=TRUE),,drop=FALSE]
  rownames(candidates) <- candidates$gene
  list(
    candidates=candidates,markers=markers,donor_detection=donor_detection,
    donor_cells=donor_cells,gnrh_cells=gnrh_cells,control_cells=control_cells,
    parameters=list(
      status_col=status_col,positive=positive,negative=negative,
      confidence_col=confidence_col,confidence_value=confidence_value,
      annotation_col=annotation_col,control_ident=control_ident,
      assay=assay,layer=layer,min_fc=min_fc,min_pct=min_pct,
      max_control_pct=max_control_pct,min_coexpr_pct=min_coexpr_pct,
      max_padj=max_padj,donor_col=donor_col,min_donor_pct=min_donor_pct,
      min_cells_donor=min_cells_donor,
      exclude_transitional_controls=exclude_transitional_controls
    )
  )
}

# ============================================================================= #
# Identify stage-specific markers in GnRH-lineage cells
# ============================================================================= #

#' Identify stage-specific markers in GnRH-lineage cells
#'
#' Identifies genes enriched in each GnRH developmental stage relative to the
#' other selected GnRH stages using one-versus-rest differential-expression
#' testing, detection specificity, statistical significance, and optional
#' donor-level reproducibility.
#'
#' @param object A Seurat object processed with GnRHcell and containing
#'   developmental-stage metadata.
#' @param stage_col Character scalar giving the metadata column containing GnRH
#'   developmental-stage labels. Default is \code{"gnrh_stage"}.
#' @param stages Character vector giving the stages to compare. Default is
#'   \code{c("identity", "migrating", "mature")}.
#' @param class_col Optional character scalar giving the GnRH detection-class
#'   metadata column. When supplied, only cells whose class belongs to
#'   \code{positive_classes} are retained. Default is \code{"gnrh_class"}.
#' @param positive_classes Character vector defining GnRH-positive classes.
#'   Default is \code{c("direct", "supported")}.
#' @param exclude_stages Character vector of stage labels excluded before
#'   analysis. Default is \code{"non-gnrh"}.
#' @param assay Optional character scalar giving the assay used for
#'   differential expression. If \code{NULL}, the default Seurat assay is used.
#' @param layer Character scalar giving the expression layer used for
#'   donor-level concordance calculations. Default is \code{"data"}.
#' @param test_use Character scalar giving the differential-expression method
#'   passed to \code{Seurat::FindMarkers()}. Default is \code{"wilcox"}.
#' @param min_cells Minimum number of selected cells required in each analyzed
#'   developmental stage. Default is \code{10L}.
#' @param min_pct Minimum fraction of cells expressing a gene passed to
#'   \code{Seurat::FindMarkers()}. Default is \code{0.10}.
#' @param min_log2fc Minimum positive average log2 fold change required for a
#'   final candidate marker. Default is \code{0.25}.
#' @param max_padj Maximum adjusted P value allowed for candidate markers.
#'   Default is \code{0.05}.
#' @param min_specificity Minimum detection-frequency difference
#'   \code{pct.1 - pct.2} required for candidate retention.
#'   Default is \code{0.05}.
#' @param donor_col Optional character scalar identifying biological donors,
#'   samples, or replicates.
#' @param min_cells_donor Minimum number of cells required in both the target
#'   stage and donor-specific reference population for a donor to be evaluable.
#'   Default is \code{3L}.
#' @param min_donor_support Minimum fraction of evaluable donors in which the
#'   target stage shows higher mean expression than its donor-specific
#'   reference. Default is \code{0.60}.
#' @param min_donors Minimum number of evaluable donors required before donor
#'   reproducibility is used as a candidate filter. Default is \code{2L}.
#' @param max_cells_per_ident Maximum number of cells sampled per identity by
#'   \code{Seurat::FindMarkers()}. Default is \code{Inf}.
#' @param exclude_pattern Optional regular expression identifying genes to
#'   exclude from the final candidate table.
#' @param seed Integer random seed used for differential-expression
#'   subsampling. Default is \code{1234L}.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers()}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{candidates}}{Filtered and ranked stage-specific candidate
#'     markers.}
#'   \item{\code{markers}}{Complete positive marker results from all
#'     stage-versus-rest comparisons before final candidate filtering.}
#'   \item{\code{donor_effects}}{Long-format donor-level expression and
#'     detection effects.}
#'   \item{\code{stage_counts}}{Number of selected cells in each analyzed
#'     developmental stage.}
#'   \item{\code{cells}}{Cell barcodes included in the analysis.}
#'   \item{\code{parameters}}{Principal marker-discovery parameters.}
#' }
#'
#' @details
#' Each requested developmental stage is compared against the union of the
#' other retained stages. Only positive markers are returned by
#' \code{Seurat::FindMarkers()}.
#'
#' Marker specificity is defined as:
#' \deqn{specificity = pct.1 - pct.2}
#' where \code{pct.1} is the detection fraction in the target stage and
#' \code{pct.2} is the detection fraction in the other GnRH stages.
#'
#' When \code{donor_col} is supplied, donors are considered evaluable only when
#' at least \code{min_cells_donor} target cells and \code{min_cells_donor}
#' reference cells are available. Mean-expression and detection-frequency
#' differences are then calculated independently within each donor.
#'
#' Donor support is the fraction of evaluable donors with a positive
#' target-versus-reference mean-expression difference. Candidates with at least
#' \code{min_donors} evaluable donors must satisfy
#' \code{min_donor_support}. Genes with insufficient donor information are
#' retained.
#'
#' The final ranking score combines fold change, detection specificity,
#' target-stage detection frequency, statistical significance, and donor
#' support:
#' \deqn{
#' max(logFC,0) \times max(specificity,0) \times
#' sqrt(pct.1) \times significance \times donor\ weight
#' }
#'
#' Cell-level differential-expression testing is intended primarily for marker
#' discovery and prioritization. Replicate-aware pseudobulk analysis is
#' preferable for formal differential-expression inference when biological
#' replicates are available.
#'
#' @examples
#' \dontrun{
#' stage_markers <- find_gnrh_stage_markers(
#'   object = wang,
#'   donor_col = "orig.ident"
#' )
#'
#' head(stage_markers$candidates)
#'
#' migrating <- subset(
#'   stage_markers$candidates,
#'   stage == "migrating"
#' )
#'
#' stage_markers <- find_gnrh_stage_markers(
#'   object = human_hypomap,
#'   donor_col = "orig.ident",
#'   max_cells_per_ident = 5000
#' )
#' }
#'
#' @seealso
#' \code{\link{gnrh_markers}},
#' \code{\link{find_gnrh_genes}},
#' \code{\link{stage_gnrh}}
#'
#' @family marker discovery
#' @export
find_gnrh_stage_markers <- function(
    object,stage_col="gnrh_stage",
    stages=c("identity","migrating","mature"),
    class_col="gnrh_class",
    positive_classes=c("direct","supported"),
    exclude_stages="non-gnrh",
    assay=NULL,layer="data",test_use="wilcox",
    min_cells=10L,min_pct=0.10,min_log2fc=0.25,max_padj=0.05,
    min_specificity=0.05,donor_col=NULL,min_cells_donor=3L,
    min_donor_support=0.60,min_donors=2L,max_cells_per_ident=Inf,
    exclude_pattern=NULL,seed=1234L,verbose=TRUE,...
) {
  if (!inherits(object,"Seurat")) stop("`object` must be a Seurat object.",call.=FALSE)
  log <- .msg(verbose)
  metadata <- object[[]]
  required <- c(stage_col,if (!is.null(class_col)) class_col,if (!is.null(donor_col)) donor_col)
  missing <- setdiff(required,colnames(metadata))
  if (length(missing)) stop("Missing metadata column(s): ",paste(missing,collapse=", "),call.=FALSE)
  assay <- assay %||% Seurat::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay `",assay,"` was not found.",call.=FALSE)
  # ------------------------------------------------------------------------- #
  # Select GnRH cells
  # ------------------------------------------------------------------------- #
  stage_values <- as.character(metadata[[stage_col]])
  keep <- !is.na(stage_values) & nzchar(stage_values) & !stage_values %in% exclude_stages
  if (!is.null(class_col)) keep <- keep & !is.na(metadata[[class_col]]) & as.character(metadata[[class_col]]) %in% positive_classes
  observed <- unique(stage_values[keep])
  stages <- intersect(stages,observed)
  if (length(stages)<2L) stop("At least two requested stages must be present. Present: ",paste(observed,collapse=", "),call.=FALSE)
  keep <- keep & stage_values %in% stages
  cells <- rownames(metadata)[keep]
  stage_factor <- factor(stage_values[keep],levels=stages)
  stage_counts <- table(stage_factor)
  too_small <- names(stage_counts)[stage_counts<min_cells]
  if (length(too_small)) stop("Stage(s) below `min_cells = ",min_cells,"`: ",paste(paste0(too_small," (n=",stage_counts[too_small],")"),collapse=", "),call.=FALSE)
  log(paste(paste0(names(stage_counts),"=",as.integer(stage_counts)),collapse=" | "))
  # ------------------------------------------------------------------------- #
  # Analysis object
  # ------------------------------------------------------------------------- #
  analysis_object <- subset(object,cells=cells)
  analysis_object$.gnrh_stage_test <- stage_factor
  SeuratObject::Idents(analysis_object) <- ".gnrh_stage_test"
  set.seed(seed)
  # ------------------------------------------------------------------------- #
  # Stage-versus-rest DE
  # ------------------------------------------------------------------------- #
  run_stage <- function(stage) {
    log(paste0("Testing ",stage," versus other GnRH stages"),type="step")
    reference <- setdiff(stages,stage)
    result <- Seurat::FindMarkers(
      object=analysis_object,ident.1=stage,ident.2=reference,assay=assay,
      test.use=test_use,only.pos=TRUE,min.pct=min_pct,
      logfc.threshold=min_log2fc,max.cells.per.ident=max_cells_per_ident,
      random.seed=seed,verbose=FALSE,...
    )
    if (!nrow(result)) return(NULL)
    result$gene <- rownames(result)
    result$stage <- stage
    result$reference <- paste(reference,collapse="+")
    fc_col <- intersect(c("avg_log2FC","avg_logFC"),colnames(result))
    if (!length(fc_col)) stop("FindMarkers result has no average logFC column.",call.=FALSE)
    result$avg_log2FC <- result[[fc_col[[1L]]]]
    result$specificity <- result$`pct.1`-result$`pct.2`
    result
  }
  marker_list <- lapply(stages,run_stage)
  marker_list <- marker_list[!vapply(marker_list,is.null,logical(1))]
  if (!length(marker_list)) stop("No stage markers were detected.",call.=FALSE)
  markers <- do.call(rbind,marker_list)
  rownames(markers) <- NULL
  # ------------------------------------------------------------------------- #
  # Donor concordance
  # ------------------------------------------------------------------------- #
  donor_effects <- data.frame()
  donor_summary <- data.frame()
  if (!is.null(donor_col)) {
    expression <- .get_expr(analysis_object,assay=assay,layer=layer)
    genes <- intersect(unique(markers$gene),rownames(expression))
    expression <- expression[genes,cells,drop=FALSE]
    donor <- as.character(metadata[cells,donor_col,drop=TRUE])
    stage_by_cell <- as.character(metadata[cells,stage_col,drop=TRUE])
    donor_rows <- lapply(stages,function(stage) {
      other_stages <- setdiff(stages,stage)
      donors <- unique(donor[!is.na(donor) & nzchar(donor)])
      evaluable <- donors[vapply(donors,function(id) {
        idx <- donor==id
        sum(idx & stage_by_cell==stage)>=min_cells_donor &&
          sum(idx & stage_by_cell %in% other_stages)>=min_cells_donor
      },logical(1))]
      if (!length(evaluable)) return(NULL)
      stage_genes <- intersect(markers$gene[markers$stage==stage],genes)
      if (!length(stage_genes)) return(NULL)
      do.call(rbind,lapply(evaluable,function(id) {
        target <- donor==id & stage_by_cell==stage
        reference <- donor==id & stage_by_cell %in% other_stages
        target_mat <- expression[stage_genes,target,drop=FALSE]
        ref_mat <- expression[stage_genes,reference,drop=FALSE]
        mean_difference <- Matrix::rowMeans(target_mat)-Matrix::rowMeans(ref_mat)
        pct_target <- Matrix::rowMeans(target_mat>0)
        pct_reference <- Matrix::rowMeans(ref_mat>0)
        data.frame(
          stage=stage,gene=stage_genes,donor=id,
          mean_difference=as.numeric(mean_difference),
          detection_difference=as.numeric(pct_target-pct_reference),
          pct_target=as.numeric(pct_target),pct_reference=as.numeric(pct_reference),
          n_target=sum(target),n_reference=sum(reference),
          stringsAsFactors=FALSE
        )
      }))
    })
    donor_rows <- donor_rows[!vapply(donor_rows,is.null,logical(1))]
    if (length(donor_rows)) donor_effects <- do.call(rbind,donor_rows)
    if (nrow(donor_effects)) {
      donor_summary <- donor_effects |>
        dplyr::group_by(.data$stage,.data$gene) |>
        dplyr::summarise(
          n_donors=dplyr::n(),
          donors_positive=sum(.data$mean_difference>0),
          donor_support=mean(.data$mean_difference>0),
          donor_detection_support=mean(.data$detection_difference>0),
          median_donor_effect=stats::median(.data$mean_difference),
          median_detection_difference=stats::median(.data$detection_difference),
          .groups="drop"
        )
      markers <- dplyr::left_join(markers,donor_summary,by=c("stage","gene"))
    }
  }
  # ------------------------------------------------------------------------- #
  # Candidate filtering
  # ------------------------------------------------------------------------- #
  donor_cols <- c("n_donors","donors_positive","donor_support","donor_detection_support","median_donor_effect","median_detection_difference")
  for (column in donor_cols) if (!column %in% colnames(markers)) markers[[column]] <- NA_real_
  candidates <- markers[
    is.finite(markers$avg_log2FC) &
      markers$avg_log2FC>=min_log2fc &
      is.finite(markers$p_val_adj) &
      markers$p_val_adj<=max_padj &
      markers$specificity>=min_specificity,
    ,drop=FALSE
  ]
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) candidates <- candidates[!grepl(exclude_pattern,candidates$gene,ignore.case=TRUE),,drop=FALSE]
  # ------------------------------------------------------------------------- #
  # Donor reproducibility
  # ------------------------------------------------------------------------- #
  candidates$donor_evaluable <- !is.na(candidates$n_donors) & candidates$n_donors>=min_donors
  candidates$donor_reproducible <- candidates$donor_evaluable & candidates$donor_support>=min_donor_support
  if (!is.null(donor_col)) candidates <- candidates[!candidates$donor_evaluable | candidates$donor_reproducible,,drop=FALSE]
  # ------------------------------------------------------------------------- #
  # Stage score
  # ------------------------------------------------------------------------- #
  significance <- pmin(-log10(pmax(candidates$p_val_adj,.Machine$double.xmin)),50)
  donor_weight <- ifelse(is.na(candidates$donor_support),1,0.5+candidates$donor_support)
  candidates$stage_score <- pmax(candidates$avg_log2FC,0)*pmax(candidates$specificity,0)*sqrt(pmax(candidates$`pct.1`,0))*significance*donor_weight
  candidates <- candidates[order(candidates$stage,-candidates$stage_score,candidates$p_val_adj),,drop=FALSE]
  # ------------------------------------------------------------------------- #
  # Return
  # ------------------------------------------------------------------------- #
  list(
    candidates=candidates,
    markers=markers,
    donor_effects=donor_effects,
    stage_counts=data.frame(stage=names(stage_counts),n_cells=as.integer(stage_counts),stringsAsFactors=FALSE),
    cells=cells,
    parameters=list(
      stages=stages,positive_classes=positive_classes,assay=assay,layer=layer,
      test_use=test_use,min_cells=min_cells,min_pct=min_pct,
      min_log2fc=min_log2fc,max_padj=max_padj,min_specificity=min_specificity,
      donor_col=donor_col,min_cells_donor=min_cells_donor,
      min_donor_support=min_donor_support,min_donors=min_donors
    )
  )
}
