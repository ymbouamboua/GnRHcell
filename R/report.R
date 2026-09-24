# ========================================================================= #
# GnRH diagnostic report
# ========================================================================= #

#' Generate a publication-ready gnrhcell diagnostic report
#'
#' Generates a compact six-panel diagnostic report summarizing GnRH detection,
#' developmental-state assignment, score separation, threshold behaviour,
#' ROC/self-discrimination and marker-program support.
#'
#' Internal ROC curves quantify diagnostic self-consistency only and must not
#' be interpreted as independent validation. When an independent truth
#' annotation is supplied, external ROC performance can be displayed instead.
#'
#' @param object Seurat object processed with `run_gnrh()` and/or
#'   `gnrh_diagnostics()`.
#' @param style Plot theme: `"test"`, `"classic"`, `"minimal"` or `"bw"`.
#' @param mode `"light"` or `"dark"`.
#' @param truth Optional metadata column containing an independent binary
#'   reference annotation.
#' @param roc_mode ROC mode: `"auto"`, `"external"`, `"internal"` or `"none"`.
#' @param positive_truth Values interpreted as positive in `truth`.
#' @param txtsize Base text size.
#' @param max_points Maximum number of cells displayed in scatter plots.
#'   All GnRH-positive cells are retained whenever possible.
#' @param seed Sampling seed.
#' @param show_class_panel Show detection composition when ROC is unavailable.
#' @param verbose Print progress messages.
#'
#' @return A patchwork object.
#' @export
gnrh_report <- function(
    object,
    style=c("test","classic","minimal","bw"),
    mode=c("light","dark"),
    truth=NULL,
    roc_mode=c("auto","external","internal","none"),
    positive_truth=c("1","TRUE","true","positive","Positive","pos","Pos"),
    txtsize=getOption("gnrhcell.base_size",14),
    max_points=100000L,
    seed=1234L,
    show_class_panel=TRUE,
    verbose=TRUE
) {

  # ----------------------------------------------------------------------- #
  # Validation
  # ----------------------------------------------------------------------- #

  style <- match.arg(style)
  mode <- match.arg(mode)
  roc_mode <- match.arg(roc_mode)

  if (!inherits(object,"Seurat"))
    stop("`object` must be a Seurat object.",call.=FALSE)

  if (!requireNamespace("patchwork",quietly=TRUE))
    stop("Package `patchwork` is required.",call.=FALSE)

  if (!requireNamespace("scales",quietly=TRUE))
    stop("Package `scales` is required.",call.=FALSE)

  if (
    length(max_points)!=1L ||
    !is.numeric(max_points) ||
    !is.finite(max_points) ||
    max_points<100
  )
    stop("`max_points` must be >= 100.",call.=FALSE)

  roc_mode <- if (roc_mode=="auto") {
    if (is.null(truth)) "internal" else "external"
  } else {
    roc_mode
  }

  if (roc_mode=="external" && is.null(truth))
    stop("`roc_mode = 'external'` requires `truth`.",call.=FALSE)

  if (roc_mode=="external" && identical(truth,"gnrh_status"))
    stop(
      "`truth = 'gnrh_status'` is not an independent reference. ",
      "Use `roc_mode = 'internal'`.",
      call.=FALSE
    )

  if (verbose)
    message("[INFO] Generating GnRH diagnostic report")

  # ----------------------------------------------------------------------- #
  # Retrieve diagnostics
  # ----------------------------------------------------------------------- #

  g <- object@misc$gnrh

  if (is.null(g$diagnostics))
    stop(
      "GnRH diagnostics are missing. Run `run_gnrh()` or ",
      "`gnrh_diagnostics()` first.",
      call.=FALSE
    )

  d <- as.data.frame(g$diagnostics)
  md <- object[[]]

  req <- c(
    "expr",
    "score",
    "status",
    "core_hits",
    "mig_hits",
    "neuro_hits"
  )

  miss <- setdiff(req,names(d))

  if (length(miss))
    stop(
      "Missing diagnostic column(s): ",
      paste(miss,collapse=", "),
      call.=FALSE
    )

  cells <- rownames(d)

  can_match <-
    length(cells)==nrow(d) &&
    !is.null(rownames(md)) &&
    all(cells %in% rownames(md))

  fields <- c(
    "gnrh_class",
    "gnrh_confident",
    "gnrh_stage",
    "gnrh_secretory",
    "gnrh_support_score",
    "gnrh_support_score_raw",
    "gnrh_alternative_score",
    "gnrh_identity_score",
    "gnrh_migration_score",
    "gnrh_neuro_score",
    "gnrh_stage_early_score",
    "gnrh_stage_migrating_score",
    "gnrh_stage_mature_score",
    "gnrh_secretory_core_hits",
    "gnrh_secretory_supportive_hits",
    "gnrh_secretory_hits",
    "gnrh_knn",
    "gnrh_direct_signal",
    "gnrh_direct_isolated",
    "gnrh_review_candidate",
    "gnrh_signal_class",
    "gnrh_transcriptomic_candidate"
  )

  if (
    !"gnrh_transcriptomic_candidate" %in% names(d) &&
    "gnrh_dropout_candidate" %in% names(md)
  ) {
    d$gnrh_transcriptomic_candidate <- if (can_match) {
      md[cells,"gnrh_dropout_candidate",drop=TRUE]
    } else if (nrow(d)==nrow(md)) {
      md$gnrh_dropout_candidate
    } else {
      NULL
    }
  }

  for (nm in intersect(fields,names(md))) {

    if (nm %in% names(d))
      next

    if (can_match) {
      d[[nm]] <- md[cells,nm,drop=TRUE]
    } else if (nrow(d)==nrow(md)) {
      d[[nm]] <- md[[nm]]
    } else {
      stop(
        "Cannot align diagnostics with Seurat metadata.",
        call.=FALSE
      )
    }
  }

  # ----------------------------------------------------------------------- #
  # Colours
  # ----------------------------------------------------------------------- #

  fg <- if (mode=="dark") "#F5F5F5" else "#1A1A1A"
  bg <- if (mode=="dark") "#111111" else "white"
  muted <- if (mode=="dark") "#BDBDBD" else "#555555"
  grid_col <- if (mode=="dark") "#353535" else "#E6E6E6"
  tile_low <- if (mode=="dark") "#252525" else "#FAFAFA"
  tile_border <- if (mode=="dark") "#111111" else "white"

  status_cols <- gnrh_palette("status")
  stage_cols <- gnrh_palette("stage")

  # ----------------------------------------------------------------------- #
  # Theme
  # ----------------------------------------------------------------------- #

  theme_report <- function(
    leg.pos="right",
    leg.dir = "horizontal",
    x.ang=0,
    axes=TRUE,
    legend.title=TRUE
  ) {

    th <- gnrh_theme(
      txtsize=txtsize,
      x.ang=x.ang,
      leg.pos=leg.pos,
      leg.dir = leg.dir,
      style=style,
      mode=mode
    ) +
      ggplot2::theme(
        plot.title.position="plot",

        plot.title=ggplot2::element_text(
          face="bold",
          size=txtsize+1,
          colour=fg,
          hjust=0,
          margin=ggplot2::margin(b=3)
        ),

        plot.subtitle=ggplot2::element_text(
          size=max(7,txtsize-1.5),
          colour=muted,
          hjust=0,
          margin=ggplot2::margin(b=5)
        ),

        axis.title=ggplot2::element_text(
          size=txtsize,
          colour=fg
        ),

        axis.text=ggplot2::element_text(
          size=max(7,txtsize-1),
          colour=fg
        ),

        legend.title=ggplot2::element_text(
          face="bold",
          size=max(7,txtsize-1),
          colour=fg
        ),

        legend.text=ggplot2::element_text(
          size=max(7,txtsize-1.5),
          colour=fg
        ),

        legend.key.height=grid::unit(0.28,"cm"),
        legend.key.width=grid::unit(0.34,"cm"),

        legend.spacing.x=grid::unit(0.08,"cm"),
        legend.spacing.y=grid::unit(0.04,"cm"),

        legend.box.spacing=grid::unit(0.12,"cm"),

        plot.margin=ggplot2::margin(
          t=5,
          r=5,
          b=5,
          l=5
        )
      )

    if (!legend.title) {
      th <- th +
        ggplot2::theme(
          legend.title=ggplot2::element_blank()
        )
    }

    if (!axes) {
      th <- th +
        ggplot2::theme(
          axis.title=ggplot2::element_blank(),
          axis.text=ggplot2::element_blank(),
          axis.ticks=ggplot2::element_blank(),
          axis.line=ggplot2::element_blank(),
          panel.border=ggplot2::element_blank()
        )
    }

    th
  }

  # ----------------------------------------------------------------------- #
  # Helpers
  # ----------------------------------------------------------------------- #

  empty <- function(title,subtitle=NULL) {

    ggplot2::ggplot() +
      ggplot2::labs(
        title=title,
        subtitle=subtitle
      ) +
      theme_report(
        leg.pos="none",
        axes=FALSE
      )
  }

  format_pct <- function(x) {

    out <- rep(NA_character_,length(x))

    out[is.finite(x) & x>=1] <-
      sprintf("%.1f%%",x[is.finite(x) & x>=1])

    out[is.finite(x) & x<1 & x>=0.01] <-
      sprintf("%.2f%%",x[is.finite(x) & x<1 & x>=0.01])

    out[is.finite(x) & x<0.01] <-
      sprintf("%.3f%%",x[is.finite(x) & x<0.01])

    out
  }

  labels_n <- function(x,lev,denominator=NULL) {

    x <- as.character(x)
    x[is.na(x)] <- "Unknown"

    tab <- table(x)

    lev <- intersect(
      lev,
      names(tab)
    )

    n <- as.integer(tab[lev])

    denom <- if (is.null(denominator)) {
      sum(tab)
    } else {
      denominator
    }

    pct <- if (denom>0) {
      100*n/denom
    } else {
      rep(NA_real_,length(n))
    }

    stats::setNames(
      sprintf(
        "%s  ·  %s",
        format(n,big.mark=","),
        format_pct(pct)
      ),
      lev
    )
  }

  compact_legend <- function(
    nrow=NULL,
    ncol=NULL,
    byrow=TRUE
  ) {

    ggplot2::guide_legend(
      nrow=nrow,
      ncol=ncol,
      byrow=byrow,
      override.aes=list(
        size=2.6,
        alpha=1
      ),
      keyheight=grid::unit(0.27,"cm"),
      keywidth=grid::unit(0.30,"cm")
    )
  }

  # ----------------------------------------------------------------------- #
  # Scatter plotting data
  # ----------------------------------------------------------------------- #

  pd <- d

  if (nrow(pd)>max_points) {

    set.seed(seed)

    pos <- which(
      as.character(pd$status)=="pos"
    )

    neg <- which(
      as.character(pd$status)!="pos"
    )

    keep_pos <- pos

    room <- max(
      0L,
      as.integer(max_points)-length(keep_pos)
    )

    keep_neg <- if (length(neg)>room) {
      sample(neg,room)
    } else {
      neg
    }

    pd <- pd[
      c(keep_neg,keep_pos),
      ,
      drop=FALSE
    ]
  }

  # ========================================================================= #
  # A. Detection landscape
  # ========================================================================= #

  status_lev <- intersect(
    names(status_cols),
    unique(as.character(d$status))
  )

  status_lab <- labels_n(
    d$status,
    status_lev
  )

  p1 <- ggplot2::ggplot(
    pd,
    ggplot2::aes(
      x=.data$expr,
      y=.data$score,
      colour=.data$status
    )
  ) +
    ggplot2::geom_point(
      alpha=0.38,
      size=1.25,
      stroke=0
    ) +
    ggplot2::scale_colour_manual(
      values=status_cols,
      breaks=status_lev,
      labels=status_lab,
      drop=FALSE
    ) +
    ggplot2::labs(
      title="GnRH detection landscape",
      x=expression(italic(GNRH1)~"expression"),
      y="Composite score",
      colour="Status"
    ) +
    theme_report(
      leg.pos="bottom"
    ) +
    ggplot2::guides(
      colour=compact_legend(
        nrow=1
      )
    )

  # ========================================================================= #
  # B. Developmental-stage landscape
  # ========================================================================= #

  if ("gnrh_stage" %in% names(d)) {

    stage <- as.character(d$gnrh_stage)
    status <- as.character(d$status)

    stage[
      is.na(stage) |
        stage=="" |
        status!="pos"
    ] <- "non-gnrh"

    pd_stage <- pd

    pd_stage$plot_stage <- as.character(
      pd_stage$gnrh_stage
    )

    pd_stage$plot_stage[
      is.na(pd_stage$plot_stage) |
        pd_stage$plot_stage=="" |
        as.character(pd_stage$status)!="pos"
    ] <- "non-gnrh"

    stage_cols2 <- c(
      "non-gnrh"=status_cols[["neg"]],
      stage_cols
    )

    stage_order <- c(
      "non-gnrh",
      "early",
      "migrating",
      "post-migratory",
      "mature",
      "transitional"
    )

    stage_lev <- intersect(
      stage_order,
      unique(stage)
    )

    n_gnrh <- sum(
      status=="pos",
      na.rm=TRUE
    )

    stage_counts <- table(stage)

    stage_lab <- stats::setNames(
      vapply(
        stage_lev,
        function(z) {

          n <- as.integer(
            stage_counts[z]
          )

          if (z=="non-gnrh") {
            sprintf(
              "non-GnRH  ·  %s",
              format(n,big.mark=",")
            )
          } else {
            sprintf(
              "%s  ·  %s",
              format(n,big.mark=","),
              format_pct(
                100*n/n_gnrh
              )
            )
          }
        },
        character(1)
      ),
      stage_lev
    )

    p2 <- ggplot2::ggplot(
      pd_stage,
      ggplot2::aes(
        x=.data$expr,
        y=.data$score,
        colour=.data$plot_stage
      )
    ) +
      ggplot2::geom_point(
        alpha=0.38,
        size=1.25,
        stroke=0
      ) +
      ggplot2::scale_colour_manual(
        values=stage_cols2,
        breaks=stage_lev,
        labels=stage_lab,
        drop=FALSE
      ) +
      ggplot2::labs(
        title="Developmental-stage landscape",
        subtitle=paste0(
          "Stage percentages among GnRH+ cells (n = ",
          format(n_gnrh,big.mark=","),
          ")"
        ),
        x=expression(italic(GNRH1)~"expression"),
        y="Composite score",
        colour=NULL
      ) +
      theme_report(
        leg.pos="bottom",
        legend.title=FALSE
      ) +
      ggplot2::guides(
        colour=compact_legend(
          nrow=2
        )
      )

  } else {

    p2 <- empty(
      "Developmental stages unavailable"
    )
  }

  # ========================================================================= #
  # C. Score separation
  # ========================================================================= #

  dd <- d[
    is.finite(d$score),
    ,
    drop=FALSE
  ]

  p3 <- ggplot2::ggplot(
    dd,
    ggplot2::aes(
      x=.data$score,
      fill=.data$status,
      colour=.data$status
    )
  ) +
    ggplot2::geom_density(
      alpha=0.22,
      linewidth=0.65,
      adjust=1,
      na.rm=TRUE
    ) +
    ggplot2::scale_fill_manual(
      values=status_cols,
      breaks=status_lev,
      labels=status_lab,
      drop=FALSE
    ) +
    ggplot2::scale_colour_manual(
      values=status_cols,
      breaks=status_lev,
      labels=status_lab,
      drop=FALSE
    ) +
    ggplot2::labs(
      title="GnRH score separation",
      x="Composite score",
      y="Density",
      fill="Status",
      colour="Status"
    ) +
    theme_report(
      leg.pos="bottom"
    ) +
    ggplot2::guides(
      colour="none",
      fill=compact_legend(
        nrow=1
      )
    )

  # ========================================================================= #
  # D. Threshold performance
  # ========================================================================= #

  threshold_fun <- function(
    score,
    ref,
    n=300L
  ) {

    score <- suppressWarnings(
      as.numeric(score)
    )

    ref <- as.logical(ref)

    ok <-
      is.finite(score) &
      !is.na(ref)

    score <- score[ok]
    ref <- ref[ok]

    if (
      length(score)<2L ||
      length(unique(ref))<2L
    )
      return(NULL)

    rng <- range(score)

    if (
      !all(is.finite(rng)) ||
      diff(rng)<=0
    )
      return(NULL)

    div <- function(a,b) {
      if (b>0) a/b else NA_real_
    }

    thr <- seq(
      rng[1],
      rng[2],
      length.out=n
    )

    out <- lapply(
      thr,
      function(t) {

        pred <- score>=t

        TP <- sum(pred & ref)
        FP <- sum(pred & !ref)
        FN <- sum(!pred & ref)
        TN <- sum(!pred & !ref)

        se <- div(
          TP,
          TP+FN
        )

        sp <- div(
          TN,
          TN+FP
        )

        pr <- div(
          TP,
          TP+FP
        )

        f1 <- if (
          is.finite(pr) &&
          is.finite(se) &&
          pr+se>0
        ) {
          2*pr*se/(pr+se)
        } else {
          NA_real_
        }

        data.frame(
          threshold=t,
          sensitivity=se,
          specificity=sp,
          precision=pr,
          F1=f1
        )
      }
    )

    do.call(
      rbind,
      out
    )
  }

  tc <- NULL
  threshold_mode <- "none"

  if (
    roc_mode=="external" &&
    !is.null(g$threshold_curve)
  ) {

    tc <- as.data.frame(
      g$threshold_curve
    )

    threshold_mode <- "external"

  } else if (roc_mode=="internal") {

    predictor <- NULL

    if ("gnrh_support_score_raw" %in% names(d)) {
      predictor <- d$gnrh_support_score_raw
    } else if ("gnrh_support_score" %in% names(d)) {
      predictor <- d$gnrh_support_score
    } else if ("gnrh_support_score_raw" %in% names(md)) {
      predictor <- if (can_match) {
        md[cells,"gnrh_support_score_raw",drop=TRUE]
      } else if (nrow(d)==nrow(md)) {
        md$gnrh_support_score_raw
      } else {
        NULL
      }
    } else if ("gnrh_support_score" %in% names(md)) {
      predictor <- if (can_match) {
        md[cells,"gnrh_support_score",drop=TRUE]
      } else if (nrow(d)==nrow(md)) {
        md$gnrh_support_score
      } else {
        NULL
      }
    }

    if (!is.null(predictor)) {

      tc <- threshold_fun(
        predictor,
        as.character(d$status)=="pos"
      )

      if (!is.null(tc))
        threshold_mode <- "internal"
    }
  }

  if (
    !is.null(tc) &&
    all(
      c(
        "threshold",
        "sensitivity",
        "specificity",
        "F1"
      ) %in% names(tc)
    )
  ) {

    valid <-
      is.finite(tc$threshold) &
      is.finite(tc$F1)

    bi <- if (any(valid)) {
      which.max(
        replace(
          tc$F1,
          !valid,
          -Inf
        )
      )
    } else {
      NA_integer_
    }

    best <- if (!is.na(bi)) {
      tc$threshold[bi]
    } else {
      NA_real_
    }

    best_f1 <- if (!is.na(bi)) {
      tc$F1[bi]
    } else {
      NA_real_
    }

    long <- rbind(
      data.frame(
        threshold=tc$threshold,
        metric="Sensitivity",
        value=tc$sensitivity
      ),
      data.frame(
        threshold=tc$threshold,
        metric="Specificity",
        value=tc$specificity
      ),
      data.frame(
        threshold=tc$threshold,
        metric="F1",
        value=tc$F1
      )
    )

    metric_cols <- c(
      Sensitivity="#0072B2",
      Specificity="#009E73",
      F1="#D55E00"
    )

    metric_lab <- c(
      Sensitivity="Sensitivity",
      Specificity="Specificity",
      F1="F1"
    )

    subtitle4 <- if (threshold_mode=="internal") {

      if (is.finite(best)) {
        sprintf(
          "Self-consistency · optimum %.3f · max F1 %.3f",
          best,
          best_f1
        )
      } else {
        "Diagnostic self-consistency"
      }

    } else {

      if (is.finite(best)) {
        sprintf(
          "Independent reference · optimum %.3f",
          best
        )
      } else {
        "Independent reference"
      }
    }

    xlab4 <- if (threshold_mode=="internal") {
      "Support-score threshold"
    } else {
      "Threshold"
    }

    p4 <- ggplot2::ggplot(
      long,
      ggplot2::aes(
        x=.data$threshold,
        y=.data$value,
        colour=.data$metric
      )
    ) +
      ggplot2::geom_line(
        linewidth=0.75,
        na.rm=TRUE
      ) +
      ggplot2::geom_vline(
        xintercept=best,
        linetype="dashed",
        colour=muted,
        linewidth=0.45,
        na.rm=TRUE
      ) +
      ggplot2::scale_colour_manual(
        values=metric_cols,
        breaks=names(metric_cols),
        labels=metric_lab
      ) +
      ggplot2::scale_y_continuous(
        limits=c(0,1),
        breaks=seq(0,1,0.25),
        expand=ggplot2::expansion(
          mult=c(0,0.025)
        )
      ) +
      ggplot2::labs(
        title="Threshold performance",
        subtitle=subtitle4,
        x=xlab4,
        y="Performance",
        colour=NULL
      ) +
      theme_report(
        leg.pos="bottom",
        legend.title=FALSE
      ) +
      ggplot2::guides(
        colour=ggplot2::guide_legend(
          nrow=1,
          override.aes=list(
            linewidth=1.2
          ),
          keyheight=grid::unit(0.25,"cm"),
          keywidth=grid::unit(0.45,"cm")
        )
      )

  } else {

    p4 <- empty(
      "Threshold curve unavailable",
      if (roc_mode=="internal") {
        "GnRH support score unavailable"
      } else {
        "Supply an independent reference annotation"
      }
    )
  }

  # ========================================================================= #
  # E. ROC / classification
  # ========================================================================= #

  roc_add <- function(
    ref,
    pred,
    label
  ) {

    ref <- as.integer(ref)

    pred <- suppressWarnings(
      as.numeric(pred)
    )

    ok <-
      !is.na(ref) &
      is.finite(pred)

    if (
      sum(ok)<2L ||
      length(unique(ref[ok]))!=2L
    )
      return(NULL)

    r <- pROC::roc(
      ref[ok],
      pred[ok],
      levels=c(0,1),
      direction="<",
      quiet=TRUE
    )

    data.frame(
      FPR=1-r$specificities,
      TPR=r$sensitivities,
      group=label,
      auc=as.numeric(
        pROC::auc(r)
      )
    )
  }

  p5 <- NULL

  if (
    roc_mode!="none" &&
    requireNamespace("pROC",quietly=TRUE)
  ) {

    rl <- list()

    if (roc_mode=="external") {

      if (!truth %in% names(md))
        stop(
          "Truth column `",
          truth,
          "` not found.",
          call.=FALSE
        )

      tv <- if (can_match) {
        md[cells,truth,drop=TRUE]
      } else if (nrow(d)==nrow(md)) {
        md[[truth]]
      } else {
        stop(
          "Cannot align external truth with diagnostics.",
          call.=FALSE
        )
      }

      rl[["External GnRH"]] <- roc_add(
        as.character(tv) %in% as.character(positive_truth),
        d$score,
        "External GnRH"
      )

    } else {

      rl[["GnRH"]] <- roc_add(
        as.character(d$status)=="pos",
        d$score,
        "GnRH"
      )

      if ("gnrh_stage" %in% names(d)) {

        stage_map <- c(
          early="gnrh_stage_early_score",
          migrating="gnrh_stage_migrating_score",
          mature="gnrh_stage_mature_score"
        )

        st <- as.character(
          d$gnrh_stage
        )

        for (nm in names(stage_map)) {

          score_col <- stage_map[[nm]]

          if (
            nm %in% st &&
            score_col %in% names(d)
          ) {
            rl[[nm]] <- roc_add(
              st==nm,
              d[[score_col]],
              nm
            )
          }
        }
      }

      if (
        all(
          c(
            "gnrh_secretory",
            "gnrh_secretory_hits"
          ) %in% names(d)
        )
      ) {

        rl[["Secretory"]] <- roc_add(
          as.character(d$gnrh_secretory)=="supported",
          d$gnrh_secretory_hits,
          "Secretory"
        )
      }
    }

    rl <- Filter(
      Negate(is.null),
      rl
    )

    if (length(rl)) {

      rd <- do.call(
        rbind,
        rl
      )

      au <- unique(
        rd[c("group","auc")]
      )

      au$label <- sprintf(
        "%s  %.3f",
        au$group,
        au$auc
      )

      rd$group <- factor(
        rd$group,
        levels=au$group
      )

      roc_cols <- c(
        "External GnRH"=status_cols[["pos"]],
        "GnRH"=status_cols[["pos"]],
        stage_cols
      )

      missing_cols <- setdiff(
        au$group,
        names(roc_cols)
      )

      if (length(missing_cols)) {

        extra <- grDevices::hcl.colors(
          length(missing_cols),
          "Dark 3"
        )

        names(extra) <- missing_cols

        roc_cols <- c(
          roc_cols,
          extra
        )
      }

      p5 <- ggplot2::ggplot(
        rd,
        ggplot2::aes(
          x=.data$FPR,
          y=.data$TPR,
          colour=.data$group
        )
      ) +
        ggplot2::geom_abline(
          slope=1,
          intercept=0,
          linetype="dashed",
          colour=muted,
          linewidth=0.4
        ) +
        ggplot2::geom_line(
          linewidth=0.85
        ) +
        ggplot2::scale_colour_manual(
          values=roc_cols,
          breaks=au$group,
          labels=stats::setNames(
            au$label,
            au$group
          )
        ) +
        ggplot2::scale_x_continuous(
          limits=c(0,1),
          breaks=seq(0,1,0.25),
          expand=c(0,0)
        ) +
        ggplot2::scale_y_continuous(
          limits=c(0,1),
          breaks=seq(0,1,0.25),
          expand=c(0,0)
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title=if (roc_mode=="external") {
            "External ROC performance"
          } else {
            "Internal score discrimination"
          },
          subtitle=if (roc_mode=="external") {
            paste0(
              "Independent reference: ",
              truth
            )
          } else {
            "Self-consistency only · AUC shown in legend"
          },
          x="False-positive rate",
          y="True-positive rate",
          colour=NULL
        ) +
        theme_report(
          leg.pos="bottom",
          legend.title=FALSE
        ) +
        ggplot2::guides(
          colour=ggplot2::guide_legend(
            ncol=2,
            byrow=TRUE,
            override.aes=list(
              linewidth=1.2
            ),
            keyheight=grid::unit(0.25,"cm"),
            keywidth=grid::unit(0.45,"cm")
          )
        )
    }
  }

  # ----------------------------------------------------------------------- #
  # Fallback status panel
  # ----------------------------------------------------------------------- #

  if (
    is.null(p5) &&
    show_class_panel
  ) {

    sd <- as.data.frame(
      table(
        factor(
          as.character(d$status),
          levels=c("neg","pos")
        )
      ),
      stringsAsFactors=FALSE
    )

    names(sd) <- c(
      "status",
      "n"
    )

    sd$pct <- 100*sd$n/sum(sd$n)

    sd$label <- paste0(
      format(
        sd$n,
        big.mark=","
      ),
      "\n",
      format_pct(sd$pct)
    )

    p5 <- ggplot2::ggplot(
      sd,
      ggplot2::aes(
        x=.data$status,
        y=.data$n,
        fill=.data$status
      )
    ) +
      ggplot2::geom_col(
        width=0.62,
        colour=tile_border,
        linewidth=0.3
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label=.data$label
        ),
        vjust=-0.25,
        size=txtsize/3.2,
        colour=fg
      ) +
      ggplot2::scale_fill_manual(
        values=status_cols,
        drop=FALSE
      ) +
      ggplot2::scale_y_continuous(
        labels=scales::label_comma(),
        expand=ggplot2::expansion(
          mult=c(0,0.16)
        )
      ) +
      ggplot2::labs(
        title="GnRH detection status",
        x=NULL,
        y="Cells"
      ) +
      theme_report(
        leg.pos="none"
      )
  }

  if (is.null(p5)) {
    p5 <- empty(
      "Detection-status panel unavailable"
    )
  }

  # ========================================================================= #
  # F. Marker-program support
  # ========================================================================= #

  m <- data.frame(
    status=rep(
      as.character(d$status),
      3
    ),
    module=rep(
      c(
        "Core",
        "Migration",
        "Neuroendocrine"
      ),
      each=nrow(d)
    ),
    hits=c(
      d$core_hits,
      d$mig_hits,
      d$neuro_hits
    )
  )

  ms <- stats::aggregate(
    hits~status+module,
    m,
    FUN=function(x)
      mean(
        x,
        na.rm=TRUE
      )
  )

  ms$module <- factor(
    ms$module,
    levels=c(
      "Core",
      "Migration",
      "Neuroendocrine"
    )
  )

  ms$status <- factor(
    ms$status,
    levels=c(
      "neg",
      "pos"
    )
  )

  # Text colour chosen according to fill intensity.
  hit_mid <- mean(
    range(
      ms$hits,
      finite=TRUE
    )
  )

  ms$text_col <- ifelse(
    ms$hits>hit_mid,
    if (mode=="dark") "#111111" else "white",
    fg
  )

  p6 <- ggplot2::ggplot(
    ms,
    ggplot2::aes(
      x=.data$module,
      y=.data$status,
      fill=.data$hits
    )
  ) +
    ggplot2::geom_tile(
      colour=tile_border,
      linewidth=0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label=sprintf(
          "%.2f",
          .data$hits
        ),
        colour=.data$text_col
      ),
      size=txtsize/3.0,
      fontface="bold",
      show.legend=FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient(
      low=tile_low,
      high=status_cols[["pos"]]
    ) +
    ggplot2::scale_x_discrete(
      labels=c(
        Core="Core",
        Migration="Migration",
        Neuroendocrine="Neuroendo."
      )
    ) +
    ggplot2::labs(
      title="Marker-program support",
      subtitle="Mean marker hits per cell",
      x=NULL,
      y=NULL,
      fill="Mean hits"
    ) +
    theme_report(
      leg.pos="bottom",
      leg.dir = "horizontal",
      x.ang=0
    ) +
    ggplot2::theme(
      axis.ticks=ggplot2::element_blank(),
      legend.key.height=grid::unit(0.65,"cm"),
      legend.key.width=grid::unit(0.20,"cm")
    ) +
    ggplot2::guides(
      fill=ggplot2::guide_colourbar(
        title.position="left",
        #title.hjust=0.5,
        barheight=grid::unit(0.5,"cm"),
        barwidth=grid::unit(5,"cm")
      )
    )

  # ========================================================================= #
  # Report header
  # ========================================================================= #

  n_pos <- sum(
    as.character(d$status)=="pos",
    na.rm=TRUE
  )

  n_conf <- if ("gnrh_confident" %in% names(d)) {

    x <- d$gnrh_confident

    if (is.logical(x)) {
      sum(x,na.rm=TRUE)
    } else {
      sum(
        as.character(x) %in%
          c(
            "TRUE",
            "true",
            "1",
            "yes",
            "Yes"
          ),
        na.rm=TRUE
      )
    }

  } else {
    NA_integer_
  }

  pct_pos <- 100*n_pos/nrow(d)

  subtitle <- paste0(
    format(
      ncol(object),
      big.mark=","
    ),
    " cells",
    "  ·  ",
    format(
      nrow(object),
      big.mark=","
    ),
    " features",
    "  ·  ",
    format(
      n_pos,
      big.mark=","
    ),
    " GnRH+ (",
    format_pct(pct_pos),
    ")",
    if (!is.na(n_conf)) {
      paste0(
        "  ·  ",
        format(
          n_conf,
          big.mark=","
        ),
        " confident"
      )
    } else {
      ""
    }
  )

  annotation_theme <- ggplot2::theme(
    plot.background=ggplot2::element_rect(
      fill=bg,
      colour=NA
    ),

    plot.title=ggplot2::element_text(
      face="bold",
      size=txtsize+5,
      colour=fg,
      hjust=0,
      margin=ggplot2::margin(b=2)
    ),

    plot.subtitle=ggplot2::element_text(
      size=txtsize,
      colour=muted,
      hjust=0,
      margin=ggplot2::margin(b=8)
    ),

    plot.tag=ggplot2::element_text(
      face="bold",
      size=txtsize+1,
      colour=fg
    )
  )

  # ========================================================================= #
  # Assemble
  # ========================================================================= #

  report <-
    (
      p1 |
        p2 |
        p3
    ) /
    (
      p4 |
        p5 |
        p6
    ) +
    patchwork::plot_layout(
      widths=c(1,1,0.92),
      heights=c(1,1),
      guides="keep"
    ) +
    patchwork::plot_annotation(
      title="gnrhcell diagnostic report",
      subtitle=subtitle,
      tag_levels="A",
      theme=annotation_theme
    )

  report
}
