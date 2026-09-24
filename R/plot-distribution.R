# ========================================================================= #
# Distribution plots
# ========================================================================= #

#' Plot gnrhcell metadata distributions
#'
#' @param object A Seurat object.
#' @param group.by Metadata column defining categories.
#' @param split.by Optional metadata column defining bars.
#' @param cols Optional named colour vector.
#' @param proportion Display within-split proportions instead of counts.
#' @param position Bar position, such as `"stack"` or `"dodge"`.
#' @param label Add value labels.
#' @param label.size Label text size.
#' @param plot.ttl Optional plot title.
#' @param txtsize Base text size.
#' @param x.ang X-axis label angle. When `NULL`, the angle is selected from
#' @param style Theme style: `"classic"`, `"minimal"`, `"bw"`, or `"test"`.
#'   the number and length of sample labels.
#' @param flip Flip coordinates.
#' @param adaptive Automatically use compact spacing, an economical legend
#'   layout, and recommended export dimensions based on the number of samples.
#' @param bar.width Bar width. When `NULL`, it is selected automatically.
#' @param bar.gap Gap between adjacent bar edges, in x-axis units. It is
#'   independent of `bar.width`; for example, `bar.width = 0.3` and
#'   `bar.gap = 0.2` produce centres 0.5 units apart. When `NULL`, an adaptive
#'   value is used.
#' @param legend.position Legend position. Use `"auto"` to place it according
#'   to the number of samples, or a standard ggplot2 legend position.
#' @return A ggplot object.
#' @export
gnrh_celldistribution <- function(
    object,
    group.by,
    split.by = NULL,
    cols = NULL,
    proportion = FALSE,
    position = "stack",
    label = TRUE,
    label.size = 3,
    plot.ttl = NULL,
    txtsize = getOption("gnrhcell.base_size", 14),
    x.ang = NULL,
    style = "classic",
    flip = FALSE,
    adaptive = TRUE,
    bar.width = NULL,
    bar.gap = NULL,
    legend.position = "auto"
) {
  .validate_seurat(object)

  style = match.arg(style)

  md <- object[[]]

  if (!group.by %in% names(md))
    stop("`group.by` not found.", call. = FALSE)

  if (!is.null(split.by) && !split.by %in% names(md))
    stop("`split.by` not found.", call. = FALSE)

  if (is.null(split.by)) {
    df <- as.data.frame(
      table(md[[group.by]]),
      stringsAsFactors = FALSE
    )
    names(df) <- c("group", "n")
    df$split <- "all"
  } else {
    df <- as.data.frame(
      table(md[[split.by]], md[[group.by]]),
      stringsAsFactors = FALSE
    )
    names(df) <- c("split", "group", "n")
  }

  if (proportion) {
    df <- df |>
      dplyr::group_by(.data$split) |>
      dplyr::mutate(value = .data$n / sum(.data$n)) |>
      dplyr::ungroup()
  } else {
    df$value <- df$n
  }

  df$label <- if (proportion)
    scales::percent(df$value, accuracy = 0.1)
  else
    format(df$n, big.mark = ",")

  groups <- unique(as.character(df$group))
  cols <- .resolve_colors(groups, cols, group.by)

  sample_labels <- unique(as.character(df$split))
  n_samples <- length(sample_labels)
  longest_label <- max(nchar(sample_labels), 1L)

  if (is.null(bar.width)) {
    bar.width <- if (isTRUE(adaptive)) {
      if (n_samples <= 4L) 0.58 else if (n_samples <= 10L) 0.68 else 0.78
    } else {
      0.75
    }
  }

  if (is.null(bar.gap)) bar.gap <- if (isTRUE(adaptive)) 0.12 else 0.25
  if (length(bar.gap) != 1L || !is.finite(bar.gap) || bar.gap < 0) {
    stop("`bar.gap` must be one non-negative number.", call. = FALSE)
  }

  split_levels <- unique(as.character(df$split))
  spacing <- bar.width + bar.gap
  split_positions <- stats::setNames(
    (seq_along(split_levels) - 1) * spacing + 1,
    split_levels
  )
  df$.split_position <- unname(split_positions[as.character(df$split)])
  outer_gap <- bar.gap / 2
  x_limits <- c(
    min(split_positions) - bar.width / 2 - outer_gap,
    max(split_positions) + bar.width / 2 + outer_gap
  )

  if (is.null(x.ang)) {
    x.ang <- if (isTRUE(flip)) {
      0
    } else if (n_samples <= 4L && longest_label <= 8L) {
      0
    } else if (n_samples <= 10L && longest_label <= 15L) {
      30
    } else {
      60
    }
  }

  if (identical(legend.position, "auto")) {
    legend.position <- if (is.null(split.by)) {
      "none"
    } else if (isTRUE(adaptive) && n_samples <= 8L) {
      "bottom"
    } else {
      "right"
    }
  }

  plot_width <- if (isTRUE(flip)) {
    max(5.5, min(9, 5 + longest_label / 10))
  } else {
    max(4.8, min(14, 3.4 + 0.62 * n_samples + longest_label / 30))
  }
  plot_height <- if (isTRUE(flip)) {
    max(4.2, min(12, 2.8 + 0.42 * n_samples))
  } else {
    if (identical(legend.position, "bottom")) 5.2 else 4.8
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$.split_position, .data$value, fill = .data$group)
  ) +
    ggplot2::geom_col(
      position = position,
      width = bar.width,
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_x_continuous(
      breaks = unname(split_positions),
      labels = names(split_positions),
      limits = x_limits,
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::labs(
      title = plot.ttl,
      x = if (is.null(split.by)) NULL else split.by,
      y = if (proportion) "Proportion" else "Cells",
      fill = group.by
    )

  if (isTRUE(proportion)) {
    p <- p + ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = function(x) ifelse(x == 1, "100", sprintf("%.2f", x)),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.015))
    )
  }

  if (label) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        position = if (position == "stack")
          ggplot2::position_stack(vjust = 0.5)
        else
          ggplot2::position_dodge(width = bar.width),
        size = label.size
      )
  }

  if (flip) p <- p + ggplot2::coord_flip()

  p <- p +
    gnrh_theme(
      txtsize = txtsize,
      x.ang = x.ang,
      leg.pos = legend.position,
      style = style
    ) +
    ggplot2::theme(
      legend.direction = if (identical(legend.position, "bottom")) "horizontal" else "vertical",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  attr(p, "recommended_size") <- c(width = plot_width, height = plot_height)
  attr(p, "n_samples") <- n_samples
  p
}


# ========================================================================= #
# Module hits
# ========================================================================= #

#' Plot GnRH module hit distributions
#'
#' @param data Data frame containing plotting columns.
#' @param x Column mapped to the x axis.
#' @param fill Column mapped to fill colour.
#' @param palette Optional named colour vector.
#' @param type Display counts or fractions.
#' @param title Optional plot title.
#' @param txtsize Base text size.
#' @return A ggplot object.
#' @export
gnrh_hits <- function(
    data,
    x,
    fill,
    palette = NULL,
    type = c("count", "fraction"),
    title = NULL,
    txtsize = 10
) {
  type <- match.arg(type)

  if (!is.data.frame(data) || !all(c(x, fill) %in% names(data)))
    stop("Invalid input data or metadata columns.", call. = FALSE)

  groups <- unique(as.character(data[[fill]]))

  if (is.null(palette))
    palette <- .resolve_colors(groups, group_by = fill)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[x]], fill = .data[[fill]])
  ) +
    ggplot2::geom_bar(
      position = if (type == "count") "stack" else "fill",
      colour = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(
      title = title,
      x = x,
      y = if (type == "count") "Cell count" else "Fraction",
      fill = fill
    )

  if (type == "fraction")
    p <- p +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format()
    )

  p + gnrh_theme(txtsize = txtsize, x.ang = 45)
}
