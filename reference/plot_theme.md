# GnRHcell ggplot2 theme

Flexible ggplot2 theme used across GnRHcell visualizations.

## Usage

``` r
plot_theme(
  style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
  txtsize = 12,
  xy.val = TRUE,
  x.ang = 0,
  hjust = NULL,
  vjust = NULL,
  xlab = TRUE,
  ylab = TRUE,
  xy.lab = TRUE,
  facet.face = "bold",
  ttl.face = "bold",
  txt.face = c("plain", "italic", "bold"),
  ttl.pos = c("center", "left", "right"),
  x.ttl = TRUE,
  y.ttl = TRUE,
  ticks = NULL,
  line = NULL,
  border = NULL,
  grid.major = NULL,
  grid.minor = NULL,
  panel.fill = "white",
  facet.bg = TRUE,
  mode = c("light", "dark"),
  leg.pos = "right",
  leg.dir = "vertical",
  leg.size = 10,
  leg.ttl = 10,
  leg.ttl.size = 10,
  leg.just = "center",
  leg.ttl.text = NULL,
  ...
)
```

## Arguments

- style:

  Theme style. One of `"classic"`, `"minimal"`, `"bw"`, `"test"`,
  `"void"`, `"dirty"`, or `"gray"`.

- txtsize:

  Base text size.

- xy.val:

  Show axis tick labels.

- x.ang:

  X-axis text angle.

- hjust, vjust:

  Horizontal and vertical justification for x-axis labels.

- xlab, ylab:

  Show x/y axis tick labels.

- xy.lab:

  Show all axis tick labels.

- facet.face:

  Facet label font face.

- ttl.face:

  Plot title font face.

- txt.face:

  Text font face.

- ttl.pos:

  Plot title position.

- x.ttl, y.ttl:

  Show x/y axis titles.

- ticks, line, border:

  Logical overrides for ticks, axis lines, and panel border.

- grid.major, grid.minor:

  Logical overrides for major/minor grid lines.

- panel.fill:

  Panel background fill.

- facet.bg:

  Show facet background.

- mode:

  Theme mode: `"light"` or `"dark"`.

- leg.pos:

  Legend position.

- leg.dir:

  Legend direction.

- leg.size:

  Legend text size.

- leg.ttl:

  Legend title size.

- leg.ttl.size:

  Legend title text size.

- leg.just:

  Legend justification.

- leg.ttl.text:

  Optional legend title text.

- ...:

  Additional arguments passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A ggplot2 theme object.
