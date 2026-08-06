# Dark ggplot2 theme

Dark theme presets for presentation-quality GnRHcell plots.

## Usage

``` r
dark_theme(
  style = c("soft", "true", "paper"),
  txtsize = 12,
  family = "Helvetica",
  axis.text = TRUE,
  axis.title = TRUE
)
```

## Arguments

- style:

  Dark theme style: `"soft"`, `"true"`, or `"paper"`.

- txtsize:

  Base text size.

- family:

  Font family.

- axis.text:

  Logical; show axis tick labels.

- axis.title:

  Logical; show axis titles.

## Value

A ggplot2 theme object.
