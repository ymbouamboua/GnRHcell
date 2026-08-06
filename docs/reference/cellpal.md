# Generate a discrete color palette

This function returns a vector of colors from predefined palettes,
thematic palettes (viridis, plasma, etc.), or a user-provided set of
base colors. It can oversample, remove grayish colors, adjust lightness
and saturation, and reverse the palette order.

## Usage

``` r
cellpal(
  n = 5,
  preset = "base",
  theme = NULL,
  base_colors = NULL,
  space = c("Lab", "rgb", "HCL"),
  oversample_factor = 1.3,
  remove_gray = TRUE,
  reverse = FALSE,
  adjust_saturation = 1,
  adjust_lightness = 1
)
```

## Arguments

- n:

  Integer. Number of colors to return.

- preset:

  Character. Name of a predefined palette. Default NULL.

- theme:

  Character. Thematic palettes: "viridis", "magma", "plasma", "inferno",
  "cividis", or any RColorBrewer palette name.

- base_colors:

  Character vector. User-defined colors. Overrides preset/theme if
  provided.

- space:

  Character. Color interpolation space: "Lab", "rgb", or "HCL".

- oversample_factor:

  Numeric. How much to oversample before picking final n colors.

- remove_gray:

  Logical. Remove grayish colors when oversampling.

- reverse:

  Logical. Reverse the color order.

- adjust_saturation:

  Numeric. Factor to adjust saturation (1 = no change).

- adjust_lightness:

  Numeric. Factor to adjust lightness (1 = no change).

## Value

Character vector of colors.
