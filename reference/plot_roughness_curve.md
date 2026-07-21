# Plot baseline-horseshoe IBM roughness components

Plot posterior pointwise summaries of either total local roughness or
the horseshoe excess roughness. Values are process standard deviations
on the original response/time scale and are evaluated at transition
midpoints.

## Usage

``` r
plot_roughness_curve(
  ibmfit,
  component = c("total", "excess"),
  n_samples = 1000,
  level = 0.95,
  log_y = TRUE,
  title = NULL,
  line_color = "purple4",
  ribbon_fill = "purple",
  ...
)
```

## Arguments

- ibmfit:

  A baseline-horseshoe `ibmfit` object.

- component:

  Either `"total"` or `"excess"`.

- n_samples:

  Maximum number of posterior draws to use.

- level:

  Pointwise credible interval level.

- log_y:

  Logical; use a log10 y-axis.

- title:

  Optional plot title.

- line_color, ribbon_fill:

  Plot colors.

- ...:

  Additional arguments, currently unused.

## Value

A `ggplot` object with its summary data in the `roughness_summary`
attribute.
