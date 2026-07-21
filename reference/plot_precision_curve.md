# Plot the locally adaptive process precision curve

Plot posterior summaries of the time-varying IBM process precision for a
locally adaptive fit. The curve is evaluated at transition midpoints, so
each plotted time value represents the interval between two adjacent
latent grid points.

## Usage

``` r
plot_precision_curve(
  ibmfit,
  n_samples = 1000,
  level = 0.95,
  inner_level = 0.8,
  log_y = TRUE,
  title = NULL,
  line_color = "purple4",
  ribbon_fill = "purple",
  alpha_outer = 0.18,
  alpha_inner = 0.3,
  ...
)
```

## Arguments

- ibmfit:

  An object of class `ibmfit` from a locally adaptive Stan or INLA fit.

- n_samples:

  Number of posterior samples to use. For INLA fits these are drawn with
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
  For Stan fits this is the maximum number of MCMC draws to keep.

- level:

  Outer pointwise credible interval level. Defaults to `0.95`.

- inner_level:

  Optional inner pointwise credible interval level. Defaults to `0.80`.
  Use `NULL` to suppress the inner interval.

- log_y:

  Logical. If `TRUE`, display the positive precision values on a log10
  y-axis.

- title:

  Optional plot title.

- line_color:

  Color for the posterior median curve.

- ribbon_fill:

  Fill color for the credible interval ribbons.

- alpha_outer:

  Alpha for the outer credible interval ribbon.

- alpha_inner:

  Alpha for the inner credible interval ribbon.

- ...:

  Additional arguments, currently unused.

## Value

A `ggplot` object. The data used to construct the plot are attached as
the `precision_summary` attribute.
