# Plot an integrated Brownian motion fit

Plot an integrated Brownian motion fit

## Usage

``` r
# S3 method for class 'ibmfit'
plot(x, truth = NULL, ...)
```

## Arguments

- x:

  An `ibmfit` object.

- truth:

  Optional list with `f` and/or `fprime` values on the latent grid.

- ...:

  Additional arguments passed to
  [`plot_curve()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/plot_curve.md).

## Value

A list containing function and derivative ggplots.
