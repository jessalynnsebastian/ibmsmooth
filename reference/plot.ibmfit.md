# Plot method for ibmfit objects

Visualize fitted results from an object of class "ibmfit". The plot
typically shows model predictions with uncertainty and, if supplied,
observed truth values for comparison.

## Usage

``` r
# S3 method for class 'ibmfit'
plot(x, truth = NULL, titles = NULL, ...)
```

## Arguments

- x:

  An object of class "ibmfit", produced by fitting the IBM smoothing
  model.

- truth:

  Optional vector, data.frame, or time series of true/observed values to
  overlay on the plot. If provided, it will be used to compare fitted
  values to ground truth.

- titles:

  Optional character vector of length one or two giving plot titles for
  the function and derivative plots.

- ...:

  Additional graphical parameters passed to underlying plotting
  functions (e.g., type, col, lwd) or other method-specific options
  recognized by plot.ibmfit.

## Value

A list with components:

- function_plot:

  A ggplot object showing the estimated function and its uncertainty.

- derivative_plot:

  A ggplot object showing the estimated derivative and its uncertainty.

## Details

This S3 method produces diagnostic and summary plots for an "ibmfit"
object, such as fitted trajectories, uncertainty intervals, and residual
summaries. When `truth` is supplied, observed values are added to the
plot for visual comparison. The exact panels and layout depend on the
contents of the `ibmfit` object.
