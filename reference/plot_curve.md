# Plot posterior function and derivative summaries

Construct ggplot objects showing posterior summaries for the latent
function and, optionally, its derivative.

## Usage

``` r
plot_curve(
  t_unique,
  f_samples,
  fprime_samples = NULL,
  dat_orig = NULL,
  truth = NULL,
  titles = NULL,
  ...
)
```

## Arguments

- t_unique:

  Numeric vector of time points on the latent grid.

- f_samples:

  Matrix of posterior samples for the latent function, with rows
  representing posterior draws and columns representing time points.

- fprime_samples:

  Optional matrix of posterior samples for the derivative, with rows
  representing posterior draws and columns representing time points.

- dat_orig:

  Optional data frame containing observed data with columns `t` and `y`.

- truth:

  Optional vector, list, or data frame containing truth values to
  overlay. Lists or data frames may contain `truth` and `deriv`.

- titles:

  Optional character vector of length one or two giving titles for the
  function and derivative plots.

- ...:

  Additional arguments, currently unused.

## Value

A list containing one or two ggplot objects.
