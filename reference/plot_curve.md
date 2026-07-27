# Plot posterior function and derivative summaries

Plot posterior function and derivative summaries

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

  Latent-grid locations.

- f_samples:

  Posterior function draws, with draws in rows.

- fprime_samples:

  Optional posterior derivative draws.

- dat_orig:

  Optional data frame with `t` and `y`.

- truth:

  Optional list with `f` and/or `fprime`.

- titles:

  Optional two plot titles.

- ...:

  Unused.

## Value

A list of ggplot objects.
