# Extract posterior curve samples

Extract posterior samples for the latent function and/or derivative over
the full latent grid.

## Usage

``` r
get_curve_samples(
  ibmfit,
  param = c("f", "fprime", "both"),
  n_samples = 1000,
  format = c("matrix", "long")
)
```

## Arguments

- ibmfit:

  An object of class `ibmfit`.

- param:

  Character. One of `"f"`, `"fprime"`, or `"both"`.

- n_samples:

  Number of samples to use for INLA fits, or maximum number of MCMC
  draws to keep for Stan fits. Use `NULL` to keep all Stan draws.

- format:

  Output format, either `"matrix"` or `"long"`.

## Value

If `format = "matrix"` and one parameter is requested, a posterior
sample matrix with rows as draws and columns as time points. If both
parameters are requested, a list of two matrices. If `format = "long"`,
a data frame with columns `draw`, `t`, `parameter`, and `value`.
