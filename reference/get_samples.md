# Extract posterior samples from an IBM fit

Extract posterior samples from an IBM fit

## Usage

``` r
get_samples(ibmfit, param = c("f", "fprime"), n_samples = NULL, ...)
```

## Arguments

- ibmfit:

  An object of class `ibmfit`.

- param:

  Character, either `"f"` or `"fprime"`.

- n_samples:

  Optional number of posterior draws to return.

- ...:

  Additional arguments passed to
  [`get_samples_inla()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_samples_inla.md)
  for INLA fits.

## Value

A matrix of posterior samples.
