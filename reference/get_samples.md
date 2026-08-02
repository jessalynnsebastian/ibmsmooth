# Extract posterior samples

Extract posterior samples

## Usage

``` r
get_samples(ibmfit, param = c("f", "fprime"), n_samples = NULL)
```

## Arguments

- ibmfit:

  An `ibmfit` object.

- param:

  Stan parameter name. Common choices are `"f"`, `"fprime"`, `"sigma"`,
  `"phi"`, `"tau"`, `"gamma"`, `"xi"`, `"xi_regularized"`, `"slab"`, and
  `"lambda_interval"`.

- n_samples:

  Maximum number of draws; `NULL` keeps all draws.

## Value

A matrix with posterior draws in rows.
