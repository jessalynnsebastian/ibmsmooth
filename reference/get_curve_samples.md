# Extract posterior curve samples

Extract posterior curve samples

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

  An `ibmfit` object.

- param:

  One of `"f"`, `"fprime"`, or `"both"`.

- n_samples:

  Maximum number of draws.

- format:

  Either `"matrix"` or `"long"`.

## Value

A matrix, list of matrices, or long data frame.
