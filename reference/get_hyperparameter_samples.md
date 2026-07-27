# Extract posterior hyperparameter samples

Extract posterior hyperparameter samples

## Usage

``` r
get_hyperparameter_samples(
  ibmfit,
  n_samples = 1000,
  format = c("long", "wide")
)
```

## Arguments

- ibmfit:

  An `ibmfit` object.

- n_samples:

  Maximum number of draws.

- format:

  Either `"long"` or `"wide"`.

## Value

A data frame.
