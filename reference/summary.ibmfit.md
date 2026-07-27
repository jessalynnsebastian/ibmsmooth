# Summarize an IBM fit

Summarize an IBM fit

## Usage

``` r
# S3 method for class 'ibmfit'
summary(object, ..., n_samples = 1000, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- object:

  An `ibmfit` object.

- ...:

  Unused.

- n_samples:

  Maximum number of draws.

- probs:

  Posterior probabilities.

## Value

An object of class `summary_ibmfit`.
