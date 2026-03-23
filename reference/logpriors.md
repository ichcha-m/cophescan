# Calculate log priors

Calculate log priors

## Usage

``` r
logpriors(
  params,
  covar = FALSE,
  alpha_mean = -10,
  alpha_sd = 0.5,
  beta_shape = 2,
  beta_scale = 2,
  gamma_shape = 2,
  gamma_scale = 2
)
```

## Arguments

- params:

  Vector of parameters: \\\alpha\\, \\\beta\\ and \\\gamma\\

- covar:

  logical: Should the covariate information be used? default: False

- alpha_mean:

  prior for the mean of alpha

- alpha_sd:

  prior for the standard deviation of alpha

- beta_shape:

  prior for the shape (gamma distribution) of beta

- beta_scale:

  prior for the scale of beta

- gamma_shape:

  prior for the shape (gamma distribution) of gamma

- gamma_scale:

  prior for the scale of gamma

## Value

log priors
