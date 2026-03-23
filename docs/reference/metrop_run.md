# Run the hierarchical mcmc model to infer priors

Run the hierarchical mcmc model to infer priors

## Usage

``` r
metrop_run(
  lbf_mat,
  nsnps,
  covar_vec,
  covar = FALSE,
  nits = 10000L,
  thin = 1L,
  alpha_mean = -10,
  alpha_sd = 0.5,
  beta_shape = 2,
  beta_scale = 2,
  gamma_shape = 2,
  gamma_scale = 2
)
```

## Arguments

- lbf_mat:

  matrix of log bayes factors: lBF.Ha and lBF.Hc

- nsnps:

  number of snps

- covar_vec:

  Vector of the covariate

- covar:

  logical: Should the covariate information be used? default: False

- nits:

  Number of iterations run in mcmc

- thin:

  thinning

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

named list of log likelihood (ll) and parameters: alpha, beta and gamma
