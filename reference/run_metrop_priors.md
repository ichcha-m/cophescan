# Run the hierarchical Metropolis Hastings model to infer priors

Run the hierarchical Metropolis Hastings model to infer priors

## Usage

``` r
run_metrop_priors(
  multi.dat,
  covar = FALSE,
  covar_vec = NULL,
  is_covar_categorical = FALSE,
  nits = 10000,
  thin = 1,
  posterior = FALSE,
  avg_pik = TRUE,
  avg_posterior = TRUE,
  pik = FALSE,
  alpha_mean = -10,
  alpha_sd = 0.5,
  beta_shape = 2,
  beta_scale = 2,
  gamma_shape = 2,
  gamma_scale = 2
)
```

## Arguments

- multi.dat:

  matrix of bf values, rows=traits, named
  columns=("lBF.Ha","lBF.Hc","nsnps")

- covar:

  whether to include covariates

- covar_vec:

  vector of covariates

- is_covar_categorical:

  only two categories supported (default=FALSE) - Experimental

- nits:

  number of iterations

- thin:

  burnin

- posterior:

  default: FALSE, estimate posterior probabilities of the hypotheses

- avg_pik:

  default: FALSE, estimate the average of the pik

- avg_posterior:

  default: FALSE, estimate the average of the posterior probabilities of
  the hypotheses

- pik:

  default: FALSE, inferred prior probabilities

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

List containing the posterior distribution of the parameters alpha,
beta, gamma (if covariate included) and the loglikelihood

if avg_posterior=TRUE matrix with average of all the posterior
probabilities of Hn, Ha and Hc

if avg_pik=TRUE matrix with average of all the priors: pn, pa and pc

data, nits and thin contain the input data, number of iterations and
burnin respectively specified for the hierarchical model
