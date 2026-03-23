# Average of posterior probabilities: Hn, Ha and Hc

Average of posterior probabilities: Hn, Ha and Hc

## Usage

``` r
average_posterior_prob(
  params,
  lbf_mat,
  nsnps,
  covar_vec,
  nits,
  thin,
  covar = FALSE
)
```

## Arguments

- params:

  Vector of parameters: \\\alpha\\, \\\beta\\ and \\\gamma\\

- lbf_mat:

  matrix of log bayes factors: lBF.Ha and lBF.Hc

- nsnps:

  number of snps

- covar_vec:

  Vector of the covariate

- nits:

  Number of iterations run in mcmc

- thin:

  thinning

- covar:

  logical: was the covariate information used? default: False

## Value

matrix with average of all the posterior probabilities: Hn, Ha and Hc
