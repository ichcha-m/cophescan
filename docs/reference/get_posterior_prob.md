# Calculation of the posterior prob of Hn, Ha and Hc

Calculation of the posterior prob of Hn, Ha and Hc

## Usage

``` r
get_posterior_prob(params, lbf_mat, nsnps, covar_vec, covar = FALSE)
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

- covar:

  logical: should the covariate information be used? default: False

## Value

posterior prob of Hn, Ha and Hc
