# Log posterior calculation

Log posterior calculation

## Usage

``` r
logpost(params, lbf_mat, nsnps, covar_vec, covar = FALSE)
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

logpost log of the posteriors
