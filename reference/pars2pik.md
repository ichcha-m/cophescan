# Conversion of parameters alpha, beta and gamma to pnk, pak and pck

Conversion of parameters alpha, beta and gamma to pnk, pak and pck

## Usage

``` r
pars2pik(params, nsnps, covar_vec, covar = FALSE)
```

## Arguments

- params:

  Vector of parameters: \\\alpha\\, \\\beta\\ and \\\gamma\\

- nsnps:

  number of snps

- covar_vec:

  Vector of the covariate

- covar:

  logical: should the covariate information be used? default: False

## Value

pik matrix of priors: pnk, pak and pck
