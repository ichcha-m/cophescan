# Average of priors: pnk, pak and pck

Average of priors: pnk, pak and pck

## Usage

``` r
average_piks(params, nsnps, covar_vec, nits, thin, covar = FALSE)
```

## Arguments

- params:

  Vector of parameters: \\\alpha\\, \\\beta\\ and \\\gamma\\

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

average pik matrix of priors: pnk, pak and pck
