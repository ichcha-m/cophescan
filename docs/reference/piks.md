# List of priors: pn, pa and pc over all iterations

List of priors: pn, pa and pc over all iterations

## Usage

``` r
piks(params, nsnps, covar_vec, covar = FALSE)
```

## Arguments

- params:

  Vector of parameters: \\\alpha\\, \\\beta\\ and \\\gamma\\

- nsnps:

  number of snps

- covar_vec:

  Vector of the covariate

- covar:

  logical: was the covariate information used? default: False

## Value

List of priors (len: iterations): pnk, pak and pck
