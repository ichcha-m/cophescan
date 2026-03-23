# Estimate the Hc.cutoff for the required FDR

Estimate the Hc.cutoff for the required FDR

## Usage

``` r
Hc.cutoff.fdr(ppHc, fdr = 0.05, return_plot = TRUE)
```

## Arguments

- ppHc:

  a vector containing the PP.Hc (the posterior probability of causal
  association) of all tests

- fdr:

  FDR default: 0.05

- return_plot:

  default: TRUE, plot the fdr estimated at the different Hc.cutoff

## Value

the Hc.cutoff value for the specified FDR, if return_plot is True
returns a plot showing the FDR calculated at different Hc thresholds
