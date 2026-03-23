# Plot region Manhattan for a trait highlighting the queried variant

Plot region Manhattan for a trait highlighting the queried variant

## Usage

``` r
plot_trait_manhat(trait.dat, querysnpid, alt.snpid = NULL)
```

## Arguments

- trait.dat:

  dataset used as input for running cophescan

- querysnpid:

  the id of the causal variant as present in trait.dat\$snp, plotted in
  red

- alt.snpid:

  the id of the other variants as a vector to be plotted, plotted in
  blue

## Value

regional manhattan plot
