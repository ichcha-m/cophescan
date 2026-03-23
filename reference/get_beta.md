# Extract beta and p-values of queried variant

Extract beta and p-values of queried variant

## Usage

``` r
get_beta(traits.dat, querysnpid, querytrait)
```

## Arguments

- traits.dat:

  list of coloc structured dataset

- querysnpid:

  vector of querysnpid

- querytrait:

  vector of querytrait names

## Value

data.frame with one column named beta_plot: indicating beta direction
(n/p) and another column named pval_plot with -log10(pval) of the
queried variant
