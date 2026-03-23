# Heatmap of multi-trait cophescan results

Heatmap of multi-trait cophescan results

## Usage

``` r
cophe_heatmap(
  multi.dat,
  querysnpid,
  query_trait_names,
  thresh_Hc = 0.5,
  thresh_Ha = 0.5,
  ...
)
```

## Arguments

- multi.dat:

  multi trait cophescan results returned from `cophe.multitrait` or
  formatted in the same way with multitrait.simplify

- querysnpid:

  query variant

- query_trait_names:

  names of phenotypes corresponding to the multi.dat results

- thresh_Hc:

  Hc threshold to be displayed

- thresh_Ha:

  Ha threshold to be displayed

- ...:

  additional arguments to be passed to
  [pheatmap](https://rdrr.io/pkg/pheatmap/man/pheatmap.html)

## Value

heatmap of posterior probabilities of the phenotypes above the set
threshold
