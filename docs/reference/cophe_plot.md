# cophe_plots showing the Ha and Hc of all traits and labelled above the specified threshold

cophe_plots showing the Ha and Hc of all traits and labelled above the
specified threshold

## Usage

``` r
cophe_plot(
  multi.dat,
  querysnpid,
  query_trait_names,
  thresh_Hc = 0.5,
  thresh_Ha = 0.5,
  beta_p = NULL,
  traits.dat = NULL
)
```

## Arguments

- multi.dat:

  multi trait cophescan results returned from cophe.multitrait or
  multitrait.simplify

- querysnpid:

  query variant (only a single variant for PheWAS plots)

- query_trait_names:

  list of phenotype names

- thresh_Hc:

  Hc threshold to be displayed

- thresh_Ha:

  Ha threshold to be displayed

- beta_p:

  data.frame (from the `get.beta` function) with four columns : 1.
  "beta_plot": indicating beta direction (p or n) 2. "pval_plot":
  -log10(pval) of the queried variant 3. "querysnp" 4. "querytrait".

- traits.dat:

  list of multi-trait coloc structured datasets

## Value

cophescan plots of Ha and Hc

## See also

[`cophe.single`](https://ichcha-m.github.io/cophescan/reference/cophe.single.md),
[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md),
[`cophe.multitrait`](https://ichcha-m.github.io/cophescan/reference/cophe.multitrait.md),
,
[`multitrait.simplify`](https://ichcha-m.github.io/cophescan/reference/multitrait.simplify.md)
