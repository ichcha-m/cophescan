# Prepare data for plotting

Prepare data for plotting

## Usage

``` r
prepare_plot_data(
  multi.dat,
  querysnpid,
  query_trait_names,
  thresh_Ha = 0.5,
  thresh_Hc = 0.5,
  hmp = FALSE,
  cophe.plot = TRUE
)
```

## Arguments

- multi.dat:

  multi trait cophescan results returned from cophe.multitrait or
  multitrait.simplify

- querysnpid:

  query variant

- query_trait_names:

  vector of names of the query traits

- thresh_Ha:

  Ha threshold to be displayed

- thresh_Hc:

  Hc threshold to be displayed

- hmp:

  return for heatmap

- cophe.plot:

  default: TRUE, return for `cophe_plot`

## Value

plot list

## See also

[`cophe_plot`](https://ichcha-m.github.io/cophescan/reference/cophe_plot.md),
[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md),
[`cophe.multitrait`](https://ichcha-m.github.io/cophescan/reference/cophe.multitrait.md),
[`multitrait.simplify`](https://ichcha-m.github.io/cophescan/reference/multitrait.simplify.md)
default NULL
