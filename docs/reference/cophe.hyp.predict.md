# Predict cophescan hypothesis for tested associations

Predict cophescan hypothesis for tested associations

## Usage

``` r
cophe.hyp.predict(
  cophe.res,
  grouping.vars = c("querysnp", "querytrait"),
  Hc.cutoff = 0.6,
  Hn.cutoff = 0.2
)
```

## Arguments

- cophe.res:

  results obtained from `cophe.single`, `cophe.susie` or
  `cophe.multitrait` or data.frame with the following columns: PP.Hn,
  PP.Hc, PP.Ha, querysnp, querytrait

- grouping.vars:

  This is important for results from `cophe.susie` where there are
  multiple signals. These will be collapsed into one call. If you want
  to return all signals set this to a single variable eg: grouping.vars
  = c('querysnp')

- Hc.cutoff:

  threshold for PP.Hc above which the associations are called Hc

- Hn.cutoff:

  threshold for PP.Hn above which the associations are called Hn

## Value

returns dataframe with posterior probabilties of Hn, Hc and Ha with the
predicted hypothesis based on the provided cut.offs.

## See also

[`cophe.single`](https://ichcha-m.github.io/cophescan/reference/cophe.single.md),
[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md),
[`cophe.multitrait`](https://ichcha-m.github.io/cophescan/reference/cophe.multitrait.md),
,
[`multitrait.simplify`](https://ichcha-m.github.io/cophescan/reference/multitrait.simplify.md)
