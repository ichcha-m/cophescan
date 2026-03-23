# extract data through Bayes factors

a dataset represented by Bayes factors from SuSIE

## Usage

``` r
cophe.bf_bf(
  sus_dat,
  cred_set,
  querysnpid,
  querytrait,
  pn = NULL,
  pa = NULL,
  pc = NULL,
  ret_pp = TRUE
)
```

## Arguments

- sus_dat:

  a list with the output of running susie

- cred_set:

  credible set extracted from susie

- querysnpid:

  Id of the query variant

- querytrait:

  Query trait name

- pn:

  prior probability that none of the SNPs/variants in the region are
  associated with the query trait

- pa:

  prior probability that a non-query variant is causally associated with
  the query trait

- pc:

  prior probability that the query variant is causally associated with
  the query trait

## Value

bayes factors of signals

## See also

[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md)
