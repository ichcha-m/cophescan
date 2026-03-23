# combine.bf

Calculate posterior probabilities for all the configurations

## Usage

``` r
combine.bf(lBF_df, pn, pa, pc)
```

## Arguments

- lBF_df:

  dataframe with log bayes factors of hypothesis Ha and Hn: column names
  should be lBF.Ha and lBF.Hc

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

named numeric vector of posterior probabilities and bayes factors

## Author

Ichcha Manipur
