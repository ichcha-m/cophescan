# cophe.susie.lbf

Calculate log bayes factors for each hypothesis (SuSIE - multiple causal
variant assumption)

## Usage

``` r
cophe.susie.lbf(
  dataset,
  querysnpid,
  querytrait,
  switch = TRUE,
  susie.args = list(),
  MAF = NULL
)
```

## Arguments

- dataset:

  a list with specifically named elements defining the query trait
  dataset to be analysed.

- querysnpid:

  Id of the query variant, (id in dataset\$snp)

- querytrait:

  Query trait name

- switch:

  Set switch=TRUE to obtain single BF when credible sets not found with
  SuSIE

- susie.args:

  a named list of additional arguments to be passed to
  [runsusie](https://rdrr.io/pkg/coloc/man/runsusie.html)

- MAF:

  Minor allele frequency vector

## Value

data frame with log bayes factors for Hn and Ha hypotheses

## See also

[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md)

## Author

Ichcha Manipur

## Examples

``` r
library(cophescan)
data(cophe_multi_trait_data)
query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
query_trait_1$LD <- cophe_multi_trait_data$LD
querysnpid <- cophe_multi_trait_data$querysnpid
res.susie.lbf <- cophe.susie.lbf(query_trait_1, querysnpid = querysnpid,
                                  querytrait='Trait_1', switch=T)
#> running max iterations: 100
#>  converged: TRUE
res.susie.lbf
#>      lBF.Ha   lBF.Hc nsnps       querysnp querytrait           hit1
#>       <num>    <num> <int>         <char>     <char>         <char>
#> 1: 15.31003 11.95277  1000 chr19-11182353    Trait_1 chr19-11182353
#>              hit2  typeBF  idx1  idx2
#>            <char>  <char> <int> <int>
#> 1: chr19-11182144 susieBF     1     1
```
