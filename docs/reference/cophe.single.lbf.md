# cophe.single.lbf

Calculate log bayes factors for each hypothesis (Single causal variant
assumption)

## Usage

``` r
cophe.single.lbf(dataset, querysnpid, querytrait, MAF = NULL)
```

## Arguments

- dataset:

  a list with specifically named elements defining the query trait
  dataset to be analysed.

- querysnpid:

  Id of the query variant, (id in dataset\$snp)

- querytrait:

  Query trait name

- MAF:

  Minor allele frequency vector

## Value

data frame with log bayes factors for Hn and Ha hypotheses

## See also

[`cophe.single`](https://ichcha-m.github.io/cophescan/reference/cophe.single.md)

## Author

Ichcha Manipur

## Examples

``` r
library(cophescan)
data(cophe_multi_trait_data)
query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
querysnpid <- cophe_multi_trait_data$querysnpid
res.single.lbf <- cophe.single.lbf(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
res.single.lbf
#>     lBF.Ha   lBF.Hc nsnps       querysnp querytrait typeBF
#> 1 15.32189 11.96576  1000 chr19-11182353    Trait_1    ABF
```
