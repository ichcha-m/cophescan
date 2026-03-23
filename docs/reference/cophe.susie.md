# run `cophe.susie` using susie to detect separate signals

Check if a variant causally associated in one trait might be causal in
another trait

## Usage

``` r
cophe.susie(
  dataset,
  querysnpid,
  querytrait,
  pa = 3.82e-05,
  pc = 0.00182,
  p1 = NULL,
  p2 = NULL,
  p12 = NULL,
  susie.args = list()
)
```

## Arguments

- dataset:

  *either* a list with specifically named elements defining the dataset
  to be analysed. (see
  [check_dataset](https://rdrr.io/pkg/coloc/man/check_dataset.html))

- querysnpid:

  Id of the query variant

- querytrait:

  Query trait name

- pa:

  prior probability that a non-query variant is causally associated with
  the query trait (cophescan prior), default 3.82e-5

- pc:

  prior probability that the query variant is causally associated with
  the query trait (cophescan prior), default 1.82e-3

- p1:

  prior probability a SNP is associated with trait 1, (coloc prior), pc
  derived by using \\pc = p12/p1+p12\\; use p1, p2, p12 only when pa and
  pc are unavailable (See vignettes)

- p2:

  prior probability a SNP is associated with trait 2, (coloc prior), pa
  derived by using \\pa = p2\\

- p12:

  prior probability a SNP is associated with both traits, (coloc prior),
  pc derived by using \\pc = p12/p1+p12\\

- susie.args:

  a named list of additional arguments to be passed to
  [runsusie](https://rdrr.io/pkg/coloc/man/runsusie.html)

## Value

a list, containing elements

- summary a data.table of posterior probabilities of each global
  hypothesis, one row per pairwise comparison of signals from the two
  traits

- results a data.table of detailed results giving the posterior
  probability for each snp to be jointly causal for both traits
  *assuming Hc is true*. Please ignore this column if the corresponding
  posterior support for H4 is not high.

- priors a vector of the priors used for the analysis

## Author

Ichcha Manipur

## Examples

``` r
library(cophescan)
data(cophe_multi_trait_data)
query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
querysnpid <- cophe_multi_trait_data$querysnpid
query_trait_1$LD <- cophe_multi_trait_data$LD
res.susie <- cophe.susie(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
#> running max iterations: 100
#>  converged: TRUE
#> SNP Priors
#> pn pa pc
#> 9.60056e-01 3.82000e-05 1.82000e-03
#> Hypothesis Priors
#> Hn Ha Hc
#> 9.60056e-01 3.81236e-02 1.82000e-03
#> Running cophe.susie...
#> Using 1000 and 999 available
#> PP.Hn PP.Ha PP.Hc
#> 2.12e-03 3.75e-01 6.23e-01
#> PP for causal query variant: 62.3%
summary(res.susie)
#>   nsnps           hit1           hit2       PP.Hn    PP.Ha     PP.Hc   lBF.Ha
#> 1  1000 chr19-11182353 chr19-11182144 0.002115712 0.375218 0.6226663 15.31003
#>     lBF.Hc       querysnp  typeBF idx1 idx2 querytrait
#> 1 11.95277 chr19-11182353 susieBF    1    1    Trait_1
```
