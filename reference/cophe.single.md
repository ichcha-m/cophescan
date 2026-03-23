# Bayesian cophescan analysis using Approximate Bayes Factors

Bayesian cophescan analysis under single causal variant assumption

## Usage

``` r
cophe.single(
  dataset,
  querysnpid,
  querytrait,
  MAF = NULL,
  pa = 3.82e-05,
  pc = 0.00182,
  p1 = NULL,
  p2 = NULL,
  p12 = NULL
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

- MAF:

  Minor allele frequency vector

- pa:

  prior probability that a non-query variant is causally associated with
  the query trait (cophescan prior), default 3.82e-5

- pc:

  prior probability that the query variant is causally associated with
  the query trait (cophescan prior), default 1.82e-3 (cophescan prior)

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

## Value

a list of two `data.frame`s:

- summary is a vector giving the number of SNPs analysed, and the
  posterior probabilities of Hn (no shared causal variant), Ha (two
  distinct causal variants) and Hc (one common causal variant)

- results is an annotated version of the input data containing log
  Approximate Bayes Factors and intermediate calculations, and the
  posterior probability SNP.PP.Hc of the SNP being causal for the shared
  signal *if* Hc is true. This is only relevant if the posterior support
  for Hc in summary is convincing.

## Details

This function calculates posterior probabilities of different causal
variant configurations under the assumption of a single causal variant
for each trait.

If regression coefficients and variances are available, it calculates
Bayes factors for association at each SNP. If only p values are
available, it uses an approximation that depends on the SNP's MAF and
ignores any uncertainty in imputation. Regression coefficients should be
used if available. Find more input data structure details in the coloc
package

## Author

Ichcha Manipur

## Examples

``` r
library(cophescan)
data(cophe_multi_trait_data)
query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
querysnpid <- cophe_multi_trait_data$querysnpid
res.single <- cophe.single(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
#> Running cophe.single...
#> SNP Priors
#> pn pa pc
#> 9.60018e-01 3.82000e-05 1.82000e-03
#> Hypothesis Priors
#> Hn Ha Hc
#> 9.60018e-01 3.81618e-02 1.82000e-03
#> PP.Hn PP.Ha PP.Hc
#> 2.09e-03 3.75e-01 6.23e-01
#> PP for causal query variant: 62.3%
summary(res.single)
#>   nsnps      PP.Hn     PP.Ha     PP.Hc   lBF.Ha   lBF.Hc       querysnp
#> 1  1000 0.00208927 0.3749638 0.6229469 15.32189 11.96576 chr19-11182353
#>   querytrait typeBF
#> 1    Trait_1    ABF
```
