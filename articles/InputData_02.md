# CoPheScan: Input data

The input dataset for a trait (querytrait) should contain the summary
data for SNPs in a genomic region around the query variant (querysnpid)
and should have the following fields:

**For a Case-control dataset**

beta: $\beta$ or effect size

varbeta: variance of $\beta$ or square of the standard error of $\beta$

snp: SNP identifier which maybe rsid or CHR_BP_REF_ALT or CHR_BP

type:‘cc’

N: sample size

**For a Quantitative dataset**

When, beta and varbeta are not available the following

beta: $\beta$ or effect size

varbeta: variance of $\beta$ or square of the standard error of $\beta$

snp: SNP identifier which maybe rsid or CHR_BP_REF_ALT or CHR_BP

type:‘quant’

N: sample size

sdY: for a quantitative trait, the population standard deviation of the
trait.

**Additional fields in case of missing beta/varbeta or sdY**

MAF: Minor allele frequency (only required when either beta/varbeta or
sdY are unavailable)

pvalues: only required when beta/varbeta are unavailable

s: fraction of samples that are cases (only for a case-control trait
when beta/varbeta are unavailable)

``` r
library(cophescan)
```

Explore the data structure of the example dataset available in the
cophescan package

``` r
data("cophe_multi_trait_data")
trait_dat = cophe_multi_trait_data$summ_stat$Trait_1
str(trait_dat)
#> List of 8
#>  $ beta   : Named num [1:1000] -0.01369 0.01666 0.09057 -0.00571 -0.05606 ...
#>   ..- attr(*, "names")= chr [1:1000] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>  $ varbeta: Named num [1:1000] 0.000516 0.000399 0.003124 0.000419 0.000473 ...
#>   ..- attr(*, "names")= chr [1:1000] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>  $ z      : Named num [1:1000] -0.603 0.834 1.62 -0.279 -2.578 ...
#>   ..- attr(*, "names")= chr [1:1000] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>  $ snp    : chr [1:1000] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>  $ MAF    : Named num [1:1000] 0.2614 0.4871 0.0318 0.4046 0.3042 ...
#>   ..- attr(*, "names")= chr [1:1000] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>  $ type   : chr "cc"
#>  $ N      : num 20000
#>  $ s      : num 0.5
```

**Additional field for `cophe.susie`**

LD: Linkage Disequilibrium matrix with row and column names being the
same as the snp field.

``` r
trait_dat$LD = cophe_multi_trait_data$LD
str(trait_dat$LD[1:10, 1:10])
#>  num [1:10, 1:10] 1 0.0267 -0.1078 -0.0627 0.1033 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:10] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
#>   ..$ : chr [1:10] "chr19-11173352" "chr19-11173626" "chr19-11173716" "chr19-11173807" ...
```

It is important to check that there is alignment of alleles for which
the beta is reported and those in the LD matrix. This can be verified
either using coloc::check_alignment or performing a diagnostic check
using the susie package
<https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html>.

**Note**

- The input summary data structure for a trait is the same format as
  required by coloc. To gain a detailed understanding of the input
  please see the coloc vignette
  <https://chr1swallace.github.io/coloc/articles/a02_data.html> and the
  test datasets provided with the coloc package.
- More information on the input can be obtained from the coloc
  documentation :
  [`?coloc::check_dataset`](https://rdrr.io/pkg/coloc/man/check_dataset.html)

------------------------------------------------------------------------
