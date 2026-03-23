# per.snp.priors

Estimate per snp priors

## Usage

``` r
per.snp.priors(
  nsnps,
  pa = 3.82e-05,
  pc = 0.00182,
  p1 = NULL,
  p2 = NULL,
  p12 = NULL
)
```

## Arguments

- nsnps:

  number of SNPs

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

priors at the query variant

## Author

Ichcha Manipur
