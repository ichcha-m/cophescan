# Prepare data for `cophe.susie`

Prepare data for `cophe.susie`

## Usage

``` r
cophe.prepare.dat.susie(dataset, querysnpid, susie.args)
```

## Arguments

- dataset:

  a list with specifically named elements defining the query trait
  dataset to be analysed.

- querysnpid:

  Id of the query variant, (id in dataset\$snp)

- susie.args:

  a named list of additional arguments to be passed to
  [runsusie](https://rdrr.io/pkg/coloc/man/runsusie.html)

## Value

a list with the output of running susie

## See also

[`cophe.susie`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md)
