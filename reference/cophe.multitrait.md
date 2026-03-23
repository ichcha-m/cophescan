# Run cophescan on multiple traits at once

Run cophescan on multiple traits at once

## Usage

``` r
cophe.multitrait(
  trait.dat,
  querysnpid,
  querytrait.names,
  LDmat = NULL,
  method = "single",
  simplify = FALSE,
  predict.hyp = TRUE,
  Hn.cutoff = 0.2,
  Hc.cutoff = 0.6,
  est.fdr.based.cutoff = FALSE,
  fdr = 0.05,
  ...
)
```

## Arguments

- trait.dat:

  Named(traits) list of coloc structured data for k traits (Total number
  of traits)

- querysnpid:

  vector of query variant ids = length(trait.dat), if the same variant

- querytrait.names:

  vector of names for the query traits, if the names of the multi.dat
  list contain the trait names please pass
  querytrait.names=names(multi.dat)

- LDmat:

  LD matrix

- method:

  either 'single' for `cophe.single` or 'susie' for `cophe.susie`

- simplify:

  if TRUE removes intermediate results from output using
  'multitrait.simplify'

- predict.hyp:

  if TRUE predicts the hypothesis based on the provided thresholds for
  pp.Hc and pp.Hn (overrides simplify) using `cophe.hyp.predict`

- Hn.cutoff:

  threshold for PP.Hn above which the associations are called Hn

- Hc.cutoff:

  threshold for PP.Hc above which the associations are called Hc

- est.fdr.based.cutoff:

  if True calculates the Hc.cutoff using 1-mean(PP.Hc)\|PP.Hc \> cutoff

- fdr:

  fdr threshold to estimate Hc.cutoff

- ...:

  additional arguments of priors for `cophe.susie` or `cophe.single`

## Value

if simplify is False returns multi-trait list of lists, each with:

- a summary data.frame of the cophescan results

- priors used

- querysnp

- querytrait

if simplify is TRUE only returns dataframe with posterior probabilties
of Hn, Hc and Ha with no intermediate results if predict.hyp is TRUE
returns a dataframe with output of simplify and the predicted hypotheses
for all associations

## Author

Ichcha Manipur
