# Package index

## CoPheScan with fixed priors

- [`per.snp.priors()`](https://ichcha-m.github.io/cophescan/reference/per.snp.priors.md)
  : per.snp.priors

- [`adjust_priors()`](https://ichcha-m.github.io/cophescan/reference/adjust_priors.md)
  : adjust_priors

- [`hypothesis.priors()`](https://ichcha-m.github.io/cophescan/reference/hypothesis.priors.md)
  : hypothesis.priors

- [`cophe.single()`](https://ichcha-m.github.io/cophescan/reference/cophe.single.md)
  : Bayesian cophescan analysis using Approximate Bayes Factors

- [`cophe.single.lbf()`](https://ichcha-m.github.io/cophescan/reference/cophe.single.lbf.md)
  : cophe.single.lbf

- [`cophe.susie()`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.md)
  :

  run `cophe.susie` using susie to detect separate signals

- [`cophe.susie.lbf()`](https://ichcha-m.github.io/cophescan/reference/cophe.susie.lbf.md)
  : cophe.susie.lbf

- [`combine.bf()`](https://ichcha-m.github.io/cophescan/reference/combine.bf.md)
  : combine.bf

- [`cophe.multitrait()`](https://ichcha-m.github.io/cophescan/reference/cophe.multitrait.md)
  : Run cophescan on multiple traits at once

- [`summary(`*`<cophe>`*`)`](https://ichcha-m.github.io/cophescan/reference/summary.cophe.md)
  : print the summary of results from cophescan single or susie

- [`multitrait.simplify()`](https://ichcha-m.github.io/cophescan/reference/multitrait.simplify.md)
  :

  Simplifying the output obtained from `cophe.multitrait`,
  `cophe.single` or `cophe.susie`

- [`logsum()`](https://ichcha-m.github.io/cophescan/reference/logsum.md)
  : logsum

## CoPheScan with hierarchical priors

- [`run_metrop_priors()`](https://ichcha-m.github.io/cophescan/reference/run_metrop_priors.md)
  : Run the hierarchical Metropolis Hastings model to infer priors
- [`average_piks()`](https://ichcha-m.github.io/cophescan/reference/average_piks.md)
  : Average of priors: pnk, pak and pck
- [`average_piks_list()`](https://ichcha-m.github.io/cophescan/reference/average_piks_list.md)
  : Average of priors: pnk, pak and pck from list (memory intensive)
- [`average_posterior_prob()`](https://ichcha-m.github.io/cophescan/reference/average_posterior_prob.md)
  : Average of posterior probabilities: Hn, Ha and Hc
- [`average_posterior_prob_list()`](https://ichcha-m.github.io/cophescan/reference/average_posterior_prob_list.md)
  : Average of posterior probabilities: Hn, Ha and Hc from list (memory
  intensive)
- [`get_posterior_prob()`](https://ichcha-m.github.io/cophescan/reference/get_posterior_prob.md)
  : Calculation of the posterior prob of Hn, Ha and Hc
- [`get_beta()`](https://ichcha-m.github.io/cophescan/reference/get_beta.md)
  : Extract beta and p-values of queried variant
- [`sample_alpha()`](https://ichcha-m.github.io/cophescan/reference/sample_alpha.md)
  : sample alpha
- [`sample_beta()`](https://ichcha-m.github.io/cophescan/reference/sample_beta.md)
  : sample beta
- [`sample_gamma()`](https://ichcha-m.github.io/cophescan/reference/sample_gamma.md)
  : sample gamma
- [`logd_alpha()`](https://ichcha-m.github.io/cophescan/reference/logd_alpha.md)
  : dnorm for alpha
- [`logd_beta()`](https://ichcha-m.github.io/cophescan/reference/logd_beta.md)
  : dgamma for beta
- [`logd_gamma()`](https://ichcha-m.github.io/cophescan/reference/logd_gamma.md)
  : dgamma for gamma
- [`loglik()`](https://ichcha-m.github.io/cophescan/reference/loglik.md)
  : Log likelihood calculation
- [`logpost()`](https://ichcha-m.github.io/cophescan/reference/logpost.md)
  : Log posterior calculation
- [`logsumexp()`](https://ichcha-m.github.io/cophescan/reference/logsumexp.md)
  : Log sum
- [`metrop_run()`](https://ichcha-m.github.io/cophescan/reference/metrop_run.md)
  : Run the hierarchical mcmc model to infer priors
- [`pars2pik()`](https://ichcha-m.github.io/cophescan/reference/pars2pik.md)
  : Conversion of parameters alpha, beta and gamma to pnk, pak and pck
- [`piks()`](https://ichcha-m.github.io/cophescan/reference/piks.md) :
  List of priors: pn, pa and pc over all iterations
- [`posterior_prob()`](https://ichcha-m.github.io/cophescan/reference/posterior_prob.md)
  : List of posterior probabilities: Hn, Ha and Hc over all iterations
- [`logpriors()`](https://ichcha-m.github.io/cophescan/reference/logpriors.md)
  : Calculate log priors
- [`pars_init()`](https://ichcha-m.github.io/cophescan/reference/pars_init.md)
  : Initiate parameters alpha, beta and gamma
- [`propose()`](https://ichcha-m.github.io/cophescan/reference/propose.md)
  : Proposal distribution
- [`target()`](https://ichcha-m.github.io/cophescan/reference/target.md)
  : Target distribution

## Predict hypothesis

- [`cophe.hyp.predict()`](https://ichcha-m.github.io/cophescan/reference/cophe.hyp.predict.md)
  : Predict cophescan hypothesis for tested associations
- [`Hc.cutoff.fdr()`](https://ichcha-m.github.io/cophescan/reference/Hc.cutoff.fdr.md)
  : Estimate the Hc.cutoff for the required FDR

## Visualization

- [`plot_trait_manhat()`](https://ichcha-m.github.io/cophescan/reference/plot_trait_manhat.md)
  : Plot region Manhattan for a trait highlighting the queried variant
- [`cophe_plot()`](https://ichcha-m.github.io/cophescan/reference/cophe_plot.md)
  : cophe_plots showing the Ha and Hc of all traits and labelled above
  the specified threshold
- [`cophe_heatmap()`](https://ichcha-m.github.io/cophescan/reference/cophe_heatmap.md)
  : Heatmap of multi-trait cophescan results
- [`prepare_plot_data()`](https://ichcha-m.github.io/cophescan/reference/prepare_plot_data.md)
  : Prepare data for plotting

## Test data

- [`cophe_multi_trait_data`](https://ichcha-m.github.io/cophescan/reference/cophe_multi_trait_data.md)
  : Simulated multi-trait data

## Package

- [`cophescan-package`](https://ichcha-m.github.io/cophescan/reference/cophescan-package.md)
  [`cophescan`](https://ichcha-m.github.io/cophescan/reference/cophescan-package.md)
  : The 'cophescan' package.
