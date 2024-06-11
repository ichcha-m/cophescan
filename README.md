# CoPheScan

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/cophescan)](https://cran.r-project.org/package=cophescan)

<span><a href="https://ichcha-m.github.io/cophescan/" class="external-link"> <img src="man/figures/logo.png" align="right" height="200" style="float:right; height:180px;"/></a>

The cophescan package implements Coloc adapted Phenome-wide Scan (CoPheScan), a Bayesian method to perform Phenome-wide association studies (PheWAS) that identifies causal associations between genetic variants and phenotypes while simultaneously accounting for confounding due to linkage disequilibrium.

See the description vignette for background and references: [Introduction to CoPheScan](https://ichcha-m.github.io/cophescan/articles/IntroductionCoPheScan_01.html)

------------------------------------------------------------------------

#### **Installation**

Install cophescan from CRAN

```r
install.packages("cophescan")
```

Install the developmental version from GitHub

``` r
if(!require("remotes"))
   install.packages("remotes") # if necessary
remotes::install_github("ichcha-m/cophescan")
```

------------------------------------------------------------------------

#### **Vignettes**

For a detailed walkthrough of cophescan browse through the vignettes in <https://ichcha-m.github.io/cophescan/>.

Vignette articles:

1\.
[Introduction to CoPheScan](https://ichcha-m.github.io/cophescan/articles/IntroductionCoPheScan_01.html)

2\.
[Input data format](https://ichcha-m.github.io/cophescan/articles/InputData_02.html)

3\.
[An example for running CoPheScan with fixed priors](https://ichcha-m.github.io/cophescan/articles/FixedPriors_03.html)

4\.
[An example for running CoPheScan with hierarchical priors](https://ichcha-m.github.io/cophescan/articles/HierarchicalPriors_04.html)

------------------------------------------------------------------------

#### **Quick start**

#### Usage

``` r
library(cophescan)
## Load the simulated summary stats data of 24 traits
data("cophe_multi_trait_data")
names(cophe_multi_trait_data)
```

##### Single trait

``` r
query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
querysnpid <- cophe_multi_trait_data$querysnpid
print(querysnpid)
query_trait_1$position <- sapply(query_trait_1$snp, function(x) as.numeric(unlist(strsplit(x, "-"))[2]))
plot_trait_manhat(query_trait_1, querysnpid)

# Run cophescan under a single causal variant assumption by providing the snpid of the query variant (querysnpid) for the query trait.
res.single <- cophe.single(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
summary(res.single)
# Run cophescan with susie (multiple variants) by providing the snpid of the query variant (querysnpid) for the query trait
query_trait_1$LD <- cophe_multi_trait_data$LD
res.susie <- cophe.susie(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
summary(res.susie)
```

##### Run multi-trait analysis

``` r
res.multi <- cophe.multitrait(cophe_multi_trait_data$summ_stat, querysnpid = querysnpid, querytrait.names = names(cophe_multi_trait_data$summ_stat), method = 'single')
```

##### Plot cophescan results

``` r
cophe.plots.res <- cophe_plot(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, querysnpid = querysnpid)

ggpubr::ggarrange(cophe.plots.res$pval, cophe.plots.res$ppHa, cophe.plots.res$ppHc, nrow=1)

# cophe.plots.hmp <- cophe_heatmap(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, querysnpid = querysnpid, color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name ="Greens")))(100))
                                    
```

##### Run hierarchical model for priors

``` r
cophe.hier.res <- run_metrop_priors(res.multi, posterior=TRUE, avg_posterior=TRUE, pik=TRUE) 
ll <- cophe.hier.res$ll
params <- cophe.hier.res$parameters

### store user parameters
old_par = par(no.readonly = TRUE)

## Plot mcmc diagnostics
par(mfrow=c(2,2))
plot(1:length(ll), ll, main="loglik",type="l", col="orange")
plot(1:ncol(params), params[1,], main="alpha",type="l", col="orange")
plot(1:ncol(params), params[2,], main="beta",type="l", col="orange")

### reset user parameters
par(old_par)
```

##### Predict

``` r
res.post.prob = cbind(cophe.hier.res$avg.posterior, cophe.hier.res$data)
res.hier.predict <- cophe.hyp.predict(as.data.frame(res.post.prob ))
tail(res.hier.predict, row.names = F)
```

------------------------------------------------------------------------

#### [NEWS: cophescan 1.4.1](https://ichcha-m.github.io/cophescan/news/index.html)

------------------------------------------------------------------------
