
## CoPheScan

<!-- badges: start -->
<!-- badges: end -->

### Installation

You can install cophescan in R with:

``` r
if(!require("remotes"))
   install.packages("remotes") # if necessary
remotes::install_github("ichcha-m/cophescan")

```

### Usage

``` r
library(cophescan)
## Load the simulated summary stats data of 30 traits
data("cophe_multi_trait_data")
attach(cophe_multi_trait_data)
names(cophe_multi_trait_data)
```

#### Single trait
```r
trait1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
causal.snpid <- cophe_multi_trait_data$causal.snpid
print(causal.snpid)
trait1$position <- sapply(trait1$snp, function(x) as.numeric(unlist(strsplit(x, "-"))[2]))
plot_trait_manhat(trait1, causal.snpid)

# Run cophescan under a single causal variant assumption by providing the snpid of the known causal variant for trait 1 = causal.snpid
res.single <- cophe.single(trait1, causal.snpid = causal.snpid)

# Run cophescan with susie (multiple variants) by providing the snpid of the known causal variant for trait 1 = causal.snpid
trait1$LD <- LD
res.susie <- cophe.susie(trait1, causal.snpid = causal.snpid)

```

#### Run multi-trait analysis
```r
res.multi <- cophe.multitrait(cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid, method = 'single')

```

#### Plot cophescan results
```r
cophe.plots.res <- cophe_plot(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid)

ggpubr::ggarrange(cophe.plots.res$pval, cophe.plots.res$ppHa, cophe.plots.res$ppHc, nrow=1)

# cophe.plots.hmp <- cophe_heatmap(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid, color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name ="Greens")))(100))
                                    

```

#### Run hierarchical model for priors
```r
cophe.hier.res <- run_metrop_priors(res.multi, posterior = T, avg_posterior=T, pik=T) 
ll <- cophe.hier.res$ll
params <- cophe.hier.res$params
par(mfrow=c(2,2))
plot(1:length(ll), ll, main="loglik",type="l", col="orange")
plot(1:ncol(params), params[1,], main="alpha",type="l", col="orange")
plot(1:ncol(params), params[2,], main="beta",type="l", col="orange")

```

