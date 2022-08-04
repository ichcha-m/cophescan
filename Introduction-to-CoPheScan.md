---
title: "Introduction to CoPheScan"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to CoPheScan}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(cophescan)
library(ggplot2)
#library(ggpubr)
```

#### Load data  


```r
data("cophe_multi_trait_data")
attach(cophe_multi_trait_data)
names(cophe_multi_trait_data)
#> [1] "summ_stat"    "causal.snpid" "LD"           "rg_vec"
```

The dataset contains the summary statistics of 30 traits with 10 traits each simulated to have high posterior probability of Hc (shared causal variant), Ha and Hn.

The query variant is chr19-11182353

#### Single trait

```r
trait1 <- cophe_multi_trait_data$summ_stat[['Trait_1']] 
causal.snpid <- cophe_multi_trait_data$causal.snpid
trait1$position <- sapply(trait1$snp, function(x) as.numeric(unlist(strsplit(x, "-"))[2]))
plot_trait_manhat(trait1, causal.snpid)
```

<img src="figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r

# Run cophescan under a single causal variant assumption by providing the snpid of the known causal variant for trait 1 = causal.snpid
res.single <- cophe.single(trait1, causal.snpid = causal.snpid)
#> [1] "Running cophe.single..."
#> [1] "SNP Priors"
#>         pn         pa         pc 
#> 0.80919091 0.00010000 0.09090909 
#> [1] "Hypothesis Priors"
#>         Hn         Ha         Hc 
#> 0.80919091 0.09990000 0.09090909 
#>    PP.Hn    PP.Ha    PP.Hc 
#> 5.49e-05 3.06e-02 9.69e-01 
#> [1] "PP for shared variant: 96.9%"

# Run cophescan with susie (multiple variants) by providing the snpid of the known causal variant for trait 1 = causal.snpid
trait1$LD <- LD
res.susie <- cophe.susie(trait1, causal.snpid = causal.snpid)
#> running max iterations: 100
#> 	converged: TRUE
#> [1] "SNP Priors"
#>         pn         pa         pc 
#> 0.80929091 0.00010000 0.09090909 
#> [1] "Hypothesis Priors"
#>         Hn         Ha         Hc 
#> 0.80929091 0.09980000 0.09090909 
#> [1] "Running cophe.susie..."
#> Using 1000 and 999 available
#>    PP.Hn    PP.Ha    PP.Hc 
#> 5.46e-05 3.06e-02 9.69e-01 
#> [1] "PP for shared variant: 96.9%"
```


#### Multi-trait analysis

```r
res.multi <- cophe.multitrait(cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid, method = 'single', simplify = T)
# res.multi <- cophe.multitrait(cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid, method = 'susie', LDmat = cophe_multi_trait_data$LD)
```

##### **Plot cophescan results**

```r
cophe.plots.res <- cophe_plot(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, causal.snpid = causal.snpid)
# ggpubr::ggarrange(cophe.plots.res$pval, cophe.plots.res$ppHa, cophe.plots.res$ppHc, ncol = 2, nrow = 2) + ggplot2::theme(text = ggplot2::element_text(size=2))
cophe.plots.res$pval
cophe.plots.res$ppHc
cophe.plots.res$ppHa
```

<img src="figure/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="48.5%" /><img src="figure/unnamed-chunk-5-2.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="48.5%" /><img src="figure/unnamed-chunk-5-3.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="48.5%" />

##### **cophescan Ternary plots** 

```r
plot_cophe_ternary(res.multi)
plot_cophe_ternary(res.multi, traits.dat = cophe_multi_trait_data$summ_stat, plot_pval = T)
```

<img src="figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="48.5%" /><img src="figure/unnamed-chunk-6-2.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="48.5%" />

#### Run hierarchical model for priors

```r
rg=F
# rg=T
cophe.hier.res <- run_metrop_priors(res.multi, posterior = T, avg_posterior=T, pik=T, avg_pik = T, rg_vec = cophe_multi_trait_data$rg_vec, rg = rg, nits = 30000) 
loglik <- cophe.hier.res$ll
parameters <- cophe.hier.res$parameters
col <- rgb(red = 0.4, green = 0.7, blue = 0.5, alpha = 0.8)
par(mfrow=c(2,2))
plot(1:length(loglik), loglik, main="loglik",type="l", col=col, ylab = "ll", xlab="")
plot(1:ncol(parameters), parameters[1,], main="alpha",type="l", col=col, ylab = "alpha", xlab="")
plot(1:ncol(parameters), parameters[2,], main="beta",type="l", col=col, ylab = "beta", xlab="")
if (rg == T)
  plot(1:ncol(parameters), parameters[3,], main="gamma",type="l", col=col, ylab = "gamma", xlab="")
 summary(cophe.hier.res$avg.pik)
#>       p0k              pak                pck        
#>  Min.   :0.4918   Min.   :0.000141   Min.   :0.3672  
#>  1st Qu.:0.4918   1st Qu.:0.000141   1st Qu.:0.3672  
#>  Median :0.4918   Median :0.000141   Median :0.3672  
#>  Mean   :0.4918   Mean   :0.000141   Mean   :0.3672  
#>  3rd Qu.:0.4918   3rd Qu.:0.000141   3rd Qu.:0.3672  
#>  Max.   :0.4918   Max.   :0.000141   Max.   :0.3672
```

<img src="figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />
