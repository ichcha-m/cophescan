
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
## Load the cophescan simulated data for 30 traits
data("cophe_multi_trait_data")
attach(cophe_multi_trait_data)
```

#### Single trait
```r
trait1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
causal.snpid <- cophe_multi_trait_data$causal.snpid

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
cophe_plot(res.multi, thresh=0.5)
```

