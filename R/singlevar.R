##' causal.priors
##'
##' Estimate priors at the known causal variant of trait 1
##' @param p1 prior probability a SNP is associated with trait 1
##' @param p12 prior probability a SNP is associated with both traits
##' @return priors at the causal variant
##' @export
##' @author Ichcha Manipur

causal.priors <- function(p1, p12){
  p12c <-  p12/(p1+p12)
  p1c <- 1- p12c
  priors=c(p1c=p1c, p12c=p12c)
  return(priors)
}

##' hypothesis.priors
##'
##' @param p2a prior probability a SNP other than the causal variant is associated with trait 2, default 1e-4
##' @param p12c prior probability the causal SNP of trait 1 is associated with both traits, default 1e-5
##' @param nsnps number of SNPs
##' @return hypotheses priors
##' @export
##' @author Ichcha Manipur
hypothesis.priors <- function(p2a, p12c, nsnps){
  Hn <- 1 - (p2a*(nsnps-1)) - p12c
  hp <- c(Hn=Hn, Ha=(p2a*(nsnps-1)), Hc=p12c)
  return(hp)
}

##' Internal function, calculate posterior probabilities for configurations, given logABFs for each SNP and prior probs and position of known causal variant of trait 1
##'
##' @title combine.bf.kc
##' @param labf log approximate bayes factors
##' @param p1c prior probability a SNP other than the causal variant is associated with trait 1
##' @param p2a prior probability a SNP other than the causal variant is associated with trait 2
##' @param p12c prior probability the causal SNP of trait 1 is associated with both traits
##' @param causalpos1 Position of trait1 causal SNP
##' @return named numeric vector of posterior probabilities and bayes factors
##' @author Ichcha Manipur
combine.bf.kc <- function(labf, p1c, p2a, p12c, causalpos1) {
  lHn.bf <- 0
  lbfak <- coloc:::logsum(labf[-causalpos1])
  lHa.bf <- (log(p2a) - log(p1c)) + lbfak
  lbfck <- labf[causalpos1]
  lHc.bf <- (log(p12c) - log(p1c)) + lbfck
  # overall bf
  bf <- c(lbfak, lbfck)
  names(bf) <- c('lbfak', 'lbfck')

  all.bf <- c(lHn.bf, lHa.bf, lHc.bf)
  my.denom.log.bf <- coloc:::logsum(all.bf)
  pp <- exp(all.bf - my.denom.log.bf)
  # pp.bf
  names(pp) <- paste("PP.H", c('n', 'a', 'c') , sep = "")
  # barplot(pp.bf)
  print(signif(pp,3))
  print(paste("PP for shared variant: ", signif(pp["PP.Hc"],3)*100 , '%', sep=''))
  return(list(pp=pp, bf=bf))
}

##' Bayesian cophescan analysis under single causal variant assumption
##'
##' This function calculates posterior probabilities of different
##' causal variant configurations under the assumption of a single
##' causal variant for each trait.
##'
##' If regression coefficients and variances are available, it
##' calculates Bayes factors for association at each SNP.  If only p
##' values are available, it uses an approximation that depends on the
##' SNP's MAF and ignores any uncertainty in imputation.  Regression
##' coefficients should be used if available. Find more input data structure details
##' in the coloc package
##'
##' @title  Bayesian cophescan analysis using Approximate Bayes Factors
##' @param dataset a list with specifically named elements defining the dataset
##'   to be analysed.
##' @param causal.snpid Id of the query variant
##' @param MAF Minor allele frequency vector
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @param p2a prior probability a SNP other that the causal variant (for a different trait) is associated with the queried trait , default \eqn{p1c = 1- p12c}
##' @param p12c prior probability that the known causal variant (for a different trait) is associated with the queried trait, default \eqn{p12c =  p12/p1+p12}
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of Hn (no shared causal variant), Ha (two distinct causal variants) and Hc (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.Hc of the SNP being causal for the shared signal *if* Hc is true. This is only relevant if the posterior support for Hc in summary is convincing.
##' }
##' @importFrom coloc check_dataset
##' @author Ichcha Manipur
##' @export
cophe.single <- function(dataset, causal.snpid, MAF=NULL, p1=1e-4, p2=1e-4, p12=1e-5,
                           p2a=NULL, p12c=NULL) {


  if(!("MAF" %in% names(dataset)) & !is.null(MAF))
    dataset$MAF <- MAF
  coloc::check_dataset(d=dataset,2)
  causalpos1 <- which(dataset$snp%in%causal.snpid)
  if(!causal.snpid %in% dataset$snp) {
    message("please check your dataset, given  causal snp not present in dataset")
    return(NULL)
  }
  print('Running cophe.single...')
  df <- coloc:::process.dataset(d=dataset, suffix="df")
  ##############################

  p2c <- p0c <- 0
  p1a <- p12a <- 0
  if (is.null(p12c)){
    cp  <-  causal.priors(p1, p12)
    p1c  <-  cp[["p1c"]]
    p12c  <-  cp[["p12c"]]
  }else{
    p1c <- 1 - p12c
  }
  if (is.null(p2a)){
    p2a <- p2
  }

  psp <- c(p1c=p1c, p2a=p2a, p12c=p12c)
  print('SNP Priors')
  print(psp)
  common.snps <- nrow(df)
  hp <- hypothesis.priors(p2a, p12c, common.snps)
  print('Hypothesis Priors')
  print(hp)

  pp.bf <- combine.bf.kc(df$lABF.df, p1c=p1c, p2a=p2a, p12c=p12c, causalpos1 = causalpos1)
  results <- c(nsnps=common.snps, pp.bf$pp, pp.bf$bf)
  output <- list(summary=results,
                 results=df,
                 priors=psp, querysnp=causal.snpid)
  attr(output, "class") <- "cophe"
  # class(output) <- c("cophe",class(output))
  return(output)
}

#' print the summary of results from cophescan single or susie
#'
#' @param cophe.res Result from either cophe.susie or cophe.single
#'
#' @return log bayes and posterior probabilities
#' @export
#'
summary.cophe <- function(cophe.res, ...){
  print(cophe.res$summary)
}
