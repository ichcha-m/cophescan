##' per.snp.priors
##'
##' Estimate per snp priors
##' @param nsnps number of SNPs
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @param pn prior probability that none of the SNPS are associated with the queried trait , default \eqn{pn = 1 - (pa*(nsnps-1)) - pc}
##' @param pa prior probability a SNP other that the causal variant (for a different trait) is associated with the queried trait , default \eqn{pa = p2}
##' @param pc prior probability that the known causal variant (for a different trait) is associated with the queried trait, default \eqn{pc =  p12/p1+p12}
##' @return priors at the causal variant
##' @export
##' @author Ichcha Manipur

per.snp.priors <- function(nsnps, p1=1e-4, p2=1e-4, p12=1e-5,
                           pa=NULL, pc=NULL){
  if (is.null(pc)){
    pc <- p12/(p1+p12)
  }
  if (is.null(pa)){
    pa <- p2
  }
  pn <- 1 - (pa*(nsnps-1)) - pc
  priors=c(pn=pn, pa=pa, pc=pc)
  return(priors)
}

##' hypothesis.priors
##'
##' @param nsnps number of SNPs
##' @param pn prior probability none of the SNPS are associated with the queried trait
##' @param pa prior probability a SNP other that the causal variant (for a different trait) is associated with the queried trait
##' @param pc prior probability that the known causal variant (for a different trait) is associated with the queried trait
##' @return hypotheses priors
##' @export
##' @author Ichcha Manipur
hypothesis.priors <- function(nsnps, pn, pa, pc){
  hp <- c(Hn=pn, Ha=(pa*(nsnps-1)), Hc=pc)
  return(hp)
}

##' Internal function, calculate posterior probabilities for configurations, given logABFs for each SNP and prior probs and position of known causal variant of trait 1
##'
##' @title combine.bf.kc
##' @param labf log approximate bayes factors
##' @param pn prior probability a SNP other than the causal variant is associated with trait 1
##' @param pa prior probability a SNP other than the causal variant is associated with trait 2
##' @param pc prior probability the causal SNP of trait 1 is associated with both traits
##' @param causalpos1 Position of trait1 causal SNP
##' @return named numeric vector of posterior probabilities and bayes factors
##' @author Ichcha Manipur
combine.bf.kc <- function(labf, pn, pa, pc, causalpos1) {
  lHn.bf <- 0
  lBF.Ha <- coloc:::logsum(labf[-causalpos1])
  lHa.bf <- (log(pa) - log(pn)) + lBF.Ha
  lBF.Hc <- labf[causalpos1]
  lHc.bf <- (log(pc) - log(pn)) + lBF.Hc
  # overall bf
  bf <- c(lBF.Ha, lBF.Hc)
  names(bf) <- c('lBF.Ha', 'lBF.Hc')

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
##' @param pa prior probability a SNP other that the causal variant (for a different trait) is associated with the queried trait , default \eqn{pn = 1- pc}
##' @param pc prior probability that the known causal variant (for a different trait) is associated with the queried trait, default \eqn{pc =  p12/p1+p12}
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of Hn (no shared causal variant), Ha (two distinct causal variants) and Hc (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.Hc of the SNP being causal for the shared signal *if* Hc is true. This is only relevant if the posterior support for Hc in summary is convincing.
##' }
##' @importFrom coloc check_dataset
##' @author Ichcha Manipur
##' @export
cophe.single <- function(dataset, causal.snpid, MAF=NULL, p1=1e-4, p2=1e-4, p12=1e-5,
                           pa=NULL, pc=NULL) {


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

    # number of snps in the region
  common.snps <- nrow(df)

  psp  <-  per.snp.priors(nsnps = common.snps, p1 = p1, p2 = p2, p12 = p12, pa = pa, pc = pc)
  print('SNP Priors')
  print(psp)

  hp <- hypothesis.priors(nsnps = common.snps, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])
  print('Hypothesis Priors')
  print(hp)

  pp.bf <- combine.bf.kc(df$lABF.df, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]], causalpos1 = causalpos1)
  results <- do.call("data.frame",c(list(nsnps=common.snps), as.list(pp.bf$pp), as.list(pp.bf$bf), querysnp=causal.snpid))
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
