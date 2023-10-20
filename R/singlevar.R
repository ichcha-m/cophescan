#' per.snp.priors
#'
#' Estimate per snp priors
#' @param nsnps number of SNPs
#' @param p1 prior probability a SNP is associated with trait 1, default 1e-4 (coloc prior)
#' @param p2 prior probability a SNP is associated with trait 2, default 1e-4 (coloc prior)
#' @param p12 prior probability a SNP is associated with both traits, default 1e-5 (coloc prior)
#' @param pa prior probability that a non-query variant is causally associated with the query trait
#' @param pc prior probability that the query variant is causally associated with the query trait
#' @return priors at the query variant
#' @export
#' @author Ichcha Manipur

per.snp.priors <- function(nsnps, p1=1e-4, p2=1e-4, p12=1e-5,
                           pa=NULL, pc=NULL){
  if (is.null(pc)){
    pc <- p12/(p1+p12)
  }
  if (is.null(pa)){
    pa <- p2
  }
  if (((pa*(nsnps-1)) + pc) >= 1){
    stop("(pa*(nsnps-1)) - pc should be < 1, Revaluate priors or see help for adjust_priors")
  }
  pn <- 1 - (pa*(nsnps-1)) - pc
  priors <- c(pn=pn, pa=pa, pc=pc)
  return(priors)
}

#' hypothesis.priors
#'
#' Estimate priors for each hypothesis
#' @param nsnps number of SNPs
#' @param pn prior probability that none of the SNPs/variants in the region are associated with the query trait
#' @param pa prior probability that a non-query variant is causally associated with the query trait
#' @param pc prior probability that the query variant is causally associated with the query trait
#' @return hypotheses priors
#' @export
#' @author Ichcha Manipur
hypothesis.priors <- function(nsnps, pn, pa, pc){
  hp <- c(Hn=pn, Ha=(pa*(nsnps
                         -1)), Hc=pc)
  return(hp)
}

#' Calculate posterior probabilities for all the configurations
#'
#' @title combine.bf
#' @param lBF_df dataframe with log bayes factors of hypothesis Ha  and Hn: column names should be lBF.Ha and lBF.Hc
#' @param pn prior probability that none of the SNPs/variants in the region are associated with the query trait
#' @param pa prior probability that a non-query variant is causally associated with the query trait
#' @param pc prior probability that the query variant is causally associated with the query trait
#' @return named numeric vector of posterior probabilities and bayes factors
#' @export
#' @author Ichcha Manipur
combine.bf <- function(lBF_df, pn, pa, pc) {
  lHn.bf <- 0
  lBF.Ha <- lBF_df$lBF.Ha
  lHa.bf <- (log(pa) - log(pn)) + lBF.Ha
  lBF.Hc <- lBF_df$lBF.Hc
  lHc.bf <- (log(pc) - log(pn)) + lBF.Hc
  # overall bf
  bf <- c(lBF.Ha, lBF.Hc)
  names(bf) <- c('lBF.Ha', 'lBF.Hc')

  all.bf <- c(lHn.bf, lHa.bf, lHc.bf)
  denom.log.bf <- logsum(all.bf)
  pp <- exp(all.bf - denom.log.bf)
  # pp.bf
  names(pp) <- paste("PP.H", c('n', 'a', 'c') , sep = "")
  # barplot(pp.bf)
  message(signif(pp,3))
  message(paste("PP for causal query variant: ", signif(pp["PP.Hc"],3)*100 , '%', sep=''))
  return(list(pp=pp, bf=bf))
}

#' Bayesian cophescan analysis under single causal variant assumption
#'
#' This function calculates posterior probabilities of different
#' causal variant configurations under the assumption of a single
#' causal variant for each trait.
#'
#' If regression coefficients and variances are available, it
#' calculates Bayes factors for association at each SNP.  If only p
#' values are available, it uses an approximation that depends on the
#' SNP's MAF and ignores any uncertainty in imputation.  Regression
#' coefficients should be used if available. Find more input data structure details
#' in the coloc package
#'
#' @title  Bayesian cophescan analysis using Approximate Bayes Factors
#' @param dataset a list with specifically named elements defining the query trait dataset
#'   to be analysed.
#' @param querysnpid Id of the query variant, (id in dataset$snp)
#' @param querytrait Query trait name
#' @param MAF Minor allele frequency vector
#' @param p1 prior probability a SNP is associated with trait 1, default 1e-4 (coloc prior)
#' @param p2 prior probability a SNP is associated with trait 2, default 1e-4 (coloc prior)
#' @param p12 prior probability a SNP is associated with both traits, default 1e-5 (coloc prior)
#' @param pa prior probability that a non-query variant is causally associated with the query trait , default \eqn{pa = p2} (cophescan prior)
#' @param pc prior probability that the query variant is causally associated with the query trait, default \eqn{pc =  p12/p1+p12} (cophescan prior)
#' @return a list of two \code{data.frame}s:
#' \itemize{
#' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of Hn (no shared causal variant), Ha (two distinct causal variants) and Hc (one common causal variant)
#' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.Hc of the SNP being causal for the shared signal *if* Hc is true. This is only relevant if the posterior support for Hc in summary is convincing.
#' }
#' @importFrom coloc check_dataset
#' @examples
#' library(cophescan)
#' data(cophe_multi_trait_data)
#' query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
#' querysnpid <- cophe_multi_trait_data$querysnpid
#' res.single <- cophe.single(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
#' summary(res.single)
#' @author Ichcha Manipur
#' @export
cophe.single <- function(dataset, querysnpid, querytrait, MAF=NULL, p1=1e-4, p2=1e-4, p12=1e-5,
                           pa=NULL, pc=NULL) {
  message('Running cophe.single...')

  lABF.df = cophe.single.lbf(dataset, querysnpid, querytrait, MAF)

  # number of snps in the region
  nsnps <- lABF.df$nsnps

  psp  <-  per.snp.priors(nsnps = nsnps, p1 = p1, p2 = p2, p12 = p12, pa = pa, pc = pc)
  message('SNP Priors')
  message(psp)

  hp <- hypothesis.priors(nsnps = nsnps, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])
  message('Hypothesis Priors')
  message(hp)

  pp.bf <- combine.bf(lABF.df, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])

  results <- do.call("data.frame",c(list(nsnps=nsnps), as.list(pp.bf$pp), as.list(pp.bf$bf), querysnp=querysnpid, querytrait=querytrait, typeBF='ABF'))
  output <- list(summary=results,
                 priors=psp, querysnp=querysnpid, querytrait=querytrait)
  attr(output, "class") <- "cophe"
  # class(output) <- c("cophe",class(output))
  return(output)
}


#' cophe.single.lbf
#'
#' Calculate log bayes factors for each hypothesis (Single causal variant assumption)
#' @param dataset a list with specifically named elements defining the query trait dataset
#'   to be analysed.
#' @param querysnpid Id of the query variant, (id in dataset$snp)
#' @param querytrait Query trait name
#' @param MAF Minor allele frequency vector
#' @seealso \code{\link{cophe.single}}
#' @return data frame with log bayes factors for Hn and Ha hypotheses
#' @export
#' @examples
#' library(cophescan)
#' data(cophe_multi_trait_data)
#' query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
#' querysnpid <- cophe_multi_trait_data$querysnpid
#' res.single.lbf <- cophe.single.lbf(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
#' res.single.lbf
#' @author Ichcha Manipur
cophe.single.lbf <- function(dataset, querysnpid, querytrait, MAF=NULL) {
  proc_data = cophe.prepare.dat.single(dataset, querysnpid, MAF)
  df=proc_data$df
  querypos = proc_data$querypos

  lBF.persnp = df$lABF.df
  lBF.Ha <- logsum(lBF.persnp[-querypos])
  lBF.Hc <- lBF.persnp[querypos]

  # overall bf
  out_bf <- data.frame(lBF.Ha = lBF.Ha, lBF.Hc = lBF.Hc, nsnps=nrow(df), querysnp=querysnpid, querytrait=querytrait, typeBF='ABF')
  return(out_bf)
}


#' Prepare data for cophe.single
#'
#' @param dataset a list with specifically named elements defining the query trait dataset
##'   to be analysed.
#' @param querysnpid Id of the query variant, (id in dataset$snp)
#'
#' @return named list containing the per snp BFs (df) and position of the query variant (querypos)
#' @keywords internal
#' @author Ichcha Manipur
cophe.prepare.dat.single <- function(dataset, querysnpid, MAF=NULL){
  if(!("MAF" %in% names(dataset)) & !is.null(MAF))
    dataset$MAF <- MAF
  coloc::check_dataset(d=dataset,2)
  querypos <- which(dataset$snp%in%querysnpid)
  if (!querysnpid %in% dataset$snp) {
    stop("Please check your dataset, queried snp not present in dataset")
  }
  df <- coloc::process.dataset(d=dataset, suffix="df")
  return(list(df=df, querypos=querypos))

}

#' adjust_priors
#'
#' adjust fixed priors when nsnps in region is high
#'
#' @param nsnps number of SNPs
#' @param p1 prior probability a SNP is associated with trait 1, default 1e-4 (coloc prior)
#' @param p2 prior probability a SNP is associated with trait 2, default 1e-4 (coloc prior)
#' @param p12 prior probability a SNP is associated with both traits, default 1e-5 (coloc prior)
#' @param pa prior probability that a non-query variant is causally associated with the query trait , default \eqn{pa = p2} (cophescan prior)
#' @param pc prior probability that the query variant is causally associated with the query trait, default \eqn{pc =  p12/p1+p12} (cophescan prior)
#'
#' @return vector of pn, pa and pc adjusted prior probabilities
#' @export
#'
adjust_priors <- function(nsnps, p1=1e-4, p2=1e-4, p12=1e-5,
                          pa=NULL, pc=NULL) {
  if (is.null(pc)){
    pc <- p12/(p1+p12)
  }
  if (is.null(pa)){
    pa <- p2
  }
  sum_priors <- (nsnps+1)*pa + pc
  pc <- pc/sum_priors
  pa <- pa/sum_priors
  pn <- pa
  return(c(pn=pn, pa=pa, pc=pc))
}



#' print the summary of results from cophescan single or susie
#'
#' @param object Result from either `cophe.susie` or `cophe.single`
#' @param ... additional arguments affecting the summary produced.
#' @return log bayes and posterior probabilities
#' @seealso \code{\link{cophe.single}}, \code{\link{cophe.susie}}
#' @export
#'
summary.cophe <- function(object, ...){
  summ =object$summary
  class(summ) = c("summary.cophe", "data.frame")
  return(summ)
  # message(cophe.res$summary)
}

##' Internal function, logsum
##' Function directly taken from coloc
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}
