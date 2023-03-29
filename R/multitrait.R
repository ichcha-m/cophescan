#' Run cophescan on multiple traits at once
#'
#' @param trait.dat  Named(traits) list of coloc structured data for k traits (Total number of traits)
#' @param query.snpid vector of query variant ids = length(trait.dat), if the same variant
#' @param querytrait.names vector of names for the query traits, if the names of
#' the multi.dat list contain the trait names please pass querytrait.names=names(multi.dat)
#' @param method either 'single' for cophe.single or 'susie' for cophe.susie
#' @param LDmat LD matrix
#' @param simplify if True removes intermediate results from output
#' @param ... additional arguments of priors for cophe.susie or cophe.single
#' @return if simplify is False returns multi-trait list of lists, each with two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of Hn (no shared causal variant), Ha (two distinct causal variants) and Hc (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.Hc of the SNP being causal for the shared signal *if* Hc is true. This is only relevant if the posterior support for Hc in summary is convincing.
##' }
##' if simplify is False only returns dataframe with posterior probabilties of Hn, Hc and Ha with no intermediate results
#' @export
##' @author Ichcha Manipur
cophe.multitrait <- function(trait.dat, query.snpid, querytrait.names, LDmat=NULL, method='single', simplify=F, ...){
  if (length(query.snpid)==1){
    query.snpid  = rep(query.snpid, length(trait.dat))
  }
  if (length(querytrait.names)==1){
    querytrait.names = rep(querytrait.names, length(trait.dat))
  }
  if (is.null(LDmat) & method == 'susie'){
      print('Please provide the LD matrix')
      return(NULL)
  }
  cophe_results <- list()
  for (idx in seq_along(trait.dat)){
    dat <- trait.dat[[idx]]
    qv <- query.snpid[idx]
    qt <- querytrait.names[idx]
    if (qv%in%dat$snp){
      if (method=='susie'){
        dat$LD <- LDmat
        rownames(dat$LD) <- colnames(dat$LD) <- dat$snp
        cophe_results[[idx]] <- cophe.susie(dat, query.snpid=qv, querytrait=qt, ...)
      } else{
        cophe_results[[idx]] <- cophe.single(dat, query.snpid=qv, querytrait=qt, ...)
      }
    }
  }
  names(cophe_results) <- names(trait.dat)
  if (simplify){
    cophe_results <- multitrait.simplify(cophe_results)
  }
  return(cophe_results)
}


#' simplify.multitrait
#' Simplifying the output obtained from cophe.multitrait
#' @param multi.dat output obtained from cophe.multitrait
#' @return  dataframe with posterior probabilties of Hn, Hc and Ha
multitrait.simplify <- function(multi.dat, cov_labels=NULL){
  pp_df <- data.frame()

  for (trait in seq_along(multi.dat)){
    dat <- multi.dat[[trait]]
    pp <-  as.data.frame((dat$summary[, c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp', 'querytrait')]))
    pp$querysnp <-dat$querysnp
    if (any(names(dat$summary)%in%'hit1')){
      rownames(pp) <- paste0(names(multi.dat)[trait], '_hit_', dat$summary$hit2)
    } else {
      rownames(pp) <-names(multi.dat)[trait]
    }
      pp_df <- rbind(pp_df, pp)
  }
  colnames(pp_df) <- c('Hn' , 'Ha', 'Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp', 'querytrait')
  return(pp_df)
}
