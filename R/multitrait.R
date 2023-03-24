#' Run cophescan on multiple traits at once
#'
#' @param trait.dat  Named(traits) list of coloc structured data for n_traits (Total number of traits)
#' @param query.snpid query variant id
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
cophe.multitrait <- function(trait.dat, query.snpid, LDmat=NULL, method='single', simplify=F, ...){
  if (is.null(names(trait.dat))){
    print('Assign names of the traits: names(trait.dat) = vector of phenotype names')
    return(NULL)
  }
  if (is.null(LDmat) & method == 'susie'){
      print('Please provide the LD matrix')
      return(NULL)
  }
  cophe_results <- list()
  for (idx in seq_along(trait.dat)){
    dat <- trait.dat[[idx]]
    if (query.snpid%in%dat$snp){
      if (method=='susie'){
        dat$LD <- LDmat
        rownames(dat$LD) <- colnames(dat$LD) <- dat$snp
        cophe_results[[idx]] <- cophe.susie(dat, query.snpid, ...)
      } else{
        cophe_results[[idx]] <- cophe.single(dat, query.snpid, ...)
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
#' @param query_trait_names vector of names of the query traits, if the names of
#' the multi.dat list contain the trait names please pass query_trait_names=names(multi.dat)
#' default NULL
#' to keep track of the querysnp-trait combinations
#' @return  dataframe with posterior probabilties of Hn, Hc and Ha
multitrait.simplify <- function(multi.dat, query_trait_names=NULL){
  pp_df <- data.frame()

  for (trait in seq_along(multi.dat)){
    dat <- multi.dat[[trait]]
    pp <-  as.data.frame((dat$summary[, c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp')]))
    pp$querysnp <-dat$querysnp
    if (any(names(dat$summary)%in%'hit1')){
      rownames(pp) <- paste0(names(multi.dat)[trait], '_hit_', dat$summary$hit2)
    } else {
      rownames(pp) <-names(multi.dat)[trait]
    }
      pp_df <- rbind(pp_df, pp)
  }
  colnames(pp_df) <- c('Hn' , 'Ha', 'Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp')
  if (!is.null(query_trait_names)){
    pp_df$querytrait <- query_trait_names
  }
  return(pp_df)
}
