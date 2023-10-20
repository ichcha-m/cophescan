#' Run cophescan on multiple traits at once
#'
#' @param trait.dat  Named(traits) list of coloc structured data for k traits (Total number of traits)
#' @param querysnpid vector of query variant ids = length(trait.dat), if the same variant
#' @param querytrait.names vector of names for the query traits, if the names of
#' the multi.dat list contain the trait names please pass querytrait.names=names(multi.dat)
#' @param method either 'single' for `cophe.single` or 'susie' for `cophe.susie`
#' @param LDmat LD matrix
#' @param simplify if TRUE removes intermediate results from output using 'multitrait.simplify'
#' @param predict.hyp if TRUE predicts the hypothesis based on the provided thresholds for pp.Hc and pp.Hn (overrides simplify) using `cophe.hyp.predict`
#' @param Hn.cutoff threshold for PP.Hc above which the associations are called Hc
#' @param Hc.cutoff threshold for PP.Hc above which the associations are called Hn
#' @param est.fdr.based.cutoff if True calculates the Hc.cutoff using 1-mean(PP.Hc)|PP.Hc > cutoff
#' @param fdr fdr threshold to estimate Hc.cutoff
#' @param ... additional arguments of priors for `cophe.susie` or `cophe.single`
#' @return if simplify is False returns multi-trait list of lists, each with:
##' \itemize{
##' \item a summary data.frame of the cophescan results
##' \item priors used
##' \item querysnp
##' \item querytrait
##' }
##' if simplify is TRUE only returns dataframe with posterior probabilties of Hn, Hc and Ha with no intermediate results
##' if predict.hyp is TRUE returns a dataframe with output of simplify and the predicted hypotheses for all associations
#' @export
##' @author Ichcha Manipur
cophe.multitrait <- function(trait.dat, querysnpid, querytrait.names, LDmat=NULL, method='single', simplify=FALSE, predict.hyp=TRUE, Hn.cutoff = 0.2, Hc.cutoff = 0.6, est.fdr.based.cutoff = FALSE, fdr = 0.05, ...){
  if (length(querysnpid)==1){
    querysnpid  = rep(querysnpid, length(trait.dat))
  }
  if (length(querytrait.names)==1){
    querytrait.names = rep(querytrait.names, length(trait.dat))
  }
  if (is.null(LDmat) & method == 'susie'){
      stop('Please provide the LD matrix')
  }
  cophe_results <- list()
  for (idx in seq_along(trait.dat)){
    dat <- trait.dat[[idx]]
    qv <- querysnpid[idx]
    qt <- querytrait.names[idx]
    if (qv%in%dat$snp){
      if (method=='susie'){
        dat$LD <- LDmat
        rownames(dat$LD) <- colnames(dat$LD) <- dat$snp
        cophe_results[[idx]] <- cophe.susie(dat, querysnpid=qv, querytrait=qt, ...)
      } else{
        cophe_results[[idx]] <- cophe.single(dat, querysnpid=qv, querytrait=qt, ...)
      }
    }
  }
  names(cophe_results) <- names(trait.dat)
  if (predict.hyp == TRUE){
    cophe_results <- cophe.hyp.predict(cophe_results)
    if (simplify){
      message("results in list simplified to data.frame when predict.hyp=TRUE")
    }
    simplify = FALSE
  }
  if (simplify){
    cophe_results <- multitrait.simplify(cophe_results)
  }
  return(cophe_results)
}


#' Simplifying the output obtained from `cophe.multitrait`, `cophe.single` or `cophe.susie`
#'
#' @param multi.dat output obtained from `cophe.multitrait`, `cophe.single` or `cophe.susie`
#' @param only_BF return only bayes factors and not posterior probabilities (default=FALSE)
#'
#' @return  dataframe with posterior probabilties of Hn, Hc and Ha
multitrait.simplify <- function(multi.dat, only_BF=FALSE){
  if (inherits(multi.dat, 'cophe')){
    multi.dat = list(multi.dat)
  }

  pp_df <- data.frame()

  for (trait in seq_along(multi.dat)){
    dat <- multi.dat[[trait]]
    if (any(grepl('PP', colnames(dat$summary)))){
      pp <-  as.data.frame((dat$summary[, c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp', 'querytrait', 'typeBF')]))
    } else{
      pp <-  as.data.frame((dat$summary[, c('lBF.Ha', 'lBF.Hc', 'nsnps',  'querysnp', 'querytrait', 'typeBF')]))
    }
    if (any(names(dat$summary)%in%'hit1')){
      pp$sus_labels <- rownames(pp) <- paste0(names(multi.dat)[trait], '_hit_', dat$summary$hit2)
    } else {
      rownames(pp) <-names(multi.dat)[trait]
      pp$sus_labels <- NA
    }
      pp_df <- rbind(pp_df, pp)
  }
  # colnames(pp_df) <- c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'lBF.Ha', 'lBF.Hc', 'querysnp', 'querytrait', 'typeBF')
  return(pp_df)
}
