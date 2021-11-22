#' cophe.multitrait
#'
#' @param trait.dat  List of coloc structured data for n_traits (Total number of traits)
#' @param causal.snpid query variant id
#' @param method either 'single' for cophe.single or 'susie' for cophe.susie
#' @param LDmat LD matrix
#' @param simplify if True removes intermediate results from output
#' @return if simplify is False returns multi-trait list of lists, each with two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of Hn (no shared causal variant), Ha (two distinct causal variants) and Hc (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.Hc of the SNP being causal for the shared signal *if* Hc is true. This is only relevant if the posterior support for Hc in summary is convincing.
##' }
##' if simplify is False only returns dataframe with posterior probabilties of Hn, Hc and Ha with no intermediate results
#' @export
##' @author Ichcha Manipur
cophe.multitrait <- function(trait.dat, causal.snpid, LDmat=NULL, method='single', simplify=F){
  if (!is.null(LDmat) & method == 'susie'){
      print('Please provide the LD matrix')
      return(NULL)
  }
  cophe_results <- list()
  for (idx in seq_along(trait.dat)){
    dat <- trait.dat[[idx]]
    if (causal.snpid%in%dat$snp){
      if (method=='susie'){
        rownames(LD) <- colnames(LD) <- dat$snp
        dat$LD <- LD
        cophe_results[[idx]] <- cophe.susie(dat, causal.snpid)
      } else{
        cophe_results[[idx]] <- cophe.single(dat, causal.snpid)
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
multitrait.simplify <- function(multi.dat){
  if (is.null(names(multi.dat))){
    options(ggrepel.max.overlaps = Inf)
    print('Please assign the names of the phenotypes to the dataset: names(multi.dat) = vector of phenotype names')
    return(NULL)
  }
  print('Simplifying...')
  print('Removing intermediate results...')
  # v=unlist(sapply(multi.dat, function(x) any(names(x$summary)%in%'hit1') ))
  pp_df <- data.frame()
  for (trait in seq_along(multi.dat)){
    dat <- multi.dat[[trait]]
    if (any(names(dat$summary)%in%'hit1')){
      pp <-  as.data.frame(t(dat$summary[, c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'bf.a', 'bf.c')]))
      pp$querysnp <-dat$querysnp
      rownames(pp) <- paste0(names(multi.dat)[trait], '_hit_', dat$summary$hit2)
      pp_df <- rbind(pp_df, pp)
    } else {
      pp <- as.data.frame(t(dat$summary[c('PP.Hn' , 'PP.Ha', 'PP.Hc', 'nsnps', 'bf.a', 'bf.c')]))
      pp$querysnp <-dat$querysnp
      rownames(pp) <-names(multi.dat)[trait]
      pp_df <- rbind(pp_df, pp)
    }
  }
  colnames(pp_df) <- c('Hn' , 'Ha', 'Hc', 'nsnps', 'querysnp', 'bf.a', 'bf.c')
  print('Done')
  return(pp_df)
}
