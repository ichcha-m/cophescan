#' Check if a variant causally associated in one trait might be causal in another trait
#'
#' @title run `cophe.susie` using susie to detect separate signals
#' @param dataset *either* a list with specifically named elements defining the dataset
#'   to be analysed. (see
#'   \link{check_dataset}), or the result of running \link{runsusie} on such a
#'   dataset
#' @param querysnpid Id of the query variant
#' @param querytrait Query trait name
#' @param pa prior probability that a non-query variant is causally associated with the query trait (cophescan prior), default 3.82e-5
#' @param pc prior probability that the query variant is causally associated with the query trait (cophescan prior), default 1.82e-3 (cophescan prior)
#' @param p1 prior probability a SNP is associated with trait 1, (coloc prior), pc derived by using \eqn{pc =  p12/p1+p12}; use p1, p2, p12 only when pa and pc are unavailable (See vignettes)
#' @param p2 prior probability a SNP is associated with trait 2,  (coloc prior), pa derived by using \eqn{pa = p2}
#' @param p12 prior probability a SNP is associated with both traits,  (coloc prior), pc derived by using \eqn{pc =  p12/p1+p12}
#' @param susie.args a named list of additional arguments to be passed to
#'   \link{runsusie}
#' @return a list, containing elements
#' * summary a data.table of posterior
#'   probabilities of each global hypothesis, one row per pairwise comparison
#'   of signals from the two traits
#' * results a data.table of detailed results giving the posterior probability
#'   for each snp to be jointly causal for both traits *assuming Hc is true*.
#'   Please ignore this column if the corresponding posterior support for H4
#'   is not high.
#' * priors a vector of the priors used for the analysis
#' @importFrom coloc runsusie
#' @examples
#' library(cophescan)
#' data(cophe_multi_trait_data)
#' query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
#' querysnpid <- cophe_multi_trait_data$querysnpid
#' query_trait_1$LD <- cophe_multi_trait_data$LD
#' res.susie <- cophe.susie(query_trait_1, querysnpid = querysnpid, querytrait='Trait_1')
#' summary(res.susie)
#' @export
#' @author Ichcha Manipur
cophe.susie=function(dataset, querysnpid, querytrait,
                     pa=3.82e-5, pc=1.82e-3,
                     p1=NULL, p2=NULL, p12=NULL,
                     susie.args=list()) {
  sus_dat = cophe.prepare.dat.susie(dataset, querysnpid, susie.args)
  cred_set = sus_dat$sets

  # number of snps in the region
  nsnps <- ncol(sus_dat$alpha) - 1

  psp  <-  per.snp.priors(nsnps = nsnps, p1 = p1, p2 = p2, p12 = p12, pa = pa, pc = pc)
  message('SNP Priors')
  message(psp)

  hp <- hypothesis.priors(nsnps = nsnps, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])
  message('Hypothesis Priors')
  message(hp)


  if (is.null(cred_set$cs) || length(cred_set$cs)==0 ){
    warning('No credible sets found with SuSIE...')
    warning('Switching to cophe.single ...')
    ret <- cophe.single(dataset, querysnpid, querytrait, pa=psp[["pa"]], pc=psp[["pc"]])
  }
  # else if (!any(sub(".*[.]","",names(unlist(sus_dat$sets$cs))) %in% querysnpid) & length(sus_dat$sets$cs) < 2){
  #   message('Queried SNP not in the susie credible set ...')
  #   message('Switching to cophe.single ...')
  #   ret <- cophe.single(dataset, querysnpid, pa=psp[["pa"]], pc=psp[["pc"]])
  # }
  else {
    if (!any(sub(".*[.]","",names(unlist(sus_dat$sets$cs))) %in% querysnpid)){
      message('Queried SNP not in the susie credible sets ...')
    }
    message('Running cophe.susie...')
    isnps=colnames(sus_dat$alpha)
    if(!length(isnps))
      return(data.table::data.table(nsnps=NA))
    message("Using ",length(isnps), " and ",ncol(sus_dat$alpha)-1," available")

    ret=cophe.bf_bf(sus_dat, cred_set, querysnpid, querytrait, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])
  }
  attr(ret, "class") <- "cophe"
  return(ret)
}


#' a dataset represented by Bayes factors from SuSIE
#' @title extract data through Bayes factors
#' @param sus_dat a list with the output of running susie
#' @param cred_set credible set extracted from susie
#' @param querysnpid Id of the query variant
#' @param querytrait Query trait name
#' @param pn prior probability that none of the SNPs/variants in the region are associated with the query trait
#' @param pa prior probability that a non-query variant is causally associated with the query trait
#' @param pc prior probability that the query variant is causally associated with the query trait
#' @seealso \code{\link{cophe.susie}}
#' @keywords internal
#' @return bayes factors of signals
cophe.bf_bf=function(sus_dat, cred_set, querysnpid, querytrait, pn=NULL, pa=NULL, pc=NULL, ret_pp=TRUE) {
  idx2=cred_set$cs_index
  sus_bf=sus_dat$lbf_variable[idx2,,drop=FALSE][,setdiff(colnames(sus_dat$lbf_variable),"")]
  if(is.vector(sus_bf))
    sus_bf=matrix(sus_bf,nrow=1,dimnames=list(NULL,names(sus_bf)))
  todo <- data.table::as.data.table(expand.grid(i=1:length(querysnpid),j=1:nrow(sus_bf)))
  # todo[,pp4:=0]
  isnps=colnames(sus_bf)
  if(!length(isnps))
    return(data.table::data.table(nsnps=NA))

  if("null" %in% colnames(sus_bf))
    sus_bf=sus_bf - matrix(sus_bf[,"null"],nrow(sus_bf),ncol(sus_bf))

  querypos <- which(colnames(sus_bf)%in%querysnpid)
  sus_bf=sus_bf[,isnps,drop=FALSE]
  results <- vector("list",nrow(todo))
  hit1=querysnpid
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, sus_bf=sus_bf[todo$j[k], ])
    nsnps <- nrow(df)

    lBF.persnp = df$sus_bf
    lBF.Ha <- logsum(lBF.persnp[-querypos])
    lBF.Hc <- lBF.persnp[querypos]
    lBF_df = data.frame(lBF.Ha = lBF.Ha, lBF.Hc = lBF.Hc)
    hit2=names(which.max(sus_bf[todo$j[k],]))
    if (ret_pp){
      pp.bf <- combine.bf(lBF_df, pn=pn, pa=pa, pc=pc)
      if(is.null(hit1)) {
        hit1="-"
        pp.bf$pp[c(1,3)]=c(0,1)
      }
      if(is.null(hit2)) {
        hit2="-"
        pp.bf$pp[c(1,2)]=c(0,1)
      }
      results[[k]]=do.call("data.frame",c(list(nsnps=nsnps, hit1=hit1, hit2=hit2), as.list(pp.bf$pp), as.list(pp.bf$bf), querysnp=querysnpid, typeBF='susieBF'))
    } else {
      results[[k]]=data.frame(lBF.Ha = lBF.Ha, lBF.Hc = lBF.Hc, nsnps=nsnps, querysnp=querysnpid, querytrait=querytrait, hit1=hit1, hit2=hit2, typeBF='susieBF')
    }
  }
  results <- data.table::as.data.table(do.call("rbind",results))
  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=paste(results$hit1,results$hit2,sep=".")
  if(any(duplicated(hits))) {
    results=results[!duplicated(hits)]
  }
  if (ret_pp){
    ret = list(summary=results,
               priors=c(pn=pn,pa=pa, pc=pc))
    ## renumber index to match
    ret$summary[,idx2:=cred_set$cs_index[idx2]]
    ret$summary[,idx1:=rep(1, length(ret$summary$idx2))]
    ret$summary[,querytrait:=querytrait]
    ret$bf=sus_bf
    ret$querysnp=querysnpid
    ret$querytrait=querytrait
    ret$typeBF='susieBF'
    return(ret)
  } else{
    return(results)
  }
}


#' cophe.susie.lbf
#'
#' Calculate log bayes factors for each hypothesis (SuSIE - multiple causal variant assumption)
#' @param dataset a list with specifically named elements defining the query trait dataset
#'   to be analysed.
#' @param querysnpid Id of the query variant, (id in dataset$snp)
# @param querytrait Query trait name
#' @param querytrait Query trait name
#' @param switch Set switch=TRUE to obtain single BF when credible sets not found with SuSIE
#' @param susie.args a named list of additional arguments to be passed to
#'   \link{runsusie}
#' @param MAF Minor allele frequency vector
#' @seealso \code{\link{cophe.susie}}
#' @return data frame with log bayes factors for Hn and Ha hypotheses
#' @export
#' @examples
#' library(cophescan)
#' data(cophe_multi_trait_data)
#' query_trait_1 <- cophe_multi_trait_data$summ_stat[['Trait_1']]
#' query_trait_1$LD <- cophe_multi_trait_data$LD
#' querysnpid <- cophe_multi_trait_data$querysnpid
#' res.susie.lbf <- cophe.susie.lbf(query_trait_1, querysnpid = querysnpid,
#'                                   querytrait='Trait_1', switch=T)
#' res.susie.lbf
#' @author Ichcha Manipur
cophe.susie.lbf <- function(dataset, querysnpid, querytrait, switch=TRUE,
                            susie.args=list(), MAF=NULL) {
  sus_dat = cophe.prepare.dat.susie(dataset, querysnpid, susie.args)
  cred_set = sus_dat$sets
  if (is.null(cred_set$cs) || length(cred_set$cs)==0){
    message('No credible sets found calculation with SuSIE...')
    if (switch){
      message('Switching to cophe.single.bf as switch=TRUE...')
      out_bf <- cophe.single.lbf(dataset, querysnpid, querytrait)
    } else{
      stop('Set switch=TRUE to obtain single BF when credible sets not found with SuSIE')
    }
  } else {
    out_bf=cophe.bf_bf(sus_dat, cred_set, querysnpid, querytrait, ret_pp=FALSE)
  }
  return(out_bf)
}


#' Prepare data for `cophe.susie`
#'
#' @param dataset a list with specifically named elements defining the query trait dataset
#'   to be analysed.
#' @param querysnpid Id of the query variant, (id in dataset$snp)
#' @param susie.args a named list of additional arguments to be passed to
#'   \link{runsusie}
#' @return a list with the output of running susie
#' @seealso \code{\link{cophe.susie}}
#' @keywords internal
#'
cophe.prepare.dat.susie <- function(dataset, querysnpid, susie.args){
  # if(!requireNamespace("susieR", quietly = TRUE)) {
  #   stop("please install susieR https://github.com/stephenslab/susieR")
  # }
  if (!"pvalues" %in% names(dataset)){
    pvalues=pnorm(-abs(dataset$beta/sqrt(dataset$varbeta)))*2
  } else {
    pvalues=dataset$pvalues
  }
  ## N Required in latest SuSIE update
  if (!"N" %in% names(dataset)){
    stop('N (sample size) of the dataset not provided')
  }

  if("susie" %in% class(dataset))
    sus_dat = dataset
  else{
    if(!querysnpid %in% dataset$snp) {
      stop("Please check your dataset, queried trait not in dataset")
    }
  sus_dat=do.call("runsusie", c(list(d=dataset,suffix=1),susie.args))
  }
  ## Remove credible sets with all variants having p-value > 0.1. Once such a set is detected,
  ## all subsequent sets are also removed
  cred_set=sus_dat$sets
  for (cs_id in seq_along(cred_set$cs)){
    if (all(pvalues[cred_set$cs[[cs_id]]] > 0.1)){
      message('Removing credible sets with snp p-vals > 0.1')
      cred_set$cs_index <- cred_set$cs_index[!cred_set$cs_index %in% cs_id:length(cred_set$cs)]
      cred_set$cs[cs_id:length(cred_set$cs)] <- NULL
      break
    }
  }
  return(sus_dat)
}
