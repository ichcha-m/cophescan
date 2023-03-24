##' Check if a variant causally associated in one trait might be causal in another trait
##'
##' @title run cophe.susie using susie to detect separate signals
##' @param dataset a list with specifically named elements defining the dataset
##'   to be analysed.
##' @param query.snpid Id of the query variant
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4 (coloc prior)
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4 (coloc prior)
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5 (coloc prior)
##' @param pn prior probability that none of the SNPs/variants in the region are associated with the query trait , default \eqn{pn = 1 - (pa*(nsnps-1)) - pc} (cophescan prior)
##' @param pa prior probability that a non-query variant is causally associated with the query trait , default \eqn{pa = p2} (cophescan prior)
##' @param pc prior probability that the query variant is causally associated with the query trait, default \eqn{pc =  p12/p1+p12} (cophescan prior)
##' @param dataset *either* an input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param susie.args a named list of additional arguments to be passed to
##'   \link{runsusie}
##' @return a list, containing elements
##' * summary a data.table of posterior
##'   probabilities of each global hypothesis, one row per pairwise comparison
##'   of signals from the two traits
##' * results a data.table of detailed results giving the posterior probability
##'   for each snp to be jointly causal for both traits *assuming Hc is true*.
##'   Please ignore this column if the corresponding posterior support for H4
##'   is not high.
##' * priors a vector of the priors used for the analysis
##' @importFrom coloc runsusie
##' @export
##' @author Ichcha Manipur
cophe.susie=function(dataset, query.snpid, p1=1e-4, p2=1e-4, p12=1e-5,
                       pa=NULL, pc=NULL,
                       susie.args=list()) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }
  if (!"pvalues" %in% names(dataset)){
    pvalues=pnorm(-abs(dataset$beta/sqrt(dataset$varbeta)))*2
  } else{
    pvalues=dataset$pvalues
  }
  ## N Required in latest SuSIE update
  if (!"N" %in% names(dataset)){
    stop('N (sample size) of the dataset not provided')
  }

  if("susie" %in% class(dataset))
    s2=dataset
  else
    s2=do.call("runsusie", c(list(d=dataset,suffix=1),susie.args))
  if(!query.snpid %in% colnames(s2$alpha)) {
    message("Please check your dataset, given trait 1 causal snp not present in trait 2 cs")
    return(NULL)
  }
  # number of snps in the region
  common.snps <- ncol(s2$alpha) - 1

  psp  <-  per.snp.priors(nsnps = common.snps, p1 = p1, p2 = p2, p12 = p12, pa = pa, pc = pc)
  print('SNP Priors')
  print(psp)

  hp <- hypothesis.priors(nsnps = common.snps, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])
  print('Hypothesis Priors')
  print(hp)

  ## Remove credible sets with all variants having p-value > 0.1. Once such a set is detected,
  ## all subsequent sets are also removed
  cs2=s2$sets
  for (cs_id in seq_along(cs2$cs)){
     if (all(pvalues[cs2$cs[[cs_id]]] > 0.1)){
       print('Removing credible sets with snp p-vals > 0.1')
       cs2$cs_index <- cs2$cs_index[!cs2$cs_index %in% cs_id:length(cs2$cs)]
       cs2$cs[cs_id:length(cs2$cs)] <- NULL
       break
     }
  }

  if (is.null(cs2$cs) || length(cs2$cs)==0 ){
    print('No credible sets found ...')
    print('Switching to cophe.single ...')
    ret <- cophe.single(dataset, query.snpid, pa=psp[["pa"]], pc=psp[["pc"]])
  }
  # else if (!any(sub(".*[.]","",names(unlist(s2$sets$cs))) %in% query.snpid) & length(s2$sets$cs) < 2){
  #   print('Queried SNP not in the susie credible set ...')
  #   print('Switching to cophe.single ...')
  #   ret <- cophe.single(dataset, query.snpid, pa=psp[["pa"]], pc=psp[["pc"]])
  # }
  else {
    if (!any(sub(".*[.]","",names(unlist(s2$sets$cs))) %in% query.snpid)){
      print('Queried SNP not in the susie credible sets ...')
    }
    print('Running cophe.susie...')
    isnps=colnames(s2$alpha)
    if(!length(isnps))
      return(data.table(nsnps=NA))
    message("Using ",length(isnps), " and ",ncol(s2$alpha)-1," available")

     #sapply(names(cs2$cs), get_idx)
    idx2=cs2$cs_index
    bf2=s2$lbf_variable[idx2,,drop=FALSE][,setdiff(colnames(s2$lbf_variable),"")]

    ret=cophe.bf_bf(bf2, query.snpid, pn=psp[["pn"]], pa=psp[["pa"]], pc=psp[["pc"]])

    ## renumber index to match
    ret$summary[,idx2:=cs2$cs_index[idx2]]
    ret$summary[,idx1:=rep(1, length(ret$summary$idx2))]
    ret$bf=bf2
    ret$querysnp=query.snpid
  }
  attr(ret, "class") <- "cophe"
  return(ret)
}


##' a dataset represented by Bayes factors
##' @title extract data through Bayes factors
##' @param bf2 named vector of BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param query.snpid Id of the query variant
##' @param pn prior probability that none of the SNPs/variants in the region are associated with the query trait
##' @param pa prior probability that a non-query variant is causally associated with the query trait
##' @param pc prior probability that the query variant is causally associated with the query trait
##' @return bayes factors of signals
cophe.bf_bf=function(bf2, query.snpid, pn, pa, pc) {

  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  todo <- as.data.table(expand.grid(i=1:length(query.snpid),j=1:nrow(bf2)))
  # todo[,pp4:=0]
  isnps=colnames(bf2)
  if(!length(isnps))
    return(data.table(nsnps=NA))

  if("null" %in% colnames(bf2))
    bf2=bf2 - matrix(bf2[,"null"],nrow(bf2),ncol(bf2))
  querypos <- which(colnames(bf2)%in%query.snpid)
  bf2=bf2[,isnps,drop=FALSE]
  results <- vector("list",nrow(todo))
  hit1=query.snpid
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf2=bf2[todo$j[k], ])
    pp.bf <- combine.bf.kc(df$bf2, pn=pn, pa=pa, pc=pc, querypos=querypos)
    common.snps <- nrow(df)
    if(is.null(hit1)) {
      hit1="-"
      pp.bf$pp[c(1,3)]=c(0,1)
    }
    hit2=names(which.max(bf2[todo$j[k],]))
    if(is.null(hit2)) {
      hit2="-"
      pp.bf$pp[c(1,2)]=c(0,1)
    }
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2), as.list(pp.bf$pp), as.list(pp.bf$bf), querysnp=query.snpid))
  }
  results <- as.data.table(do.call("rbind",results))
  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=paste(results$hit1,results$hit2,sep=".")
  if(any(duplicated(hits))) {
    results=results[!duplicated(hits)]
  }
  list(summary=results,
       priors=c(pn=pn,pa=pa, pc=pc))
}

