##' Check if a variant causally associated in one trait might be causal in another trait
##'
##' @title run cophe.susie using susie to detect separate signals
##' @param dataset2 a list with specifically named elements defining the dataset
##'   to be analysed.
##' @param causal.snpid Id of the query variant
##' @param p2a prior probability a SNP other than the causal variant is associated with trait 2, default 1e-4
##' @param p12c prior probability the causal SNP of trait 1 is associated with both traits, default 1e-5
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @param dataset2 *either* an input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param susie.args a named list of additional arguments to be passed to
##'   \link{runsusie}
##' @param ... other arguments passed to \link{cophe.bf_bf}, in particular prior
##'   values for causal association with one trait (p1, p2) or both (p12)
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
cophe.susie=function(dataset2, causal.snpid, p1=1e-4, p2=1e-4, p12=1e-5,
                       p2a=NULL, p12c=NULL,
                       susie.args=list(),  ...) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }

  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=do.call("runsusie", c(list(d=dataset2,suffix=1),susie.args))
  if(!causal.snpid %in% colnames(s2$alpha)) {
    message("please check your dataset, given trait 1 causal snp not present in trait 2 cs")
    return(NULL)
  }

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
  common.snps <- ncol(s2$alpha) - 1
  hp <- hypothesis.priors(p2a, p12c, common.snps)
  print('Hypothesis Priors')
  print(hp)

  cs2=s2$sets

  if (is.null(cs2$cs) || length(cs2$cs)==0 ){
    print('No credible sets found ...')
    print('Switching to cophe.single ...')
    ret <- cophe.single(dataset2, causal.snpid)
  } else if (!any(sub(".*[.]","",names(unlist(s2$sets$cs))) %in% causal.snpid) & length(s2$sets$cs) < 2){
    print('Queried SNP not in the susie credible set ...')
    print('Switching to cophe.single ...')
    ret <- cophe.single(dataset2, causal.snpid)
  } else {
    if (!any(sub(".*[.]","",names(unlist(s2$sets$cs))) %in% causal.snpid)){
      print('Queried SNP not in the susie credible sets ...')
    }
    print('Running cophe.susie...')
    isnps=colnames(s2$alpha)
    if(!length(isnps))
      return(data.table(nsnps=NA))
    message("Using ",length(isnps), " and ",ncol(s2$alpha)-1," available")

    idx2=cs2$cs_index #sapply(names(cs2$cs), get_idx)

    bf2=s2$lbf_variable[idx2,,drop=FALSE][,setdiff(colnames(s2$lbf_variable),"")]

    ret=cophe.bf_bf(bf2,causal.snpid, p2a, p12c,p1=p1,p2=p2,p12=p12)
    ## renumber index to match
    ret$summary[,idx2:=cs2$cs_index[idx2]]
    ret$summary[,idx1:=rep(1, length(ret$summary$idx2))]
    ret$bf=bf2
  }
  return(ret)
}


##' a dataset represented by Bayes factors
##' @title extract data through Bayes factors
##' @param bf2 named vector of BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param causal.snpid Id of the query variant
##' @param p2a prior probability a SNP other than the causal variant is associated with trait 2, default 1e-4
##' @param p12c prior probability the causal SNP of trait 1 is associated with both traits, default 1e-5
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @return bayes factors of signals
cophe.bf_bf=function(bf2, causal.snpid, p2a, p12c,p1, p2, p12) {

  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  todo <- as.data.table(expand.grid(i=1:length(causal.snpid),j=1:nrow(bf2)))
  # todo[,pp4:=0]
  isnps=colnames(bf2)
  if(!length(isnps))
    return(data.table(nsnps=NA))

  if("null" %in% colnames(bf2))
    bf2=bf2 - matrix(bf2[,"null"],nrow(bf2),ncol(bf2))
  causalpos1 <- which(colnames(bf2)%in%causal.snpid)
  bf2=bf2[,isnps,drop=FALSE]
  results <- vector("list",nrow(todo))
  hit1=causal.snpid
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf2=bf2[todo$j[k], ])
    pp.bf <- combine.bf.kc(df$bf2, p2a, p12c, causalpos1)
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
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2), as.list(pp.bf$pp), as.list(pp.bf$bf), querysnp=causal.snpid))
  }
  results <- as.data.table(do.call("rbind",results))
  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=paste(results$hit1,results$hit2,sep=".")
  if(any(duplicated(hits))) {
    results=results[!duplicated(hits)]
  }
  list(summary=results,
       priors=c(p1=p1,p2=p2,p12=p12, p2a=p2a, p12c=p12c))
}

