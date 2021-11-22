#' run_metrop_priors
#'
#' @param multi.dat matrix of bf values
#' @param rg whether to run rg
#' @param rg_vec vector of genetic correlations
#' @param nits number of iterations
#' @param thin burnin
#'
#' @return list containing posterior of the parameters
#' @export
#' @import matrixStats
#' @examples
run_metrop_priors <- function(multi.dat, rg=FALSE, rg_vec=NULL, nits=10000, thin=1000){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  if (is.null(rg_vec)){
    rg_vec <- 0
  }
  # data=data.frame(lbfak=pp_df[,"bf.a"],
  #                 lbfck=pp_df[,"bf.c"],
  #                 nsnps=pp_df[,"nsnps"],
  #                 rg=rg_vec)
  N=length(data[[1]]) # number of observations

  ## pars2pik and logpriors both have arguments rg=FALSE/TRUE which determine which kind of model is running.
  ## other functions (loglik, posterior, target) use ... to pass the rg on to those functions

  sample_alpha=function(n=1)
    rnorm(n, mean=-10, sd=0.5)
  sample_beta=function(n=1)
    rgamma(n, 0.5,0.5)
  sample_gamma=function(n=1)
    rgamma(n, 0.5,0.5)
  logd_alpha=function(a,...)
    dnorm(a, mean=-10, sd=0.5,log=TRUE)
  logd_beta=function(b)
    dgamma(b, 0.5,0.5,log=TRUE)
  logd_gamma=function(g)
    dgamma(g, 0.5,0.5,log=TRUE)
  ## logpriors_norg=function(pars)
  ##   logd_alpha(pars[1]) + logd_beta(pars[2])
  logpriors=function(pars, rg=FALSE)
    logd_alpha(pars[1]) + logd_beta(pars[2]) + if(rg) { logd_gamma(pars[3]) } else { 0 }


  target=function(pars,...)
    loglik(pars, data[, c('lbfak', 'lbfck')], data$nsnps, rg_vec, rg=RG) + logpriors(pars, ...)

  ## check no error
  # loglik(c(1,1,1))
  # target(c(1,1,1))

  ## proposal
  propose=function(curval,propsd=0.5)
    curval + rnorm(length(curval),mean=0,sd=propsd)


  ## init
  RG=rg
  alpha = sample_alpha()
  beta = sample_beta()
  if(RG) {
    gamma=sample_gamma()
    pars=c(alpha,beta,gamma)
  } else {
    pars=c(alpha,beta)
  }

  ## storage
  nits=if(RG) { 100000 # in simulated data, need 100,000 for convergence when RG==TRUE
  } else { 10000 }
  thin=if(RG) { 10 } else { 1 } # so store 10,000 samples whatever
  params=matrix(0, length(pars), nits/thin)
  ll=numeric(nits/thin)
  rownames(params)=c("alpha","beta","gamma")[1:length(pars)]

  ## start sampling
  L0=loglik(pars, data[, c('lbfak', 'lbfck')], data$nsnps, rg_vec, rg=RG)
  T0=L0 + logpriors(pars, rg=RG)
  for(i in 1:nits) {
    ## message()
    if(nits %% thin == 0) {
      params[,i/thin]=c(pars)
      ll[i/thin]=T0
    }
    newpars=propose(c(pars))
    accept=exp(target(newpars,rg=RG) - T0)
    ## print(accept)
    if(runif(1) < accept) {
      pars=newpars
      L0=loglik(pars,rg=RG)
      T0=L0 + logpriors(pars,rg=RG)
    }
  }

  ab_est=rowMeans(params[ , -c(1:(ncol(params)/2))])
  head(pars2pik(ab_est, data$nsnps, rg_vec, rg = RG), 1)

  ## calculate posteriors
  post=lapply(1:ncol(params), function(i)
    posterior(params[,i], data[, c('lbfak', 'lbfck')], data$nsnps, rg_vec, rg=RG))
  piks=lapply(1:ncol(params), function(i)
    pars2pik(params[,i], data$nsnps, rg_vec, rg = RG))

  ## average posterior
  st=(nits/thin/2+1)
  en=nits/thin
  avpost=post[[ st ]]
  for(i in (st+1):en)
    avpost=avpost+post[[ i ]]
  avpost=avpost/(en - st + 1)
  ret = list(ll=ll, piks=piks)
  return(ret)
}
