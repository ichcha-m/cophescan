#' Run the hierarchical metropolis hastings model to infer priors
#'
#' @param multi.dat matrix of bf values, rows=traits, named columns=("lBF.Ha","lBF.Hc","nsnps")
#' @param rg whether to run rg
#' @param rg_vec vector of genetic correlations
#' @param nits number of iterations
#' @param thin burnin
#' @param posterior default: F, estimate posterior probabilities of the hypotheses
#' @param avg_pik default: F, estimate the average of the pik
#' @param avg_posterior default: F, estimate the average of the posterior probabilities of the hypotheses
#' @param pik default: F, inferred prior probabilities
#' @param chains default: 1, number of chains
#' @param cores default: 1, number of cores
#' @param alpha_mean prior for the mean of  alpha
#' @param alpha_sd prior for the standard deviation of  alpha
#' @param beta_shape prior for the shape (gamma distibution) of beta
#' @param beta_scale prior for the scale of beta
#' @param gamma_shape prior for the shape (gamma distibution) of gamma
#' @param gamma_scale prior for the scale of gamma
#' @return matrix with average of all the posterior probabilities: Hn, Ha and Hc
#' @return list containing posterior of the parameters
#' @export
run_metrop_priors <- function(multi.dat, rg=FALSE, rg_vec=NULL, nits=10000,
                              thin=1, posterior=F, avg_pik=T, avg_posterior=T,
                              pik=F, chains=1, cores=1, alpha_mean =-10,
                               alpha_sd=0.5,  beta_shape=2,  beta_scale=2,
                               gamma_shape=2,  gamma_scale=2){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  if (rg) {
    if (!(length(rg_vec) == nrow(pp_df) )){
      stop('Length of rg_vec should be equal to the number of traits (length(multi.dat))')
    } else if (is.null(rg_vec)){
        stop('rg set to true but rg_vec not supplied')
    }
  }
  if (!rg & is.null(rg_vec)){
      message('rg_vec supplied but rg set to false, setting rg_vec to 1')
      rg_vec <- rep(1, nrow(pp_df))
  }

  lbf_mat <- as.matrix(pp_df[, c('lBF.Ha', 'lBF.Hc')])
  nsnps <- pp_df[, "nsnps"]

  ## Run hierarchical model
  res.metrop <- metrop_run(lbf_mat = lbf_mat, nsnps = nsnps,
                           rg_vec = rg_vec, rg = rg, nits = nits, thin = thin,  alpha_mean =alpha_mean,
                           alpha_sd=alpha_sd,  beta_shape=beta_shape,  beta_scale=beta_scale,
                           gamma_shape=gamma_shape,  gamma_scale=gamma_scale)

  if (rg){
    rownames(res.metrop$parameters) <- c("alpha","beta","gamma")
  }else{
    rownames(res.metrop$parameters) <- c("alpha","beta")}
  ## piks from parameters
  if (pik){
    pik <- piks(res.metrop$parameters, nsnps = nsnps, rg_vec = rg_vec, rg=rg)
    res.metrop$pik <- pik
  }

  if (posterior){
    posterior <- posterior_prob(res.metrop$parameters, lbf_mat, nsnps = nsnps,
                                rg_vec = rg_vec, rg=rg)
    res.metrop$posterior <- posterior
  }

  if (avg_posterior){
    avg.posterior <- average_posterior_prob(res.metrop$parameters, lbf_mat,
                                            nsnps = nsnps, rg_vec = rg_vec,
                                            nits = nits, thin = thin, rg = rg)
    colnames(avg.posterior) <- c('Hn', 'Ha', 'Hc')
    res.metrop$avg.posterior <- avg.posterior
  }

  if (avg_pik){
    avg.pik <- average_piks(res.metrop$parameters, nsnps = nsnps,
                            rg_vec = rg_vec, nits = nits, thin = thin, rg = rg)
    colnames(avg.pik) <- c('p0k', 'pak', 'pck')
    res.metrop$avg.pik <- avg.pik
  }

  return(res.metrop)
}

#' Calculate posterior probabilities from priors, given logABFs for each SNP
#' and priors inferred from the hierarchical model
#'
#' @title combine.bf.kc.hier
#' @param pik_vec named vector (p0k, pak, pck)
#' @param lbfk_vec named log bayes factor vector (lBF.Ha, lBF.Hc)
#' @return named numeric vector of posterior probabilities and bayes factors
#' @author Ichcha Manipur
#' @export
combine.bf.kc.hier <- function(pik_vec, lbfk_vec) {

  lHn.bf <- 0
  lHa.bf <- (log(pik_vec$pak) - log(pik_vec$p0k)) +  lbfk_vec$lBF.Ha
  lHc.bf <- (log(pik_vec$pck) - log(pik_vec$p0k)) + lbfk_vec$lBF.Hc

  all.bf <- c(lHn.bf, lHa.bf, lHc.bf)
  my.denom.log.bf <- coloc:::logsum(all.bf)
  pp <- exp(all.bf - my.denom.log.bf)
  # pp.bf
  names(pp) <- paste("PP.H", c('n', 'a', 'c') , sep = "")
  print(signif(pp,3))
  print(paste("Priors: ", pik_vec))
  print(paste("PP for shared variant: ", signif(pp["PP.Hc"],3)*100 , '%',
              sep = ''))
  return(list(pp = pp, bf = lbfk_vec, pik = pik_vec))
}
