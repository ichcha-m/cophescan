#' Run the hierarchical Metropolis Hastings model to infer priors
#'
#' @param multi.dat matrix of bf values, rows=traits, named columns=("lBF.Ha","lBF.Hc","nsnps")
#' @param covar whether to include covariates
#' @param covar_vec vector of covariates
#' @param is_covar_categorical only two categories supported (default=FALSE) - Experimental
#' @param nits number of iterations
#' @param thin burnin
#' @param posterior default: FALSE, estimate posterior probabilities of the hypotheses
#' @param avg_pik default: FALSE, estimate the average of the pik
#' @param avg_posterior default: FALSE, estimate the average of the posterior probabilities of the hypotheses
#' @param pik default: FALSE, inferred prior probabilities
#' @param alpha_mean prior for the mean of  alpha
#' @param alpha_sd prior for the standard deviation of  alpha
#' @param beta_shape prior for the shape (gamma distibution) of beta
#' @param beta_scale prior for the scale of beta
#' @param gamma_shape prior for the shape (gamma distibution) of gamma
#' @param gamma_scale prior for the scale of gamma
#' @return List containing the posterior distribution of the parameters alpha, beta, gamma (if covariate included) and the loglikelihood
#' @return if avg_posterior=TRUE matrix with average of all the posterior probabilities of Hn, Ha and Hc
#' @return if avg_pik=TRUE matrix with average of all the priors: pn, pa and pc
#' @return data, nits and thin contain the input data, number of iterations and burnin respectively specified for the hierarchical model
#' @export
run_metrop_priors <- function(multi.dat, covar=FALSE, covar_vec=NULL, is_covar_categorical=FALSE, nits=10000,
                              thin=1, posterior=FALSE, avg_pik=TRUE, avg_posterior=TRUE,
                              pik=FALSE, alpha_mean =-10,
                               alpha_sd=0.5,  beta_shape=2,  beta_scale=2,
                               gamma_shape=2,  gamma_scale=2){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  if ('hit1'%in% colnames(pp_df)){
    pp_df$sus_labels <- rownames(pp_df) <- paste0(pp_df$querytrait, '_', pp_df$hit1, '_hit_', pp_df$hit2)
    pp_df$sus_labels[which(pp_df$typeBF=="ABF")] <- NA
  }
  if (covar) {
    if (!(length(covar_vec) == nrow(pp_df) )){
      stop('Length of covar_vec should be equal to the number of traits (length(multi.dat))')
    } else if (is.null(covar_vec)){
        stop('covar set to true but covar_vec not supplied')
    } else if (!is.vector(covar_vec)){
      stop('covar_vec has to be a vector')
    }
  }
  if (covar & is.character(covar_vec)){
    if (is_covar_categorical == FALSE){
      stop('covar vector of type character but is_covar_categorical set to FALSE. Note: only 2 categories supported')
    } else if (length(unique(covar_vec))>2){
      stop('Only two categories supported')
    } else {
    covar_vec_categ = covar_vec
    covar_vec = as.numeric(as.factor(covar_vec ))
    }

  }
  if (!covar & is.null(covar_vec)){
      message('covar_vec supplied but covar set to false, setting covar_vec to 1. Set covar to True if covariate to be included')
      covar_vec <- rep(1, nrow(pp_df))
  }

  lbf_mat <- as.matrix(pp_df[, c('lBF.Ha', 'lBF.Hc')])
  nsnps <- as.vector(pp_df$nsnps)
  message(paste0('Running the cophescan hierarchical model... \nNumber of iterations: ', nits,'\nBurn-in: ', thin, '\nNumber of input QV/QT pairs: ', nrow(pp_df), '\nCovariate included: ', covar ))

  ## Run hierarchical model
  res.metrop <- metrop_run(lbf_mat = lbf_mat, nsnps = nsnps,
                           covar_vec = covar_vec, covar = covar, nits = nits, thin = thin,  alpha_mean =alpha_mean,
                           alpha_sd=alpha_sd,  beta_shape=beta_shape,  beta_scale=beta_scale,
                           gamma_shape=gamma_shape,  gamma_scale=gamma_scale)

  if (covar){
    rownames(res.metrop$parameters) <- c("alpha","beta","gamma")
  }else{
    rownames(res.metrop$parameters) <- c("alpha","beta")}
  ## piks from parameters
  if (pik){
    pik <- piks(res.metrop$parameters, nsnps = nsnps, covar_vec = covar_vec, covar=covar)
    res.metrop$pik <- pik
  }

  if (posterior){
    posterior <- posterior_prob(res.metrop$parameters, lbf_mat, nsnps = nsnps,
                                covar_vec = covar_vec, covar=covar)
    res.metrop$posterior <- posterior
  }

  if (avg_posterior){
    avg.posterior <- average_posterior_prob(res.metrop$parameters, lbf_mat,
                                            nsnps = nsnps, covar_vec = covar_vec,
                                            nits = nits, thin = thin, covar = covar)
    colnames(avg.posterior) <- c('PP.Hn', 'PP.Ha', 'PP.Hc')
    res.metrop$avg.posterior <- avg.posterior
  }

  if (avg_pik){
    avg.pik <- average_piks(res.metrop$parameters, nsnps = nsnps,
                            covar_vec = covar_vec, nits = nits, thin = thin, covar = covar)
    colnames(avg.pik) <- c('pnk', 'pak', 'pck')
    res.metrop$avg.pik <- avg.pik
  }
  res.metrop$data = pp_df %>% dplyr::select_if(!grepl('PP.Hn|PP.Ha|PP.Hc', colnames(pp_df)))
  res.metrop$nits = nits
  res.metrop$thin = thin

  if (covar){
    res.metrop$covar_vec = covar_vec
    if (is_covar_categorical){
      res.metrop$covar_vec_categ = covar_vec_categ
    }
  }
  return(res.metrop)
}


