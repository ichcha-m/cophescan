#' run_metrop_priors
#'
#' @param multi.dat matrix of bf values
#' @param rg whether to run rg
#' @param rg_vec vector of genetic correlations
#' @param nits number of iterations
#' @param thin burnin
#' @param posterior default: F, estimate posterior probabilities of the hypotheses
#' @param avg_posterior default: F, estimate the average of the posterior probabilities of the hypotheses
#' @param pik default: F, inferred prior probabilities
#'
#' @return list containing posterior of the parameters
#' @export
run_metrop_priors <- function(multi.dat, rg=FALSE, rg_vec=NULL, nits=10000, thin=1, posterior=F, avg_posterior=F, pik=T){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  if (rg)  {
    if (!(length(rg_vec) == nrow(pp_df) )){
      stop('Length of rg_vec should be equal to the number of traits (length(multi.dat))')
    }
  } else if (is.null(rg_vec)){
    rg_vec <- 0
  }

  lbf_mat <- as.matrix(pp_df[, c('lbfak', 'lbfck')])
  nsnps <- pp_df[, "nsnps"]
  res.metrop <- metrop_run(lbf_mat = lbf_mat, nsnps = nsnps, rg_vec = rg_vec, rg = rg, nits = nits, thin = thin)
  if (rg){
  rownames(res.metrop$params) <- c("alpha","beta","gamma")
  }else{
  rownames(res.metrop$params) <- c("alpha","beta")}
  if (pik){
    pik <- piks(res.metrop$params, nsnps = nsnps, rg_vec = rg_vec, rg=rg)
    res.metrop$pik <- pik
  }

  if (posterior){
    posterior <- post_prob(res.metrop$params, lbf_mat, nsnps = nsnps, rg_vec = rg_vec, rg=rg)
    res.metrop$posterior <- posterior
  }

  if (avg_posterior){
    avg.posterior <- average_post(posterior, nits, thin)
    colnames(avg.posterior) <- c('Hn', 'Hc', 'Ha')
    res.metrop$avg.posterior <- avg.posterior
  }

  return(res.metrop)
}

