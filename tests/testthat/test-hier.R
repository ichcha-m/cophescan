## check if density values from R and Rcpp are the same
log_alpha_R  = dnorm(-9 , -10, 0.5, TRUE)
log_alpha_Rcpp = logd_alpha(-9 , -10, 0.5)
log_beta_R  = dgamma(2 , shape=2, rate=2, log = TRUE)
log_beta_Rcpp = logd_beta(2 , 2, 0.5)
log_gamma_R  = dgamma(2 , shape=2, rate=2, log = TRUE)
log_gamma_Rcpp = logd_beta(2 , 2, 0.5)

testthat::expect_equal(log_alpha_Rcpp, log_alpha_R)
testthat::expect_equal(log_beta_Rcpp, log_beta_R)
testthat::expect_equal(log_gamma_Rcpp, log_gamma_R)

## check pik generation
pars2pik_R = function(pars, nsnps, rg_vec, rg=FALSE) {
  ## gamma=pars[3]
  alpha=pars[1]
  beta=pars[2]
  if(rg)
    gamma=pars[3]
  pik=cbind(1,
            exp(alpha),
            exp(alpha + beta + if(rg) { gamma*rg_vec } else { 0 }))
  denom=diag(pik %*% rbind(1, nsnps - 1, 1))
  pik/denom
}

pars = c(-9, 2, 2)
nsnps = 1000
rg_vec = 0.9

pik_R = pars2pik_R(pars, nsnps, rg_vec, rg=FALSE)
pik_Rcpp = pars2pik(pars, nsnps, rg_vec, rg=FALSE)
pik_R_cov = pars2pik_R(pars, nsnps, rg_vec, rg=TRUE)
pik_Rcpp_cov = pars2pik(pars, nsnps, rg_vec, rg=TRUE)


testthat::expect_equal(pik_Rcpp, pik_R)
testthat::expect_equal(pik_Rcpp_cov, pik_R_cov)

## check log-likelihood
loglik_Rf=function(pars, lbf_mat, nsnps, rg_vec, rg=FALSE) {
  pik=pars2pik_R(pars, nsnps, rg_vec, rg=rg)
  tmp=cbind(log(pik[,1]),
            log(pik[,2]) + lbf_mat$lBF.Ha,
            log(pik[,3]) + lbf_mat$lBF.Hc)
  ll_k=apply(tmp,1, coloc:::logsum)
  sum(ll_k)
  return(ll_k)
}

lbf_mat = data.frame(lBF.Ha=3, lBF.Hc=10)
loglik_R = loglik_Rf(pars, lbf_mat, nsnps, rg_vec, rg=FALSE)
loglik_Rcpp = loglik(pars, as.matrix(lbf_mat), nsnps, rg_vec, rg=FALSE)
loglik_R_cov = loglik_Rf(pars, lbf_mat, nsnps, rg_vec, rg=TRUE)
loglik_Rcpp_cov = loglik(pars, as.matrix(lbf_mat), nsnps, rg_vec, rg=TRUE)

testthat::expect_equal(loglik_Rcpp, loglik_R)
testthat::expect_equal(loglik_Rcpp_cov, loglik_R_cov)
