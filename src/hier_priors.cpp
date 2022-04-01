// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;

//' Conversion of parameters α, β and γ to p0k, pak and pck
//'
//' @param params Vector of parameters: α, β and γ
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: should the covariate information be used? default: False
//' @return pik matrix of priors: p0k, pak and pck
// [[Rcpp::export]]
arma::mat pars2pik(arma::vec params, NumericVector nsnps, NumericVector rg_vec, bool rg=false) {
  double alpha = params[0];
  double beta = params[1];
  double gamma;
  if (rg==false){
    gamma=0;
  } else {
    gamma=params[2];
  }
  int k = nsnps.length();
  arma::mat pik = arma::ones(k , 3);
  pik.col(1) =  arma::vec (rep( exp(alpha), k )) ;
  pik.col(2) = arma::vec (exp(alpha + beta + (gamma*rg_vec)));

  arma::vec v = nsnps - 1;
  arma::mat m2 = arma::ones(3, k);
  m2.row(1) = arma::trans(v);
  arma::vec denom = diagvec(pik * m2);
  pik.each_col() /= denom ;
  return pik;
}

//' Log sum
//' @param x vector of log scale values to be added
// [[Rcpp::export]]
double logsumexp(arma::rowvec x) {
  return(log(sum(exp(x - x.max()))) + x.max());
}

//' Log posterior calculation
//'
//' @param params Vector of parameters: α, β and γ
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: should the covariate inflormation be used? default: False
//' @return logpost flog of the posteriors
// [[Rcpp::export]]
arma::mat logpost(arma::vec params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false) {
  arma::mat pik = pars2pik(params, nsnps, rg_vec=rg_vec, rg=rg);
  int k = nsnps.length();
  arma::mat logpost = arma::ones(k , 3);
  arma::mat logpik = log(pik);

  logpost.col(0) =  logpik.col(0);
  logpost.col(1) =  logpik.col(1) + lbf_mat.col(0);
  logpost.col(2) =  logpik.col(2) + lbf_mat.col(1);

  return logpost;
}

//' Log likelihood calculation
//'
//' @param params Vector of parameters: α, β and γ
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: should the covariate inflormation be used? default: False
//' @return logpost flog of the posteriors
// [[Rcpp::export]]
double loglik(arma::vec params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false) {
  arma::mat logpos = logpost(params, lbf_mat, nsnps, rg_vec, rg=rg);
  int k = nsnps.length();
  arma::vec ll_k(k);
  for (int i = 0; i < k; i++){
    arma::rowvec y = logpos.row(i);
    ll_k(i) = logsumexp(y);
  }
  double loglik = accu(ll_k);
  return loglik;
}

//' Calculation of the posterior prob of Hn, Ha and Hc
//'
//' @param params Vector of parameters: α, β and γ
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: should the covariate inflormation be used? default: False
//' @return posterior prob of Hn, Ha and Hc
// [[Rcpp::export]]
arma::mat get_posterior_prob(arma::vec params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false) {
  arma::mat pik = pars2pik(params, nsnps, rg_vec=rg_vec, rg=rg);
  int k = nsnps.length();
  arma::mat logpost = arma::ones(k , 3);
  arma::mat logpik = log(pik);

  logpost.col(0) =   logpik.col(0) - logpik.col(0);
  logpost.col(1) =  (logpik.col(1) - logpik.col(0)) + lbf_mat.col(0);
  logpost.col(2) =  (logpik.col(2) - logpik.col(0)) + lbf_mat.col(1);

  arma::vec denom(k);
  for (int i = 0; i < k; i++){
    arma::rowvec y = logpost.row(i);
    denom(i) = logsumexp(y);
  }
  logpost.each_col() -= denom ;
  arma::mat post_prob = exp(logpost);
  return post_prob;
}

//' sample alpha
//' @param alpha_mean prior for the mean of  alpha
//' @param alpha_sd prior for the standard deviation of  alpha
// [[Rcpp::export]]
arma::vec sample_alpha(double alpha_mean=-10, double alpha_sd=-0.5){
  return rnorm(1, alpha_mean, alpha_sd);
}

//' sample beta
//' @param beta_shape prior for the shape (gamma distibution) of beta
//' @param beta_scale prior for the scale of beta
// [[Rcpp::export]]
arma::vec sample_beta(double beta_shape=2, double beta_scale=2){
  // scale set to 2, to correspond to the R function where we set 0.5 for the rate (scale = 1/0.5)
  return rgamma(1, beta_shape, beta_scale);
}

//' sample gamma
//' @param gamma_shape prior for the shape (gamma distibution) of gamma
//' @param gamma_scale prior for the scale of gamma
// [[Rcpp::export]]
arma::vec sample_gamma( double gamma_shape=2, double gamma_scale=2){
  // scale converted from  required rate(0.5)
  return rgamma(gamma_shape, gamma_scale);
}

//' dnorm for alpha
//' @param a current alpha
//' @param alpha_mean prior for the mean of  alpha
//' @param alpha_sd prior for the standard deviation of  alpha
// [[Rcpp::export]]
double logd_alpha(double a, double alpha_mean=-10, double alpha_sd=0.5){
  return R::dnorm(a, alpha_mean, alpha_sd, true);
}

//' dgamma for beta
//' @param b current beta
//' @param beta_shape prior for the shape (gamma distibution) of beta
//' @param beta_scale prior for the scale of beta
// [[Rcpp::export]]
double logd_beta(double b, double beta_shape=2, double beta_scale=2){
  // scale set to 2, to correspond to the R function where we set 0.5 for the rate (scale = 1/0.5)
  // Also for logd_gamma
  return R::dgamma(b, beta_shape, beta_scale, true);
}

//' dgamma for gamma
//' @param g current gamma
//' @param gamma_shape prior for the shape (gamma distibution) of gamma
//' @param gamma_scale prior for the scale of gamma
// [[Rcpp::export]]
double logd_gamma(double g, double gamma_shape=2, double gamma_scale=2){
  return R::dgamma(g, gamma_shape, gamma_scale, true);
}

// [[Rcpp::export]]
double logpriors(arma::vec params, bool rg=false, double alpha_mean =-10, double alpha_sd=0.5, double beta_shape=2, double beta_scale=2,
                 double gamma_shape=2, double gamma_scale=2) {
  double loggamma;
  if (rg == false){
    loggamma = 0;
  } else{
    loggamma = logd_gamma(params[2], gamma_shape, gamma_scale);
  }
  double logprior = logd_alpha(params[0],  alpha_mean, alpha_sd) + logd_beta(params[1], beta_shape, beta_scale) + loggamma;
  return logprior;
}

// [[Rcpp::export]]
double target(arma::vec params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false) {
  // Target distribution
  double target = loglik(params, lbf_mat, nsnps, rg_vec, rg) + logpriors(params, rg);
  return target;
}

// [[Rcpp::export]]
arma::vec propose(arma::vec params, double propsd=0.5){
  arma::vec propose = params + arma::vec(rnorm(params.n_elem, 0, propsd));
  return propose;
}

// [[Rcpp::export]]
arma::vec pars_init(bool rg=false, double alpha_mean =-10, double alpha_sd=0.5, double beta_shape=2, double beta_scale=2,
                    double gamma_shape=2, double gamma_scale=2){
  double alpha = arma::as_scalar(sample_alpha(alpha_mean, alpha_sd));
  double beta = arma::as_scalar(sample_beta(beta_shape, beta_scale));
  arma::vec params;

  if(rg==true) {
    double gamma = arma::as_scalar(sample_gamma(gamma_shape, gamma_scale));
    params = {alpha, beta, gamma};
  } else{
    params = {alpha, beta};
  }
  return params;
}

//' Run the hierarchical mcmc model to infer priors
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param nits Number of iterations run in mcmc
//' @param thin thinning
//' @param rg logical: Should the covariate inflormation be used? default: False
//' @param alpha_mean prior for the mean of  alpha
//' @param alpha_sd prior for the standard deviation of  alpha
//' @param beta_shape prior for the shape (gamma distibution) of beta
//' @param beta_scale prior for the scale of beta
//' @param gamma_shape prior for the shape (gamma distibution) of gamma
//' @param gamma_scale prior for the scale of gamma
//' @return matrix with average of all the posterior probabilities: Hn, Ha and Hc
// [[Rcpp::export]]
List metrop_run(arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false, int nits=10000,
                int thin=1, double alpha_mean =-10, double alpha_sd=0.5, double beta_shape=2, double beta_scale=2,
                double gamma_shape=2, double gamma_scale=2){
  arma::vec pars = pars_init(rg=rg, alpha_mean, alpha_sd, beta_shape, beta_scale,
                            gamma_shape, gamma_scale);
  arma::mat params = arma::zeros(pars.n_elem, nits/thin);
  arma::vec ll(nits/thin);
  double L0;
  double T0;
  double accept;
  arma::vec newpars;
  L0 = loglik(pars, lbf_mat, nsnps, rg_vec,  rg=rg);
  T0 = L0 + logpriors(pars, rg=rg);
  for(int i = 0; i < nits; i++){
    if(i % thin == 0) {
      params.col(i/thin) = pars;
      ll(i/thin) = T0;
    }
    newpars=propose(pars);
    accept = exp(target(newpars, lbf_mat, nsnps, rg_vec,  rg=rg) - T0);
    // Rcout <<  accept << "\n";
    if (R::runif(0, 1) < accept) {
      pars = newpars;
      L0 = loglik(pars, lbf_mat, nsnps, rg_vec,  rg=rg);
      T0 = L0 + logpriors(pars, rg=rg);
    }
  }
  return Rcpp::List::create(Rcpp::Named("ll") = ll,
                            Rcpp::Named("parameters") = params);
}

//' Average of posterior probabilities: Hn, Ha and Hc
//'
//' @param params Vector of parameters: α, β and γ
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: was the covariate inflormation  used? default: False
//' @return params List of posterior probabilties (len: iterations): Hn, Ha and Hc
// [[Rcpp::export]]
List posterior_prob(arma::mat params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, bool rg=false){
  double k = params.n_cols;
  List post(k);
  for (int i = 0; i < k; i++){
    post[i] = get_posterior_prob(params.col(i),  lbf_mat, nsnps, rg_vec,  rg=rg);
    colnames(post[i]) = CharacterVector::create("Hn", "Ha", "Hc");
  }
  return post;
}

//' Average of posterior probabilities: Hn, Ha and Hc
//'
//' @param params Vector of parameters: α, β and γ
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param rg logical: was the covariate inflormation  used? default: False
//' @return List of priors (len: iterations): p0k, pak and pck
// [[Rcpp::export]]
List piks(arma::mat params, NumericVector nsnps, NumericVector rg_vec, bool rg=false){
  double k = params.n_cols;
  List piks(k);
  for (int i = 0; i < k; i++){
    piks[i] = pars2pik(params.col(i), nsnps, rg_vec,  rg=rg);
    colnames(piks[i]) = CharacterVector::create("p0k", "pak", "pck");
  }
  return piks;
}

//' Average of posterior probabilities: Hn, Ha and Hc
//'
//' @param params Vector of parameters: α, β and γ
//' @param lbf_mat matrix of log bayes factors: lBF.Ha and lBF.Hc
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param nits Number of iterations run in mcmc
//' @param thin thinning
//' @param rg logical: was the covariate inflormation  used? default: False
//' @return matrix with average of all the posterior probabilities: Hn, Ha and Hc
// [[Rcpp::export]]
arma::mat average_posterior_prob(arma::mat params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec, int nits, int thin, bool rg=false){
  double st=(nits/thin/2+1);
  double en=nits/thin;
  List posterior_prob_mat = posterior_prob(params, lbf_mat, nsnps, rg_vec, rg=rg);
  arma::mat avpost = posterior_prob_mat[st];
  for (int i = (st+1); i < en; i++){
    arma::mat post = posterior_prob_mat[i];
    avpost = avpost + post;
  }
  avpost = avpost/(en - st + 1);
  return avpost;
}

//' Average of priors: p0k, pak and pck
//'
//' @param params Vector of parameters: α, β and γ
//' @param nsnps number of snps
//' @param rg_vec Vector of the covariate
//' @param nits Number of iterations run in mcmc
//' @param thin thinning
//' @param rg logical: was the covariate inflormation  used? default: False
//' @return average pik matrix of priors: p0k, pak and pck
// [[Rcpp::export]]
arma::mat average_piks(arma::mat params, NumericVector nsnps, NumericVector rg_vec, int nits, int thin, bool rg=false){
  double st=(nits/thin/2+1);
  double en=nits/thin;
  List piks_list;
  piks_list = piks(params, nsnps, rg_vec, rg=rg);
  arma::mat avpik = piks_list[st];
  for (int i = (st+1); i < en; i++){
    arma::mat piks_mat = piks_list[i];
    avpik = avpik + piks_mat;
  }
  avpik=avpik/(en - st + 1);
  return avpik;
}
