// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]
arma::mat pars2pik_norg(arma::vec pars, NumericVector nsnps) {
  // double alpha=pars[0];
  // double beta=pars[1];

  int k = nsnps.length();
  arma::mat pik (1, 3);
  pik.row(0) = {1, exp(pars[0]), exp(pars[0] + pars[1])};

  arma::vec v = nsnps - 1;
  arma::mat m2 = arma::ones(3, k);
  m2.row(1) = arma::trans(v);
  arma::vec denom = diagvec(pik * m2);
  pik.each_col() /= denom ;
  return pik;
}

// [[Rcpp::export]]
arma::mat pars2pik_rg(arma::vec pars, NumericVector nsnps, NumericVector rg_vec) {
  // double alpha=pars[0];
  // double beta=pars[1];
  // double gamma=pars[2];

  int k = nsnps.length();
  arma::mat pik = arma::ones(k , 3);
  pik.col(1) =  arma::vec (rep( exp(pars[0]), k )) ;
  pik.col(2) = arma::vec (exp(pars[0] + pars[1] + (pars[2]*rg_vec)));

  arma::vec v = nsnps - 1;
  arma::mat m2 = arma::ones(3, k);
  m2.row(1) = arma::trans(v);
  arma::vec denom = diagvec(pik * m2);
  pik.each_col() /= denom ;
  return pik;
}

// [[Rcpp::export]]
arma::mat pars2pik(arma::vec pars, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
  arma::mat pik;
  if (rg==false){
    pik = pars2pik_norg(pars, nsnps);
  } else if ( rg_vec.length() == nsnps.length()) {
    pik = pars2pik_rg(pars, nsnps, rg_vec);
  } else{
    stop("Required rg vector when rg=TRUE.\n");
  }
  return pik;
}

// [[Rcpp::export]]
double logsumexp(arma::rowvec x) {
  return(log(sum(exp(x - x.max()))) + x.max());
}

// [[Rcpp::export]]
arma::mat logpost(arma::vec pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
  arma::mat pik = pars2pik(pars, nsnps, rg_vec=rg_vec, rg=rg);
  int k = nsnps.length();
  arma::mat logpost = arma::ones(k , 3);
  arma::mat logpik = log(pik);
  if (rg == false){
    logpost.col(0).fill(arma::as_scalar(logpik.col(0)));
    logpost.col(1) =  arma::as_scalar(logpik.col(1)) + lbf_mat.col(0);
    logpost.col(2) =  arma::as_scalar(logpik.col(2)) + lbf_mat.col(1);
  } else {
    logpost.col(0) =  logpik.col(0);
    logpost.col(1) =  logpik.col(1) + lbf_mat.col(0);
    logpost.col(2) =  logpik.col(2) + lbf_mat.col(1);
  }
  return logpost;
}

// [[Rcpp::export]]
double loglik(arma::vec pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
  arma::mat logpos = logpost(pars, lbf_mat, nsnps, rg_vec, rg=rg);
  int k = nsnps.length();
  arma::vec ll_k(k);
  for (int i = 0; i < k; i++){
    arma::rowvec y = logpos.row(i);
    ll_k(i) = logsumexp(y);
  }
  double loglik = accu(ll_k);
  return loglik;
}

// [[Rcpp::export]]
arma::mat posterior(arma::vec pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
  arma::mat logpos = logpost(pars, lbf_mat, nsnps, rg_vec, rg=rg);
  int k = nsnps.length();
  arma::vec denom(k);
  for (int i = 0; i < k; i++){
    arma::rowvec y = logpos.row(i);
    denom(i) = logsumexp(y);
  }
  logpos.each_col() -= denom ;
  arma::mat post = exp(logpos);
  return post;
}

// [[Rcpp::export]]
arma::vec sample_alpha(int n=1){
  return(rnorm(n, -10, 0.5));
}

// [[Rcpp::export]]
arma::vec sample_beta(int n=1){
  return(rgamma(n, 0.5, 2));
}

// [[Rcpp::export]]
arma::vec sample_gamma(int n=1){
  return(rgamma(n, 0.5, 2));
}

// [[Rcpp::export]]
double logd_alpha(double a, double mean=-10, double sd=0.5, bool log=true){
  return(R::dnorm(a, mean, sd, log));
}

// [[Rcpp::export]]
double logd_beta(double b){
  // scale set to 2, to correspond to the R function where we set 0.5 for the rate (scale = 1/0.5)
  // Also for logd_gamma
  return R::dgamma(b, 0.5, 2, true);
}

// [[Rcpp::export]]
double logd_gamma(double g){
  return R::dgamma(g, 0.5, 2, true);
}

// [[Rcpp::export]]
double logpriors(arma::vec pars, bool rg=false) {
  double loggamma;
  if (rg == false){
    loggamma = 0;
  } else{
    loggamma = logd_gamma(pars[2]);
  }
  double logprior = logd_alpha(pars[0]) + logd_beta(pars[1]) + loggamma;
  return logprior;
}

// [[Rcpp::export]]
double target(arma::vec pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
  double target = loglik(pars, lbf_mat, nsnps, rg_vec, rg) + logpriors(pars, rg);
  return target;
}

// [[Rcpp::export]]
arma::vec propose(arma::vec pars, double propsd=0.5){
  arma::vec propose = pars + arma::vec(rnorm(pars.n_elem, 0, propsd));
  return propose;
}

// [[Rcpp::export]]
arma::vec pars_init(bool rg=false){
  double alpha = arma::as_scalar(sample_alpha());
  double beta = arma::as_scalar(sample_beta());
  arma::vec pars;

  if(rg==true) {
    double gamma = arma::as_scalar(sample_gamma());
    pars = {alpha, beta, gamma};
  } else{
    pars = {alpha, beta};
  }
  return pars;
}

// [[Rcpp::export]]
arma::vec metrop_init(bool rg=false){
  double alpha = arma::as_scalar(sample_alpha());
  double beta = arma::as_scalar(sample_beta());
  arma::vec pars;

  if(rg==true) {
    double gamma = arma::as_scalar(sample_gamma());
    pars = {alpha, beta, gamma};
  } else{
    pars = {alpha, beta};
  }
  return pars;
}

// [[Rcpp::export]]
List metrop_run(arma::mat lbf_mat, NumericVector nsnps,
                NumericVector rg_vec=0, bool rg=false, int nits=10000,
                int thin=1){
  arma::vec pars = pars_init(rg=rg);
  arma::mat params = arma::zeros(pars.n_elem, nits/thin);
  arma::vec ll(nits/thin);
  arma::vec L0(nits);
  arma::vec T0(nits);
  L0[0] = loglik(pars, lbf_mat, nsnps, rg_vec,  rg=rg);
  T0[0] = L0[0] + logpriors(pars, rg=rg);
  for(int i = 0; i < nits; i++){
    if(i % thin == 0) {
      params.col(i) = pars;
      ll[i] = T0[i];
      }

    arma::vec newpars=propose(pars);
    double accept = exp(target(pars, lbf_mat, nsnps, rg_vec,  rg=rg) - T0[i]);
    if (R::runif(0, 1) < accept) {
      pars = newpars;
      L0[i+1] = loglik(pars, lbf_mat, nsnps, rg_vec,  rg=rg);
      T0[i+1] = L0[i+1] + logpriors(pars, rg=rg);
      }
  }
  return Rcpp::List::create(Rcpp::Named("ll") = ll,
                            Rcpp::Named("params") = params);
}

// [[Rcpp::export]]
List post_prob(arma::mat params, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false){
  List post;
  double k = params.n_cols;
  for (int i = 0; i < k; i++){
    post(i) = posterior(params.col(i),  lbf_mat, nsnps, rg_vec,  rg=rg);
  }
  return post;
}

// [[Rcpp::export]]
List piks(arma::mat params, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false){
  List piks;
  double k = params.n_cols;
  for (int i = 0; i < k; i++){
    piks(i) = pars2pik(params.col(i), nsnps, rg_vec,  rg=rg);
  }
  return piks;
}
