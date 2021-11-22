// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]
arma::mat pars2pik_norg(NumericVector pars, NumericVector nsnps) {
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
arma::mat pars2pik_rg(NumericVector pars, NumericVector nsnps, NumericVector rg_vec) {
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
arma::mat pars2pik(NumericVector pars, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
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
arma::mat logpost(NumericVector pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
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
double loglik(NumericVector pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
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
arma::mat posterior(NumericVector pars, arma::mat lbf_mat, NumericVector nsnps, NumericVector rg_vec=0, bool rg=false) {
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

