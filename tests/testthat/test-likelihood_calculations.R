library(cophescan)
library(coloc)

## some fake data
# data(coloc_test_data)
# fm_causal=head(finemap.abf(coloc_test_data$D2),-1)
# fm_null=head(finemap.abf(coloc_test_data$N1),-1)
# fm_null=head(finemap.abf(coloc_test_data$D1),-1)

data("cophe_multi_trait_data")
fm_causal = head(finemap.abf(cophe_multi_trait_data$summ_stat$Trait_1),-1)
fm_null = head(finemap.abf(cophe_multi_trait_data$summ_stat$Trait_21),-1)
idx.causal=which.max(fm_causal$SNP.PP)
bf=matrix(c(rep(c(sum(fm_causal[-idx.causal,4][1:99]),fm_causal[idx.causal,4], 100),10), # Hc
            rep(c(sum(fm_causal[-500,4]), fm_causal[500,4], nrow(fm_causal)),10), # Ha
            rep(c(sum(fm_null[-500,4]), fm_null[500,4], nrow(fm_null)),80)), # Ha
          100,3,byrow=TRUE,
          dimnames=list(NULL,c("lBF.Ha","lBF.Hc","nsnps")))
lbf=as.data.frame(bf)
alpha=-10
beta=4
lbf_mat=as.matrix(lbf)[1,,drop=FALSE]
nsnps=nrow(fm_causal)

## SANITY CHECKS
## TRUE Hc
lbf_mat=as.matrix(lbf)[1,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
# TRUE Ha
lbf_mat=as.matrix(lbf)[11,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
# TRUE Hn
lbf_mat=as.matrix(lbf)[21,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)

## generate some key quantities in R
my_pars2pik=function(alpha,beta,nsnps) {
  pnk=1/( 1 + (nsnps-1) * exp(alpha) + exp(alpha+beta) )
  pak=exp(alpha) * pnk
  pck=exp(alpha+beta) * pnk
  cbind(pnk=pnk,pak=pak,pck=pck)
}
my_pik2ph=function(pik,nsnps) {
  ph=pik * matrix(c(1,nsnps-1,1),nrow(pik),ncol(pik),byrow=TRUE)
  ph/rowSums(ph)
}
mypik=my_pars2pik(alpha,beta,nsnps)
# piks=my_pars2pik( result$parameters["alpha",], result$parameters["beta",], nsnps)
# ph=piks * matrix(c(1,nsnps,1), nrow=nrow(piks), ncol=ncol(piks), byrow=TRUE)
# ph=ph/matrix(rowSums(ph),nrow(ph), ncol(ph))

## compare them to cpp
test_that("pars2pik", {
  result=pars2pik(c(alpha,beta),nsnps=nsnps,covar_vec=1,covar=FALSE)
  expect_equal(result %*% c(1,nsnps-1,1),matrix(1,1,1))
  expect_equal(as.vector(mypik %*% c(1,nsnps-1,1)),1)
  expect_equal(as.vector(mypik), as.vector(result))
})

test_that("logsumexp", {
  x=1:100
  expect_equal(log(sum(x)), logsumexp(log(x)))
})

test_that("loglik", {
  my_result=log(mypik[1,1] + exp(lbf_mat[1,1])*mypik[1,2] + exp(lbf_mat[1,2])*mypik[1,3])
  names(my_result)=NULL
  result=loglik(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(my_result, result)
})

test_that("logpost", {
  my_logpost=log(mypik) + c(0, lbf_mat[1,1], lbf_mat[1,2])
  dimnames(my_logpost)=NULL
  result=logpost(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(result, my_logpost)
})

test_that("get_posterior_prob", {
  my_pp=log(mypik) - log(rep(mypik[1],3)) + c(0, lbf_mat[1,1], lbf_mat[1,2])
  denom=logsumexp(my_pp)
  my_pp=my_pp - denom
  dimnames(my_pp)=NULL
  ## result_orig=get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  result=get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(exp(my_pp), result)
  expect_equal(rowSums(result),1)
})
## test_that("posterior_prob", {}) - just calls get_posterior_prob on a loop and returns a list of results

test_that("logpriors", {
  my_result= dnorm(alpha, -10, 0.5, log=TRUE) + dgamma(beta, 2, 0.5, log=TRUE)
  result=logpriors(c(alpha,beta), covar=FALSE)
  expect_equal(my_result, result)
})

test_that("target", {
  my_result=logpriors(c(alpha, beta), covar=FALSE) + loglik(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  result=target(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(my_result, result)
})
## test_that("piks", {})
## test_that("average_posterior_prob", {})
## test_that("average_piks", {})

