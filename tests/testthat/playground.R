library(testthat)
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

## library(cophescan)
## library(coloc)

devtools::load_all("~/RP/coloc")
devtools::load_all("~/RP/cophescan")
data(coloc_test_data)
fm_causal=head(finemap.abf(coloc_test_data$D2),-1)
fm_null=head(finemap.abf(coloc_test_data$N1),-1)
plot_dataset(coloc_test_data$D4)

idx.causal=which.max(fm_causal$SNP.PP)
bf=matrix(c(rep(c(sum(fm_causal[-idx.causal,4]),fm_causal[idx.causal,4], nrow(fm_causal)),10), # Hc
            rep(c(sum(fm_causal[-500,4]), fm_causal[500,4], nrow(fm_causal)),10), # Ha
            rep(c(sum(fm_null[-500,4]), fm_null[500,4], nrow(fm_null)),80)), # Ha
          100,3,byrow=TRUE,
          dimnames=list(NULL,c("lbfak","lbfck","nsnps")))

lbf=as.data.frame(bf)
alpha=-10
beta=4
lbf_mat=as.matrix(lbf)[1,,drop=FALSE]
nsnps=nrow(fm_causal)

## CHECKS
devtools::load_all()
## TRUE Hc
f=function(lbfc) {
  lbf_mat=as.matrix(lbf)[4,,drop=FALSE]
  lbf_mat[1,2]=lbfc
  get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
}
lbfc=1:10000
res=sapply(lbfc,f)
plot(lbfc, res[3,])
lines(lbfc,res[2,],col="blue")
lines(lbfc,res[1,],col="grey")


  lbf_mat=as.matrix(lbf)[1,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
# TRUE Ha
lbf_mat=as.matrix(lbf)[11,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
# TRUE Hn
lbf_mat=as.matrix(lbf)[21,,drop=FALSE]
get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)


## look at some real data
library(magrittr)
proc_file=function(f) {
  (load(f))
  wh=which(rownames(bf_df)==bf_df$querysnp)
  lbf_mat=cbind(logsumexp(bf_df$lABF[-wh]), bf_df$lABF[wh], nrow(bf_df))
}
files=list.files("~/Downloads/bayesfac_test/ukbb_ferkingstad",full=TRUE)
ukbb_ferk=lapply(files, proc_file) %>% do.call("rbind",.) %>% as.data.frame()
files=list.files("~/Downloads/bayesfac_test/ukbb_finn",full=TRUE)
ukbb_finn=lapply(files, proc_file) %>% do.call("rbind",.) %>% as.data.frame()
colnames(ukbb_ferk)=colnames(ukbb_finn)=c("lbfak","lbfck","nsnps")
ukbb=rbind(ukbb_finn,ukbb_ferk)

## TRUE Hc
alpha=-10.5
beta=4
f=function(lbfc) {
  lbf_mat=as.matrix(ukbb)[8,,drop=FALSE]
  lbf_mat[1,2]=lbfc
  get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
}
lbfc=1:100
res=sapply(lbfc,f)
plot(lbfc, res[3,],type="l")
lines(lbfc,res[2,],col="blue")
lines(lbfc,res[1,],col="grey")


result_ukbb=run_metrop_priors(ukbb)
colMeans(result_ukbb$avg.posterior)
result_ukbb$avg.posterior

hist(bf_df$lABF)

## post=matrix(NA, 100, 3)
## j=1
## for(i in sample(2000:10000)) {
##   post=get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
##   j=j+1
## }

## pp=posterior_prob(matrix(c(alpha,beta),nrow=1), lbf_mat, nsnps, 1, FALSE)

result=run_metrop_priors(lbf)
colMeans(result$avg.posterior)



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


## result1=run_metrop_priors(lbf,beta_scale=1,beta_shape=1)
## colMeans(result1$avg.posterior)
## result30=run_metrop_priors(lbf,beta_scale=30,beta_shape=1/5)
## colMeans(result30$avg.posterior)

## resultn=run_metrop_priors(lbf,beta_shape=0,beta_scale=2)
## colMeans(resultn$avg.posterior)
## piks=my_pars2pik( resultn$parameters["alpha",], resultn$parameters["beta",], nsnps)

## resultn2=run_metrop_priors(lbf,beta_shape=5,beta_scale=3)
## colMeans(resultn2$avg.posterior)
piks=my_pars2pik( result$parameters["alpha",], result$parameters["beta",], nsnps)
ph=piks * matrix(c(1,nsnps,1), nrow=nrow(piks), ncol=ncol(piks), byrow=TRUE)
ph=ph/matrix(rowSums(ph),nrow(ph), ncol(ph))

par(mfrow=c(3,3))
## hist(piks[,"pnk"])
## hist(piks[,"pak"])
## hist(piks[,"pck"])
plot(1:nrow(ph), ph[,"pnk"],type="l"); abline(h=mean(ph[,"pnk"]),col="red")
plot(1:nrow(ph), ph[,"pak"],type="l"); abline(h=mean(ph[,"pak"]),col="red")
plot(1:nrow(ph), ph[,"pck"],type="l"); abline(h=mean(ph[,"pck"]),col="red")
plot(1:nrow(piks), piks[,"pnk"],type="l"); abline(h=mean(piks[,"pnk"]),col="red")
plot(1:nrow(piks), piks[,"pak"],type="l"); abline(h=mean(piks[,"pak"]),col="red")
plot(1:nrow(piks), piks[,"pck"],type="l"); abline(h=mean(piks[,"pck"]),col="red")
plot(1:ncol(result$parameters), result$parameters["alpha",],type="l");
plot(1:ncol(result$parameters), result$parameters["beta",],type="l");

## result=run_metrop_priors(lbf)
## colMeans(result_orig$avg.posterior)

## result2=run_metrop_priors(lbf,alpha_sd=sqrt(0.5), beta_scale=0.5)
## colMeans(result2$avg.posterior)

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


test_that("pars2pik", {
  result=pars2pik(c(alpha,beta),nsnps=nsnps,rg_vec=1,rg=FALSE)
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

## test_that("logpost", { # TODO
##   dnorm(alpha, -10, 0.5, log=TRUE) + dbeta(beta, 2, 2, log=TRUE) + my_loglik(lbf)
##   result=logpost(c(alpha,beta), as.matrix(lbf)[1,,drop=FALSE], nsnps, 1, FALSE)
## })

test_that("logpost, get_posterior_prob", {
  logpost=log(mypik) - log(rep(mypik[1],3)) + c(0, lbf_mat[1,1], lbf_mat[1,2])
  denom=logsumexp(logpost)
  logpost=logpost - denom
  dimnames(logpost)=NULL
  ## result_orig=get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  result=get_posterior_prob(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(exp(logpost), result)
  expect_equal(rowSums(result),1)
})
## test_that("posterior_prob", {}) - just calls get_posterior_prob on a loop

test_that("logpriors", {
  my_result= dnorm(alpha, -10, 0.5, log=TRUE) + dgamma(beta, 2, 0.5, log=TRUE)
  result=logpriors(c(alpha,beta), rg=FALSE)
  expect_equal(my_result, result)
})

test_that("target", {
  my_result=logpriors(c(alpha, beta), rg=FALSE) + loglik(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  result=target(c(alpha,beta), lbf_mat, nsnps, 1, FALSE)
  expect_equal(my_result, result)
})
## test_that("piks", {})
## test_that("average_posterior_prob", {})
## test_that("average_piks", {})
