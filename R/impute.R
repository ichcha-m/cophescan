#' cophe.impute
#' Imputation of untyped SNPs with summary statistics using the SSIMP algorithm
#' @param dataset dataset conforming to coloc, Use coloc::check_dataset() for
#' @param r2.thr tags are selected with max r2 <= r2.thr
#' @param LD LD matrix including all SNPs to be imputed
#' @param MAF Minimum allele frequency
#' @param nsample number of samples from which the LD matrix is derived
#' @return Imputed dataset
#' @export
#'
cophe.impute=function(dataset, LD, MAF, nsample, r2.thr=1){
  x <- df
  lambda=2/sqrt(nsample)
  LDl=(1-lambda) * LD + lambda * diag(nrow(LD))

  drop=which(x$varbeta==0 | is.infinite(x$beta) | is.na(x$beta))
  if(length(drop)) {
    x$snp=x$snp[-drop]
    x$beta=x$beta[-drop]
    x$varbeta=x$varbeta[-drop]
  }
  # str(x)

  ## find snps to impute - do not exceed boundary of snps typed (eqtl windows)
  # + 5kb window either side
  bp_snp=sub(".*-","",x$snp) %>% as.numeric()
  bp_LD=sub(".*-","",colnames(LD)) %>% as.numeric()
  LD_snps_window=colnames(LD)[bp_LD >= (min(bp_snp)-5e+4) & bp_LD <= (max(bp_snp)+5e+4)]
  snps_toimp=setdiff(LD_snps_window,x$snp)
  message("snps to impute: ",length(snps_toimp)," out of ",ncol(LD))
  if(length(snps_toimp)==0)
    return(x)

  ## find tags )
  LD2=LD[tags,tags]^2
  diag(LD2)=0
  drop=which(LD2>r2.thr, arr.ind=TRUE)
  tagdrop=rep(FALSE,length(tags))
  while(nrow(drop)) {
    idrop=sample(drop[,"row"],1)
    tagdrop[idrop]=TRUE
    drop=drop[ drop[,"row"]!=idrop & drop[,"col"]!=idrop , ]
  }
  tags=tags[!tagdrop]
  message("tags selected with max r2 <= ",r2.thr," : ",length(tags))

  ## ssimp algorithm
  cl=LDl[snps_toimp,tags]
  x$Z=(x$beta/sqrt(x$varbeta))
  Z_tags=x$Z[ match(tags,x$snp) ]
  ## A=solve(LDl[tags,tags], Z_tags)
  ## Z_imp=cl %*% matrix(A,ncol=1)
  Cl=solve(LDl[tags,tags])
  Z_imp=cl %*% Cl %*% matrix(Z_tags,ncol=1)
  names(Z_imp) <- snps_toimp
  se_imp=1/sqrt(2*MAF[names(Z_imp)]*(1-MAF[names(Z_imp)])*x$N)
  beta_imp=Z_imp*se_imp

  ## prediction quality
  r2pred=diag(cl %*% Cl %*% t(cl))
  ## append to dataset
  x$r2pred=c(rep(NA,length(x$snp)),r2pred)
  x$snp=c(x$snp,rownames(Z_imp))
  x$varbeta=c(x$varbeta, se_imp^2)
  x$beta=c(x$beta, beta_imp)
  x$Z=c(x$Z,Z_imp)
  x$imputed=x$snp %in% snps_toimp
  x$LD = LD
  ## save
  return(x)
}
