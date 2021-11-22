#' The 'cophescan' package.
#'
#' @description Coloc adapted Phenome-wide Scans
#'
#' @docType package
#' @name cophescan-package
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @aliases cophescan
#' @useDynLib cophescan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import RcppArmadillo
#' @import ggplot2
#' @importFrom ggplot2 aes
#' @import ggrepel
#' @import pheatmap
#' @importFrom graphics abline axis box par
#' @importFrom methods as is new slot
#' @importFrom stats as.dist as.formula coef coefficients complete.cases cor cutree glm integrate lm optimize pchisq pf prcomp qnorm sd var vcov hclust
#' @importFrom utils combn
#' @importFrom coloc runsusie
#' @import data.table
#' @importFrom viridis viridis
#' @importFrom graphics layout legend matplot mtext rect text title hist
#' @importFrom grDevices colorRampPalette palette rgb
#' @importFrom stats pnorm uniroot
#' @references

#'utils::globalVariables(c(".","dfsane","dmvnorm","Hn","Ha","Hc","hit1","hit2","lABF.df2","lABF.h3","lbf1","lbf2","nsnps","snp","snp1","snp2","varbeta","z"))
NULL
