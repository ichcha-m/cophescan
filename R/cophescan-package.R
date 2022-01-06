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
#' @importFrom ggplot2 ggplot aes scale_color_viridis_c geom_point labs scale_shape_manual aes_string element_blank element_line element_rect element_text ylab xlab theme ggtitle unit ylim xlim
#' @importFrom ggtern ggtern Tarrowlab Larrowlab Rarrowlab theme_arrowsmall theme_nomask theme_linedraw
#' @import ggrepel
#' @import pheatmap
#' @import matrixStats
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
#' @import magrittr
#' @references

utils::globalVariables(c(".","dfsane","dmvnorm","Hn","Ha","Hc","hit1","hit2","lbfak","lbfck","nsnps","snp","snp1","snp2","varbeta","z"))
NULL
