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
#' @importFrom ggplot2 ggplot aes scale_color_viridis_c geom_point labs scale_shape_manual aes_string element_blank element_line element_rect element_text ylab xlab theme ggtitle unit ylim xlim
#' @import ggrepel
#' @import pheatmap
#' @import matrixStats
#' @importFrom graphics abline axis box par
#' @importFrom methods as is new slot
#' @importFrom stats as.dist as.formula coef coefficients complete.cases cor cutree glm integrate lm optimize pchisq pf prcomp qnorm sd var vcov hclust
#' @importFrom utils combn
#' @importFrom coloc runsusie
#' @importFrom viridis viridis
#' @importFrom graphics layout legend matplot mtext rect text title hist
#' @importFrom grDevices colorRampPalette palette rgb
#' @importFrom stats pnorm uniroot
#' @import magrittr
## usethis namespace: start
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table .SD
#' @importFrom data.table .BY
#' @importFrom data.table .N
#' @importFrom data.table .I
#' @importFrom data.table .GRP
#' @importFrom data.table .NGRP
#' @importFrom data.table .EACHI
## usethis namespace: end
NULL
utils::globalVariables(c(".","pval_plot","beta_plot","Hn","Ha","Hc","ppHn","ppHa","ppHc","hit1","hit2","lBF.Ha","lBF.Hc","nsnps","snp","idx1","i","varbeta","z","x","y","points","j","df","querysnpid","L1","L2", "PP.Hn", "PP.Ha", "PP.Hc", "grp"))
NULL
