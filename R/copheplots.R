#' Prepare data for plotting
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait or multitrait.simplify
#' @param querysnpid query variant
#' @param query_trait_names vector of names of the query traits
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param hmp return for heatmap
#' @param cophe.plot default: TRUE, return for `cophe_plot`
#' @seealso \code{\link{cophe_plot}}, \code{\link{cophe.susie}}, \code{\link{cophe.multitrait}},  \code{\link{multitrait.simplify}}
#' default NULL
#'
#' @return plot list
prepare_plot_data <- function(multi.dat, querysnpid, query_trait_names, thresh_Ha=0.5, thresh_Hc=0.5, hmp=FALSE, cophe.plot=TRUE){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat, query_trait_names)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  pp_dfcat1 <- pp_df[which(pp_df$PP.Hc>thresh_Hc), ]
  pp_dfcat1 <- pp_dfcat1[order(pp_dfcat1$PP.Hc, decreasing = TRUE), ]
  pp_dfcat2 <- pp_df[which(pp_df$PP.Ha>thresh_Ha), ]
  pp_dfcat2 <- pp_dfcat2[order(pp_dfcat2$PP.Ha, decreasing = TRUE), ]
  pp_dfcat <- rbind(pp_dfcat1, pp_dfcat2)

  if (cophe.plot){
    pp_df$x <- 1:nrow(pp_df)
    pp_df$L1 <- pp_df$querytrait
    pp_df$L1[!pp_df$querytrait%in%pp_dfcat1$querytrait] <- NA
    pp_df$L2 <-pp_df$querytrait
    pp_df$L2[!pp_df$querytrait%in%pp_dfcat2$querytrait] <- NA
    pp_df$ppHa <- as.factor(rep('ppHa', nrow(pp_df)))
    pp_df$ppHc <- as.factor(rep('ppHc', nrow(pp_df)))
    return(pp_df)
  }

  if (hmp){
    return(pp_dfcat)
  }
}

#' Extract beta and p-values of queried variant
#' @param traits.dat list of coloc structured dataset
#' @param querysnpid vector of querysnpid
#' @param querytrait vector of querytrait names
#'
#' @return data.frame with one column named beta_plot: indicating beta direction (n/p) and another column named pval_plot with -log10(pval) of the queried variant
get_beta <- function(traits.dat, querysnpid, querytrait){
  pval_plot <- sapply(traits.dat, function(d) -(pnorm(-abs(d$beta[querysnpid])/sqrt(d$varbeta[querysnpid]), log.p = TRUE) +
                                                  log(2))/log(10))
  names(pval_plot) <- names(traits.dat)

  betap <- sapply(traits.dat, function(d) d$beta[querysnpid])
  beta_plot <- as.integer(betap > 0)
  beta_plot[beta_plot==0] <- 'n'
  beta_plot[beta_plot==1] <- 'p'
  # names(beta_plot) <- names(traits.dat)
  ret <- data.frame(beta_plot=beta_plot, pval_plot=pval_plot, querysnp=querysnpid, querytrait=querytrait)
  return(ret)
}

#' cophe_plots showing the Ha and Hc of all traits and labelled above the specified threshold
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait or multitrait.simplify
#' @param querysnpid query variant (only a single variant for PheWAS plots)
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param query_trait_names list of phenotype names
#' @param beta_p data.frame (from the `get.beta` function) with four columns : 1. "beta_plot": indicating beta direction (p or n) 2. "beta_plot": -log10(pval) of the queried variant 3. "querysnp" 4. "querytrait".
#' @param traits.dat list of multi-trait coloc structured datasets
#' @param group_pheno Vector with additional grouping of phenotypes
#' @seealso \code{\link{cophe.single}}, \code{\link{cophe.susie}}, \code{\link{cophe.multitrait}}, , \code{\link{multitrait.simplify}}
#' @return cophescan plots of Ha and Hc
#' @export
#'
cophe_plot <- function(multi.dat, querysnpid, query_trait_names, thresh_Hc=0.5, thresh_Ha=0.5,  beta_p=NULL, traits.dat=NULL, group_pheno=NULL){
  if (length(querysnpid)!=1){
    stop("length(querysnpid)!=1, Pass traits corresponding only to a single queryvariant")
  }
  if (is.null(traits.dat) & is.null(beta_p)){
    warning('Trait summary stat data or beta_p data.frame required for pval PheWAS plot')
  }
  if (!is.null(beta_p) & !any(colnames(beta_p)%in%c("beta_plot", "pval_plot", "querysnpid", "querytrait"))){
    stop("Check the colnames of the beta_p data.frame or pass NULL if pval PheWAS plot not required")
  }

  plot_list <- list()
  pp_df <- prepare_plot_data(multi.dat,querysnpid = querysnpid, query_trait_names=query_trait_names, thresh_Hc=thresh_Hc, thresh_Ha=thresh_Ha, cophe.plot = TRUE, hmp=FALSE)
  L1 <- pp_df$L1
  L2 <- pp_df$L2
  g1 <- suppressWarnings(ggplot(aes(x=x, y=PP.Hc, label=L1), data=pp_df) +
                           geom_point(col='royalblue',alpha=0.7, size=5, aes(shape=ppHc)) +
                           scale_shape_manual('    ', values = c('ppHc'=18)) +
                           ylab("ppHc") +
                           xlab("Phenotypes") +
                           theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11),axis.ticks.x = element_blank(), axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.7, "lines"), max.iter  = 100000, max.overlaps = Inf)+
                           ylim(0, 1))

  g2 <- suppressWarnings(ggplot(aes(x=x, y=PP.Ha, label=L2), data=pp_df) +
                           geom_point(col='forestgreen',alpha=0.8, size=5, aes(shape=ppHa)) +
                           scale_shape_manual('    ', values = c('ppHa'=18))+
                           ylab("ppHa") + xlab("Phenotypes") +
                           theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11), axis.ticks.x = element_blank(), axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.8, "lines"), max.iter  = 100000, max.overlaps = Inf)+   ylim(0, 1))
  plot_list[['ppHc']] <- g1
  plot_list[['ppHa']] <- g2

  if ((!is.null(traits.dat)) | (!is.null(beta_p))){
    if (is.null(beta_p)){
      beta_p <- get_beta(traits.dat, querysnpid, query_trait_names)
    }
    beta_p <- beta_p[order(match(beta_p$querytrait, pp_df$querytrait)),]

    L1 <- beta_p$querytrait
    L1[beta_p$pval_plot<4] <- NA
    pval_plot <- beta_p$pval_plot
    beta_plot <- beta_p$beta_plot
    df_p <- data.frame(x=seq_len(nrow(beta_p)), y=pval_plot, label=L1, beta=as.factor(as.character(beta_plot)))
    g3 <- ggplot(aes(x=x, y=y, label=L1), data=df_p) +
      geom_point(col='maroon',alpha=0.7, aes(shape=beta), size=5) +
      scale_shape_manual('beta', values = c('n'="\u25BC", 'p'="\u25B2"))+
      ylab("-log10(pval)") + xlab(label="Phenotypes") +
      theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size=11),axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.8, "lines"), max.iter  = 100000, max.overlaps = Inf) + ggtitle(querysnpid)
    plot_list[['pval']] <- g3
  }
  return(plot_list)
}


#' Heatmap of multi-trait cophescan results
#'
#' @param multi.dat multi trait cophescan results returned from `cophe.multitrait` or formatted in the same way with multitrait.simplify
#' @param querysnpid query variant
#' @param query_trait_names names of phenotypes corresponding to the multi.dat results
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param ... additional arguments to be passed to
#'   \link[pheatmap]{pheatmap}
#' @return heatmap of posterior probabilities of the phentypes above the set threshold
#' @export
#'
cophe_heatmap <- function(multi.dat, querysnpid, query_trait_names, thresh_Hc=0.5, thresh_Ha=0.5, ...){
  pp_dfcat <-  prepare_plot_data(multi.dat, querysnpid = querysnpid, query_trait_names = query_trait_names, thresh_Hc=thresh_Hc, thresh_Ha=thresh_Ha, cophe.plot = FALSE, hmp=TRUE)
  hmp <- pheatmap::pheatmap(pp_dfcat[, c('Hn', 'Ha', 'Hc') ], ...)
  return(hmp)
}

#' Plot region Manhattan for a trait highlighting the queried variant
#'
#' @param trait.dat dataset used as input for running cophescan
#' @param querysnpid the id of the causal variant as present in trait.dat$snp, plotted in red
#' @param alt.snpid the id of the other variants as a vector to be plotted, plotted in blue
#'
#' @return regional manhattan plot
#' @export
#'
plot_trait_manhat <- function(trait.dat, querysnpid, alt.snpid=NULL){
  if (is.null(trait.dat$position)){
    stop("no snp position element given")
  } else {
    x <- trait.dat$position
  }
  queryidx <- which(trait.dat$snp%in%querysnpid)
  y <- -(pnorm(-abs(trait.dat$beta)/sqrt(trait.dat$varbeta), log.p = TRUE) +
           log(2))/log(10)
  plot(x, y, xlab = "Position", ylab = "-log10(p)", pch = 16,
       col = "grey", sub=querysnpid)
  points(x[queryidx], y[queryidx], col="red", pch=16, cex=1.2)
  if (!is.null(alt.snpid)){
    altidx <- which(trait.dat$snp%in%alt.snpid)
    points(x[altidx], y[altidx], col="blue", pch=16, cex=1.2)
  }
}

