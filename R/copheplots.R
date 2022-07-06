#' Prepare data for plotting
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param causal.snpid query variant
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param hmp return for heatmap
#' @param cophe.plot return for cophe.plots
#' @param query_trait_names vector of names of the query traits, if the names of
#' the multi.dat list contain the trait names please pass query_trait_names=names(multi.dat)
#' default NULL
#'
#' @return plot list
prepare_plot_data <- function(multi.dat, causal.snpid, thresh_Ha=0.5, thresh_Hc=0.5, hmp=F, cophe.plot=T, query_trait_names=NULL){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat, query_trait_names)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  pp_dfcat1 <- pp_df[which(pp_df$Hc>thresh_Hc), ]
  pp_dfcat1 <- pp_dfcat1[order(pp_dfcat1$Hc, decreasing = T), ]
  pp_dfcat2 <- pp_df[which(pp_df$Ha>thresh_Ha), ]
  pp_dfcat2 <- pp_dfcat2[order(pp_dfcat2$Ha, decreasing = T), ]
  pp_dfcat <- rbind(pp_dfcat1, pp_dfcat2)

  if (cophe.plot){
    pp_df$x <- 1:nrow(pp_df)
    pp_df$L1 <-rownames(pp_df)
    pp_df$L1[!rownames(pp_df)%in%rownames(pp_dfcat1)] <- NA
    pp_df$L2 <-rownames(pp_df)
    pp_df$L2[!rownames(pp_df)%in%rownames(pp_dfcat2)] <- NA
    pp_df$ppHa <- as.factor(rep('ppHa', nrow(pp_df)))
    pp_df$ppHc <- as.factor(rep('ppHc', nrow(pp_df)))
    return(pp_df)
  }

  if (hmp){
    return(pp_dfcat)
  }
}

#' Extract beta and p-values of queried variant
#' @param traits.dat list of cmulti-trait oloc structured datasets
#' @param causal.snpid causal.snpid
#'
#' @return data.frame with one column indicating beta direction and another column with -log10(pval) of the queried variant
get_beta <- function(traits.dat, causal.snpid){
  pval_plot <- sapply(traits.dat, function(d) -(pnorm(-abs(trait.dat$beta)/sqrt(trait.dat$varbeta), log.p = TRUE) +
                                                  log(2))/log(10))
  names(pval_plot) <- names(traits.dat)

  betap <- sapply(traits.dat, function(d) d$beta[causal.snpid])
  beta_plot <- as.integer(betap > 0)
  beta_plot[beta_plot==0] <- 'n'
  beta_plot[beta_plot==1] <- 'p'
  names(beta_plot) <- names(traits.dat)
  return(data.frame(beta_plot=beta_plot, pval_plot=pval_plot))
}

#' cophe_plots showing the Ha and Hc of all traits and labelled above the specified threshold
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait (Takes row names as the phenotypes if you want an alternative naming pass a vector to pheno_names)
#' @param causal.snpid query variant
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param traits.dat list of multi-trait oloc structured datasets
#' @param pheno_names list of phenotype names
#' @param group_pheno Vector with additional grouping of phenotypes
#' @return cophescan plots of Ha and Hc
#' @export
#'
cophe_plot <- function(multi.dat, causal.snpid, thresh_Hc=0.5, thresh_Ha=0.5, traits.dat=NULL, pheno_names=NULL, group_pheno=NULL, beta_p=NULL){
  if (is.null(traits.dat)){
    print('Trait summary stat data required for pval PheWAS plot')
  }
  options(ggrepel.max.overlaps = Inf)
  plot_list <- list()
  pp_df <- prepare_plot_data(multi.dat,causal.snpid = causal.snpid, thresh_Hc=thresh_Hc, thresh_Ha=thresh_Ha, cophe.plot = T, hmp=F, query_trait_names=pheno_names)
  L1 <- pp_df$L1
  L2 <- pp_df$L2
  g1 <- suppressWarnings(ggplot(aes(x=x, y=Hc, label=L1), data=pp_df) +
                           geom_point(col='royalblue',alpha=0.7, size=5, aes(shape=ppHc)) +
                           scale_shape_manual('    ', values = c('ppHc'=18)) +
                           ylab("ppHc") +
                           xlab("Phenotypes") +
                           theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11),axis.ticks.x = element_blank(), axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.7, "lines"), max.iter  = 100000)+
                           ylim(0, 1))

  g2 <- suppressWarnings(ggplot(aes(x=x, y=Ha, label=L2), data=pp_df) +
                           geom_point(col='forestgreen',alpha=0.8, size=5, aes(shape=ppHa)) +
                           scale_shape_manual('    ', values = c('ppHa'=18))+
                           ylab("ppHa") + xlab("Phenotypes") +
                           theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11), axis.ticks.x = element_blank(), axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.8, "lines"), max.iter  = 100000)+   ylim(0, 1))
  plot_list[['ppHc']] <- g1
  plot_list[['ppHa']] <- g2

  if ((!is.null(traits.dat)) | (!is.null(beta_p))){
    if (is.null(beta_p)){
      beta_p <- get_beta(traits.dat, causal.snpid)
    }
    L1 <- rownames(pp_df)
    L1[beta_p$pval_plot<4] <- NA
      pval_plot <- beta_p$pval_plot
      beta_plot <- beta_p$beta_plot
    df_p <- data.frame(x=1:nrow(pp_df), y=pval_plot, label=L1, beta=as.factor(as.character(beta_plot)))
    g3 <- ggplot(aes(x=x, y=y, label=L1), data=df_p) +
      geom_point(col='maroon',alpha=0.7, aes(shape=beta), size=5) +
      scale_shape_manual('beta', values = c('n'="\u25BC", 'p'="\u25B2"))+
      ylab("-log10(pval)") + xlab(label="Phenotypes") +
      theme(axis.title=element_text(size=11, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size=11),axis.text.y = element_text(size=11), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size=3,box.padding = unit(0.8, "lines"), max.iter  = 100000) + ggtitle(causal.snpid)
    plot_list[['pval']] <- g3
  }
  return(plot_list)
}


#' Heatmap of multi-trait cophescan results
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param ... additional arguments to be passed to
#'   \link{pheatmap}
#' @return heatmap of posterior probabilities of the phentypes above the set threshold
#' @export
#'
cophe_heatmap <- function(multi.dat, thresh_Hc=0.5, thresh_Ha=0.5, ...){
  pp_dfcat <-  prepare_plot_data(multi.dat, causal.snpid = causal.snpid, thresh_Hc=thresh_Hc, thresh_Ha=thresh_Ha, cophe.plot = F, hmp=T)
  hmp <- pheatmap::pheatmap(pp_dfcat[, c('Hn', 'Ha', 'Hc') ], ...)
  return(hmp)
}

#' Plot region Manhattan for a trait highlighting the queried variant
#'
#' @param trait.dat dataset used as input for running cophescan
#' @param causal.snpid the id of the causal variant as present in trait.dat$snp
#'
#' @return region manhattan plot
#' @export
#'
plot_trait_manhat <- function(trait.dat, causal.snpid){
  x <- trait.dat$position
  cvidx <- which(trait.dat$snp%in%causal.snpid)
  y <- -(pnorm(-abs(trait.dat$beta)/sqrt(trait.dat$varbeta), log.p = TRUE) +
           log(2))/log(10)
  plot(x, y, xlab = "Position", ylab = "-log10(p)", pch = 16,
       col = "grey", sub=causal.snpid)
  points(x[cvidx], y[cvidx], col="red", pch=16)
}

#' Ternary plots of multi-trait cophescan output
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param traits.dat list of multi-trait coloc structured datasets
#' @param plot_pval Map p_value of Queried variant
#' @param thresh_Hc Hc threshold to be displayed
#' @param thresh_Ha Ha threshold to be displayed
#' @return ternary plot of the posterior probabilities of Ha, Hc and Hn
#' @export
#'
plot_cophe_ternary <- function(multi.dat, traits.dat=NULL, plot_pval=F, thresh_Hc=0.5, thresh_Ha=0.5){
  pp_df <- prepare_plot_data(multi.dat, thresh_Hc=thresh_Hc, thresh_Ha=thresh_Ha, cophe.plot = T, hmp=F)

  if (is.null(traits.dat) & plot_pval){
    stop("Please provide traits.dat to plot p_val")
  } else if (!is.null(traits.dat) & plot_pval) {
    beta_p <- get_beta(traits.dat, causal.snpid)
    betaPch <- ifelse(beta_p$beta_plot=="n", "\u25BC", "\u25B2")
    pp_df<- cbind(pp_df, beta_p)

    trn = ggtern(data=pp_df, aes_string(x="Hn",y="Hc",z="Ha"))
    trn = trn+geom_point(aes(color=pval_plot, shape=beta_plot),fill='black', size=4)+ ggtern::theme_linedraw() + scale_color_viridis_c(alpha = 0.5, n.breaks = 6) + theme_arrowsmall()+ labs(color="-log10(pval)",shape="beta", size=10) + Tarrowlab("")+ Larrowlab("")+ Rarrowlab("") +  theme_nomask()+  scale_shape_manual(values= c("\u25BC", "\u25B2"))
  } else {
    trn = ggtern(data=pp_df, aes_string(x="Hn",y="Hc",z="Ha"), aes(label=L1))
    trn = trn+geom_point(aes(color=Hc),fill='black',  size=4.5)+ theme_linedraw() + scale_color_viridis_c(alpha = 0.65, n.breaks = 6, option = "viridis") + theme_arrowsmall()+ labs(color="Hc", size=10)  + Tarrowlab("")+ Larrowlab("")+ Rarrowlab("")+  theme_nomask()+ theme(text = element_text(size=12)) #
    # + geom_label(mapping = aes(label=L1),position=position_nudge_tern(x=-0.2, y=0.3, z=0.2))
    #geom_text(position = pn, aes(label=L1),size=5)
  }

  return(trn)
}
