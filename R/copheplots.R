#' Prepare data for plotting
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh_Ha Ha threshold to be displayed
#' @param thresh_Hc Hc threshold to be displayed
#' @param hmp return for heatmap
#' @param cophe.plot return for cophe.plots
#' @param causal.snpid query variant
#'
#' @return plot list
prepare_plot_data <- function(multi.dat, causal.snpid, thresh_Ha=0.5, thresh_Hc=0.5, hmp=F, cophe.plot=T){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
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
#' @return
get_beta <- function(traits.dat, causal.snpid){
  pval_plot <- sapply(traits.dat, function(d) - log10(pnorm(-abs(d$beta[causal.snpid]/sqrt(d$varbeta[causal.snpid]))) * 2))
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
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above this threshold for Hc and Ha will be labelled
#' @return cophescan plots of Ha and Hc
#' @export
#'
cophe_plot <- function(multi.dat, causal.snpid, thresh_Hc=0.5, thresh_Ha=0.5, traits.dat=NULL){
  if (is.null(traits.dat)){
      print('Trait summary stat data required for pval PheWAS plot')
  }
  options(ggrepel.max.overlaps = Inf)
  plot_list <- list()
  pp_df <- prepare_plot_data(multi.dat, thresh_Hc, thresh_Ha, cophe.plot = T, hmp=F)
  L1 <- pp_df$L1
  L2 <- pp_df$L2
  g1 <- suppressWarnings(ggplot(aes(x=x, y=Hc, label=L1), data=pp_df) +
                           geom_point(col='royalblue', size=5, aes(shape=ppHc)) +
                           scale_shape_manual('    ', values = c('ppHc'=18)) +
                           ylab("Hc") +
                           xlab("Phenotypes") +
                           theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11),axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+
                           ylim(0, 1))

  g2 <- suppressWarnings(ggplot(aes(x=x, y=Ha, label=L2), data=pp_df) +
                           geom_point(col='forestgreen', size=5, aes(shape=ppHa)) +
                           scale_shape_manual('    ', values = c('ppHa'=18))+
                           ylab("Ha") + xlab("Phenotypes") +
                           theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+   ylim(0, 1))
  plot_list[['ppHc']] <- g1
  plot_list[['ppHa']] <- g2

  if (!is.null(traits.dat)){
    beta_p <- get_beta(traits.dat, causal.snpid)
    pval_plot <- beta_p$pval_plot
    beta_plot <- beta_p$beta_plot

    df_p <- data.frame(x=1:nrow(pp_df), y=pval_plot, label=L1, beta=as.factor(as.character(beta_plot)))
    g3 <- ggplot(aes(x=x, y=y, label=L1), data=df_p) +
      geom_point(col='maroon', aes(shape=beta), size=5) +
      scale_shape_manual('beta', values = c('n'="\u25BC", 'p'="\u25B2"))+
      ylab("-log10(pval)") + xlab(label="Phenotypes") +
      theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size=11),axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000) + ggtitle(causal.snpid)
    plot_list[['pval']] <- g3
  }
  return(plot_list)
}


#' Heatmap of multi-trait cophescan results
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above this threshold for Hc and Ha will be displayed
#'
#' @return heatmap of posterior probabilities of the phentypes above the set threshold
#' @export
#'
cophe_heatmap <- function(multi.dat, thresh_Hc=0.5, thresh_Ha=0.5, cluster_columns = F, cluster_rows = F, cell_width = 20, cell_height = 20, ...){
  pp_dfcat <-  prepare_plot_data(multi.dat, thresh_Hc, thresh_Ha, cophe.plot = F, hmp=T)
  hmp <- pheatmap::pheatmap(pp_dfcat[, c('Hn', 'Ha', 'Hc') ], cluster_cols = cluster_columns, cluster_rows = cluster_rows, cellwidth = cell_width, cellheight = cell_height, ...)
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
    y <- -log10(pnorm(-abs(trait.dat$beta/sqrt(trait.dat$varbeta))) * 2)
    plot(x, y, xlab = "Position", ylab = "-log10(p)", pch = 16,
         col = "grey", main=causal.snpid)
    points(x[cvidx], y[cvidx], col="red", pch=16)
}

#' Ternary plots of multi-trait cophescan output
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#'
#' @return ternary plot of the posterior probabilities of Ha, Hc and Hn
#' @export
#'
plot_cophe_ternary <- function(multi.dat, traits.dat=NULL, plot_pval=F, thresh_Hc=0.5, thresh_Ha=0.5, ...){
  pp_df <- prepare_plot_data(multi.dat, thresh_Hc, thresh_Ha, cophe.plot = T, hmp=F)

    if (is.null(traits.dat) & plot_pval){
    stop("Please provide traits.dat to plot p_val")
  } else if (!is.null(traits.dat) & plot_pval) {
    beta_p <- get_beta(traits.dat, causal.snpid)
    betaPch <- ifelse(beta_p$beta_plot=="n", "\u25BC", "\u25B2")
    pp_df<- cbind(pp_df, beta_p)

    trn = ggtern(data=pp_df, aes_string(x="Hn",y="Hc",z="Ha"))
    trn = trn+geom_point(aes(color=pval_plot, shape=beta_plot),fill='black', size=6) + scale_color_viridis_c(alpha = 0.8, n.breaks = 6) + theme_arrowsmall()+ labs(color="-log10(pval)",shape="beta", size=10) + Tarrowlab("")+ Larrowlab("")+ Rarrowlab("") + theme_linedraw()+  theme_nomask()+  scale_shape_manual(values= c("\u25BC", "\u25B2"))
  } else {
    trn = ggtern(data=pp_df, aes_string(x="Hn",y="Hc",z="Ha"))
    trn = trn+geom_point(aes(color=Hc),fill='black',  size=6) + scale_color_viridis_c(alpha = 0.8, n.breaks = 6, option = "plasma") + theme_arrowsmall()+ labs(color="Hc", size=10)  + theme_linedraw()+ Tarrowlab("")+ Larrowlab("")+ Rarrowlab("")+  theme_nomask()
  }

  return(trn)
}
