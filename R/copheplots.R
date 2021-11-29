#' multi_plot
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above the threshold for Hc and Ha  will be labelled
#' @param heatmap Plot heatmap
#' @param cophe.plots Plot cophe.plots. Similar to PheWAS plots with Hc and Ha on the y axis
#' @import ggplot2
#' @return plot list
#' @export
multi_plot <- function(multi.dat, thresh=0.5, pltheatmap=F, cophe.plots=T, addplot=NULL, adddat=NULL, causal.snpid=NULL){
  if (!is.data.frame(multi.dat)){
    pp_df <- multitrait.simplify(multi.dat)
    if (is.null(multi.dat)){return(NULL)}
  } else {
    pp_df <- multi.dat
  }

  options(ggrepel.max.overlaps = Inf)

  pp_dfcat1 <- pp_df[which(pp_df$Hc>thresh), ]
  pp_dfcat1 <- pp_dfcat1[order(pp_dfcat1$Hc, decreasing = T), ]
  pp_dfcat2 <- pp_df[which(pp_df$Ha>thresh), ]
  pp_dfcat2 <- pp_dfcat2[order(pp_dfcat2$Ha, decreasing = T), ]
  pp_dfcat <- rbind(pp_dfcat1, pp_dfcat2)
  plot_list <- list()
  if (pltheatmap){
    hmp <- pheatmap::pheatmap(pp_dfcat[, c('Hn', 'Ha', 'Hc') ], cluster_cols = F, cluster_rows = F, cellwidth = 20, cellheight = 20)
    plot_list[['hmp']] <- hmp
  }
  pp_df$x <- 1:nrow(pp_df)
  if (cophe.plots){
    L1 <-rownames(pp_df)
    L1[!rownames(pp_df)%in%rownames(pp_dfcat1)] <- NA
    pp_df$ppHa <- as.factor(rep('ppHa', nrow(pp_df)))
    pp_df$ppHc <- as.factor(rep('ppHc', nrow(pp_df)))
    g1 <- suppressWarnings(ggplot2::ggplot(ggplot2::aes(x=x, y=Hc, label=L1), data=pp_df) +
      geom_point(col='royalblue', size=5, aes(shape=ppHc)) +
      scale_shape_manual('    ', values = c('ppHc'=18)) +
      ylab("Hc") +
      xlab("Phenotypes") +
      theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11),axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+
      ylim(0, 1))
    pp_df$x <- 1:nrow(pp_df)
    L2 <-rownames(pp_df)
    L2[!rownames(pp_df)%in%rownames(pp_dfcat2)] <- NA
    g2 <- suppressWarnings(ggplot2::ggplot(ggplot2::aes(x=x, y=Ha, label=L2), data=pp_df) +
      ggplot2::geom_point(col='forestgreen', size=5, aes(shape=ppHa)) +
        scale_shape_manual('    ', values = c('ppHa'=18))+
      ylab("Ha") + xlab("Phenotypes") +
      theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), legend.text = element_text(size=11), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+   ylim(0, 1))
    pp_df$x <- 1:nrow(pp_df)
    plot_list[['ppHc']] <- g1
    plot_list[['ppHa']] <- g2
  }
  if (!is.null(addplot)){
    if (is.null(adddat)){
      stop('Trait summary stat data required when addplot not NULL')
    } else {
      if ('pval' %in% addplot){
        pval_plot <- sapply(adddat, function(d) - log10(pnorm(-abs(d$beta[causal.snpid]/sqrt(d$varbeta[causal.snpid]))) * 2))
        names(pval_plot) <- row.names(pp_df)
      }
        betap <- sapply(adddat, function(d) d$beta[causal.snpid])
        beta_plot <- as.integer(betap > 0)
        beta_plot[beta_plot==0] <- 'n'
        beta_plot[beta_plot==1] <- 'p'
        names(beta_plot) <- row.names(pp_df)
    }
    df_p <- data.frame(x=1:nrow(pp_df), y=pval_plot, label=L1, beta=as.factor(as.character(beta_plot)))
    g4 <- ggplot2::ggplot(ggplot2::aes(x=x, y=y, label=L1), data=df_p) +
      geom_point(col='maroon', aes(shape=beta), size=5) +
      scale_shape_manual('beta', values = c('n'="\u25BC", 'p'="\u25B2"))+
      ylab("-log10(pval)") + xlab(label="Phenotypes") +
      theme(axis.title=element_text(size=15, face = 'bold'), legend.title = element_text(size=11), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size=11),axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000) + ggtitle(causal.snpid)
    plot_list[['pval']] <- g4
  }
  return(plot_list)
}


#' cophe_plot
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above this threshold for Hc and Ha will be labelled
#'
#' @return cophescan plots of Ha and Hc
#' @export
#'
cophe_plot <- function(multi.dat, thresh=0.5, addplot=NULL, adddat=NULL){
  if (!is.null(addplot) & (is.null(adddat))){
      stop('Trait summary stat data required when addplot required')
  }
  ret_plot <- multi_plot(multi.dat, thresh, cophe.plots = T, pltheatmap = F, addplot=addplot, adddat=adddat)
  return(ret_plot)
}


#' Title
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above this threshold for Hc and Ha will be displayed
#'
#' @return heatmap of posterior probabilities of the phentypes above the set threshold
#' @export
#'
cophe_heatmap <- function(multi.dat, thresh=0.5){
  ret_plot <- multi_plot(multi.dat, thresh, cophe.plots = F, pltheatmap = T)
  return(ret_plot)
}

