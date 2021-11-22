#' multi_plot
#'
#' @param multi.dat multi trait cophescan results returned from cophe.multitrait
#' @param thresh Traits with posterior probabilities above the threshold for Hc and Ha  will be labelled
#' @param heatmap Plot heatmap
#' @param cophe.plots Plot cophe.plots. Similar to PheWAS plots with Hc and Ha on the y axis
#' @import ggplot2
#' @return plot list
#' @export
multi_plot <- function(multi.dat, thresh=0.5, pltheatmap=F, cophe.plots=T, addplot=NULL, adddat=NULL){
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
    hmp <- pheatmap::pheatmap(pp_dfcat, cluster_cols = F, cluster_rows = F, cellwidth = 20, cellheight = 20)
    plot_list[['hmp']] <- hmp
  }
  pp_df$x <- 1:nrow(pp_df)
  if (cophe.plots){
    L1 <-rownames(pp_df)
    L1[!rownames(pp_df)%in%rownames(pp_dfcat1)] <- NA
    g1 <- suppressWarnings(ggplot2::ggplot(ggplot2::aes(x=x, y=Hc, label=L1), data=pp_df) +
      geom_point(col='royalblue', size=5, shape=18) +
      xlab(label = element_blank()) +
      ylab("Hc") +
      xlab("Phenotypes") +
      theme(axis.title=element_text(size=18, face = 'bold'), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+
      ylim(0, 1))
    pp_df$x <- 1:nrow(pp_df)
    L2 <-rownames(pp_df)
    L2[!rownames(pp_df)%in%rownames(pp_dfcat2)] <- NA
    g2 <- suppressWarnings(ggplot2::ggplot(ggplot2::aes(x=x, y=Ha, label=L2), data=pp_df) +
      ggplot2::geom_point(col='forestgreen', size=5, shape=18) +
      xlab(label = element_blank()) +
      ylab("Ha") + xlab("Phenotypes") +
      theme(axis.title=element_text(size=18, face = 'bold'), panel.background = element_rect(fill = 'white'), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=12), axis.line = element_line(color='grey22')) + ggrepel::geom_label_repel(size = 4,box.padding = unit(0.8, "lines"), max.iter  = 100000)+   ylim(0, 1))
    pp_df$x <- 1:nrow(pp_df)
    plot_list[['ppHc']] <- g1
    plot_list[['ppHa']] <- g2
  }
  if (!is.null(addplot)){
    if (is.null(adddat)){
      stop('Trait summary stat data required when addplot not NULL')
    } else {
      if ('pval' %in% addplot){
        pval_plot <- sapply(adddat, function(d) -log10(pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2))
        rownames(pval_plot) <- row.names(pp_df)
      }
      if ('beta' %in% addplot){
        beta_plot <- sapply(adddat, function(d) d$beta)
        rownames(beta_plot) <- row.names(pp_df)
      }
    }
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

