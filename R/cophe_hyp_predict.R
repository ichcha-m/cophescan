#' Predict cophescan hypothesis for tested associations
#'
#' @param cophe.res results obtained from `cophe.single`, `cophe.susie` or `cophe.multitrait` or data.frame with the following columns: PP.Hn, PP.Hc, PP.Ha, querysnp, querytrait
#' @param grouping.vars This is important for results from `cophe.susie` where there are multiple signals. These will be collapsed into one call. If you want to return all signals set this to a single variable eg: grouping.vars = c('querysnp')
#' @param Hc.cutoff threshold for PP.Hc above which the associations are called Hc
#' @param Hn.cutoff threshold for PP.Hn above which the associations are called Hn
#' @seealso \code{\link{cophe.single}}, \code{\link{cophe.susie}}, \code{\link{cophe.multitrait}}, , \code{\link{multitrait.simplify}}
#' @return returns dataframe with posterior probabilties of Hn, Hc and Ha with the predicted hypothesis based on the provided cut.offs.
#' @export
#' @importFrom dplyr group_by filter
cophe.hyp.predict <- function(cophe.res, grouping.vars = c('querysnp', 'querytrait'), Hc.cutoff= 0.6, Hn.cutoff= 0.2){
  message(paste0('Hc.cutoff = ', Hc.cutoff))
  message(paste0('Hn.cutoff = ', Hn.cutoff))
  if (!is.data.frame(cophe.res)){
    df = multitrait.simplify(cophe.res)
  } else {
    df = cophe.res
  }

  check_columns = c('PP.Hn', 'PP.Hc', 'PP.Ha', 'querysnp', 'querytrait', grouping.vars)
  if (any(!check_columns %in% (colnames(df)))){
    out = paste0(check_columns[which(!check_columns %in% (colnames(df)))], collapse=', ')
    stop( out, ' not found in the columns')
  }


  df$grp = apply( df[ , grouping.vars ] , 1 , paste , collapse = "_" )
  # df$cophe.hyp.call <- sapply(seq_along(df$grp), function(o) names(which.max(df[o,c('Hn', 'Ha', 'Hc')])))
  df$cophe.hyp.call <- 'Ha'
  df$cophe.hyp.call[df$PP.Hc>=Hc.cutoff]='Hc'
  df$cophe.hyp.call[df$cophe.hyp.call!='Hc'&df$PP.Hn>Hn.cutoff]='Hn'
  df$cophe.hyp.call[!df$cophe.hyp.call%in%c('Hn', 'Hc')]='Ha'

  dfc=df[df$cophe.hyp.call=='Hc',]
  if (nrow(dfc) > 0) {
    dfc=dfc %>% group_by(grp) %>% filter(PP.Hc == max(PP.Hc))
    dfc = dfc[!duplicated(dfc$grp), ]
  }
  dfa=df[df$cophe.hyp.call=='Ha',]
  if (length(which(dfa$grp%in%dfc$grp))>0)
    dfa = dfa[-which(dfa$grp%in%dfc$grp), ]
  if (nrow(dfa > 0)){
    dfa=dfa %>% group_by(grp) %>% filter(PP.Ha == min(PP.Ha))
    dfa = dfa[!duplicated(dfa$grp), ]
  }
  dfn=df[df$cophe.hyp.call=='Hn',]
  if (nrow(dfn > 0)){
    dfn=dfn %>% group_by(grp) %>% filter(PP.Hn == max(PP.Hn))
    dfn = dfn[!duplicated(dfn$grp), ]
  }
  if (length(which(dfn$grp%in%dfa$grp))>0){
    dfn = dfn[-which(dfn$grp%in%dfa$grp), ]
  }
  if (length(which(dfn$grp%in%dfc$grp))>0){
    dfn = dfn[-which(dfn$grp%in%dfc$grp), ]
  }
  # counts_df=data.frame(Hn= nrow(dfn), Ha= nrow(dfa), Hc= nrow(dfc) )
  df_out <- rbind(dfn, dfa, dfc)
  # return(list(counts_df=counts_df, hyp_cophe.hyp.call_df=df_out))
  return(as.data.frame(df_out))

}


#' Estimate the Hc.cutoff for the required FDR
#'
#' @param ppHc a vector containing the PP.Hc (the posterior probability of causal association) of all tests
#' @param fdr FDR default: 0.05
#' @param return_plot default: TRUE, plot the fdr estimated at the different Hc.cutoff
#'
#' @return the Hc.cutoff value for the specified FDR, if return_plot is True returns a plot showing the FDR calculated at different Hc thresholds
#' @export
#'
Hc.cutoff.fdr = function(ppHc, fdr = 0.05, return_plot=TRUE){
  res =  data.frame(Hc.cutoff = seq(0, 1, by=0.1))
  res$estimated_fdr = sapply(res$Hc.cutoff, function(x) mean(1-ppHc[ppHc>=x]))
  if (return_plot == TRUE) {
    yticks = seq(0, 1, by=fdr)
    plot(res$Hc.cutoff,  res$estimated_fdr, ylim = c(0,1), yaxt = "n", xlab = 'Hc.cutoff', ylab = 'estimated FDR')
    axis(side=2, at = yticks)
    abline(h=0.05, lty=2)
  }
  Hc.cutoff = min(res$Hc.cutoff[res$estimated_fdr < fdr], na.rm = TRUE)
  return(Hc.cutoff)
}
