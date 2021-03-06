% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copheplots.R
\name{prepare_plot_data}
\alias{prepare_plot_data}
\title{Prepare data for plotting}
\usage{
prepare_plot_data(
  multi.dat,
  causal.snpid,
  thresh_Ha = 0.5,
  thresh_Hc = 0.5,
  hmp = F,
  cophe.plot = T,
  query_trait_names = NULL
)
}
\arguments{
\item{multi.dat}{multi trait cophescan results returned from cophe.multitrait}

\item{causal.snpid}{query variant}

\item{thresh_Ha}{Ha threshold to be displayed}

\item{thresh_Hc}{Hc threshold to be displayed}

\item{hmp}{return for heatmap}

\item{cophe.plot}{return for cophe.plots}

\item{query_trait_names}{vector of names of the query traits, if the names of
the multi.dat list contain the trait names please pass query_trait_names=names(multi.dat)
default NULL}
}
\value{
plot list
}
\description{
Prepare data for plotting
}
