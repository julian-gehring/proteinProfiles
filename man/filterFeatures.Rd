\name{filterFeatures}

\alias{filterFeatures}

\title{filterFeatures}

\description{Filter out protein exceeding a certain fraction of missing
  data points.
}

\usage{

filterFeatures(x, theta, plot=TRUE, verbose=FALSE)

}

\arguments{
  
  \item{x}{List containing the data to analyze, as returned from the
    \sQuote{readRatios} function.}
  
  \item{theta}{Threshold $\\Theta$ specifying the fraction of data points
    in order to keep a protein.}
  
  \item{plot}{Logical whether to plot a diagnostic plot showing
    distribution of the fraction of missing data points associated with
    the proteins (default: TRUE).}

  \item{verbose}{Logical, whether to print the number of proteins before
    and after filtering (default: FALSE).}

}

\value{
  List with the same structure as \sQuote{x}, with elements not matching
  the filter criteria removed.
}

\examples{

example(ips_sample)

y <- filterFeatures(x, 0.2, verbose=TRUE)

}