\name{filterFeatures}

\alias{filterFeatures}

\title{filterFeatures}

\description{Filter out proteins exceeding a certain fraction of missing
  data points.
}

\usage{

filterFeatures(values, maxNAfraction, verbose=FALSE, plot=FALSE, ...)

}

\arguments{
  
  \item{values}{Numeric matrix containing the data to analyze.}
  
  \item{maxNAfraction}{Numeric threshold specifying the maximum fraction
	of data points that can be missing (\sQuote{NA}) to still keep the
	protein.}

  \item{verbose}{Logical indicating whether to print the number of proteins before
    and after filtering (default: FALSE).}

  \item{plot}{Logical indicationg whether to plot a diagnostic plot showing the
    distribution of the fraction of missing data points associated with
    the proteins (default: FALSE).}

  \item{...}{Arguments passed to the '\sQuote{plot}' method.}

}

\value{
  Numeric matrix with the same structure as \sQuote{values},
  with elements not matching the filter criteria removed.
}

\examples{

data(ips_sample)

ratios_filtered <- filterFeatures(ratios, 0.2, verbose=TRUE)

}