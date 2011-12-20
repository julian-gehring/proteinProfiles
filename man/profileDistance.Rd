\name{profileDistance}

\alias{profileDistance}

\title{profileDistance}

\description{Compute mean euclidean distance between a protein group of
  interest and assess its significance by comparing it to the distances
  of randomly selected groups of proteins.}

\usage{

profileDistance(x, index, nSample=10000, plot=TRUE, main="")

}

\arguments{
  
  \item{x}{List containing the data to analyze, as returned from the
    \sQuote{readRatios} or \sQuote{filterRatios} function.}
  
  \item{index}{Vector of integers with the indices for the protein group
    of interest in the data set \sQuote{x}, as returned from
    \sQuote{grepAnnotation} or an equivalent function.}
  
  \item{nSample}{Integer specifying the number of randomly selected
    groups of proteins with the same number of elements as specified by
    \sQuote{index}, serving as a reference for the significance
    assessment (default: 10000).}
  
  \item{plot}{Logical whether to plot the distribution of distances
    (default: TRUE).}

  \item{main}{Character string, specifying the title for the plot
    (default: \"\").}
  
}

\value{
  List with elements:
  \describe{
    \item{$d0$}{Numeric with mean distance $d_{0}$ of profiles for the proteins of interest.}
    \item{$1$}{Vector of numerics with mean distances $d_{1}$ of randomly selected groups of proteins.}
    \item{$p$}{Numeric with p-value $p$.}
  }
}

\examples{

file <- system.file("extdata", "ips_sample.csv", package="proteinProfiles")
x <- readProteinData(file)
index <- grepAnnotation(x, "^28S")

z <- profileDistance(x, index)

}