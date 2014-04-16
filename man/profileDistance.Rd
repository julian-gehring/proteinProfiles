\name{profileDistance}

\alias{profileDistance}
\alias{plotProfileDistance}

\title{profileDistance}

\description{Compute distance between a protein group of
  interest and assess its significance by comparing it to the distances
  of randomly selected groups of proteins.}

\usage{

profileDistance(values, index, nSample=1000, seed)

plotProfileDistance(z, ...)

}

\arguments{
  
  \item{values}{Numeric matrix containing the protein data to analyze.}
  
  \item{index}{Protein identifiers for the protein group
    of interest in the data set \sQuote{values}, as returned from
    \sQuote{grepAnnotation} or an equivalent function.}
  
  \item{nSample}{Integer specifying the number of randomly selected
    groups of proteins with the same number of elements as specified by
    \sQuote{index}, serving as a reference for the significance
    assessment (default: 1000).}

  \item{seed}{Random seed, relevant for the permutation step. For
	details, see the \sQuote{seed} function.}

  \item{z}{Return object from \sQuote{profileDistance.}}

  \item{...}{Optional arguments passed to the \sQuote{plot} function.}
  
}

\value{
  List with elements:
  \describe{
    \item{d0}{Numeric with mean distance \eqn{$d_{0}$}{d_0} of profiles for the proteins of interest.}
    \item{d1}{Vector of numerics with mean distances \eqn{$d_{1}$}{d_1} of randomly selected groups of proteins.}
    \item{p}{Numeric with p-value \eqn{$p$}{p}.}
  }
}

\examples{

data(ips_sample)

index_28S <- grepAnnotation(annotation, pattern="^28S",
column="Protein.Name")

z <- profileDistance(ratios, index_28S)

z$p.value

plotProfileDistance(z)
}