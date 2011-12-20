\name{proteinProfiles-package}

\docType{package}

\alias{proteinProfiles}
\alias{proteinProfiles-package}

\title{proteinProfiles package}

\description{
  The \pkg{GSRI} package estimates the number of differentially
  expressed genes in gene sets, utilizing the concept of the Gene Set
  Regulation Index (GSRI).
}

\details{


  The package is published under the GPL-3 license.
}

\author{
  Julian Gehring, Wolfgang Huber

  Maintainer: Julian Gehring <julian.gehring@embl.de>
}

\seealso{
  filterFeatures, grepAnnotation, profileDistance, readProteinData
}

\examples{
file <- system.file("extdata", "ips_sample.csv", package="proteinProfiles")
x <- readProteinData(file)
index <- grepAnnotation(x, "^28S")
z <- profileDistance(x, index)
}

\keyword{package}
\keyword{htest}
\keyword{distribution}