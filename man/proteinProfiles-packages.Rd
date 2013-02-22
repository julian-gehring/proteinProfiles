\name{proteinProfiles-package}

\docType{package}

\alias{proteinProfiles}
\alias{proteinProfiles-package}

\title{proteinProfiles package}

\description{
  Significance assessment for distance measures of time-course protein profiles
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