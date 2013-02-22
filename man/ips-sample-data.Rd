\name{ips_sample-data}

\docType{data}

\alias{ips_sample}
\alias{ips_sample-data}

\alias{ratios}
\alias{annotation}

\title{IPS sample data}

\description{
  IPS sample data, separated in protein ratios and annotation.
}
\format{
  \describe{
    \item{ratios}{Matrix with protein profiles.  Rows correspond to
	  proteins, columns to samples.}
    \item{annotation}{Data frame with annotation columns, with rows
  matching those of \sQuote{ratios.}}
  }
}


\usage{data(ips_sample)}

\examples{
  data(ips_sample)
  x =list(data=ratios, annotation=annotation)
}


\keyword{datasets}

