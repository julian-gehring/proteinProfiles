\name{readProteinData}

\alias{readProteinData}

\title{readProteinData}

\description{Read protein profiles and annotation from a text file.}

\usage{
readProteinData(file, sep="\t", dataPattern="^Log2")
}

\arguments{
  \item{file}{Path to the =.csv= file containing data and annotation.}
  \item{sep}{Separator for columns in =.csv= file (default: \sQuote{\\t}
    for tab-delimited files).}
  \item{dataPattern}{Character string specifying the pattern to identify
    'data' columns containing measurement data in (contrast to
    annotation) (default: "^Log2").  'dataPattern' is matched against
    the header of the imported file.}
}

\value{
  List with elements:
  \describe{
    \item{data}{Matrix with $log_{2}$ ratios of all replicates (rows:
      proteins, columns: replicates).}
    \item{annotation}{Data frame with annotation corresponding to the
      proteins in \sQuote{data} (rows: proteins, columns: variables).}
  }
}


\examples{

file <- system.file("extdata", "ips_sample.csv", package="proteinProfiles")

x <- readProteinData(file)

}