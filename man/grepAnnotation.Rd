\name{annotation}

\alias{grepAnnotation}

\title{grepAnnotation}

\description{Find protein group of interest based on pattern
  matching against the annotation data.
}

\usage{

grepAnnotation(anno, pattern, column, ...)

}

\arguments{
  
  \item{anno}{Data frame with annotation to analyze.  The row names have to
	match  rownames in the matrix containing the protein abundance
	data.  Different columns correspond to different annotation columns,
	e.g. from different sources.}
  
  \item{pattern}{Character string representing a (sub)string or regular
    expression matched against the annotation.  For details, see
    \sQuote{grep}.}

  \item{column}{Column of the annotation data \sQuote{pattern} is matched
    against.}

  \item{...}{Optional arguments, passed to the \sQuote{grep} function.}

}

\value{
  Protein identifiers (rownames of \sQuote{anno}) specifying the protein
  group of interest in the data set.
}

\seealso{
  grep, filterRatios, proteinProfiles
}

\examples{

data(ips_sample)

index_28S <- grepAnnotation(annotation, pattern="^28S",
column="Protein.Name")

index_ribosome <- grepAnnotation(annotation, "Ribosome", "KEGG")

}