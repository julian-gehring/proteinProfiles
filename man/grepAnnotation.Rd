\name{annotation}

\alias{grepAnnotation}
\alias{geneSetAnnotation}
\alias{geneSetCollectionAnnotation}

\title{grepAnnotation}

\description{Find indices for protein group of interest based on pattern
  matching against the annotation data, or matching of gene sets.
}

\usage{

grepAnnotation(x, pattern, column="Protein.Names", ignore.case=FALSE)

geneSetAnnotation(x, geneSet, column="Gene.Names")

geneSetCollectionAnnotation(x, geneSetCollection, column="Gene.Names")

}

\arguments{
  
  \item{x}{List containing the data and annotation to analyze, as
    returned from the \sQuote{readProteinData} or \sQuote{filterRatios}
    function.}
  
  \item{pattern}{Character string representing a (sub)string or regular
    expression matched against the annotation.  For details, see
    \sQuote{grep}.}

  \item{geneSet}{Gene set of class \sQuote{GeneSet} from the
    \sQuote{GSEABase} package.}

  \item{geneSetCollection}{Collection of gene set of class
    \sQuote{GeneSetCollection} from the \sQuote{GSEABase} package.}
  
  \item{column}{Column of the annotation data \sQuote{pattern} is matched
    against (default: "Protein.Names" matching against the protein
    names for \sQuote{grepAnnotation}); or "Gene.Names" matching against
    the gene names for \sQuote{geneSetAnnotation} and
    \sQuote{geneSetCollectionAnnotation}.}

  \item{ignore.case}{Logical whether case sensitivity should be ignored
    in string matching (default: FALSE).  For details, see
    sQuote{ignore.case} for \sQuote{grep}.}

}

\value{
  Vector of integers specifying the indices for the protein group of
  interest in the data set.
}

\seealso{
  grep, readProteinData, filterRatios, GSEABase
}

\examples{

file <- system.file("extdata", "ips_sample.csv", package="proteinProfiles")
x <- readProteinData(file)

index1 <- grepAnnotation(x, "^28S")


library(GSEABase)

c3sub <- system.file("extdata", "c3.all.v3.0.sym.4.gmt", package="proteinProfiles")
gsc <- getGmt(c3sub)
gs <- gsc[[1]]

index2 <- geneSetAnnotation(x, gs)

index3 <- geneSetCollectionAnnotation(x, gsc)

}