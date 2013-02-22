## compute the p-value associated with the distance in n-dim space ##

## functions ##
readProteinData <- function(file, sep="\t", dataPattern="^Log2") {

  raw <- read.csv(file,
                  sep="\t", header=TRUE, stringsAsFactors=FALSE)

  ldx <- grepl(dataPattern, colnames(raw))
  data <- as.matrix(raw[ ,ldx])
  rownames(data) <- raw$id
  annotation <- data.frame(raw[ ,!ldx])

  return(list(data=data, annotation=annotation))
}


filterFeatures <- function(x, theta, plot=TRUE, verbose=FALSE) {

  nNa <- apply(x$data, 1, function(x) sum(is.na(x)))/ncol(x$data)

  if(plot) {
    plot(ecdf(nNa), xlim=c(0, 1),
         xlab="Fraction of missing data points per feature",
         ylab="Cumulative density",
         main="Diagnostic plot for filtering of missing data points")
    abline(v=theta)
  }

  ldx <- nNa <= theta  ## remove features with more than 'theta' missing
  x$data <- x$data[ldx, ]
  x$annotation <- x$annotation[ldx, ]
  x$theta <- theta

  if(verbose) {
    print(paste("Before filtering:", length(nNa)))
    print(paste("After filtering:", nrow(x$data)))
  }

  return(x)
}


grepAnnotation <- function(x, pattern, column="Protein.Names", ignore.case=FALSE) {

  index <- grep(pattern, x$annotation[[column]], ignore.case=ignore.case)

  return(index)
}


geneSetAnnotation <- function(x, geneSet, column="Gene.Names") {

  index <- which(toupper(x$annotation[[column]]) %in% toupper(geneIds(geneSet)))

  return(index)
}


geneSetCollectionAnnotation <- function(x, geneSetCollection, column="Gene.Names") {

  l <- lapply(geneSetCollection, geneSetAnnotation, x=x, column=column)

  return(l)
}


profileDistance <- function(x, index, nSample=10000,
                            plot=TRUE, main="") {

  d0 <- .distance(x$data[index, ])
  d1 <- replicate(nSample, .distance(x$data[sample.int(nrow(x$data), length(index)), ]))
  d1 <- d1[!is.nan(d1)]
  
  p <- sum(d1 <= d0)/length(d1)

  if(plot) {
    r <- range(d1)
    plot(ecdf(d1), do.points=FALSE,
         xlim=c(min(r[1], d0), max(r[2], d0)),
         xlab="Distance", ylab="Cumulative Density", main=main)
    abline(v=d0, col=2)
    axis(4, p, paste("p =", round(p, 2)))
    text(d0, 1, paste(" d0 =", round(d0, 2)), adj=c(0, 0))
  }

  return(list(d0=d0, d1=d1, p.value=p))
}


.distance <- function(x) mean(dist(x), na.rm=TRUE)





#gs.c3 <- getGmt("../genesets/c3.all.v3.0.symbols.gmt", collectionType=BroadCollection(category="c3"), geneIdType=SymbolIdentifier())

#gs.c2 <- getGmt("../genesets/c2.cp.v3.0.symbols.gmt", collectionType=BroadCollection(category="c2"), geneIdType=SymbolIdentifier())


#sapply(gs.c2, foo, gn)

#foo <- function(gs, geneNames) {
#  sum(geneNames %in% geneIds(gs))
#}
