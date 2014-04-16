## compute the p-value associated with the distance in n-dim space ##

## functions ##

filterFeatures <- function(values, maxNAfraction, verbose=FALSE, plot=FALSE, ...) {

  fNA <- apply(values, 1, function(x) {sum(is.na(x))})/ncol(values)

  if(plot) {
    plot(ecdf(fNA), xlim=c(0, 1),
         xlab="Fraction of missing data points per feature",
         ylab="Cumulative density",
         main="Diagnostic plot for filtering of missing data points",
         ...)
    abline(v=maxNAfraction)
  }

  ## remove features with more than 'max_na_fraction' missing
  ldx <- fNA <= maxNAfraction
  values <- values[ldx, ,drop=FALSE]

  if(verbose) {
    print(paste("Before filtering:", length(fNA)))
    print(paste("After filtering:", nrow(values)))
  }

  if(nrow(values) == 0)
    warning("No values left.")

  return(values)
}


grepAnnotation <- function(anno, pattern, column, ...) {

  ## check for one match in colunm names
  if(sum(column %in% colnames(anno)) != 1)
    stop(sprintf("Column name '%s' must have exactly one match the annotation.", column))

  sub_anno = anno[[column]]
  ids = rownames(anno)
  index <- ids[grep(pattern, sub_anno, ...)]

  if(length(index) == 0)
    warning(sprintf("No matches for pattern '%s' in the annotation.", pattern))

  return(index)
}


profileDistance <- function(values, index, nSample=1000,
                            seed) {

  ## set a user defined seed
  if(!missing(seed))
    set.seed(seed)

  idx = which(rownames(values) %in% index)
  if(length(idx) == 0)
    stop("Indices do not match any protein identifier.")

  d0 <- .distance(values[idx, ])
  d1 <- replicate(nSample, .distance(values[sample.int(nrow(values), length(idx)), ]))
  d1 <- d1[!is.nan(d1)]
  if(length(d1) == 0)
    stop("Permutation only yields invalid distances.")
  
  p <- sum(d1 <= d0)/length(d1)

  return(list(d0=d0, d1=d1, p.value=p))
}


plotProfileDistance <- function(z, ...) {

  r <- range(z$d1)
  plot(ecdf(z$d1), do.points=FALSE,
       xlim=c(min(r[1], z$d0), max(r[2], z$d0)),
       xlab="Distance", ylab="Cumulative Density", ...)
  abline(v=z$d0, col=2)
  axis(4, z$p, paste("p =", round(z$p, 2)))
  text(z$d0, 1, paste(" d0 =", round(z$d0, 2)), adj=c(0, 0))
}


.distance <- function(values) mean(dist(values), na.rm=TRUE)
