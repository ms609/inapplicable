PrepareData <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  at$names <- NULL
  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  data <- unlist(data, FALSE, FALSE)
  ret <- tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inapp.level <- which(at$levels == "-")
  if (!any(inapp.level)) warning ("    No inapplicable tokens detected in data.
    Use a hyphen ('-', not a dash) to denote taxa in which a transformation series is inapplicable.
    Does typing attributes(mydataname, 'contrast') return a column labelled '-'?")
  attr(ret, 'inapp.level') <- 2^(inapp.level - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  attr(ret, 'unique.tokens') <- apply(ret, 1, function(x) quick.min(x, inapp.level))
  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
  attr(ret, 'split.sizes') <- apply(ret, 1, function(x) vapply(applicable.tokens, function (y) sum(x==y), integer(1)))
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- '*phyDat'
  ret
}

quick.min <- function (x, inapp.level) {
  if (length(inapp.level)) return(sum(2^(c(0:(inapp.level-2), inapp.level:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}

AsBinary <- function(x) {
  # Adapted from code posted to R mailing list by Spencer Graves
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for (i in 1:ndigits){#i <- 1
		Base.b[, i] <- (x %% 2)
		x <- (x %/% 2)
	}
	if(N == 1) Base.b[1, ] else Base.b
}