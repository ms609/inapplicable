prepare.data <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  at$names <- NULL
  powers.of.2 = 2L^c(0L:(nLevel - 1L))
  tmp = cont %*% powers.of.2
  tmp = as.integer(tmp)
  data = unlist(data, FALSE, FALSE)
  ret = tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inapp.level <- which(at$levels == "-")
  attr(ret, 'inapp.level') <- 2^(inapp.level - 1)
  attr(ret, 'dim') <- c(nChar, nTip)
  
  attr(ret, 'min.steps') = apply(ret, 1, function (x) {
    vals = as.binary(unique(x))
    if (ncol(vals) >= inapp.level) vals[,inapp.level] = 0
    possibles <- matrix(FALSE, nrow(vals), ncol(vals))
    definites <- rep(FALSE, ncol(vals))
    for (i in 1:nrow(vals)) {
      V = vals[i,]
      if (any(V & definites)) next
      if (sum(V) == 1) {
        definites = definites | V
        next
      }
      possibles[i,] = V
    }
    while (any(as.logical(possibles))) {
      for (p in 1:i) if (any(possibles[p,] & definites)) possibles[p,] = FALSE
      sums = colSums(possibles)
      definites[which (sums == max(sums))[1]] = TRUE
    }
    return (max(0, sum(definites) - 1))
  })
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- '*phyDat'
  ret
}

as.binary <- function(x) {
  # Adapted from code by Spencer Graves, on R mailing list
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for (i in 1:ndigits){#i <- 1
		Base.b[, i] <- (x %% 2)
		x <- (x %/% 2)
	}
	if(N ==1) Base.b[1, ] else Base.b
}