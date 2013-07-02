prepare.data <- function (data) {
  at <- attributes(data)
  cont <- at$contrast
  nLevel <- length(at$level)
  nTip <- length(data)
  nChar <- attr(data, 'nr')
  ret <- mclapply(data, function (x) {colSums(2^(0:(nLevel-1))*t(cont[x,]))})
  ret <- mclapply(seq_along(data), function(x) { # for each node
    if (x <= nTip) matrix(as.logical(intToBits(ret[[x]])), nrow=32)[-((nLevel+1):32),]
  })
  attributes(ret) <- at
  names(ret) <- names(data)
  class(ret) <- '*phyDat'
  ret
}