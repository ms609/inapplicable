prepare.data <- function (data) {
## Prepares a phyDat object for analysis with modified Fitch algorithm
# ARGUMENTS:
#   «data», a phydat object
# RETURN:
#   A list, with each list item representing a tip and housing a matrix.
#     Each column of the matrix corresponds to a transformation series 'pattern'
#     Each row of the matrix corresponds to a token
#     Each cell contains a logical denoting whether the character could take the value indicated by the row
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