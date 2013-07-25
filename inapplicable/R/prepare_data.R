prepare.data <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  at$names <- NULL
  tmp = cont %*% 2L^c(0L:(nLevel - 1L))
  tmp = as.integer(tmp)
  data = unlist(data, FALSE, FALSE)
  ret = tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  attr(ret, 'inapp.level') <- 2^(which(at$levels == "-") - 1)
  attr(ret, 'dim') <- c(nChar, nTip)
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- '*phyDat'
  ret
}
