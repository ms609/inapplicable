#' TITLE GOES HERE
#'
#' \code{FUNCTIONNAME} does something useful
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin R. Smith
#' @export
MinSteps <- function (x, inapp.power2) {
  if (is.null(x)) return (NULL)
  vals = AsBinary(unique(as.integer(x)))
  number.of.vals = ncol(vals)
  if (is.null(number.of.vals)) return (0) 
  if (number.of.vals >= inapp.power2) vals[,inapp.power2] = 0
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
    if (max(sums)) definites[which (sums == max(sums))[1]] = TRUE
  }
  return (max(0, sum(definites) - 1))
}
#' TITLE GOES HERE
#'
#' \code{FUNCTIONNAME} does something useful
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin R. Smith
#' @export
MinStepsInapp <- function (dat, inapp.level) {
  # This function gives unconventional results; see text.
  apply(dat, 1, function (x) {
    vals = AsBinary(unique(x))
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
    inapps = sum(x == 2^(inapp.level-1)) # Two inapps are the minimum necessary to imply an additional origin of the character
    return (max(0, sum(definites) - 1 - (inapps/2)))
  })
}