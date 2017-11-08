#' Parsimony score of random postorder tree
#' 
#' @param nTip number of tips
#' @template morphyObjParam
#'
#' @return the parsimony score of a random tree, for the given Morphy object.
#'
#' @export
RandomTreeScore <- function (nTip, morphyObj) {  
  # Return:
  .Call('BUILD_POSTORDER', as.integer(nTip), morphyObj)
}