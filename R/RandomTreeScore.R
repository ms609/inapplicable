#' Parsimony score of random postorder tree
#' 
#' @param nTip number of tips (minimum 3)
#' @template morphyObjParam
#'
#' @return the parsimony score of a random tree, for the given Morphy object.
#'
#' @export
RandomTreeScore <- function (nTip, morphyObj) {  
  if (nTip < 3) {
    warning("nTip < 3 not implemented, as there's only one unrooted topology.")
    return (0)
  }
  # Return:
  .Call('BUILD_POSTORDER', as.integer(nTip), morphyObj)
}