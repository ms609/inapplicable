#' Random postorder tree
#' 
#' @param nTip number of tips
#' @return a list of three integer vectors: 
#'          First entry: parentOf: For each node, numbered in postorder, the number of its parent node.
#'          Second entry: leftChild: For each internal node, numbered in postorder, the number of its left 
#'                   child node or tip.
#'          Third entry: rightChild: For each internal node, numbered in postorder, the number of its right
#'                   child node or tip.
#' @export
RandomPostorder <- function (nTip, morphyObj) {  
  # Return:
  .Call('BUILD_POSTORDER', as.integer(nTip), morphyObj)
}