#' @title Calculate parsimony score with inapplicable data
#'
#' @description Uses code modified from the Morphy library to calculate a parsimony score 
#' in datasets that contain inapplicable data
#'
#' @template treeParam 
#' @template datasetParam
#' 
#' @examples
#' data(inapplicable.datasets)
#' tree <- TreeSearch::RandomTree(inapplicable.phyData[[1]])
#' result <- InapplicableFitch(tree, inapplicable.phyData[[1]])
#' 
#' @return This function returns the elements from a list containing:
#'    \itemize{
#' \item     The total parsimony score
#' \item     The parsimony score associated with each character 
#' \item     A matrix comprising character reconstructions for each node after the final pass
#'   }
#' The elements to return are specified by the parameter \code{detail}.  
#' If a single element is requested (default) then just that element will be returned
#' If multiple elements are requested then these will be returned in a list.
#' 
#' @seealso \code{\link{MorphyDat}}
#' @seealso \code{\link{TreeSearch}}
#' 
#' @author Martin R. Smith (using C code adapted from MorphyLib, author Martin Brazeau)
#' @importFrom phangorn phyDat
#' @importFrom TreeSearch Renumber RenumberTips
#' @export
InapplicableFitch <- function (tree, dataset) {
  if (class(dataset) != 'phyDat') stop('Invalid data type ', class(dataset), '; should be phyDat.')
  tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
  morphyObj <- LoadMorphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  MorphyLength(tree, morphyObj)
}

#' Calculate parsimony score with inapplicable data
#' 
#' @template labelledTreeParam
#' @template morphyObjParam
#'
#' @return The length of the tree (after weighting)
#'
#' @seealso LoadMorphy
#'
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeSearch Postorder
#' @export
MorphyLength <- function (tree, morphyObj) {
  nTaxa <- mpl_get_numtaxa(morphyObj)
  if (nTaxa < 1) stop("Error: ", mpl_translate_error(nTaxa))
  if (nTaxa != length(tree$tip.label)) stop ("Number of taxa in morphy object (", nTaxa, ") not equal to number of tips in tree")
  treeOrder <- attr(tree, 'order')
  if (is.null(treeOrder) || treeOrder != "postorder") tree <- TreeSearch::Postorder(tree)
  tree.edge <- tree$edge
  parent <- tree.edge[, 1]
  child <- tree.edge[, 2]
  maxNode <- nTaxa + mpl_get_num_internal_nodes(morphyObj)
  rootNode <- nTaxa + 1
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(1:maxNode, child )]
  # parentOf[rootNode] <- maxNode + 1 # Root node's parent is a dummy node
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]
  
  ret <- .Call('MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj)
  return(ret)
}


