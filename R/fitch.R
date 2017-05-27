#' @title Calculate parsimony score with inapplicable data
#'
#' @description Uses code modified from the Morphy library to calculate a parsimony score 
#' in datasets that contain inapplicable data
#'
#' @param tree A tree of class \code{phylo}
#' @param morphyData A \code{phyDat} or \code{morphyDat} object, perhaps generated with 
#'  \code{\link{phyDat}} or \code{\link{MorphyDat}}
#' @param detail Leave as 1 to just return parsimony score, or specify c(1, 2, 3) for additional detail (see below)
#' 
#' @examples
#' data(SigSut)
#' taxa <- names(SigSut.phy)
#' tree <- rtree(length(taxa), tip.label=taxa, br=NULL)
#' result <- InapplicableFitch(tree, SigSut.phy)
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
#' @author Martin Smith (using C code adapted from MorphyLib, author Martin Brazeau)
#' importFrom phangorn phyDat
#' importFrom ape tiplabels nodelabels read.tree is.rooted .PlotPhyloEnv rtree
#' importFrom graphics plot locator par text
#' importFrom stats runif reorder na.omit
#' @export
InapplicableFitch <- function (tree, dataset, ...) {
  if (class(dataset) != 'phyDat') stop('Invalid data type ', class(phyData), '; should be phyDat.')
  morphyObj <- LoadMorphy(dataset)
  result <- MorphyLength(tree, morphyObj)
  morphyObj <- UnloadMorphy(morphyObj)
  result
}

#' @title Calculate parsimony score with inapplicable data
#' 
#' @param tree Phylogenetic tree (of class \code{phylo})
#' @param morphyObj Morphy object containing character data, constructed using \code{\link{LoadMorphy}}
#' @return The length of the tree (after weighting)
#'
#' @seealso LoadMorphy
#'
#' @author Martin R. Smith
MorphyLength <- function (tree, morphyObj) {
  nTaxa <- mpl_get_numtaxa(morphyObj)
  if (nTaxa < 1) stop("Error: ", mpl_translate_error(nTaxa))
  if (nTaxa != length(tree$tip.label)) stop ("Number of taxa in morphy object (", nTaxa, ") not equal to number of tips in tree")
  treeOrder <- attr(tree, 'order')
  if (is.null(treeOrder) || treeOrder == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[ ,1]
  child <- tree.edge[, 2]
  maxNode <- parent[1] #max(parent)
  rootNode <- nTaxa + 1
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(1:maxNode, child )]
  # parentOf[rootNode] <- maxNode + 1 # Root node's parent is a dummy node
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]
  
  ret <- .Call('MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj, PACKAGE='inapplicable')
  return(ret)
}


