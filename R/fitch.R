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
#' importFrom ape tiplabels nodelabels read.tree is.rooted .PlotPhyloEnv
#' importFrom parallel makeCluster clusterCall clusterEvalQ setDefaultCluster
#' importFrom graphics plot locator par text
#' importFrom stats runif reorder na.omit
#' @export

InapplicableFitch <- function (tree, morphyData, detail=1, ...) {
  if (class(morphyData) == 'phyDat') morphyData <- MorphyDat(morphyData)
  if (class(morphyData) != 'morphyDat') stop('Invalid data type ', class(morphyData), '; try InapplicableFitch(tree, data <- MorphyData(valid.phyDat.object)).')
  treeOrder <- attr(tree, 'order')
  if (is.null(treeOrder) || treeOrder == "cladewise") tree <- reorder(tree, "postorder")
  at <- attributes(morphyData)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  tree.edge <- tree$edge
  parent <- tree.edge[,1]
  child <- tree.edge[,2]
  tipLabel <- tree$tip.label
  maxNode <- parent[1] #max(parent)
  nTip <- length(tipLabel)
  inappLevel <- at$inapp.level  
  inappChars <- at$inapp.chars
  parentOf <- parent[match(1:maxNode, child )]
  parentOf[nTip + 1] <- nTip + 1 # Root node "is its own parent"
  allNodes <- (nTip + 1L):maxNode
  childOf <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  ret <- .Call("MORPHYFITCH", as.integer(t(morphyData[tipLabel, ])), as.integer(nChar), 
               as.integer(nTip), as.integer(parent), as.integer(child),
               as.integer(parentOf), as.integer(childOf), as.double(weight), 
               as.integer(inappLevel), as.integer(inappChars), PACKAGE='inapplicable')
  
  if (length(detail) == 1) return (ret[[detail]])
  return (ret[detail])
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
#' @export
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
  parentOf[rootNode] <- maxNode + 1 # Root node's parent is a dummy node
  childOf <- child[c(match(allNodes, parent), rootNode, length(parent) + 1L - match(allNodes, rev(parent)), rootNode)]
  
  ret <- .Call('MORPHYLENGTH', as.integer(childOf -1L), as.integer(parentOf -1L), morphyObj, 
               PACKAGE='inapplicable')
  return(ret)
}
