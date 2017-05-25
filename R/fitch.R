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
#' @export

InapplicableFitch <- function (tree, morphyData, detail=1, ...) {
  # Data
  if (class(morphyData) == 'phyDat') morphyData <- MorphyDat(morphyData)
  if (class(morphyData) != 'morphyDat') stop('Invalid data type ', class(morphyData), '; try InapplicableFitch(tree, data <- MorphyData(valid.phyDat.object)).')
  at <- attributes(morphyData)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
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
  
  ret <- .Call("MORPHYFITCH", t(morphyData[tipLabel, ]), as.integer(nChar), as.integer(nTip), 
               as.integer(parent), as.integer(child), as.integer(parentOf), as.integer(childOf), 
               as.double(weight), as.integer(inappLevel), as.integer(inappChars), PACKAGE='inapplicable')
  
  if (length(detail) == 1) return (ret[[detail]])
  return (ret[detail])
}
