
#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{pscore}, the tree's parsimony score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function;
#' @template morphyObjParam
#' @param Rearrange a rearrangement function that returns a tree: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};
#' @param  min.score trees longer than \code{min.score}, probably the score of the starting tree,
#'     will be discarded;
#' @param  return.single returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
#'   \item{iter}{iteration number of calling function, for reporting to user only;
#' @template clusterParam
#' @template verbosityParam
#' 
#' @return{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}
#'  as described for argument \code{tree}, and with tip labels ordered to match morphyObj.}
#' @author Martin R. Smith
#' @seealso
#'   \itemize{
#'     \item \code{\link{RootedNNI}}
#'     \item \code{\link{RootedSPR}}
#'     \item \code{\link{RootedTBR}}
#'   }
#' 
#' @importFrom parallel clusterCall
#' @importFrom TreeSearch Renumber RenumberTips
#' @export
RearrangeTree <- function (tree, morphyObj, Rearrange, min.score=NULL, return.single=TRUE,
                           iter='?', cluster=NULL, verbosity=0) {
  if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  tipOrder <- tree$tip.label
  if (is.null(cluster)) {
    rearrangedTree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), tipOrder)
    trees <- list(rearrangedTree)
    min.score <- MorphyLength(rearrangedTree, morphyObj)
    best.trees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'pscore') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    warning("Not tested; likely to fail.")
    candidates <- clusterCall(cluster, Rearrange, tree)
    candidates <- lapply(candidates, TreeSearch::RenumberTips, tipOrder)
    scores <- vapply(candidates, MorphyLength, 1, morphyObj) # ~3x faster to do this in serial in r233.
    min.score <- min(scores)
    best.trees <- scores == min.score
    trees <- candidates[best.trees]
  }
  if (best.score < min.score) {
    if (verbosity > 3) cat("\n    . Iteration", iter, '- Min score', min.score, ">", best.score)
  } else if (best.score == min.score) {
    hits <- hits + sum(best.trees)
    if (verbosity > 2) cat("\n    - Iteration", iter, "- Best score", min.score, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (verbosity > 1) cat("\n    * Iteration", iter, "- New best score", min.score, "found on", hits, "trees")
  }
  if (length(return.single) && return.single) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'pscore') <- min.score
  trees
}
