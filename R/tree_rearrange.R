
#' Rearrange phylogenetic tree
#' @details \code{MorphyRearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{score}, the tree's optimality score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function.
#' @template morphyObjParam
#' @param Rearrange a rearrangement function that returns a tree: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}}.
#' @param  minScore trees longer than \code{minScore}, probably the score of the starting tree,
#'     will be discarded.
#' @param  returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
#'   \item{iter}{iteration number of calling function, for reporting to user only.
#' @template clusterParam
#' @template verbosityParam
#' 
#' @return{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{score}
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
MorphyRearrangeTree <- function (tree, morphyObj, Rearrange, minScore=NULL, returnSingle=TRUE,
                           iter='?', cluster=NULL, verbosity=0L) {
  if (is.null(attr(tree, 'score'))) bestScore <- 1e+07 else bestScore <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    rearrangedTree <- Rearrange(tree)
    #rearrangedTree <- RenumberTips(Rearrange(tree), tipOrder)
    trees <- list(rearrangedTree)
    minScore <- MorphyTreeLength(rearrangedTree, morphyObj)
    best.trees <- c(TRUE)
  } else {
    stop("Cluster not implemented.")
    # candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'score') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    # scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    # candidates <- lapply(seq_along(cl), function (x) Rearrange(tree)) # TODO don't pick the same tree twice
    # warning("Not tested; likely to fail.")
    # 
    # scores <- parLapply(cluster, seq_along(cluster), function (i) MorphyTreeLength(candidates[[i]], morphyObj[[i]])) # ~3x faster to do this in serial in r233.
    # minScore <- min(scores)
    # best.trees <- scores == minScore
    # trees <- candidates[best.trees]
  }
  if (bestScore < minScore) {
    if (verbosity > 3L) cat("\n    . Iteration", iter, '- Min score', minScore, ">", bestScore)
  } else if (bestScore == minScore) {
    hits <- hits + sum(best.trees)
    if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", minScore, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", minScore, "found on", hits, "trees")
  }
  if (length(returnSingle) && returnSingle) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- minScore
  trees
}

#' @describeIn MorphyRearrangeTree optimised version that requires parent and child vectors to be extracted from a tree
MorphyRearrange <- function (parent, child, morphyObj, Rearrange, minScore=NULL, returnSingle=TRUE,
                           iter='?', cluster=NULL, verbosity=0L) {
  if (is.null(attr(tree, 'score'))) bestScore <- 1e+07 else bestScore <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    rearrangedTree <- Rearrange(tree)
    #rearrangedTree <- RenumberTips(Rearrange(tree), tipOrder)
    trees <- list(rearrangedTree)
    minScore <- MorphyTreeLength(rearrangedTree, morphyObj)
    best.trees <- c(TRUE)
  } else {
    stop("Cluster not implemented.")
    # candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'score') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    # scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    # candidates <- lapply(seq_along(cl), function (x) Rearrange(tree)) # TODO don't pick the same tree twice
    # warning("Not tested; likely to fail.")
    # 
    # scores <- parLapply(cluster, seq_along(cluster), function (i) MorphyTreeLength(candidates[[i]], morphyObj[[i]])) # ~3x faster to do this in serial in r233.
    # minScore <- min(scores)
    # best.trees <- scores == minScore
    # trees <- candidates[best.trees]
  }
  if (bestScore < minScore) {
    if (verbosity > 3L) cat("\n    . Iteration", iter, '- Min score', minScore, ">", bestScore)
  } else if (bestScore == minScore) {
    hits <- hits + sum(best.trees)
    if (verbosity > 2L) cat("\n    - Iteration", iter, "- Best score", minScore, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (verbosity > 1L) cat("\n    * Iteration", iter, "- New best score", minScore, "found on", hits, "trees")
  }
  if (length(returnSingle) && returnSingle) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- minScore
  trees
}
