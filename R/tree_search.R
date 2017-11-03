#' Parsimony Ratchet
#'
#' \code{RatchetSearch} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @template datasetParam
#' @param keepAll Set to \code{TRUE} to report all MPTs encountered during the search, perhaps to analyze consensus
#' @param maxIt   maximum ratchet iterations to perform;
#' @param maxIter maximum rearrangements to perform on each bootstrap or ratchet iteration;
#' @param maxHits maximum times to hit best score before terminating a tree search within a pratchet iteration;
#' @param k stop when k ratchet iterations have found the same best score;
#' @template verbosityParam
#' @param rearrangements (list of) function(s) to use when rearranging trees
#'        e.g. \code{list(TreeSearch::RootedTBR, TreeSearch::RootedNNI)}
#' @param \dots other arguments to pass to subsequent functions.
#' @param nSearch Number of Ratchet searches to conduct (for RatchetConsensus)
#' 
#' @return This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.
#'
#' @references Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.}
#'  Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
#'
#' @author Martin R. Smith
#' 
#' Adapted from \code{\link[phangorn]{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
#' 
#' @seealso \code{\link[phangorn]{pratchet}}
#' @seealso \code{\link{BasicSearch}}
#' @seealso \code{\link{SectorialSearch}}
#' 
#' @examples{
#' data('inapplicable.datasets')
#' my.phyDat <- inapplicable.phyData[[1]]
#' RatchetSearch(tree=TreeSearch::RandomTree(my.phyDat, root=names(my.phyDat)[1]), 
#'         dataset=my.phyDat, maxIt=1, maxIter=50)
#' }
#' @keywords  tree 
#' @importFrom TreeSearch Renumber RenumberTips
#' @export
RatchetSearch <- function 
(tree, dataset, keepAll=FALSE, maxIt=100, maxIter=5000, 
  maxHits=40, k=10, verbosity=1, rearrangements=list(TreeSearch::RootedTBR, TreeSearch::RootedSPR, 
  TreeSearch::RootedNNI), ...) {
  if (class(dataset) != 'phyDat') stop("dataset must be of class phyDat, not", class(dataset))
  morphyObj <- LoadMorphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
  eps <- 1e-08
  if (is.null(attr(tree, "score"))) {
    attr(tree, "score") <- MorphyLength(tree, morphyObj, ...)
  }
  best.pars <- attr(tree, "score")
  if (verbosity > 0) cat("* Initial score:", best.pars)
  if (keepAll) forest <- vector('list', maxIter)

  if (class(rearrangements) == 'function') rearrangements <- list(rearrangements)
  kmax <- 0 
  for (i in 1:maxIt) {
    if (verbosity > 0) cat ("\n* Ratchet iteration", i, "- Running NNI on bootstrapped dataset. ")
    candidate <- inapplicable::MorphyBootstrapTree(tree=tree, morphyObj=morphyObj, maxIter=maxIter, maxHits=maxHits,
                               verbosity=verbosity-1, ...)
    
    for (Rearrange in rearrangements) {
      if (verbosity > 0) cat ("\n - Rearranging new candidate tree...")
      candidate <- DoTreeSearch(candidate, morphyObj, Rearrange=Rearrange, 
                                verbosity=verbosity-1, maxIter=maxIter, maxHits=maxHits, ...)
    }
    cand.pars <- attr(candidate, 'score')
    if((cand.pars+eps) < best.pars) {
      if (keepAll) {
        forest <- vector('list', maxIter)
        forest[[i]] <- candidate
      }
      tree <- candidate
      best.pars <- cand.pars
      kmax <- 1
    } else {
      if (best.pars+eps > cand.pars) { # i.e. best == cand, allowing for floating point error
        kmax <- kmax + 1
        candidate$tip.label <- names(dataset)
        tree <- candidate
        if (keepAll) forest[[i]] <- candidate
      }
    }
    if (verbosity > 0) cat("\n* Best score after", i, "/", maxIt, "Ratchet iterations:", 
                           best.pars, "( hit", kmax, "/", k, ")")
    if (kmax >= k) break()
  } # for
  if (verbosity > 0)
    cat ("\nCompleted parsimony ratchet with score", best.pars, "\n")
    
  if (keepAll) {
    forest <- forest[!vapply(forest, is.null, logical(1))]
    class(forest) <- 'multiPhylo'
    ret <- unique(forest)
    cat('Found', length(ret), 'unique MPTs.')
  } else {
    ret <- tree
    ret$tip.label <- names(dataset)
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}

#' @describeIn RatchetSearch returns a list of optimal trees produced by nSearch Ratchet searches
#' @export
RatchetConsensus <- function (tree, dataset, maxIt=5000, maxIter=500, maxHits=20, k=10, verbosity=0, 
  rearrangements=list(TreeSearch::NNI), nSearch=10, ...) {
  trees <- lapply(1:nSearch, function (x) inapplicable::RatchetSearch(tree, dataset, maxIt, maxIter, maxHits, 
                                                  k=1, verbosity, rearrangements, ...))
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from ', nSearch, 'searches.')
  return (trees)
}

#' Bootstrap tree search with inapplicable data
#' 
#' @template labelledTreeParam
#' @template morphyObjParam
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits maximum number of hits to accomplish in tree search
#' @template verbosityParam
#' @param \dots further parameters to send to DoTreeSearch
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @importFrom TreeSearch RootedNNI
#' @export
MorphyBootstrapTree <- function (tree, morphyObj, maxIter, maxHits, verbosity=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  v <- rep(eachChar, startWeights)
  BS <- tabulate(sample(v, replace=TRUE), length(startWeights))
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, BS[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)#
  attr(tree, 'score') <- NULL
  res <- DoTreeSearch(tree, morphyObj, Rearrange=TreeSearch::RootedNNI, maxIter=maxIter, maxHits=maxHits, verbosity=verbosity-1, ...)
  attr(res, 'score') <- NULL
  attr(res, 'hits') <- NULL
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  mpl_apply_tipdata(morphyObj)
  res
}

#' DoTreeSearch
#'
#' Performs a tree search
#' 
#' Does the hard work of searching for a most parsimonious tree.
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by RatchetSearch and Sectorial functions
#'
#' @template labelledTreeParam
#' @template morphyObjParam
#' @param Rearrange Function to use to rearrange trees; example: \code{TreeSearch::\link[TreeSearch]{RootedTBR}}.
#' @param maxIter maximum iterations to conduct.
#' @param maxHits stop search after this many hits.
#' @param forestSize how many trees to hold.
#' @template clusterParam
#' @template verbosityParam
#' @param \dots additional variables to pass to \code{\link{MorphyRearrangeTree}}.
#'
#' @author Martin R. Smith
#' 
#' @keywords internal
#' @export

DoTreeSearch <- function 
(tree, morphyObj, Rearrange, maxIter=100, maxHits=20, forestSize=1, cluster=NULL, 
 verbosity=1, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 0
  if (!is.null(forestSize) && length(forestSize)) {
    if (forestSize > 1) {
      forest <- empty.forest <- vector('list', forestSize)
      forest[[1]] <- tree
    } else {
      forestSize <- 1 
    }
  }
  if (is.null(attr(tree, 'score'))) attr(tree, 'score') <- MorphyLength(tree, morphyObj)
  best.score <- attr(tree, 'score')
  if (verbosity > 0) cat("\n  - Performing search.  Initial score:", best.score)
  return.single <- !(forestSize > 1)
  
  for (iter in 1:maxIter) {
    trees <- inapplicable::MorphyRearrangeTree(tree, morphyObj, Rearrange, min.score=best.score, 
                           return.single=return.single, iter=iter, cluster=cluster,
                           verbosity=verbosity, ...)
    iter.score <- attr(trees, 'score')
    if (length(forestSize) && forestSize > 1) {
      hits <- attr(trees, 'hits')
      if (iter.score == best.score) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      } else if (iter.score < best.score) {
        best.score <- iter.score
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iter.score <= best.score) {
        best.score <- iter.score
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxHits) break
  }
  if (verbosity > 0) cat("\n  - Final score", attr(tree, 'score'), "found", attr(tree, 'hits'), "times after", iter, "rearrangements\n")  
  if (forestSize > 1) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- best.score
    return (unique(forest))
  } else tree
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; 
#'        edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @param Rearrange Function used to rearrange trees; default: \code{\link{RootedTBR}}
#' @param maxIter the maximum number of iterations to perform before abandoning the search;
#' @param maxHits the maximum times to hit the best score before abandoning the search;
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}};
#' @template verbosityParam
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return{
#' This function returns a tree, with an attribute \code{score} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @seealso
#' \itemize{
#' \item \code{\link{InapplicableFitch}}, calculates parsimony score, supports inapplicable tokens;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{RatchetSearch}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' library('ape'); library('phangorn')
#' data('inapplicable.datasets')
#' my.phyDat <- inapplicable.phyData[[1]]
#' outgroup <- names(my.phyDat)[1]
#' njtree <- ape::root(ape::nj(phangorn::dist.hamming(my.phyDat)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL
#' njtree <- ape::root(njtree, outgroup, resolve.root=TRUE)
#'
#' \dontrun{
#' TreeSearch(njtree, my.phyDat, maxIter=20, Rearrange=TreeSearch::NNI)
#' TreeSearch(njtree, my.phyDat, maxIter=20, Rearrange=TreeSearch::RootedSPR)
#' }
#' 
#' @keywords  tree 
#' 
#' @importFrom TreeSearch Renumber RenumberTips
#' @export
BasicSearch <- function 
(tree, dataset, Rearrange=TreeSearch::RootedTBR, maxIter=100, maxHits=20, forestSize=1, verbosity=1, ...) {
  # Initialize morphy object
  if (class(dataset) != 'phyDat') stop ("dataset must be of class phyDat, not ", class(dataset))
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
  morphyObj <- LoadMorphy(dataset)
  # TODO this is the place to initialise a cluster, if that's worth doing.
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  ret <- DoTreeSearch(tree, morphyObj, Rearrange, maxIter, maxHits, forestSize, cluster=NULL, 
                      verbosity, ...)
  return (ret)
}

### #' Sectorial Search
### #'
### #' \code{SectorialSearch} performs a sectorial search on a tree, preserving the position of the root.
### #'
### #' \code{InapplicableSectorial} performs a sectorial search on the tree specified. A sectorial search 
### #' detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
### #' to the main tree (Goloboff, 1999).
### #' The improvement to local \var{score} hopefully (but not necessarily) improves the overall \var{score}.
### #' As such, the output of \code{InapplicableSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
### #' rearrangements and only retained if the ultimate parsimony score is better than 
### #' that of the original tree.
### #' 
### #' \code{SectorialSearch} is a basic recipe that runs \code{InapplicableSectorial} followed by a few rounds
### #' of tree rearrangement, returning a tree whose \var{score} is no worse than that of \code{start.tree}.
### #' 
### #' @param tree a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;
### #' @template datasetParam
### #' @template concavityParam
### #' @param maxIter maximum number of rearrangements to perform on each sectorial iteration;
### #' @template clusterParam
### #' @template verbosityParam
### #' @param rearrangements method to use when rearranging subtrees: NNI, SPR or TBR;
### #' @param \dots other arguments to pass to subsequent functions.
### #' 
### #' @return a rooted tree of class \code{phylo}.
### #' 
### #' @references Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.} Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
### #' 
### #' @author Martin R. Smith
### #' 
### #' @seealso \code{\link{BasicSearch}}
### #' @seealso \code{\link{RatchetSearch}}
### #' 
### #' @examples
### #' require('ape')
### #' data('SigSut')
### #' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
### #' njtree <- ape::root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
### #' njtree$edge.length <- NULL; njtree<-ape::root(njtree, outgroup, resolve.root=TRUE)
### #' InapplicableSectorial(njtree, SigSut.phy, maxIt=1, maxIter=50, largest.sector=7)
### #' \dontrun{SectorialSearch(njtree, SigSut.phy) # Will be time-consuming }
### #' 
### #' 
### #' @keywords  tree 
### #' @importFrom TreeSearch Renumber RenumberTips RootedNNI RootedSPR RootedTBR
### #' @export
### SectorialSearch <- function
### (tree, dataset, concavity = NULL, rearrangements='NNI', maxIter=2000, cluster=NULL, verbosity=3, ...) {
###   best.score <- attr(tree, 'score')
###   tree <- TreeSearch::RenumberTips(TreeSearch::Renumber(tree), names(dataset))
###   if (length(best.score) == 0) best.score <- InapplicableFitch(tree, dataset, ...)[[1]]
###   sect <- InapplicableSectorial(tree, dataset, cluster=cluster,
###     verbosity=verbosity-1, maxIt=30, maxIter=maxIter, maxHits=15, smallest.sector=6, 
###     largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
###   sect <- TreeSearch(sect, dataset, Rearrange=TreeSearch::RootedNNI, maxIter=maxIter, maxHits=30, cluster=cluster, verbosity=verbosity)
###   sect <- TreeSearch(sect, dataset, Rearrange=TreeSearch::RootedTBR, maxIter=maxIter, maxHits=20, cluster=cluster, verbosity=verbosity)
###   sect <- TreeSearch(sect, dataset, Rearrange=TreeSearch::RootedSPR, maxIter=maxIter, maxHits=50, cluster=cluster, verbosity=verbosity)
###   sect <- TreeSearch(sect, dataset, Rearrange=TreeSearch::RootedNNI, maxIter=maxIter, maxHits=60, cluster=cluster, verbosity=verbosity)
###   if (attr(sect, 'score') <= best.score) {
###     return (sect)
###   } else return (tree)
### }
