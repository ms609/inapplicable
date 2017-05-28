#' TITLE GOES HERE
#'
#' \code{FUNCTIONNAME} does something useful
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin Smith
#' @export
InapplicableSectorial <- function (tree, data, maxit=100, 
    maxiter=500, k=5, verbosity=0, smallest.sector=4, largest.sector=1e+06, rearrangements="NNI", criterion=NULL, ...) {
  if (class(data) != 'phyDat') stop("data must be a phyDat object.")
  if (is.null(tree)) stop("a starting tree must be provided")
  if (verbosity >= 0) cat('InapplicableSectorial search: optimizing sectors of', smallest.sector, 'to', floor(largest.sector), 'tips')
  
  SectorData <- function (X, tips) {
    at <- attributes(X)
    dec <- X[tips, ]
    nBits <- floor(log2(max(X))) + 1L
    bin <- array(FALSE, dim=c(nrow(dec), ncol(dec), nBits))
    for (i in 0:(nBits-1)) {
      bin[, , nBits-i] <- as.logical(dec %% 2)
      dec <- (dec %/% 2)
    }
    state.union <- apply(bin, c(1,3), all)
    parsimony.informative <- !as.logical(rowSums(state.union))
    if (!any(parsimony.informative)) return (NULL)
    X <- X[tips, parsimony.informative]
    informative.chars <- sum(parsimony.informative)
    SECTOR_ROOT <- rep(2^nBits-1, informative.chars)
    X <- cbind(X, SECTOR_ROOT)
    attr(X, 'nr') <- informative.chars
    attr(X, 'inapp.level') <- at$inapp.level
    inapp.power2 <- log2(at$inapp.level) + 1
    #attr(X, 'min.steps') <- apply(X, 1, function(x) min.steps(x, inapp.power2))
    attr(X, 'levels') <- at$levels
    attr(X, 'weight') <- at$weight[parsimony.informative]
    class(X) <- 'morphyDat'
    warning("#TODO this is not yet tested")
    X
  }
  
  eps <- 1e-08
  kmax <- 1
  for (i in 1:maxit) {
    edge1 <- tree$edge[,1]
    nodes <- unique(edge1)[-1]
    node.lengths <- sapply(GetDescendants(tree, nodes), length) # (10x quicker than DoDescendants)
    candidate.nodes <- nodes[node.lengths >= smallest.sector & node.lengths <= largest.sector]
    if (verbosity >= 0) cat ("\n - Iteration", i, "- attempting sectorial search on node ")
    repeat {
      sector <- sample(candidate.nodes, 1)
      crown <- ExtractClade(tree, sector)
      crown.tips <- crown$tip.label
      sector.size <- length(crown.tips)
      cat(sector, 'size', sector.size, '...')
      crown.data <- SectorData(data, crown.tips)
      if (!is.null(crown.data)) break else cat('unsuitable (no data); trying')
      candidate.nodes <- candidate.nodes[-which(candidate.nodes==sector)]
      if (length(candidate.nodes == 0)) stop('No selectable sectors contain parsimony information! Either "largest.sector" is close to "smallest.sector" or your dataset is short of parsimony information.')
    } 
    if (verbosity >= 0) cat(' Sector OK.')
    crown <- root(AddTip(crown, 0, 'SECTOR_ROOT'), length(crown$tip.label) + 1, resolve.root=TRUE) ## TODO use Root or add ape::root to includeFrom in NAMESPACE
    initial.p <- InapplicableFitch(crown, crown.data, ...)
    attr(crown, 'pscore') <- initial.p
    if (verbosity >= 0) cat("\n - Running", rearrangements, "search on sector", sector)
    candidate <- TreeSearch(crown, crown.data, 'SECTOR_ROOT', method=rearrangements, criterion=criterion, verbosity=verbosity-1, maxiter=maxiter, ...)
    candidate.p <- attr(candidate, 'pscore')
    
    if((candidate.p + eps) < initial.p) {
      kmax <- kmax + 1
      stump <- DropTip(tree, GetDescendants(tree, sector)[[1]], subtree=TRUE)
      stump.edge <- 1:nrow(stump$edge)
      stump$root.edge <- 1
      crown <- DropTip(candidate, 'SECTOR_ROOT')
      tree <- CollapseSingles((BindTree(stump, crown, where=which(stump$tip.label==paste('[', sector.size, '_tips]', sep="")), position=0)))
      if (verbosity > 0) cat(' : improved local pscore, updated tree')
    } else if (verbosity > 0) cat (' : no improvement to local pscore')
    if (kmax == k) break()
  } # for
  if (verbosity >= 0)
    cat ("\nCompleted sectorial rearrangements.\n")
  attr(tree, 'pscore') <- NULL
  attr(tree, 'hits') <- NULL
  tree
}  # InapplicableSectorial

#' @name InapplicablePratchet
#' @aliases InapplicablePratchet
#'  Parsimony ratchet
#' @description This function uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#' \usage{
#' InapplicablePratchet(tree, data, concavity = NULL, all = FALSE, outgroup = NULL, maxit = 100, 
#'   maxiter = 5000, maxhits = 40, k = 10, verbosity = 0, rearrangements = "NNI", ...)
#' }
#' \arguments{
#'   \item{tree}{An object of class \code{phyDat} denoting the topology to begin the search from;}
#'   \item{data}{A matrix opf class \code{morphyDat} format (or \code{\link{phyDat}} format,
#'     which the function will itself pass through \code{MorphyDat});}
#'   \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
#'     see \code{\link{InapplicableParsimony}};}
#'   \item{all}{Set to TRUE to report all MPTs encountered during the search, perhaps to analyze consensus}
#'   \item{outgroup}{a vector specifying all tips in the outgroup; if unspecified then identical trees with different roots will be considered unique;}
#'   \item{maxit}{maximum ratchet iterations to perform;}
#'   \item{maxiter}{maximum rearrangements to perform on each bootstrap or ratchet iteration;}
#'   \item{maxhits}{maximum times to hit best score before terminating a tree search within a pratchet iteration;}
#'   \item{k}{stop when k ratchet iterations have found the same best score;}
#'   \item{verbosity}{larger numbers provides more verbose feedback to the user;}
#'   \item{rearrangements}{method to use when rearranging trees: NNI, SPR or TBR;}
#'   \item{\dots}{other arguments to pass to subsequent functions.}
#' }
#' @return This function returns a tree modified by parsimony ratchet iteration, retaining the position of the root.
#' @references Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
#' \author{
#' Martin R. Smith
#' 
#' Adapted from \code{\link{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
#' }
#' \seealso{
#' \itemize{
#' \item \code{\link{pratchet}}
#' \item \code{\link{TreeSearch}}
#' \item \code{\link{SectorialSearch}}
#' }
#' }
#' @examples{
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#' InapplicablePratchet(njtree, SigSut.phy, outgroup, maxit=1, maxiter=50)
#' }
#' @keywords  tree 

#' @export
InapplicablePratchet <- function (tree, dataset, keepAll=FALSE, outgroup=NULL, maxit=100, maxiter=5000, maxhits=40, k=10, verbosity=0, rearrangements="NNI", criterion=NULL, ...) {
  if (class(dataset) != 'phyDat') stop("dataset must be of class phyDat, not", class(dataset))
  morphyObj <- LoadMorphy(dataset)
  eps <- 1e-0
  if (is.null(attr(tree, "pscore"))) {
    attr(tree, "pscore") <- MorphyLength(tree, morphyObj, ...)
  }
  best.pars <- attr(tree, "pscore")
  if (verbosity >= 0) cat("* Initial pscore:", best.pars)
  if (keepAll) forest <- vector('list', maxiter)

  kmax <- 0 
  for (i in 1:maxit) {
    if (verbosity >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bstree <- BootstrapTree(tree=tree, mDat=dataset, maxiter=maxiter, maxhits=maxhits, criterion=criterion, verbosity=verbosity-1, ...)
    
    if (verbosity >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new cannew candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- TreeSearch(bstree,    dataset, criterion=criterion, method='TBR', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, dataset, criterion=criterion, method='SPR', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, dataset, criterion=criterion, method='NNI', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "TBR only") {            
      candidate <- TreeSearch(bstree,    dataset, criterion=criterion, method='TBR', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR") {                  
      candidate <- TreeSearch(bstree,    dataset, criterion=criterion, method='SPR', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, dataset, criterion=criterion, method='NNI', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR only") {             
      candidate <- TreeSearch(bstree,    dataset, criterion=criterion, method='SPR', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
    } else {                                               
      candidate <- TreeSearch(bstree,    dataset, criterion=criterion, method='NNI', verbosity=verbosity, maxiter=maxiter, maxhits=maxhits, ...)
    }
    #if(class(result)=="phylo") m <- 1
    #else m = length(result)
    #if(m > 0) trees[2 : (1+m)] = result[1:m]
    #pscores <- sapply(trees, function(dataset) attr(dataset, "pscore"))
    cand.pars <- attr(candidate, 'pscore')
    if((cand.pars+eps) < best.pars) {
      if (keepAll) {
        forest <- vector('list', maxiter)
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      }
      tree <- candidate
      best.pars <- cand.pars
      kmax <- 1
    } else {
      if (best.pars+eps > cand.pars) { # i.e. best == cand, allowing for floating point error
        kmax <- kmax + 1
        tree <- candidate
        if (keepAll) forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      }
    }
    if (verbosity >= 0) cat("\n* Best pscore after", i, "/", maxit, "pratchet iterations:", best.pars, "( hit", kmax, "/", k, ")")
    if (kmax >= k) break()
  } # for
  if (verbosity >= 0)
    cat ("\nCompleted parsimony ratchet with pscore", best.pars, "\n")
    
  if (keepAll) {
    forest <- forest[!vapply(forest, is.null, logical(1))]
    class(forest) <- 'multiPhylo'
    ret <- unique(forest)
    cat('Found', length(ret), 'unique MPTs.')
    if (is.null(outgroup)) warning('"outgroup" not specified, so some "unique" trees may have same topology but distinct roots.')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}
#' TITLE GOES HERE
#'
#' \code{FUNCTIONNAME} does something useful
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin Smith
#' @export
PratchetConsensus <- function (tree, data, maxit=5000, maxiter=500, maxhits=20, k=10, verbosity=0, rearrangements="NNI", criterion=NULL, nSearch=10, ...) {
  trees <- lapply(1:nSearch, function (x) InapplicablePratchet(tree, data, maxit, maxiter, maxhits, k=1, verbosity, rearrangements, criterion=criterion, ...))
  scores <- vapply(trees, function (x) attr(x, 'pscore'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from ', nSearch, 'searches.')
  return (trees)
}

#' Bootstrap tree search with inapplicable data
#' 
#' @return A tree that is optimal under a random sampling of the original characters
#' @export
BootstrapTree <- function (tree, morphyObj, maxiter, maxhits, criterion=criterion, verbosity=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  v <- rep(eachChar, startWeights)
  BS <- tabulate(sample(v, replace=TRUE), length(startWeights))
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, BS[i], morphyObj), integer(1))
  res <- DoTreeSearch(tree, morphyObj, method='NNI', criterion=criterion, maxiter=maxiter, maxhits=maxhits, verbosity=verbosity-1, ...)
  attr(res, 'pscore') <- NULL
  attr(res, 'hits') <- NULL
  vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  res
}

#' DoTreeSearch
#'
#' Performs a tree search
#' 
#' Does the hard work of searching for a most parsimonious tree.
#' End-users are expected to access this function through its wrapper, TreeSearch
#' It is also called directly by Ratchet and Sectorial functions
#' TODO remove export line
#'
#' @author Martin R. Smith
#' @export

DoTreeSearch <- function (tree, morphyObj, method='NNI', maxiter=100, maxhits=20, forest.size=1, cluster=NULL, verbosity=1, criterion=NULL, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (!is.null(forest.size) && length(forest.size)) {
    if (forest.size > 1) {
      forest <- empty.forest <- vector('list', forest.size); forest[[1]] <- tree
    } else {
      forest.size <-1 
    }
  }
  if (is.null(attr(tree, 'pscore'))) attr(tree, 'pscore') <- MorphyLength(tree, morphyObj)
  best.pscore <- attr(tree, 'pscore')
  if (verbosity > 0) cat("\n  - Performing", method, "search.  Initial pscore:", best.pscore)
  rearrange.func <- switch(method, 'TBR' = TBR, 'SPR' = SPR, 'NNI' = NNI)
  return.single <- !(forest.size > 1)
  
  for (iter in 1:maxiter) {
    trees <- RearrangeTree(tree, morphyObj, rearrange.func, min.score=best.pscore, return.single=return.single, iter=iter, cluster=cluster, criterion=criterion, verbosity=verbosity)
    iter.pscore <- attr(trees, 'pscore')
    if (length(forest.size) && forest.size > 1) {
      hits <- attr(trees, 'hits')
      if (iter.pscore == best.pscore) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'pscore') <- iter.pscore
        attr(tree, 'hits') <- hits
      } else if (iter.pscore < best.pscore) {
        best.pscore <- iter.pscore
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'pscore') <- iter.pscore
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iter.pscore <= best.pscore) {
        best.pscore <- iter.pscore
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxhits) break
  }
  if (verbosity > 0) cat("\n  - Final score", attr(tree, 'pscore'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
  if (forest.size > 1) {
    if (hits < forest.size) forest <- forest[-((hits+1):forest.size)]
    attr(forest, 'hits') <- hits
    attr(forest, 'pscore') <- best.pscore
    return (unique(forest))
  } else tree
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#' 
#' @usage TreeSearch(tree, data, outgroup, concavity = NULL, method = "NNI", maxiter = 100, 
#'   maxhits = 20, forest.size = 1, cluster = NULL, verbosity = 1, ...)
#' 
#' 
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @param data a matrix, ideally of class \code{morphyDat}, generated by \code{\link{MorphyData}}, 
#'   alternatively as a \code{\link{phyDat}} object.  May contain inapplicable data;
#' @param outgroup a vector listing the taxa in the outgroup;
#' @param concavity concavity constant for implied weighting (not currently implemented!); 
#'        see \code{\link{InapplicableParsimony}};
#' @param method rearrangements to perform; one of \kbd{NNI}, \kbd{SPR}, or \kbd{TBR};
#' @param maxiter the maximum number of iterations to perform before abandoning the search;
#' @param maxhits the maximum times to hit the best pscore before abandoning the search;
#' @param forest.size the maximum number of trees to return - useful in concert with \code{\link{consensus}};
#' @param cluster a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;
#' @param verbosity higher values provide more verbose user feedback in stdout;
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return{
#' This function returns a tree, with an attribute \code{pscore} conveying its parsimony score.
#' Note that the parsimony score will be inherited from the tree's attributes, which is only valid if it 
#' was generated using the same \code{data} that is passed here.
#' }
#' @author Martin R. Smith
#'
#' @seealso{
#' \itemize{
#' \item \code{\link{InapplicableParsimony}}, calculates parsimony score, supports inapplicable tokens;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{SectorialSearch}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{InapplicablePratchet}}, alternative heuristic, useful to escape local optima.
#' }}
#' @examples{
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#' \dontrun{
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='NNI')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='SPR')
#' TreeSearch(njtree, SigSut.phy, outgroup, maxiter=20, method='TBR')}
#' }
#' 
#' @keywords  tree 
#' 
#' @export
TreeSearch <- function 
(tree, dataset, method='NNI', maxiter=100, maxhits=20, forest.size=1, 
 cluster=NULL, verbosity=1, criterion=NULL, ...) {
  # Initialize morphy object
  if (class(dataset) != 'phyDat') stop ("dataset must be of class phyDat, not", class(dataset))
  morphyObj <- LoadMorphy(dataset)
  ret <- DoTreeSearch(tree, morphyObj, method, maxiter, maxhits, forest.size, cluster, 
                      verbosity, criterion, ...)
  morphyObj <- UnloadMorphy(morphyObj)
  return (ret)
}

#' @name SectorialSearch
#' @aliases SectorialSearch InapplicableSectorial 
#' 
#' \code{SectorialSearch} performs a sectorial search on a tree, preserving the position of the root.
#' \usage{
#' SectorialSearch(tree, data, outgroup, concavity = NULL, rearrangements = "NNI",
#'   maxiter = 2000, cluster = NULL, verbosity = 3)
#' InapplicableSectorial(tree, data, outgroup = NULL, concavity = NULL, maxit = 100, maxiter = 500, k = 5,
#'   verbosity = 0, smallest.sector = 4, largest.sector = 1e+06, rearrangements = "NNI", ...)
#' }
#' \arguments{
#'   \item{tree}{a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;}
#'   \item{data}{a data matrix in \code{morphyDat} format, perhaps created with \code{\link{MorphyData}}
#'               (\code{phyDat} format also accepted);}
#'   \item{outgroup}{a vector listing the taxa that form the outgroup;}
#'   \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
#'     see \code{\link{InapplicableParsimony}};}
#'   \item{maxit}{maximum number of sectorial iterations to perform;}
#'   \item{maxiter}{maximum number of rearrangements to perform on each sectorial iteration;}
#'   \item{cluster}{a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;}
#'   \item{k}{stop when \code{k} searches have improved their sectorial score;}
#'   \item{verbosity}{integer determining how verbose the reporting to stdout will be;}
#'   \item{smallest.sector}{sectors with fewer than \code{smallest.sector} taxa will not be selected; \kbd{4} is the smallest sensible value;}
#'   \item{largest.sector}{sectors with more than \code{largest.sector} taxa will not be selected;}
#'   \item{rearrangements}{method to use when rearranging subtrees: NNI, SPR or TBR;}
#'   \item{\dots}{other arguments to pass to subsequent functions.}
#' }
#' @details{
#' /code{InapplicableSectorial} performs a sectorial search on the tree specified. A sectorial search 
#' detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
#' to the main tree (Goloboff, 1999).
#' The improvement to local \var{pscore} hopefully (but not necessarily) improves the overall \var{pscore}.
#' As such, the output of \code{InapplicableSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
#' rearrangements and only retained if the ultimate parsimony score is better than 
#' that of the original tree.
#' 
#' \code {SectorialSearch} is a basic recipe that runs \code{InapplicableSectorial} followed by a few rounds
#' of tree rearrangement, returning a tree whose \var{pscore} is no worse than that of \code{start.tree}.
#' }
#' @return a rooted tree of class \code{phylo}.
#' 
#' \references{
#' Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.} Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
#' }
#' \author{
#' Martin R. Smith
#' }
#' 
#' \seealso{
#' \itemize{ 
#' \item \code{\link{TreeSearch}}
#' \item \code{\link{InapplicablePratchet}}
#' }
#' }
#' 
#' @examples{
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#' InapplicableSectorial(njtree, SigSut.phy, outgroup, maxit=1, maxiter=50, largest.sector=7)
#' \dontrun {SectorialSearch(njtree, SigSut.phy, outgroup, 'SPR') # Will be time-consuming}
#' 
#' ## SectorialSearch is currently defined as
#' function (start.tree, data, outgroup, rearrangements='NNI') {
#'   best.score <- attr(start.tree, 'pscore')
#'   if (length(best.score) == 0) best.score <- InapplicableParsimony(start.tree, data)
#'   sect <- InapplicableSectorial(start.tree, data, outgroup=outgroup, verbosity=0, maxit=30, maxiter=200, maxhits=15, smallest.sector=6, largest.sector=length(start.tree$edge[,2])*0.25, rearrangements=rearrangements)
#'   sect <- TreeSearch(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=20, verbosity=3)
#'   sect <- TreeSearch(sect, data, outgroup, method='TBR', maxiter=2000, maxhits=25, verbosity=3)
#'   sect <- TreeSearch(sect, data, outgroup, method='SPR', maxiter=2000, maxhits=50, verbosity=3)
#'   sect <- TreeSearch(sect, data, outgroup, method='NNI', maxiter=2000, maxhits=50, verbosity=3)
#'   if (attr(sect, 'pscore') <= best.score) {
#'     return (sect)
#'   } else return (SetOutgroup(start.tree, outgroup))
#' }
#' }
#' @keywords  tree 
#' @export
SectorialSearch <- function (tree, data, concavity = NULL, rearrangements='NNI', maxiter=2000, cluster=NULL, verbosity=3) {
  best.score <- attr(tree, 'pscore')
  if (length(best.score) == 0) best.score <- InapplicableFitch(tree, data, ...)[[1]]
  sect <- InapplicableSectorial(tree, data, cluster=cluster,
    verbosity=verbosity-1, maxit=30, maxiter=maxiter, maxhits=15, smallest.sector=6, 
    largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
  sect <- TreeSearch(sect, data, method='NNI', maxiter=maxiter, maxhits=30, cluster=cluster, verbosity=verbosity)
  sect <- TreeSearch(sect, data, method='TBR', maxiter=maxiter, maxhits=20, cluster=cluster, verbosity=verbosity)
  sect <- TreeSearch(sect, data, method='SPR', maxiter=maxiter, maxhits=50, cluster=cluster, verbosity=verbosity)
  sect <- TreeSearch(sect, data, method='NNI', maxiter=maxiter, maxhits=60, cluster=cluster, verbosity=verbosity)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (tree)
}
