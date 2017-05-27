#' Reorder pruning
#' Modified from phangorn:::reorderPruning
ReorderPruning <- function (x) {
  edge <- x$edge
  parents <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parents[!match(parents, child, 0)][1])
  n_edge <- length(parents)
  max_edge <- max(edge)
  neworder <- .C("phangorn_reorder", parents, child, as.integer(n_edge), 
      as.integer(max_edge), integer(n_edge), as.integer(root - 1L), PACKAGE = "inapplicable")[[5]]
  x$edge <- edge[neworder, ]
  x$edge.length <- x$edge.length[neworder]
  attr(x, "order") <- "pruningwise"
  x
}

#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' @usage #' RearrangeTree(tree, data, rearrange, min.score = NULL, concavity = NULL, return.single = TRUE,
#'  iter = "<unknown>", cluster = NULL, verbosity = 0)
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes
#'     \code{pscore}, the tree's parsimony score, and 
#'     \code{hits}, the number of times the best score has been hit in the calling function;
#'   
#' @param dataset a data matrix in \code{morphyDat} format, perhaps created with \code{\link{MorphyData}};
#' @param Rearrange a rearrangement function: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};
#' @param  min.score trees longer than \code{min.score}, probably the score of the starting tree,
#'     will be discarded;
#' @param concavity concavity constant for implied weighting (not currently implemented!); 
#'     see \code{\link{InapplicableParsimony}};
#' @param  return.single returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
#'   \item{iter}{iteration number of calling function, for reporting to user only;
#' @param  cluster a cluster, prepared with \code{\link{PrepareCluster}}, to accelerate 
#'     searches on multicore machines;
#' @param verbosity determines how much information to output to screen.
#' 
#' @return{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}
#'  as described for argument \code{tree}.}
#' @author Martin R. Smith
#' @seealso
#'   \itemize{
#'     \item \code{\link{RootedNNI}}
#'     \item \code{\link{RootedSPR}}
#'     \item \code{\link{RootedTBR}}
#'   }
#' 
#' @examples
#' data('SigSut')
#' random.tree <- rtree(34, tip.label=names(SigSut.data), br=NULL)
#' RearrangeTree(random.tree, SigSut.preparedata, RootedNNI)
#' 
#' @importFrom parallel clusterCall
#' @export
RearrangeTree <- function (tree, morphyObj, Rearrange, min.score=NULL, concavity=NULL, return.single=TRUE, iter='<unknown>', cluster=NULL, criterion=NULL, verbosity=0) {
  if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    trees <- list(re.tree<-Rearrange(tree))
    min.score <- MorphyLength(re.tree, morphyObj)
    best.trees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'pscore') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    candidates <- clusterCall(cluster, Rearrange, tree)
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
#' @name RootedNNI
#' @aliases RootedNNI
#' @aliases RootedSPR
#' @aliases RootedTBR
#' Rearrange a rooted tree
#' @description This function performs a rearrangement iteration on a tree, retaining the position of the root.
#' \usage{
#' RootedNNI(tree)
#' RootedSPR(tree)
#' RootedTBR(tree)
#' }
#' \arguments{
#'   \item{tree}{An object of class \code{\link{phylo}}, with all nodes resolved (bifurcating).}
#' }
#' \details{
#' A single \acronym{NNI}, \acronym{SPR} or \acronym{TBR} rearrangement is performed, subject to the constraint that 
#' no taxon may be moved to the opposite side of the root node.
#' Branch lengths are not (yet) supported.
#' }
#' @return This function returns a tree, in \code{phylo} format.
#'
#' @author Martin R. Smith
#' 
#' \code{RootedNNI} is abridged from the \pkg{phangorn} function \code{nnin}
#' 
#' @seealso{
#' \itemize{
#' \item \code{\link{SetOutgroup}}, set the outgroup of the phylogenetic tree
#' \item \code{\link{NNI}}, unrooted \acronym{NNI} and \acronym{SPR}
#' \item \code{\link{TBR}}, unrooted \acronym{TBR}
#' }
#' }
#' @examples{
#'   tree <- read.tree(text='(((a,b),c),(d,(e,f)));')
#'   tree <- SetOutgroup(tree, c('e', 'f'))
#'   plot(tree)
#'   dev.new()
#'   plot(RootedNNI(tree))
#'   plot(RootedSPR(tree))
#'   plot(RootedTBR(tree))
#' }
#' 
#' @export
RootedNNI <- function (tree) {
  edge <- matrix(tree$edge, ncol = 2)
  parent <- edge[, 1]
  child <- edge[, 2]
  safe.child <- child
  safe.child[which(parent == as.integer(parent[!match(parent, child, 0)][1]))] <- -1 # Don't want to switch across the root
  k <- min(parent) - 1
  sampleable <- length(na.omit(match(safe.child, parent)))
  n <- sample(sampleable, 1)
  ind <- which(safe.child > k)[n] # Internal nodes
  p1 <- parent[ind]
  p2 <- child[ind]
  ind1 <- which(parent == p1)
  ind1 <- ind1[ind1 != ind][1]
  ind2 <- which(parent == p2)[sample(2,1)]
  e1 <- child[ind1]
  e2 <- child[ind2]
  tree$edge[ind1, 2] <- e2
  tree$edge[ind2, 2] <- e1
  if(!is.null(tree$edge.length)) {
    warning("Edge lengths have been deleted")
    tree$edge.length <- NULL
  }
  tree <- Renumber(ReorderPruning(tree))  
}

#' @export
RootedSPR <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- SetOutgroup(tree, outgroup).")
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- nTips + 1L # Assumes fully-resolved bifurcating tree
  root.children <- child[parent==root]
  left.nodes <- DoDescendants(parent, child, nTips, root.children[1L])
  right.nodes <- !left.nodes
  left.nodes[root.children] <- right.nodes[root.children] <- right.nodes[root] <- FALSE
  size <- c(sum(left.nodes), sum(right.nodes))
  moves <- (size-2L) * (size-1L) / 2
  moves[size < 3] <- 0
  if (!max(moves)) return (tree)
  choose.right <- runif(1, min=0, max=sum(moves)) > moves[1]
  pruning.candidates <- if (choose.right) which(right.nodes) else which(left.nodes)
  subtree.base <- child[parent==root.children[choose.right + 1L]]
  subtree.basal.tip <- subtree.base < root
  if (any(subtree.basal.tip)) pruning.candidates <- pruning.candidates[-match(subtree.base[!subtree.basal.tip], pruning.candidates)]
  prune.node <- sample(pruning.candidates, 1)
  moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
  moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
  dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
  graft.candidates <- c(root.children[choose.right + 1L], pruning.candidates)
  graft.candidates <- c(graft.candidates[!graft.candidates %in% dont.graft.here])
  graft.node <- sample(graft.candidates, 1)
  graft.edge <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  
  leading.edge <- match(prune.parent, child)
  prune.edge <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge <- match(prune.parent, parent.duplicate)
  edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  
  nEdge <- length(child)
  reordered.edge <- .C('order_edges', as.integer(edge[,1]), as.integer(edge[,2]), as.integer(nTips-1L), as.integer(nEdge), PACKAGE='inapplicable')
  numbered.edge <- .C('number_nodes', as.integer(reordered.edge[[1]]), as.integer(reordered.edge[[2]]), as.integer(root), as.integer(nEdge), PACKAGE='inapplicable')
  tree$edge <- matrix(c(numbered.edge[[1]], numbered.edge[[2]]), ncol=2)
  tree
}

#' @export
RootedTBR <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- SetOutgroup(tree, outgroup).")
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- 1 + (nTips <- dim(edge)[1] - tree$Nnode + 1L)
  root.children <- child[parent==root]
  size <- c(left.size <- sum(DoDescendants(parent, child, nTips, root.children[1L])) + 1L,
            (nTips * 2L - 1L) - left.size - 1L)
  if (min(size) == 1L) {
    outgroup <- tree$tip.label[root.children[size==1L]]
    tree <- TBR(tree)
    return(Root(tree, outgroup))
  } else {
    tip.label <- tree$tip.label
    moves <- (size-3L) * (size-2L) / 2
    subtree.root <- root.children[1L + (runif(1, min=0, max=sum(moves)) > moves[1L])]
    in.crown <- DoDescendants(parent, child, nTips, subtree.root)
    in.crown[subtree.root] <- TRUE
    crown.edges <- parent %in% which(in.crown)
    in.stump <- !in.crown
    in.stump[root] <- FALSE
    stump.edges <- parent %in% which(in.stump)
    stump <- KeepEdges(edge, tip.label, nTips, stump.edges) # faster than DropTip
    crown <- ExtractClade(tree, subtree.root) # faster than KeepEdges
    new.crown <- TBR(crown)
    new.crown$root.edge <- 1
    return (stump + new.crown)
  }
}

#' @export
NNI <- function (tree) {
  n      <- sample(tree$Nnode - 1L, 1L)
  edge   <- tree$edge
  parent <- edge[, 1L]
  child  <- edge[, 2L]
  k      <- min(parent) - 1L
  ind    <- which(child > k)[n]
  p1     <- parent[ind]
  p2     <- child[ind]
  ind1   <- which(parent == p1)
  ind1   <- ind1[ind1 != ind][1L]
  ind2   <- which(parent == p2)[sample(2L,1L)]
  tree$edge[c(ind1, ind2), 2L] <- child[c(ind2, ind1)]
  Renumber(ReorderPruning(tree))
}

#' @export
SPR <- function(tree) {
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  edge  <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  nEdge <- length(child)
  root  <- nTips + 1L
  if (nTips < 4L) stop ('must be >3 tips for SPR rearrangement!')
  pruning.candidates <- seq(nEdge + 1L)[-root]
  repeat {
    prune.node <- sample(pruning.candidates, 1)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- sample(graft.node, 1)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2L] <- edge[graft.edge, 2L]
    edge[graft.edge, 2L] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  reordered.edge <- .C('order_edges', as.integer(edge[,1]), as.integer(edge[,2]), as.integer(nTips-1L), as.integer(nEdge), PACKAGE='inapplicable')
  numbered.edge <- .C('number_nodes', as.integer(reordered.edge[[1]]), as.integer(reordered.edge[[2]]), as.integer(root), as.integer(nEdge), PACKAGE='inapplicable')
  tree$edge <- matrix(c(numbered.edge[[1]], numbered.edge[[2]]), ncol=2)
  tree
}

#' @name TBR
#' 
#'  Tree bisection and reconnection
#' @description This function performs a single random \acronym{TBR} iteration.
#' \usage{
#' TBR(tree, edge.to.break = NULL)
#' }
#' @param tree a fully resolved tree in \code{\link{phyDat}} format;
#' @param edge.to.break the index of an edge to bisect, generated randomly if not specified.
#' 
#' @details Branch lengths are not (yet) supported.
#' 
#' @return This function returns a tree in \code{phyDat} format that has undergone one \acronym{TBR} iteration.
#' @references The \acronym{TBR} algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' 
#' @author Martin R. Smith
#' 
#' @seealso RootedTBR useful when the position of the root node should be retained.
#' 
#' @examples{
#' tree <- rtree(20, br=NULL)
#' TBR(tree)
#' }
#' @export
TBR <- function(tree, edge.to.break=NULL) {
# Improvement targets: Root; ExtractClade; DropTip
  nTips <- tree$Nnode + 1
  if (nTips < 3) return (tree)
  tree.edge <- tree$edge
  tree.parent <- tree.edge[,1]
  tree.child <- tree.edge[,2]
  if (nTips == 3) return (Root(tree, sample(tree.child[tree.parent==max(tree.parent)], 1L)))
  all.nodes <- 1:(2*(nTips-1))
  root <- nTips + 1
  if (is.null(edge.to.break)) edge.to.break <- sample(2L:nrow(tree.edge), 1L) # Only include one root edge
  subtree.root <- tree.child[edge.to.break]
  #cat("\n - ", edge.to.break, subtree.root)
  stump <- if (subtree.root <= nTips) {
    DropTipNoSubtree(tree, subtree.root, root.edge=1)
  } else {
    in.crown <- DoDescendants(tree.parent, tree.child, nTips, subtree.root, just.tips=TRUE)
    DropTipNoSubtree (tree, which(in.crown), root.edge=1)
  }
  stump.len <- dim(stump$edge)[1]
  crown <- ExtractClade(tree, subtree.root) # ~ 2x faster than DropTip
  crown.edge <- crown$edge
  crown.len <- dim(crown.edge)[1]  
  if (crown.len > 1) {
    if (crown.len == 2) {
      rerooted.crown <- crown
    } else {
      crown.parent <- crown.edge[,1]
      crown.child <- crown.edge[,2]
      crown.nNode <- crown$Nnode
      crown.tips <- crown.nNode + 1L
      crown.root <- min(crown.parent)
      new.root.candidates <- crown.child[-1] # Include existing root once only
      new.root.node <- sample(new.root.candidates, 1L)
      if (new.root.node <= crown.tips) new.outgroup <- new.root.node else new.outgroup <- which(DoDescendants(crown.parent, crown.child, crown.tips, new.root.node, just.tips=TRUE))
      rerooted.crown <- Root(crown, new.outgroup)
    }
    rerooted.crown$root.edge <- 1L
    if (stump.len > 1) {
      bind.location <- sample(seq_len(stump.len), 1L)
      ret <- BindTree(stump, rerooted.crown, position=1, where=bind.location)
    } else {
      ret <- AddTip(rerooted.crown, crown.root, stump$tip.label)
    }
  } else {
    bind.location <- stump$edge[sample(2L:stump.len, 1L), 2L]
    ret <- AddTip(stump, bind.location, crown$tip.label)
  }
  Renumber(ret)
}

#' Generate random tree topology from dataset
#' 
#' @param dataset A dataset in \code{phyDat} format
#' 
#' @author Martin R. Smith 
#' @importFrom ape rtree
#' @export
RandomTree <- function (dataset) rtree(length(dataset), tip.label=names(dataset), br=NULL)
