rearrange.tree <- function (tree, data, rearrange, concavity=NULL, return.single=TRUE, iter='<unknown>', cluster=NULL, trace=0) {
  if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    trees <- list(rearrange(tree))
    min.score <- parsimony.inapp(trees[[1]], data, concavity)
    best.trees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'pscore') <- parsimony.inapp(ret, cl.data, k); ret}, rearrange, tree, concavity)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    candidates <- clusterCall(cluster, rearrange, tree)
    scores <- vapply(candidates, parsimony.inapp, 1, data, concavity) # ~3x faster to do this in serial in r233.
    min.score <- min(scores)
    best.trees <- scores == min.score
    trees <- candidates[best.trees]
  }
  if (best.score < min.score) {
    if (trace > 3) cat("\n    . Iteration", iter, '- Min score', min.score, ">", best.score)
  } else if (best.score == min.score) {
    hits <- hits + sum(best.trees)
    if (trace > 2) cat("\n    - Iteration", iter, "- Best score", min.score, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (trace > 1) cat("\n    * Iteration", iter, "- New best score", min.score, "found on", hits, "trees")
  }
  if (return.single) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'pscore') <- min.score
  trees
}

rooted.nni <- function (tree) {
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
  tree <- renumber(phangorn:::reorderPruning(tree))  
}

rooted.spr <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- set.outgroup(tree, outgroup).")
  nTips <- length(tree$tip.label)
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- nTips + 1L # Assumes fully-resolved bifurcating tree
  root.children <- child[parent==root]
  left.nodes <- do.descendants(parent, child, nTips, root.children[1L])
  right.nodes <- !left.nodes
  left.nodes[root.children] <- right.nodes[root.children] <- right.nodes[root] <- FALSE
  size <- c(sum(left.nodes), sum(right.nodes))
  moves <- (size-2L) * (size-1L) / 2
  moves[size < 3] <- 0
  choose.right <- runif(1, min=0, max=sum(moves)) > moves[1]
  candidate.nodes <- if (choose.right) which(right.nodes) else which(left.nodes)
#  singletons <- candidate.nodes[1:2] < root #1:2 from Descendants
#  if (any(singletons)) candidate.nodes <- candidate.nodes[-which(!singletons)]
  prune.node <- sample(candidate.nodes, 1)
  prune.tips <- do.descendants(edge1, edge2, nTips, prune.node, just.tips=TRUE, include.ancestor=TRUE)
  pruning <- extract.clade.robust(tree, prune.node); pruning$root.edge <- 1
  tree$tip.label[prune.tips] <- 'PRUNED_TIP'
  affected.nodes <- c(prune.parent <- parent[child==prune.node], child[parent==prune.parent], which(do.descendants(edge1, edge2, nTips, prune.node)))
  candidate.nodes <- c(root.children[chosen.subtree], candidate.nodes[!candidate.nodes %in% affected.nodes])
  tree <- bind.tree(tree, pruning, where=graft.site <- sample(candidate.nodes, 1), position=1)
  tree <- drop.tip(tree, 'PRUNED_TIP')
  tree <- renumber(phangorn:::reorderPruning(tree))  
  tree
}

rooted.tbr <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- set.outgroup(tree, outgroup).")
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- 1 + (nTips <- dim(edge)[1] - tree$Nnode + 1L)
  root.children <- child[parent==root]
  size <- c(left.size <- sum(do.descendants(parent, child, nTips, root.children[1L])) + 1L,
            (nTips * 2L - 1L) - left.size - 1L)
  if (min(size) == 1L) {
    outgroup <- tree$tip.label[root.children[size==1L]]
    tree <- tbr(tree)
    return(root.robust(tree, outgroup))
  } else {
    moves <- (size-3L) * (size-2L) / 2
    subtree.root <- root.children[1L + (runif(1, min=0, max=sum(moves)) > moves[1L])]
    subtree.tips <- which(do.descendants(parent, child, nTips, subtree.root, just.tips=TRUE))

    stump <- drop.tip.fast(tree, subtree.tips, subtree=FALSE)
    stump$root.edge <- 1
    crown <- extract.clade.robust(tree, subtree.root)
    new.crown <- tbr(crown)
    new.crown$root.edge <- 1

    return (stump + new.crown)
  }
}

lb <- function () {nodelabels(); edgelabels(); tiplabels(adj=c(2, 0.5))}

tbr <- function(tree, edge.to.break=NULL) {
  tree.edge <- tree$edge
  tree.parent <- tree.edge[,1]
  tree.child <- tree.edge[,2]
  nTips <- tree$Nnode + 1
  all.nodes <- 1:(2*(nTips-1))
  root <- nTips + 1
  edge.to.avoid <- which(tree.parent==root) # This isn't ideal: it always avoids BOTH root edges.  But using one of the root edges can cause an error, and I can't tell whether there's a way to know (a) when this is the case; (b) which one.
  if (is.null(edge.to.break)) edge.to.break <- sample(setdiff(1:nrow(tree.edge), edge.to.avoid), 1)
  subtree.root <- tree.child[edge.to.break]
  subtree.tips <- descendants(tree, subtree.root, TRUE)
  #unrooted.tree <- unroot(tree)
  #tree.edge   <- unrooted.tree$edge
  #tree.parent <- tree.edge[,1]
  #tree.child  <- tree.edge[,2]
  #nTips <- length(tree.child) - unrooted.tree$Nnode + 1
  #if (is.null(edge.to.break)) edge.to.break <- sample(seq_along(tree.parent), 1)
  #subtree.root <- tree.child[edge.to.break]
  #subtree.tips <- descendants(unrooted.tree, subtree.root, TRUE)
  if (!length(subtree.tips)) subtree.tips <- subtree.root
  stump <- drop.tip.fast(tree, subtree.tips, subtree=FALSE)
  stump$root.edge <- 1
  crown <- extract.clade.robust(tree, subtree.root)
  
  if (dim(crown$edge)[1] > 1) {
    crown.tips <- crown$Nnode + 1L
    new.root.location <- crown$edge[sample(2L:nrow(crown$edge), 1L), 2L]
    nrp <- "NEW_ROOT_PLACEHOLDER"

    tip.added <- add.tip(crown, where=new.root.location, nrp)
    rooted <- root.robust(tip.added, nrp)
    rooted.crown <- drop.tip.fast(root.robust(add.tip(crown, where=new.root.location, nrp), nrp), nrp)

    rooted.crown <- drop.tip.fast(root.robust(add.tip(crown, where=new.root.location, nrp), nrp), nrp)
    bind.location <- sample(1L:nrow(stump$edge), 1L)
    rooted.crown$root.edge <- 1L
    ret <- bind.tree.fast(stump, rooted.crown, position=1, where=bind.location) #### DELETED collapse.singles as outer func
  } else {
    bind.location <- stump$edge[sample(seq_len(nrow(stump$edge)), 1L), 2L]
    ret <- add.tip(stump, bind.location, crown$tip.label)
  }
  renumber(ret)
}

set.outgroup <- root.robust