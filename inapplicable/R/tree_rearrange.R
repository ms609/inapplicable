rearrange.tree <- function (tree, data, rearrange, return.single=TRUE, iter='«unknown»', cores=4, trace=0) {
## Perform one «method» rearrangement on «tree»
# ARGUMENTS:
#   «tree», a rooted bifurcating phylogenetic tree with the desired outgroup, with the attributes:
#           : pscore, tree's parsimony score
#           : hits, the number of times the best score has been hit in the calling function
#   «data», output from optimize.data
#   «rearrange», a function to produce rearrangements
#   «return.single», returns all trees if FALSE or a randomly selected tree if TRUE
#   «iter», number of the iteration that called this function, for screen output only
#   «cores», number of separate rearrangements to evaluate
# RETURN:
#   A the most parsimonious of the «cores» trees generated, with attributes 'hits' and 'pscore'
  if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  candidates <- mclapply(1:cores, function (i) {rearrange(tree)})
  scores <- as.integer(mclapply(candidates, function(cand) {parsimony.inapp(cand, data)}))
  min.score <- min(scores)
  best.trees <- scores == min.score
  trees <- candidates[best.trees]
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
  ## abridged from phangorn::treeManipulation.R::nnin
  # ARGUMENTS:
  #   «tree», a rooted phylogenetic tree with the desired outgroup, perhaps generated using set.outgroup(tree, outgroup).
  # RETURN: 
  #   the input tree modified by a random NNI iteration, retaining the position of the root
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
    warning("Edge lengths have been deleted")  # I don't need this information so won't bother retaining it.
    tree$edge.length <- NULL
    #nnin has a line akin to: tree1$edge.length[c(ind1, ind2)] <- tree$edge.length[c(ind2, ind1)]
  }
  
  tree <- renumber(phangorn:::reorderPruning(tree))  
}

rooted.spr <- function(tree) {
## Returns a tree modified by a random SPR iteration, retaining the position of the root
# ARGUMENTS:
#   «tree», a rooted phylogenetic tree (phyDat format)
# RETURN:
#   a tree modified by a random SPR iteration, retaining the position of the root
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- set.outgroup(tree, outgroup).")
  nTips <- length(tree$tip.label)
  root <- nTips + 1L
  root.children <- Children(tree, root)
  size <- c(length(Descendants(tree, root.children[1], "all")),
            length(Descendants(tree, root.children[2], "all")))
  if (min(size) == 1) {
    spr.tree <- phangorn:::kSPR(tree, k=1)
    spr.tree <- renumber(spr.tree)
    tree <- root.robust(spr.tree, tree$tip.label[root.children[which(size==1)]])
    return (tree)
  }
  moves <- (size-2L) * (size-1L)
  chosen.subtree <- 1L + (runif(1, min=0, max=sum(moves)) >= moves[1])
  outgroup.root <- root.children[3L-chosen.subtree]
  outgroup.tips <- Descendants(tree, outgroup.root, type='tips')[[1]]
  subtree.root <- root.children[chosen.subtree]
  subtree.tips <- setdiff(1:nTips, outgroup.tips) 
  stump <- drop.tip(tree, subtree.tips, subtree=FALSE)
  stump$root.edge <- 1
  crown <- phangorn:::kSPR(extract.clade.robust(tree, subtree.root), k=1)
  crown$root.edge <- 1
  root.robust(stump + crown, tree$tip.label[outgroup.tips])
}

rooted.tbr <- function(tree) {
## Returns all trees produced by TBR at a random edge, retaining the position of the root
# ARGUMENTS:
#   «tree», a rooted and resolved phylogenetic tree (phyDat format, without branch length information)
# RETURN:
#   a tree topology produced by one TBR rearrangement
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- set.outgroup(tree, outgroup).")
  root <- 1 + (nTips <- dim(tree$edge)[1] - tree$Nnode + 1)
  root.children <- Children(tree, root)
  size <- c(length(Descendants(tree, root.children[1], "all")),
            length(Descendants(tree, root.children[2], "all")))
  if (min(size) == 1) {
    return(root.robust(tbr(tree), tree$tip.label[root.children[which(size==1)]]))
  } else {
    moves <- (size-2) * (size-1)
    subtree.root <- root.children[1 + (runif(1, min=0, max=sum(moves)) > moves[1])]
    subtree.tips <- Descendants(tree, subtree.root, type='tips')[[1]]

    stump <- drop.tip.fast(tree, subtree.tips, subtree=FALSE)
    stump$root.edge <- 1
    crown <- extract.clade.robust(tree, subtree.root)
    new.crown <- tbr(crown)
    new.crown$root.edge <- 1

    return (stump + new.crown)
  }
}

tbr <- function(tree, edge.to.break=NULL) {
## Perform a TBR operation on a tree
# ARGUMENTS:
#   «tree», a phylo object with nodes in postorder
#   «edge.to.break» (optional), an edge to bisect, generated randomly if not specified
# RETURN:
#   A tree topology generated after one TBR operation
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

set.outgroup <- function (tree, outgroup) {
  root.robust(tree, outgroup)
  #root(root(tree, setdiff(tree$tip, outgroup), resolve.root=TRUE), outgroup, resolve.root=TRUE)
}