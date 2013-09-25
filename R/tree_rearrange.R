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
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- nTips + 1L # Assumes fully-resolved bifurcating tree
  root.children <- child[parent==root]
  left.nodes <- do.descendants(parent, child, nTips, root.children[1L])
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
  moving.subnodes <- c(prune.node, which(do.descendants(parent, child, nTips, prune.node)))
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
  edge[,1] <- numbered.edge[[1]]
  edge[,2] <- numbered.edge[[2]]
  tree$edge <- edge
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
    tip.label <- tree$tip.label
    moves <- (size-3L) * (size-2L) / 2
    subtree.root <- root.children[1L + (runif(1, min=0, max=sum(moves)) > moves[1L])]
    in.crown <- do.descendants(parent, child, nTips, subtree.root)
    in.crown[subtree.root] <- TRUE
    crown.edges <- parent %in% which(in.crown)
    in.stump <- !in.crown
    in.stump[root] <- FALSE
    stump.edges <- parent %in% which(in.stump)
    stump <- keep.edges(edge, tip.label, nTips, stump.edges) # faster than drop.tip.fast
    crown <- extract.clade.robust(tree, subtree.root) # faster than keep.edges
    new.crown <- tbr(crown)
    new.crown$root.edge <- 1
    return (stump + new.crown)
  }
}

lb <- function () {nodelabels(); edgelabels(); tiplabels(adj=c(2, 0.5))}

tbr <- function(tree, edge.to.break=NULL) {
# Improvement targets: root.robust; extract.clade.robust; drop.tip.fast
  nTips <- tree$Nnode + 1
  if (nTips < 3) return (tree)
  tree.edge <- tree$edge
  tree.parent <- tree.edge[,1]
  tree.child <- tree.edge[,2]
  if (nTips == 3) return (root.robust(tree, sample(tree.child[tree.parent==max(tree.parent)], 1L)))
  all.nodes <- 1:(2*(nTips-1))
  root <- nTips + 1
  if (is.null(edge.to.break)) edge.to.break <- sample(2L:nrow(tree.edge), 1L) # Only include one root edge
  subtree.root <- tree.child[edge.to.break]
  #cat("\n - ", edge.to.break, subtree.root)
  stump <- if (subtree.root <= nTips) {
    drop.tip.no.subtree(tree, subtree.root, root.edge=1)
  } else {
    in.crown <- do.descendants(tree.parent, tree.child, nTips, subtree.root, just.tips=TRUE)
    drop.tip.no.subtree (tree, which(in.crown), root.edge=1)
  }
  stump.len <- dim(stump$edge)[1]
  crown <- extract.clade.robust(tree, subtree.root) # ~ 2x faster than drop.tip.fast
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
      if (new.root.node <= crown.tips) new.outgroup <- new.root.node else new.outgroup <- which(do.descendants(crown.parent, crown.child, crown.tips, new.root.node, just.tips=TRUE))
      rerooted.crown <- root.robust(crown, new.outgroup)
    }
    rerooted.crown$root.edge <- 1L
    if (stump.len > 1) {
      bind.location <- sample(seq_len(stump.len), 1L)
      ret <- bind.tree.fast(stump, rerooted.crown, position=1, where=bind.location)
    } else {
      ret <- add.tip(rerooted.crown, crown.root, stump$tip.label)
    }
  } else {
    bind.location <- stump$edge[sample(2L:stump.len, 1L), 2L]
    ret <- add.tip(stump, bind.location, crown$tip.label)
  }
  renumber(ret)
}

set.outgroup <- root.robust