parsimony.inapp <- function (tree, data, inherit.ancestral = FALSE, target = NULL, criterion=NULL) {
  if (is.null(criterion)) return (fitch.inapp(tree, data)[[1]]) 
  return (fitch.info(tree, data))
}

fitch.inapp <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try fitch.inapp(tree, data <- prepare.data(valid.phyDat.object)).')
  at <- attributes(data)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[,1]
  child <- tree.edge[,2]
  tip.label <- tree$tip.label
  nEdge <- length(parent)
  maxNode <- parent[1] #max(parent)
  nTip <- length(tip.label)
  nNode <- maxNode - nTip
  inapp <- at$inapp.level  
  parentof <- parent[match((nTip + 2L):maxNode, child )]
  allNodes <- (nTip + 1L):maxNode
  childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  ret <- .Call("FITCHINAPP", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof), as.integer(nEdge), as.integer(nNode), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # 
  
  return (ret)
}

fitch.info <- function (tree, data) {
    # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try fitch.inapp(tree, data <- prepare.data(valid.phyDat.object)).')
  at <- attributes(data)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[,1]
  child <- tree.edge[,2]
  tip.label <- tree$tip.label
  nEdge <- length(parent)
  maxNode <- parent[1] #max(parent)
  nTip <- length(tip.label)
  nNode <- maxNode - nTip
  inapp <- at$inapp.level  
  parentof <- parent[match((nTip + 2L):maxNode, child )]
  allNodes <- (nTip + 1L):maxNode
  childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  nLevel <- length(at$level)
  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
  inapp.level <- which(at$levels == "-")
  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
  
  fitch <- .Call("FITCHINAPP", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof), as.integer(nEdge), as.integer(nNode), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # 
  
  steps <- fitch[[2]]
  nodes.inapp <- fitch[[5]]
  splits <- attr(data, 'split.sizes')
  info <- vapply(1:attr(data, 'nr'), function (TS) {
    inapp.nodes <- nodes.inapp[TS,allNodes] < 0
    if (any(inapp.nodes)) {
      nodecount <- vapply(1:nNode, function(node) {
        if (inapp.nodes[node]) do.descendants(parent, child, nTip, node + nTip, include.ancestor=FALSE) else logical(nNode + nTip)
      }, logical(nNode + nTip))
      for (i in order(-colSums(nodecount))) {
        nodecount[apply(as.matrix(nodecount[,nodecount[allNodes, i]]), 1, any), i] <- FALSE
      }
      nodecount <- nodecount[1:nTip,]
      number.of.origins <- sum(colSums(nodecount) > 1)
      my.splits <- splits[,TS]
      my.splits <- my.splits[my.splits > 1]
      if (length(my.splits) < 2) return (1)
      n.inapplicables <- sum(nodes.inapp[TS,1:nTip] < 0)
      # need to multiply by unrooted(n.inapplicables) #TODO
      vapply(1:number.of.origins, function(n.origins) {
        ways.to.add.next.tip(n.inapplicables, 0, 0, 0, 0, 0, 0, n.origins, 0, my.splits[1], TRUE)
        * possibilities(n.inapplicables, my.splits[1], my.splits[-1], steps[TS], n.origins, number.of.origins)
      }, double(1))
      
      
      number.of.trees <- unrooted(n.inapplicables <- sum(fitch[[5]][TS,1:nTip] < 0)) * 
      ways.to.add.next.tip(n.inapplicables, 0, 0, steps[TS] + 1, 0, 0, number.of.origins, 0, my.splits[1])
    } else {
      return (proportion.of.trees.consistent(splits[,TS], steps[TS]))
    }
  }, double(1))
  log2(prod(info^attr(data, 'weight'))) # prod then log is faster than log then sum
}

#fitch.info <- function (tree, data) {
#  if (class(data) == 'phyDat') data <- prepare.data(data)
#  if (class(data) != '*phyDat') stop('Invalid data type; try fitch.inapp(tree, data <- prepare.data(valid.phyDat.object)).')
#  at <- attributes(data)
#  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
#  weight <- at$weight
#  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
#  tree.edge <- tree$edge
#  parent <- tree.edge[,1]
#  child <- tree.edge[,2]
#  tip.label <- tree$tip.label
#  nEdge <- length(parent)
#  nTip <- length(tip.label)
#  nNode <- nTip - 1
#  maxNode <- nNode + nTip
#  inapp <- at$inapp.level
#  tips <- seq(nTip)
#  nodes <- nTip + seq(nNode)
#  parentof <- parent[match((nTip + 2L):maxNode, child )] # Exclude the root, which has no parent
#  allNodes <- (root <- nTip + 1L):maxNode
#  childEdge <- c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))
#  childof <- child[childEdge]
#  if (any(is.na(data[, tip.label]))) stop("Tree's tip labels could not all be found in data matrix")
#  ret <- .Call("FITCHTRANS", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof),  as.integer(childEdge), as.integer(nEdge), as.integer(nNode), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # 
#  #extras <- apply(ret[[4]], 1, function (x) {
#  #  unx <- unique(x)
#  #  unx <- unx[unx>0 & unx != inapp]
#  #  uniques <- vapply(unx, function (u) sum(x==u), integer(1))
#  #  sum(uniques - 1)
#  #})
#  #il <- vapply(1:nChar, function (i) {
#  #  information.loss(colSums(as.split(data[i, ], attr(data, 'inapp.level'))), extras[i])
#  #}, double(1))
#  #-sum(il)
#  bric <- apply(ret[[5]], 1, brooks.info.content)
#  sum(bric)
#}