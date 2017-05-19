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
  n.edge <- length(parent)
  max.node <- parent[1] #max(parent)
  n.tip <- length(tip.label)
  n.node <- max.node - n.tip
  inapp <- at$inapp.level  
  parentof <- parent[match((n.tip + 2L):max.node, child )]
  allNodes <- (n.tip + 1L):max.node
  childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  ret <- .Call("FITCHINAPP", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof), as.integer(n.edge), as.integer(n.node), as.double(weight), as.integer(max.node), as.integer(n.tip), as.integer(inapp), PACKAGE='inapplicable') # 
  
  return (ret)
}

FitchInfo <- function (tree, data) {
    # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try FitchInfo(tree, data <- prepare.data(valid.phyDat.object)).')
  at <- attributes(data)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[,1]
  child <- tree.edge[,2]
  tip.label <- tree$tip.label
  n.edge <- length(parent)
  max.node <- parent[1] #max(parent)
  n.tip <- length(tip.label)
  n.node <- max.node - n.tip
  inapp <- at$inapp.level  
  parentof <- parent[match((n.tip + 2L):max.node, child )]
  allNodes <- (n.tip + 1L):max.node
  childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  nLevel <- length(at$level)
  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
  inapp.level <- which(at$levels == "-")
  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
  
  fitch <- .Call("FITCHINAPP", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof), as.integer(n.edge), as.integer(n.node), as.double(weight), as.integer(max.node), as.integer(n.tip), as.integer(inapp), PACKAGE='inapplicable') # 
  
  steps <- fitch[[2]]
  nodes.inapp <- fitch[[5]]
  splits <- attr(data, 'split.sizes')
  info <- vapply(1:attr(data, 'nr'), function (TS) {
    inapp.nodes <- nodes.inapp[TS,allNodes] < 0
    if (any(inapp.nodes)) {
      nodecount <- vapply(1:n.node, function(node) {
        if (inapp.nodes[node]) do.descendants(parent, child, n.tip, node + n.tip, include.ancestor=FALSE) else logical(n.node + n.tip)
      }, logical(n.node + n.tip))
      for (i in order(-colSums(nodecount))) {
        nodecount[apply(as.matrix(nodecount[,nodecount[allNodes, i]]), 1, any), i] <- FALSE
      }
      nodecount <- nodecount[1:n.tip,]
      number.of.origins <- sum(colSums(nodecount) > 1)
      my.splits <- splits[,TS]
      my.splits <- my.splits[my.splits > 1]
      if (length(my.splits) < 2) return (1)
      n.inapplicables <- sum(nodes.inapp[TS,1:n.tip] < 0)
      # need to multiply by unrooted(n.inapplicables) #TODO
      vapply(1:number.of.origins, function(n.origins) (
        ways.to.add.next.tips(n.inapplicables,0,0,0, 0,0,0,0, 0,n.origins,0,my.splits[1], TRUE)
        * possibilities(n.inapplicables, my.splits[1], my.splits[-1], steps[TS], n.origins, 0, number.of.origins)
      ), double(1))
      
      
      number.of.trees <- unrooted(n.inapplicables <- sum(fitch[[5]][TS,1:n.tip] < 0)) * 
      ways.to.add.next.tip(n.inapplicables, 0, 0, steps[TS] + 1, 0, 0, number.of.origins, 0, my.splits[1])
    } else {
      return (ProportionOfTreesConsistent(splits[,TS], steps[TS]))
    }
  }, double(1))
  log2(prod(info^attr(data, 'weight'))) # prod then log is faster than log then sum
}

FitchInfoFast <- function (tree, data) {
    # Data
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') stop('Invalid data type; try FitchInfoFast(tree, data <- PrepareDataFitch(valid.phyDat.object)).')
  at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  info <- at$info.amounts
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[,1]
  child <- tree.edge[,2]
  tip.label <- tree$tip.label
  n.edge <- length(parent)
  max.node <- parent[1] #max(parent)
  n.tip <- length(tip.label)
  n.node <- max.node - n.tip
  inapp <- at$inapp.level
  parent.of <- parent[match((n.tip + 2L):max.node, child )]
  allNodes <- (n.tip + 1L):max.node
  child.of <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  fitch <- .Call("FITCH", data[, tree$tip.label], as.integer(n.char), 
        as.integer(parent), as.integer(child), as.integer(n.edge), 
        as.double(weight), as.integer(max.node), as.integer(n.tip), package='phangorn')
#  
#  nLevel <- length(at$level)
#  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
#  inapp.level <- which(at$levels == "-")
#  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
 
#  fitch <- .Call("FITCHINAPP", data[, tip.label], as.integer(n.char), as.integer(parent), as.integer(child), as.integer(parent.of), as.integer(child.of), as.integer(n.edge), as.integer(n.node), as.double(weight), as.integer(max.node), as.integer(n.tip), as.integer(inapp), PACKAGE='inapplicable') 

  steps <- fitch[[2]]
  # Return a negative value because algorithms assume that smaller numbers are better
  return(-sum(info[(steps - 1) * n.char + seq_len(n.char)] * weight))
}
