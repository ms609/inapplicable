parsimony.inapp <- function (tree, data, inherit.ancestral = FALSE, target = NULL, criterion=NULL) {
  if (is.null(criterion)) return (fitch.inapp(tree, data, target=target, inherit.ancestral=inherit.ancestral)[[1]]) 
  return (tree.information(tree, data))
}

fitch.inapp <- function (tree, data, inherit.ancestral = FALSE, target = NULL) {
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
  
  return (ret[1:3])
}