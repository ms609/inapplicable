InapplicableParsimony <- function (tree, data, inherit.ancestral = FALSE, target = NULL, criterion=NULL) {
  if (is.null(criterion)) return (InapplicableFitch(tree, data)[[1]]) 
  return (FitchInfo(tree, data))
}

InapplicableFitch <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- PrepareData(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try InapplicableFitch(tree, data <- PrepareData(valid.phyDat.object)).')
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
