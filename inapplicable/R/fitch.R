parsimony.inapp <- function (tree, data) return (fitch.inapp(tree, data)[[1]])

fitch.inapp <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try data <- prepare.data(valid.phyDat.object).')

  nChar <- attr(data, 'nr') # strictly, transformation series patterns; these'll be upweighted later
  weight <- attr(data, 'weight')
  # Tree
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") tree <- reorder(tree, "postorder")
  node <- tree$edge[,1]
  edge <- tree$edge[,2]
  tip.label <- tree$tip.label
  maxNode <- max(node) # m
  nTip <- length(tip.label) # q
  inapp <- attr(data, 'inapp.level')
  nNode <- tree$Nnode
  
  ret <- .Call("FITCHI", data[, tip.label], as.integer(nChar), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(rep(inapp, nChar)), PACKAGE='inapplicable')

  return (list(ret[[1]], ret[[2]]))
}