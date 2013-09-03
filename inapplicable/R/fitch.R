parsimony.inapp <- function (tree, data, concavity = NULL) {
  if (is.null(concavity)) return (fitch.inapp(tree, data)[[1]]) 
  # Implied weights
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try data <- prepare.data(valid.phyDat.object).')
  fit <- fitch.inapp(tree, data)
  e <- fit[[2]]
  inapp.level <- attr(data, 'inapp.level')
  inapp.power2 <- log2(inapp.level) + 1
  fit3 <- fit[[3]]
  fit3.inapp <- fit3 == inapp.level
  inapp.present <- rowSums(fit3.inapp) > 3 # two tips and one node causes no problems.
  if (any(inapp.present)) {
    nTips <- length(tree$tip.label)
    nChar <- length(e)
    edge <- tree$edge
    parent <- edge[,1]
    child <- edge[,2]
    children <- function (node) {child[parent==node]}
    root.node <- nTips + 1L # Assumes fully-resolved bifurcating tree
  }
  for (i in which(inapp.present)) {
    this.line  <- fit3[i,]
    nEntries <- length(this.line)
    this.inapp <- fit3.inapp[i,]
    lists <- vector('list', nTips - 1)
    current.list <- 1
    lists[[current.list]] <- matrix(inapp.level, ncol=1, nrow=nTips)
    traverse.order <- root.node:nEntries
    for (node in traverse.order) {
      if (this.inapp[node] && any(lists[[current.list]] != inapp.level)) {
        current.list <- current.list + 1
        lists[[current.list]] <- matrix(inapp.level, ncol=1, nrow=nTips)
      }
      node.children <- children(node)
      child.tips <- node.children[node.children <= nTips]
      lists[[current.list]][child.tips,] <- data[i, child.tips]
    }
    i.min <- unlist(sapply(lists, function (x) {min.steps(x, inapp.power2)}))
    e[i] <- e[i] - sum(i.min)
  }
  e[!inapp.present] <- e[!inapp.present] - attr(data, 'min.steps')[!inapp.present]
  weighted.fit <- attr(data, 'weight') * e / (concavity + e) # Corresponds to 1 - f = e / (e + k).  f = k / (e + k)
  return (sum(weighted.fit))
}

fitch.inapp <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try data <- prepare.data(valid.phyDat.object).')
  nChar  <- attr(data, 'nr') # strictly, transformation series patterns; these'll be upweighted later
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

  return (list(ret[[1]], ret[[2]], ret[[3]]))
}