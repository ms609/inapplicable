parsimony.inapp <- function (tree, data) return (fitch.inapp(tree, data)[[1]])

fitch.inapp <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try data <- prepare.data(valid.phyDat.object).')
  weight <- attr(data, 'weight')
  nChar <- attr(data, 'nr') # strictly, transformation series patterns; these'll be upweighted later
  # Tree
  node <- tree$edge[,1]
  edge <- tree$edge[,2]
  tip.label <- tree$tip.label
  nNode <- tree$Nnode
  nTip <- nNode + 1
  children <- sapply(1:nNode, function (x) {edge[node==nTip+x]})
  root.node.index <- min(node)
  # Workspace: matrix of possible states
  na.matrix <- matrix(NA, ncol=nChar, nrow=32)
  charlist <<- mclapply(1:(length(tip.label)+nTip-1), function (tl) {if (tl <= nTip) data[[tip.label[tl]]] else na.matrix} )
  local.pscore <<- rep(0, nChar)
  nLevel <- length(attr(data, 'levels'))
  left.child   <- c(rep(0,nTip), children[1,])
  right.child  <- c(rep(0,nTip), children[2,])
  # Downpass function
  fitch.downpass <- function (oNode) {
    oLeft <- left.child[oNode]
    oRight <- right.child[oNode]
    left.state <- charlist[[oLeft]]
    right.state <- charlist[[oRight]]
    if (is.na(left.state[1])) left.state <- fitch.downpass(oLeft)
    if (is.na(right.state[1])) right.state <- fitch.downpass(oRight)
    fitch.result <- fitch.combine(left.state, right.state, nLevel)
    
    extra.pscore <- fitch.result[[2]]
    local.pscore <<- local.pscore + extra.pscore
    
    charlist[[oNode]] <<- fitch.result[[1]]
    return(fitch.result[[1]])
  }
  fitch.downpass(root.node.index)
  pvec <- weight * local.pscore
  return (list(sum(pvec), pvec))
}

prepare.data.c <- function (data) {
  data <- phangorn:::prepareDataFitch(data)
  d <- attributes(data)
  data <- as.integer(data)
  attributes(data) <- d
  attr(data, 'inapp.level') <- 2^(which(attr(data, 'levels') == "-")-1)
  class(data) <- '*phyDat'
  data
}

fitch.inapp.c <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data.c(data)
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
  
  ret <- .Call("FITCHI", data[, tip.label], as.integer(nChar), as.integer(node), as.integer(edge), as.integer(length(edge)), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(rep(inapp, nChar)))

  return (list(ret[[1]], ret[[2]]))
}

fitch.combine <- function (a, b, inapplicable.token.index) {
  shared <- a & b
  unions <- as.logical(colSums(shared))
  al<-a[-inapplicable.token.index,]; bl<-b[-inapplicable.token.index,]
  one.child.inapp.only <- !colSums(al) | !colSums(bl)
  intersects <- !unions
  
  ret <- matrix(FALSE, nrow(a), nCol <- ncol(a))
  ret[,unions] <- shared[,unions]
  ret[,intersects] <- a[,intersects] | b[,intersects]
  ret[inapplicable.token.index, !one.child.inapp.only] <- FALSE
  
  col.score <- rep(TRUE, nCol)
  col.score[unions] <- FALSE
  col.score[one.child.inapp.only] <- FALSE
  
  return (list(ret, col.score))
}

fitch.combine.single <- function (a, b, inapplicable.token.index) {
## data$levels[inapplicable.token.index] must be the inapplicable token
# REQUIRE
#   "a", states that the left node could take - a vector comprising TRUE or FALSE statements
#        for tokens 1:(inapplicable.token.index-1), with a[inapplicable.token.index] representing the inapplicable token
#   "b", as "a", for the right node
#   "inapplicable.token.index" <- length("data"$level)
# RETURN
#   A list comprising: [[1]], state at node, in format of "a"; [[2]], unweighted addition to pscore
  shared <- a & b
  if (any(shared)) return (list(shared, 0))
  if (!any(a[-inapplicable.token.index]) || !any(b[-inapplicable.token.index])) return (list(a|b, 0))
  ret <- a|b
  ret[inapplicable.token.index] <- FALSE
  return (list(ret, 1))
}