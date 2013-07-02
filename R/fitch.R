parsimony.inapp <- function (tree, data) {
## Determine true parsimony score for a rooted tree
# ARGUMENTS: 
#   «tree», a rooted tree
#   «data», the output of optimize.data()
#   «target» (optional), a best score; abandon search if score will exceed «target»
# RETURN:
#   if «target» is exceeded: a non-integer larger than «target»
#   if «target» is not exceeded: the parsimony score of «tree» 
  return (fitch.inapp(tree, data)[[1]])
}

fitch.inapp <- function (tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try «data» <- prepare.data(«valid.phyDat.object»).')
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

fitch.combine <- function (a, b, nLevel) {
  shared <- a & b
  unions <- as.logical(colSums(shared))
  al<-a[-nLevel,]; bl<-b[-nLevel,]
  one.child.inapp.only <- !colSums(al) | !colSums(bl)
  intersects <- !unions
  
  ret <- matrix(FALSE, nrow(a), nCol <- ncol(a))
  ret[,unions] <- shared[,unions]
  ret[,intersects] <- a[,intersects] | b[,intersects]
  ret[nLevel, !one.child.inapp.only] <- FALSE
  
  col.score <- rep(TRUE, nCol)
  col.score[unions] <- FALSE
  col.score[one.child.inapp.only] <- FALSE
  
  return (list(ret, col.score))
}

fitch.combine.single <- function (a, b, nLevel) {
## data$levels[nLevel] must be the inapplicable token
# REQUIRE
#   «a», states that the left node could take - a vector comprising TRUE or FALSE statements
#        for tokens 1:(nLevel-1), with a[nLevel] representing the inapplicable token
#   «b», as «a», for the right node
#   «nLevel» <- length(«data»$level)
# RETURN
#   A list comprising: [[1]], state at node, in format of «a»; [[2]], unweighted addition to pscore
  shared <- a & b
  if (any(shared)) return (list(shared, 0))
  if (!any(a[-nLevel]) || !any(b[-nLevel])) return (list(a|b, 0))
  ret <- a|b
  ret[nLevel] <- FALSE
  return (list(ret, 1))
}