parsimony.inapp <- function (tree, data, concavity = NULL, target = NULL) {
  return (fitch.inapp(tree, data, target=target)[[1]]) 
#  
#  if (is.null(concavity)) return (fitch.inapp(tree, data, target=target)[[1]]) 
#  # Implied weights
#  if (class(data) == 'phyDat') data <- prepare.data(data)
#  if (class(data) != '*phyDat') stop('Invalid data type; try data <- prepare.data(valid.phyDat.object).')
#  fit <- fitch.inapp(tree, data)
#  e <- fit[[2]]
#  inapp.level <- attr(data, 'inapp.level')
#  if (!is.null(inapp.level)) {
#    inapp.power2 <- log2(inapp.level) + 1
#    fit3 <- fit[[3]]  
#    fit3.inapp <- fit3 == inapp.level
#    inapp.present <- rowSums(fit3.inapp) > 3 # two tips and one node causes no problems.
#  } else inapp.present <- FALSE
#  min.step <- attr(data, 'min.steps')
#  if (any(inapp.present)) {
#    nTips <- length(tree$tip.label)
#    nChar <- length(e)
#    nEntries <- ncol(fit3)
#    edge <- tree$edge
#    parent <- edge[,1]
#    child <- edge[,2]
#    children <- function (node) {child[parent==node]}
#    root.node <- nTips + 1L # Assumes fully-resolved bifurcating tree
#    
#    for (i in which(inapp.present)) {
#      this.line <- fit3[i,]
#      this.inapp <- fit3.inapp[i,]
#      lists <- matrix(inapp.level, nTips - 2, nTips)
#      current.list <- 1
#      traverse.order <- root.node
#      start.new.list <- rep(FALSE, nEntries)
#      while (length(traverse.order)) {
#        node = traverse.order[1]; traverse.order <- traverse.order[-1]
#        if (start.new.list[node]) current.list <- current.list + 1
#        node.children <- children(node)
#        a.tip <- node.children <= nTips
#        child.tips <- node.children[a.tip]
#        lists[current.list,child.tips] <- this.line[child.tips]
#        child.nodes <- node.children[!a.tip]
#        if (this.inapp[node]) {
#          start.new.list[child.nodes] = TRUE
#          traverse.order <- c(traverse.order, child.nodes)
#        } else traverse.order <- c(child.nodes, traverse.order)
#      }
#      lists <- lists[apply(lists, 1, function(x) {any(x!=inapp.level)}),] # Probably redundant but helps debugging
#      lists #<- lists[apply(lists, 1, function(x) {any(x!=inapp.level)}),apply(lists, 2, function(x) {any(x!=inapp.level)})] # Probably redundant but helps debugging
#      #apply(lists, 1, unique)  # TODO delete
#      if (is.null(dim(lists))) {
#        min.step[i] <- sum(min.steps(lists, inapp.power2))
#      } else {
#      # unlist(apply(lists, 1, function (x) {min.steps(x, inapp.power2)})) # TODO delete when debugged
#        min.step[i] <- sum(unlist(apply(lists, 1, function (x) {min.steps(x, inapp.power2)})))
#      }
#    }
#  }
#  e <- e - min.step
#  weighted.fit <- attr(data, 'weight') * e / (concavity + e) # Corresponds to 1 - f = e / (e + k).  f = k / (e + k)
#  return (sum(weighted.fit))
}

fitch.inapp <- function (tree, data, target = NULL, inherit.ancestral = TRUE) {
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
  
  if (inherit.ancestral) {
    ret <- .Call("FITCHDOWNIA", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(nEdge), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # Return: (1), pscore; (2), pars; (3), DAT; (4), pvec; (5), need_up
    parentof <- parent[match((nTip + 2L):maxNode, child )]
    allNodes <- (nTip + 1L):maxNode
    childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
    ups <- .Call("FITCHUPIA", as.integer(ret[[3]]), as.integer(nChar), as.integer(parentof), as.integer(childof), as.integer(nNode), as.double(weight[need.uppass]), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable')
    ret[[1]] <- ret[[1]] + ups[[1]]
    ret[[2]][need.uppass] <- ret[[2]][need.uppass] + ups[[2]]
    ret[[3]][need.uppass] <- ups[[3]]
  } else {
    ret <- .Call("FITCHDOWN", data[, tip.label], as.integer(nChar), as.integer(parent), as.integer(child), as.integer(nEdge), as.double(weight), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # Return: (1), pscore; (2), pars; (3), DAT; (4), pvec; (5), need_up
    if (any(need.uppass <- as.logical(ret[[5]])) && (is.null(target) || ret[[1]] <= target)) {
      parentof <- parent[match((nTip + 2L):maxNode, child )]
      allNodes <- (nTip + 1L):maxNode
      childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
      ups <- .Call("FITCHUP", as.integer(ret[[3]][need.uppass,]), as.integer(sum(need.uppass)), as.integer(parentof), as.integer(childof), as.integer(nNode), as.double(weight[need.uppass]), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable')
      ret[[1]] <- ret[[1]] + ups[[1]]
      ret[[2]][need.uppass] <- ret[[2]][need.uppass] + ups[[2]]
      ret[[3]][need.uppass] <- ups[[3]]
    }
  }
  return (ret[1:3])
}