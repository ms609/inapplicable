visualize.character <- visualise.character <- visualize.char <- visualise.char <- function (tree, data, char.no, plot.fun = plot) {
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop('Invalid data type; try fitch.inapp(tree, data <- prepare.data(valid.phyDat.object)).')
  at <- attributes(data)
  char.dat <- data[char.no,]
  char.index <- at$index[char.no]
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
  allNodes <- (nTip + 1L):maxNode
  parentof <- parent[match((nTip + 2L):maxNode, child )]
  childof <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
  
  down <- .Call("FITCHDOWN", char.dat[tip.label], as.integer(1), as.integer(parent), as.integer(child), as.integer(nEdge), as.double(1), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # Return: (1), pscore; (2), pars; (3), DAT; (4), pvec; (5), need_up
  up <- .Call("FITCHUP", as.integer(down[[3]]), as.integer(1), as.integer(parentof), as.integer(childof), as.integer(nNode), as.double(1), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable')
  
  plot.fun(tree)
  
  downpass.states <- down[[3]]
  down.scorers <- down[[4]]
  down.no.change <- non.scorer(parent, child, nTip, down.scorers)
  tiplabels(possible.tokens(at$levels, downpass.states[1,seq_along(tree$tip.label)]), adj=c(0.3,0.5), bg='white', col='#000088', cex=0.75)
  nodelabels(possible.tokens(at$levels, downpass.states[1,nTip + seq(tree$Nnode)]), adj=rep(1.25,2), frame='no', bg='white', font=ifelse(down.no.change, 1,2) , col=ifelse(down.no.change, '#880000cc', '#cc3333'), cex=ifelse(down.no.change, 0.6,1))
  
  uppass.states <- up[[3]]
  uppass.scorers <- up[[4]]
  up.no.change <- non.scorer(parent, child, nTip, uppass.scorers)  
  nodelabels(possible.tokens(at$levels, uppass.states[1,nTip + seq(tree$Nnode)]), adj=c(1.25,-0.75), frame='no', bg='white', font=ifelse(up.no.change, 1,2) , col=ifelse(up.no.change, '#008800cc', '#008800'), cex=ifelse(up.no.change, 0.6,1))
  text(1,1,paste0('downpass +', down[[2]], '; uppass +', up[[2]], '; total +', down[[2]] + up[[2]]), pos=4, cex=0.8)
}

non.scorer <- function (parent, child, nTip, states) {
  tips <- 1:nTip
  nodes <- (nTip + 1):length(states)
  sapply(nodes, function(n) {
    children <- child[parent==n]
    return (states[n] == states[children[1]] ||states[n] == states[children[2]])   
  })
}

possible.tokens <- function (lvls, number) {
  nTokens <- length(lvls)
  nNumber <- length(number)
  output <- function (x) {paste0(x, collapse='')}
  if (nNumber == 1) {    
    if (number == 2^nTokens - 1) return('?')
    which.levels <- rep(FALSE, nTokens)
    binary <- as.binary(number)
    which.levels[seq_along(binary)] <- binary
    return (output(lvls[as.logical(which.levels)]))
  }
  which.levels <- matrix(FALSE, nNumber, nTokens)
  binary <- as.binary(number)
  which.levels[,seq_along(binary[1,])] <- as.logical(binary)
  apply(which.levels, 1, function(x) {if (all(x)) return ('?') else y <- x; y[lvls=='-']<-TRUE; if (all(y)) return ('+') else return (output(lvls[x]))})
}
