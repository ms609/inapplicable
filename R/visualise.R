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
  nTip <- length(tip.label)
  nNode <- nTip - 1
  maxNode <- nNode + nTip
  inapp <- at$inapp.level
  tips <- seq(nTip)
  nodes <- nTip + seq(nNode)
  parentof <- parent[match((nTip + 2L):maxNode, child )] # Exclude the root, which has no parent
  childof <- child [c(match(nodes, parent), length(parent) + 1L - match(nodes, rev(parent)))]
  
  down <- .Call("FITCHDOWN", char.dat[tip.label], as.integer(1), as.integer(parent), as.integer(child), as.integer(nEdge), as.double(1), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable') # Return: (1), pscore; (2), pars; (3), DAT; (4), pvec; (5), need_up
  up <- .Call("FITCHUP", as.integer(down[[3]]), as.integer(1), as.integer(parentof), as.integer(childof), as.integer(nNode), as.double(1), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable')
  
  plot.fun(tree)
  text(1,1,paste0('TS', paste(which(at$index == char.no), collapse=', '), ': downpass +', down[[2]], if (down[[5]]) paste0('; uppass +', up[[2]]) else paste0('; uppass skipped (+', up[[2]], ')'), '; total +', down[[2]] + up[[2]]), pos=4, cex=0.8)
    
  downpass.states <- down[[3]]
  down.scorers <- down[[4]]
  down.change <- sapply(nodes, function(n) {
      children <- child[parent==n]
      return (down.scorers[n] != down.scorers[children[1]] + down.scorers[children[2]])   
    })
  #nodelabels(down.scorers[nodes])
  tipcols = c('#fafafa', '#fafafa', '#fafabb', '#ffbbbb', '#bbffbb', '#bbbbff', '#bbbbff', '#bbffbb', '#ffbbbb', '#bbddff', '#ffbbdd')
  names(tipcols) <- c(NA, max(downpass.states[1,]), max(downpass.states[1,])-inapp, 2^(0:7))
  tipcols[as.character(inapp)] <- '#999999'
  tipcols <- rev(tipcols)
  bgcols <- tipcols[as.character(downpass.states[1,tips])]
  bgcols[is.na(bgcols)] <- '#ffffbb'
  tiplabels(possible.tokens(at$levels, downpass.states[1,tips]), adj=c(0.3,0.5), bg=bgcols, col='#000088', cex=0.85)
  nodelabels(possible.tokens(at$levels, downpass.states[1,nodes]), adj=rep(1.25,2), frame='no', bg=ifelse(down.change, 'white', 'yellow'), font=ifelse(down.change, 2, 1) , col=ifelse(down.change, '#cc3333', '#880000cc'), cex=ifelse(down.change,1,0.6))
  
  uppass.states <- up[[3]]
  up.scorers <- up[[4]]
  up.change <- as.logical(up.scorers[nodes])
  nodelabels(possible.tokens(at$levels, uppass.states[1,nodes]), adj=c(1.25,-0.75), frame='no', bg='white', font=ifelse(up.change, 2, 1), col=ifelse(up.change, '#008800', '#008800cc'), cex=ifelse(up.change, 1,0.6))
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
  apply(which.levels, 1, function(x) {if (all(x)) return ('?') else y <- x; y[lvls=='-']<-TRUE; if (nTokens > 4 && all(y)) return ('+') else return (output(lvls[x]))})
}
