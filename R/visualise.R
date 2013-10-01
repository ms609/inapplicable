visualize.character <- visualise.character <- visualize.char <- visualise.char <- function (tree, data, char.no, plot.fun = plot) {
  plot.fun(tree)
  optimized.states <- fitch.inapp(tree, data)[[3]]
  dat.at <- attributes(data)
  char.index <- dat.at$index[char.no]
  is.change <- fitch.switch(tree, optimized.states[char.index,])
  tiplabels(possible.tokens(dat.at$levels, optimized.states[char.index,seq_along(tree$tip.label)]), adj=c(1,1), bg='white', frame='no')
  nodelabels(possible.tokens(dat.at$levels, optimized.states[char.index,length(tree$tip.label) + seq(tree$Nnode)]), adj=rep(1.25,2), frame='no', bg='white', font=ifelse(is.change, 2, 1), col=ifelse(is.change, 'red', '#000066'))
}

fitch.switch <- function (tree, states) {
  binary <- as.binary(states)
  e <- tree$edge
  parent <- e[,1]
  child <- e[,2]
  tips <- seq_along(tree$tip.label)
  nodes <- (length(tree$tip.label) + 1):length(states)
  sapply(nodes, function(n) {
    children <- child[parent==n]
    if (isTRUE(all.equal(binary[children[1],], binary[children[2],]))) return (FALSE)
    if (isTRUE(all.equal(binary[children[1],] | binary[children[2],], as.logical(binary[n,])))) return (TRUE)
    return (FALSE)
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
