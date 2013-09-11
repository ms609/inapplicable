
tree <- best; data <- dat; charno <- 33;
visualize.character <- visualise.character <- visualize.char <- visualise.char <- function (tree, data, charno, plotfun = plot) {
  plotfun(tree)
  optimized.states <- fitch.inapp(tree, data)[[3]]
  dat.at <- attributes(data)
  char.index <- dat.at$index[charno]
  tiplabels(possible.tokens(dat.at$levels, optimized.states[char.index,seq_along(tree$tip.label)]), adj=c(1,1), bg='white', frame='no', font=2)
  nodelabels(possible.tokens(dat.at$levels, optimized.states[char.index,length(tree$tip.label) + seq(tree$Nnode)]), adj=c(1.5,1.5), frame='no', bg='white', font=2)
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
