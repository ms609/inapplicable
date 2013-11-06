#visualize.char(best, data, i <- 7, colplot); text(2,2,paste('i =', i))
#visualize.char(best, data, i <- i+1, colplot); text(2,2,paste('i =', i))

visualize.character <- visualise.character <- visualize.char <- visualise.char <- function (tree, data, char.no, plot.fun = plot) {
  plot.fun(tree)
  fitch.out <- fitch.inapp(tree, data)
  optimized.states <- fitch.out[[3]]
  dat.at <- attributes(data)
  char.index <- dat.at$index[char.no]
  is.change <- fitch.switch(tree, optimized.states[char.index,], dat.at$inapp.level)
  tiplabels(possible.tokens(dat.at$levels, optimized.states[char.index,seq_along(tree$tip.label)]), adj=c(1,1), bg='white', frame='no')
  nodelabels(possible.tokens(dat.at$levels, optimized.states[char.index,length(tree$tip.label) + seq(tree$Nnode)]), adj=rep(1.25,2), frame='no', bg='white', font=ifelse(is.change != 1, 2, 1), col=is.change)
  text(1,1,paste0('+', fitch.out[[2]][char.no], ' of which ', fitch.out[[4]][char.no], ' on uppass'), pos=4, cex=0.8)
}

fitch.switch <- function (tree, states, inapp) {
  binary <- as.binary(states)
  inapp <- log2(inapp) + 1L
  e <- tree$edge
  parent <- e[,1]
  child <- e[,2]
  tips <- seq_along(tree$tip.label)
  nodes <- (length(tree$tip.label) + 1):length(states)
  sapply(nodes, function(n) {
    children <- child[parent==n]
    #if (isTRUE(all.equal(binary[children[1],] | binary[children[2],], as.logical(binary[n,])))) return (TRUE)
    c1 <- binary[children[1],]
    c2 <- binary[children[2],]
    if (isTRUE(all.equal(c1, c2))) return (1)
    if (!any(c1 & c2)) return ('red')
    c1[inapp] <- FALSE
    c2[inapp] <- FALSE
    if ((any(c1 | c2)) && !any(c1 & c2)) return ('#33aa33')
    return (1)
  })
}
#visualize.char(best, data, i <- 7, colplot); text(2,2,paste('i =', i))

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
