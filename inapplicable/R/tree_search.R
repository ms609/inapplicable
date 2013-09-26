sectorial.inapp <- function (start.tree, data, outgroup=NULL, concavity=NULL, maxit=100, 
    maxiter=500, k=5, trace=0, smallest.sector=4, largest.sector=1e+06, rearrangements="NNI", ...) {
  if (class(data) == 'phyDat') data <- prepare.data(data)
  if (class(data) != '*phyDat') stop("data must be a phyDat object, or the output of prepare.data(phyDat object).")
  if (is.null(start.tree)) stop("a start.tree must be provided")
  tree <- start.tree
  if (trace >= 0) cat('Sectorial search: optimizing sectors of', smallest.sector, 'to', floor(largest.sector), 'tips')
  
  sector.data <- function (X, tips) {
    at <- attributes(X)
    dec <- X[,tips]
    nBits <- floor(log2(max(X))) + 1L
    bin <- array(FALSE, dim=c(nrow(dec), ncol(dec), nBits))  ## TODO compare with as.binary
    for (i in 0:(nBits-1)) {
      bin[, , nBits-i] <- as.logical(dec %% 2)
      dec <- (dec %/% 2)
    }
    state.union <- apply(bin, c(1,3), all)
    parsimony.informative <- !as.logical(rowSums(state.union))
    if (!any(parsimony.informative)) return (NULL)
    X <- X[parsimony.informative, tips]
    informative.chars <- sum(parsimony.informative)
    SECTOR_ROOT <- rep(2^nBits-1, informative.chars)
    X <- cbind(X, SECTOR_ROOT)
    attr(X, 'nr') <- informative.chars
    attr(X, 'inapp.level') <- at$inapp.level
    inapp.power2 <- log2(at$inapp.level) + 1
    attr(X, 'min.steps') <- apply(X, 1, function(x) min.steps(x, inapp.power2))
    attr(X, 'levels') <- at$levels
    attr(X, 'weight') <- at$weight[parsimony.informative]
    class(X) <- '*phyDat'
    X
  }
  
  eps <- 1e-08
  kmax <- 1
  for (i in 1:maxit) {
    edge1 <- tree$edge[,1]
    nodes <- unique(edge1)[-1]
    node.lengths <- sapply(Descendants(tree, nodes), length)
    candidate.nodes <- nodes[node.lengths >= smallest.sector & node.lengths <= largest.sector]
    if (trace >= 0) cat ("\n - Iteration", i, "- attempting sectorial search on node ")
    repeat {
      sector <- sample(candidate.nodes, 1)
      crown <- extract.clade.robust(tree, sector)
      crown.tips <- crown$tip.label
      sector.size <- length(crown.tips)
      cat(sector, 'size', sector.size, '...')
      crown.data <- sector.data(data, crown.tips)
      if (!is.null(crown.data)) break else cat('unsuitable (no data); trying')
      candidate.nodes <- candidate.nodes[-which(candidate.nodes==sector)]
      if (length(candidate.nodes == 0)) stop('No selectable sectors contain parsimony information! Either "largest.sector" is close to "smallest.sector" or your dataset is short of parsimony information.')
    } 
    if (trace >= 0) cat(' Sector OK.')
    crown <- root(add.tip(crown, 0, 'SECTOR_ROOT'), length(crown$tip.label) + 1, resolve.root=TRUE)
    initial.p <- parsimony.inapp(crown, crown.data, concavity)
    attr(crown, 'pscore') <- initial.p
    if (trace >= 0) cat("\n - Running", rearrangements, "search on sector", sector)
    candidate <- tree.search(crown, crown.data, 'SECTOR_ROOT', concavity, method=rearrangements, concavity=concavity, trace=trace-1, maxiter=maxiter, ...)
    candidate.p <- attr(candidate, 'pscore')
    
    if((candidate.p + eps) < initial.p) {
      kmax <- kmax + 1
      stump <- drop.tip.fast(tree, Descendants(tree, sector)[[1]], subtree=TRUE)
      stump.edge <- 1:nrow(stump$edge)
      stump$root.edge <- 1
      crown <- drop.tip.fast(candidate, 'SECTOR_ROOT')
      tree <- collapse.singles.fast((bind.tree.fast(stump, crown, where=which(stump$tip.label==paste('[', sector.size, '_tips]', sep="")), position=0)))
      if (trace > 0) cat(' : improved local pscore, updated tree')
    } else if (trace > 0) cat (' : no improvement to local pscore')
    if (kmax == k) break()
  } # for
  if (trace >= 0)
    cat ("\nCompleted sectorial rearrangements.\n")
  if (!is.null(outgroup)) tree <- set.outgroup(tree, outgroup)
  attr(tree, 'pscore') <- NULL
  attr(tree, 'hits') <- NULL
  tree
}  # sectorial.inapp

pratchet.inapp <- function (start.tree, data, outgroup=NULL, concavity=NULL, maxit=5000, maxiter=500, maxhits=20, k=10, trace=0, rearrangements="NNI", ...) {
  if (class(data) == 'phyDat') data <- prepare.data(data)
  tree <- start.tree; start.tree <- NULL
  if (class(data) != '*phyDat') stop("data must be a phyDat object, or the output of prepare.data(phyDat object).")
  eps <- 1e-08
  if (is.null(attr(tree, "pscore"))) attr(tree, "pscore") <- parsimony.inapp(tree, data, concavity)
  mp <- attr(tree, "pscore")
  if (trace >= 0) cat("* Initial pscore:", mp)

  kmax <- 1
  for (i in 1:maxit) {
    if (trace >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bstree <- bootstrap.inapp(phy=tree, x=data, outgroup=outgroup, concavity=concavity, maxiter=maxiter, trace=trace-1, ...) #15% of function time in r233
    
    if (trace >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- tree.search(bstree,    data, outgroup, concavity, method='TBR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- tree.search(candidate, data, outgroup, concavity, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- tree.search(candidate, data, outgroup, concavity, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "TBR only") {
      candidate <- tree.search(bstree,    data, outgroup, concavity, method='TBR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR") {
      candidate <- tree.search(bstree,    data, outgroup, concavity, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- tree.search(candidate, data, outgroup, concavity, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR only") {
      candidate <- tree.search(bstree,    data, outgroup, concavity, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else {                             
      candidate <- tree.search(bstree,    data, outgroup, concavity, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    }
    #if(class(result)=="phylo") m <- 1
    #else m = length(result)
    #if(m > 0) trees[2 : (1+m)] = result[1:m]
    #pscores <- sapply(trees, function(data) attr(data, "pscore"))
    mp1 <- attr(candidate, 'pscore')
    if((mp1+eps) < mp) {
      kmax <- 1
      tree <- candidate
      mp   <- mp1
    } else {
      if (mp+eps > mp1) kmax <- kmax + 1
      if ((sample(2,1) == 1) & (mp1 < (mp+eps))) tree <- candidate
    }
    if (trace >= 0) cat("\n* Best pscore after", i, "/", maxit, "pratchet iterations:", mp, "( hit", kmax, "/", k, ")")
    if (kmax == k) break()
  } # for
  if (trace >= 0)
    cat ("\nCompleted parsimony ratchet with pscore", mp, "\n")
    
  attr(tree, 'hits') <- NULL
  tree
}

bootstrap.inapp <- function (phy, x, outgroup, concavity, maxiter, trace=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(x)
  weight <- at$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE),length(weight)) 
  keep <- BS > 0
  ind <- which(keep)
  x <- x[ind,]
  attr(x, 'weight') <- BS[ind]
  attr(x, 'min.steps') <- at$min.steps[keep]
  attr(x, 'nr') <- length(ind)
  attr(x, 'inapp.level') <- at$inapp.level
  attr(phy, 'pscore') <- NULL
  class(x) <- '*phyDat'
  res <- tree.search(phy, x, outgroup, concavity, method='NNI', maxiter, trace=trace-1, ...)
  attr(res, 'pscore') <- NULL
  attr(res, 'hits') <- NULL
  res
}

tree.search <- function (start.tree, data, outgroup, concavity=NULL, method='NNI', maxiter=100, maxhits=20, forest.size=1, cluster=NULL, trace=1, ...) {
  start.tree$edge.length <- NULL # Edge lengths are not supported
  tree <- set.outgroup(start.tree, outgroup)
  attr(tree, 'hits') <- 1
  if (forest.size > 1) {forest <- empty.forest <- vector('list', forest.size); forest[[1]] <- tree}
  if (is.null(attr(tree, 'pscore'))) attr(tree, 'pscore') <- parsimony.inapp(tree, data, concavity)
  best.pscore <- attr(tree, 'pscore')
  if (trace > 0) cat("\n  - Performing", method, "search.  Initial pscore:", best.pscore)
  rearrange.func <- switch(method, 'TBR' = rooted.tbr, 'SPR' = rooted.spr, 'NNI' = rooted.nni)
  for (iter in 1:maxiter) {
    trees <- rearrange.tree(tree, data, rearrange.func, min.score=best.pscore, concavity=concavity, return.single=forest.size==1, iter=iter, cluster=cluster, trace=trace)
    iter.pscore <- attr(trees, 'pscore')
    if (forest.size > 1) {
      hits <- attr(trees, 'hits')
      if (iter.pscore == best.pscore) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'pscore') <- iter.pscore
        attr(tree, 'hits') <- hits
      } else if (iter.pscore < best.pscore) {
        best.pscore <- iter.pscore
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'pscore') <- iter.pscore
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iter.pscore <= best.pscore) {
        best.pscore <- iter.pscore
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxhits) break
  }
  if (trace > 0) cat("\n  - Final score", attr(tree, 'pscore'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
  if (forest.size > 1) {
    if (hits < forest.size) forest <- forest[-((hits+1):forest.size)]
    attr(forest, 'hits') <- hits
    attr(forest, 'pscore') <- best.pscore
    return (unique(forest))
  } else tree
}

sectorial.search <- function (start.tree, data, outgroup, concavity = NULL, rearrangements='NNI', maxiter=2000, trace=3, cluster=NULL) {
  best.score <- attr(start.tree, 'pscore')
  if (length(best.score) == 0) best.score <- parsimony.inapp(start.tree, data, concavity)
  if (length(outgroup) == 0) warning('"outgroup" parameter not specified')
  sect <- sectorial.inapp(start.tree, data, outgroup=outgroup, concavity=concavity, cluster=cluster,
    trace=trace-1, maxit=30, maxiter=maxiter, maxhits=15, smallest.sector=6, 
    largest.sector=length(start.tree$edge[,2L])*0.25, rearrangements=rearrangements)
  sect <- tree.search(sect, data, outgroup, method='NNI', concavity=concavity, maxiter=maxiter, maxhits=30, cluster=cluster, trace=trace)
  sect <- tree.search(sect, data, outgroup, method='TBR', concavity=concavity, maxiter=maxiter, maxhits=20, cluster=cluster, trace=trace)
  sect <- tree.search(sect, data, outgroup, method='SPR', concavity=concavity, maxiter=maxiter, maxhits=50, cluster=cluster, trace=trace)
  sect <- tree.search(sect, data, outgroup, method='NNI', concavity=concavity, maxiter=maxiter, maxhits=60, cluster=cluster, trace=trace)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (set.outgroup(start.tree, outgroup))
}
