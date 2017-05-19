InapplicableSectorial <- function (tree, data, maxit=100, 
    maxiter=500, k=5, trace=0, smallest.sector=4, largest.sector=1e+06, rearrangements="NNI", criterion=NULL, ...) {
  if (class(data) == 'phyDat') data <- PrepareData(data)
  if (class(data) != '*phyDat') stop("data must be a phyDat object, or the output of PrepareData(phyDat object).")
  if (is.null(tree)) stop("a starting tree must be provided")
  if (trace >= 0) cat('InapplicableSectorial search: optimizing sectors of', smallest.sector, 'to', floor(largest.sector), 'tips')
  
  sector.data <- function (X, tips) {
    at <- attributes(X)
    dec <- X[,tips]
    nBits <- floor(log2(max(X))) + 1L
    bin <- array(FALSE, dim=c(nrow(dec), ncol(dec), nBits))
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
    #attr(X, 'min.steps') <- apply(X, 1, function(x) min.steps(x, inapp.power2))
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
    node.lengths <- sapply(Descendants(tree, nodes), length) # (10x quicker than DoDescendants)
    candidate.nodes <- nodes[node.lengths >= smallest.sector & node.lengths <= largest.sector]
    if (trace >= 0) cat ("\n - Iteration", i, "- attempting sectorial search on node ")
    repeat {
      sector <- sample(candidate.nodes, 1)
      crown <- ExtractClade(tree, sector)
      crown.tips <- crown$tip.label
      sector.size <- length(crown.tips)
      cat(sector, 'size', sector.size, '...')
      crown.data <- sector.data(data, crown.tips)
      if (!is.null(crown.data)) break else cat('unsuitable (no data); trying')
      candidate.nodes <- candidate.nodes[-which(candidate.nodes==sector)]
      if (length(candidate.nodes == 0)) stop('No selectable sectors contain parsimony information! Either "largest.sector" is close to "smallest.sector" or your dataset is short of parsimony information.')
    } 
    if (trace >= 0) cat(' Sector OK.')
    crown <- root(AddTip(crown, 0, 'SECTOR_ROOT'), length(crown$tip.label) + 1, resolve.root=TRUE)
    initial.p <- MorphyParsimony(crown, crown.data, ...)
    attr(crown, 'pscore') <- initial.p
    if (trace >= 0) cat("\n - Running", rearrangements, "search on sector", sector)
    candidate <- TreeSearch(crown, crown.data, 'SECTOR_ROOT', method=rearrangements, criterion=criterion, trace=trace-1, maxiter=maxiter, ...)
    candidate.p <- attr(candidate, 'pscore')
    
    if((candidate.p + eps) < initial.p) {
      kmax <- kmax + 1
      stump <- DropTip(tree, Descendants(tree, sector)[[1]], subtree=TRUE)
      stump.edge <- 1:nrow(stump$edge)
      stump$root.edge <- 1
      crown <- DropTip(candidate, 'SECTOR_ROOT')
      tree <- CollapseSingles((BindTree(stump, crown, where=which(stump$tip.label==paste('[', sector.size, '_tips]', sep="")), position=0)))
      if (trace > 0) cat(' : improved local pscore, updated tree')
    } else if (trace > 0) cat (' : no improvement to local pscore')
    if (kmax == k) break()
  } # for
  if (trace >= 0)
    cat ("\nCompleted sectorial rearrangements.\n")
  attr(tree, 'pscore') <- NULL
  attr(tree, 'hits') <- NULL
  tree
}  # InapplicableSectorial

InapplicablePratchet <- function (tree, data, all=FALSE, outgroup=NULL, maxit=100, maxiter=5000, maxhits=40, k=10, trace=0, rearrangements="NNI", criterion=NULL, ...) {
  if (class(data) == 'phyDat') data <- PrepareData(data)
  if (class(data) != '*phyDat') stop("data must be a phyDat object, or the output of PrepareData(phyDat object).")
  eps <- 1e-08
  if (is.null(attr(tree, "pscore"))) attr(tree, "pscore") <- MorphyParsimony(tree, data, ...)
  best.pars <- attr(tree, "pscore")
  if (trace >= 0) cat("* Initial pscore:", best.pars)
  if (all) forest <- vector('list', maxiter)

  kmax <- 0
  for (i in 1:maxit) {
    if (trace >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bstree <- bootstrap.inapp(phy=tree, x=data, maxiter=maxiter, maxhits=maxhits, criterion=criterion, trace=trace-1, ...)
    
    if (trace >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new cannew candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- TreeSearch(bstree,    data, criterion=criterion, method='TBR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, data, criterion=criterion, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, data, criterion=criterion, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "TBR only") {            
      candidate <- TreeSearch(bstree,    data, criterion=criterion, method='TBR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR") {                  
      candidate <- TreeSearch(bstree,    data, criterion=criterion, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
      candidate <- TreeSearch(candidate, data, criterion=criterion, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else if (rearrangements == "SPR only") {             
      candidate <- TreeSearch(bstree,    data, criterion=criterion, method='SPR', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    } else {                                               
      candidate <- TreeSearch(bstree,    data, criterion=criterion, method='NNI', trace=trace, maxiter=maxiter, maxhits=maxhits, ...)
    }
    #if(class(result)=="phylo") m <- 1
    #else m = length(result)
    #if(m > 0) trees[2 : (1+m)] = result[1:m]
    #pscores <- sapply(trees, function(data) attr(data, "pscore"))
    cand.pars <- attr(candidate, 'pscore')
    if((cand.pars+eps) < best.pars) {
      if (all) {
        forest <- vector('list', maxiter)
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      }
      tree <- candidate
      best.pars <- cand.pars
      kmax <- 1
    } else {
      if (best.pars+eps > cand.pars) { # i.e. best == cand, allowing for floating point error
        kmax <- kmax + 1
        tree <- candidate
        if (all) forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      }
    }
    if (trace >= 0) cat("\n* Best pscore after", i, "/", maxit, "pratchet iterations:", best.pars, "( hit", kmax, "/", k, ")")
    if (kmax >= k) break()
  } # for
  if (trace >= 0)
    cat ("\nCompleted parsimony ratchet with pscore", best.pars, "\n")
    
  if (all) {
    forest <- forest[!vapply(forest, is.null, logical(1))]
    class(forest) <- 'multiPhylo'
    ret <- unique(forest)
    cat('Found', length(ret), 'unique MPTs.')
    if (is.null(outgroup)) warning('"outgroup" not specified, so some "unique" trees may have same topology but distinct roots.')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}

PratchetConsensus <- function (tree, data, maxit=5000, maxiter=500, maxhits=20, k=10, trace=0, rearrangements="NNI", criterion=NULL, nSearch=10, ...) {
  trees <- lapply(1:nSearch, function (x) InapplicablePratchet(tree, data, maxit, maxiter, maxhits, k=1, trace, rearrangements, criterion=criterion, ...))
  scores <- vapply(trees, function (x) attr(x, 'pscore'), double(1))
  trees <- unique(trees[scores == min(scores)])
  cat ("Found", length(trees), 'unique trees from ', nSearch, 'searches.')
  return (trees)
}

bootstrap.inapp <- function (phy, x, maxiter, maxhits, criterion=criterion, trace=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(x)
  weight <- at$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
  keep <- BS > 0
  ind <- which(keep)
  x <- x[ind, ]
  attr(x, 'weight') <- BS[ind]
  attr(x, 'min.steps') <- at$min.steps[keep]
  attr(x, 'unique.tokens') <- at$unique.tokens[keep]
  attr(x, 'nr') <- length(ind)
  attr(x, 'inapp.level') <- at$inapp.level
  attr(phy, 'pscore') <- NULL
  class(x) <- '*phyDat'
  res <- TreeSearch(phy, x, method='NNI', criterion=criterion, maxiter=maxiter, maxhits=maxhits, trace=trace-1, ...)
  attr(res, 'pscore') <- NULL
  attr(res, 'hits') <- NULL
  res
}

TreeSearch <- function (tree, data, method='NNI', maxiter=100, maxhits=20, forest.size=1, cluster=NULL, trace=1, criterion=NULL, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (!is.null(forest.size) && length(forest.size)) {
    if (forest.size > 1) {
      forest <- empty.forest <- vector('list', forest.size); forest[[1]] <- tree
    } else {
      forest.size <-1 
    }
  }
  if (is.null(attr(tree, 'pscore'))) attr(tree, 'pscore') <- MorphyParsimony(tree, data)
  best.pscore <- attr(tree, 'pscore')
  if (trace > 0) cat("\n  - Performing", method, "search.  Initial pscore:", best.pscore)
  rearrange.func <- switch(method, 'TBR' = TBR, 'SPR' = SPR, 'NNI' = QuickNNI)
  return.single <- !(forest.size > 1)
  
  for (iter in 1:maxiter) {
    trees <- RearrangeTree(tree, data, rearrange.func, min.score=best.pscore, return.single=return.single, iter=iter, cluster=cluster, criterion=criterion, trace=trace)
    iter.pscore <- attr(trees, 'pscore')
    if (length(forest.size) && forest.size > 1) {
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

SectorialSearch <- function (tree, data, concavity = NULL, rearrangements='NNI', maxiter=2000, cluster=NULL, trace=3) {
  best.score <- attr(tree, 'pscore')
  if (length(best.score) == 0) best.score <- MorphyParsimony(tree, data, ...)
  sect <- InapplicableSectorial(tree, data, cluster=cluster,
    trace=trace-1, maxit=30, maxiter=maxiter, maxhits=15, smallest.sector=6, 
    largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
  sect <- TreeSearch(sect, data, method='NNI', maxiter=maxiter, maxhits=30, cluster=cluster, trace=trace)
  sect <- TreeSearch(sect, data, method='TBR', maxiter=maxiter, maxhits=20, cluster=cluster, trace=trace)
  sect <- TreeSearch(sect, data, method='SPR', maxiter=maxiter, maxhits=50, cluster=cluster, trace=trace)
  sect <- TreeSearch(sect, data, method='NNI', maxiter=maxiter, maxhits=60, cluster=cluster, trace=trace)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (tree)
}
