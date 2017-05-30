library(inapplicable); library(ape)
data(SigSut)
mDat <- morphy <- dataset  <- MorphyDat(phy <- SigSut.phy)
tree <- rtree(nrow(mDat), tip.label=rownames(mDat), br=NULL)
ReorderPruning <- inapplicable:::ReorderPruning
Renumber <- inapplicable:::Renumber
BSI <- inapplicable:::BootstrapInapp
#tree <- read.tree(text='((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));')
#mp <- mDat <- morphyData <- StringToMorphy((data_string<-'1??--??--100'), tree$tip)


 maxiter=200; maxhits=20; forest.size=1; trace=99; method='SPR'; cluster=NULL; Rearrange=TBR

 ### B OOTSTRAP I NNAP
  at <- attributes(mDat)
  weight <- at$weight
  inappls <- at$inapp.chars
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
  keep <- BS > 0
  mDat <- mDat[, keep, drop=FALSE]
  attr(mDat, 'weight') <- BS[keep]
  attr(mDat, 'min.steps') <- at$min.steps[keep]
  attr(mDat, 'unique.tokens') <- at$unique.tokens[keep]
  attr(mDat, 'nr') <- sum(keep)
  attr(mDat, 'inapp.level') <- at$inapp.level
  attr(mDat, 'inapp.chars') <- sum(BS[1:inappls])
  attr(tree, 'pscore') <- NULL
  class(mDat) <- 'morphyDat'
  ### TREE SEARCH
    tree <- tree; dataset <- mDat
    tree$edge.length <- NULL # Edge lengths are not supported
    attr(tree, 'hits') <- 1
    attr(tree, 'pscore')  ##<- InapplicableFitch(tree, dataset)
    #### INAPPLICABLE FITCH
      morphyData <- dataset
      
      if (class(morphyData) == 'phyDat') morphyData <- MorphyDat(morphyData)
      if (class(morphyData) != 'morphyDat') stop('Invalid data type ', class(morphyData), '; try InapplicableFitch(tree, data <- MorphyData(valid.phyDat.object)).')
      treeOrder <- attr(tree, 'order')
      if (is.null(treeOrder) || treeOrder == "cladewise") tree <- reorder(tree, "postorder")
      at <- attributes(morphyData)
      nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
      weight <- at$weight
      tree.edge <- tree$edge
      parent <- tree.edge[,1]
      child <- tree.edge[,2]
      tipLabel <- tree$tip.label
      maxNode <- parent[1] #max(parent)
      nTip <- length(tipLabel)
      inappLevel <- at$inapp.level  
      inappChars <- at$inapp.chars
      parentOf <- parent[match(1:maxNode, child )]
      parentOf[nTip + 1] <- nTip + 1 # Root node "is its own parent"
      allNodes <- (nTip + 1L):maxNode
      childOf <- child [c(match(allNodes, parent), length(parent) + 1L - match(allNodes, rev(parent)))]
      # was
      t(morphyData[tipLabel, ])
      charData <- matrix(as.integer(morphyData[tipLabel, ]), nrow=nChar, ncol=nTip, byrow=TRUE)
      list(charData, as.integer(nChar), 
                   as.integer(nTip), as.integer(parent), as.integer(child),
                   as.integer(parentOf), as.integer(childOf), as.double(weight), 
                   as.integer(inappLevel), as.integer(inappChars), PACKAGE='inapplicable')

                   
                   ret <- .Call("MORPHYFITCH", charData, as.integer(nChar), 
                   as.integer(nTip), as.integer(parent), as.integer(child),
                   as.integer(parentOf), as.integer(childOf), as.double(weight), 
                   as.integer(inappLevel), as.integer(inappChars), PACKAGE='inapplicable')
      
      if (length(detail) == 1) return (ret[[detail]])
      return (ret[detail])
    
    
    
    
    
    ### END INAPP FITHC
    
    best.pscore <- attr(tree, 'pscore')
    if (trace > 0) cat("\n  - Performing", method, "search.  Initial pscore:", best.pscore)
    rearrange.func <- switch(method, 'TBR' = TBR, 'SPR' = SPR, 'NNI' = NNI)
    return.single <- !(forest.size > 1)
    
    iter <- 1 ### for (iter in 1:maxiter) {
    ##trees <- RearrangeTree(tree, dataset, rearrange.func, min.score=best.pscore, return.single=return.single, iter=iter, cluster=cluster, criterion=criterion, trace=trace)
    
    ### REARRRANGE TREE
      if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
      if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
      if (is.null(cluster)) {
        trees <- list(re.tree<-Rearrange(tree))
        min.score <- InapplicableFitch(re.tree, dataset)
        best.trees <- c(TRUE)
      } else {
        #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'pscore') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
        #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
        candidates <- clusterCall(cluster, Rearrange, tree)
        scores <- vapply(candidates, InapplicableFitch, 1, dataset) # ~3x faster to do this in serial in r233.
        min.score <- min(scores)
        best.trees <- scores == min.score
        trees <- candidates[best.trees]
      }
      if (best.score < min.score) {
        if (trace > 3) cat("\n    . Iteration", iter, '- Min score', min.score, ">", best.score)
      } else if (best.score == min.score) {
        hits <- hits + sum(best.trees)
        if (trace > 2) cat("\n    - Iteration", iter, "- Best score", min.score, "hit", hits, "times")
      } else {
        hits <- sum(best.trees)
        if (trace > 1) cat("\n    * Iteration", iter, "- New best score", min.score, "found on", hits, "trees")
      }
      if (length(return.single) && return.single) trees <- sample(trees, 1)[[1]]
      attr(trees, 'hits') <- hits
      attr(trees, 'pscore') <- min.score
      trees
    
    
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
    ### /end iter
    if (trace > 0) cat("\n  - Final score", attr(tree, 'pscore'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
    if (forest.size > 1) {
      if (hits < forest.size) forest <- forest[-((hits+1):forest.size)]
      attr(forest, 'hits') <- hits
      attr(forest, 'pscore') <- best.pscore
      return (unique(forest))
    } else tree
  
  
  res <- TreeSearch(tree, mDat, method='NNI', criterion=criterion, maxiter=maxiter, maxhits=maxhits, trace=trace-1, ...)
  attr(res, 'pscore') <- NULL
  attr(res, 'hits') <- NULL
  res


  # rp <- ReorderPruning(tree)
  x <- tree

  edge <- x$edge
  parents <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parents[!match(parents, child, 0)][1])
  n_edge <- length(parents)
  max_edge <- max(edge)
  neworder <- .C("phangorn_reorder", parents, child, as.integer(n_edge), 
      as.integer(max_edge), integer(n_edge), as.integer(root - 1L), PACKAGE = "inapplicable")[[5]]
  x$edge <- edge[neworder, ]
  x$edge.length <- x$edge.length[neworder]
  attr(x, "order") <- "pruningwise"
  rp <- x


Renumber(rp)

NNI(x)
SPR(x)
TBR(x)
  
  
  








 

 NNI(tree)
 
 
 
 
 
 
 
 
# RearrangeTree(tree, dataset, rearrange.func, min.score=best.pscore, return.single=return.single, iter=iter, cluster=cluster, criterion=criterion, trace=trace)
 if (is.null(attr(tree, 'pscore'))) best.score <- 1e+07 else best.score <- attr(tree, 'pscore')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    trees <- list(re.tree<-Rearrange(tree))
    min.score <- InapplicableFitch(re.tree, dataset)
    best.trees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'pscore') <- InapplicableFitch(ret, cl.data, k); ret}, rearrange, tree, concavity)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    candidates <- clusterCall(cluster, rearrange, tree)
    scores <- vapply(candidates, InapplicableFitch, 1, dataset) # ~3x faster to do this in serial in r233.
    min.score <- min(scores)
    best.trees <- scores == min.score
    trees <- candidates[best.trees]
  }
  if (best.score < min.score) {
    if (trace > 3) cat("\n    . Iteration", iter, '- Min score', min.score, ">", best.score)
  } else if (best.score == min.score) {
    hits <- hits + sum(best.trees)
    if (trace > 2) cat("\n    - Iteration", iter, "- Best score", min.score, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (trace > 1) cat("\n    * Iteration", iter, "- New best score", min.score, "found on", hits, "trees")
  }
  if (length(return.single) && return.single) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'pscore') <- min.score
  trees   