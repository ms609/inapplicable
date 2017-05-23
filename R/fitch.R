CallMorphy <- function (n_char, n_taxa, desc, ancs, rawmatrix) {
  .Call('RMorphy', as.integer(n_char), as.integer(n_taxa), as.integer(desc - 1L),
                  as.integer(ancs - 1L), as.character(rawmatrix))[[1]]
}

InapplicableFitch <- function (tree, morphyData) {
  # Data
  if (class(morphyData) == 'phyDat') morphyData <- MorphyDat(morphyData)
  if (class(morphyData) != 'morphyDat') stop('Invalid data type ', class(morphyData), '; try InapplicableFitch(tree, data <- MorphyData(valid.phyDat.object)).')
  at <- attributes(morphyData)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
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
  
  ret <- .Call("MORPHYFITCH", t(morphyData[tipLabel, ]), as.integer(nChar), as.integer(nTip), 
               as.integer(parent), as.integer(child), as.integer(parentOf), as.integer(childOf), 
               as.double(weight), as.integer(inappLevel), as.integer(inappChars))#, PACKAGE='inapplicable')
  
  return (ret)
}


MorphyParsimony <- function (tree, data) {
  tipNames <- tree$tip.label
  if (class(data) == 'phyDat') {
    at <- attributes(data)
    tipNames <- tipNames[tipNames %in% at$names]
    weight <- at$weight
    n.char <- sum(weight) # At present, MorphyLib doesn't allow us to leverage this efficient 
                          # solution, so we'll have to be inefficient.
    contrast <- at$contrast
    dimContrast <- dim(contrast)
    nRowContrast <- dimContrast[1]
    contrast <- matrix(as.logical(contrast), nRowContrast, dimContrast[2])
    levels <- at$levels
    matrixAsString <- paste0(c(vapply(tipNames, function (name) paste0(vapply(
      rep(as.integer(data[[name]]), weight), function (x) {
        if (x == nRowContrast) return('?')
        tokens <- levels[contrast[x ,]]
        if (length(tokens) > 1) return (paste0(c('{', tokens, '}'), collapse=''))
        return (tokens)
    }, character(1)), collapse=''), character(1)), ';'), collapse='')
  } else if (class(data) == 'matrix') {
    tipNames <- tipNames[tipNames %in% rownames(data)]
    dim.dat <- dim(data)
    n.char <- dim.dat[2]
    nTip  <- dim.dat[1]
    matrixAsString <- paste0(c(unlist(data), ';'), collapse='')
    weight <- rep(1, n.char)
  }  else if (class(data) == 'list') {
    tipNames <- tipNames[tipNames %in% names(data)]
    dim.dat <- dim(dat.matrix)
    n.char <- dim.dat[1]
    nTip  <- dim.dat[2]
    dat.matrix <- vapply(dat, as.vector, dat[[1]])
    matrixAsString <- paste0(c(t(dat.matrix), ';'), collapse='') ##CHECK: we probably don't need the t()
    weight <- rep(1, n.char)
  } else {
    stop ("Unrecognized data format. Try phyDat format; see ?phyDat.")
  }
  nTip <- length(tipNames)
  if (nTip < 3) stop ("Sorry, couldn't find enough tips with available data.")
  
  edge <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]

  maxNode <- nTip * 2 - 1
  root.node <- nTip + 1
  dummy.root.node <- maxNode + 1
  
  if (maxNode != max(parent)) stop ("Tree must be binary")
  if (root.node != min(parent)) stop ("Root node miscalculated")

  preorder <- root.node:maxNode

  ancestor <- function (x) parent[child==x]
  descendant <- function (x) child[parent==x]
  ancestors <- as.integer(c(vapply(1:nTip, ancestor, double(1)), 0, vapply((nTip + 2):(nTip * 2 - 1), ancestor, double(1))))
  ancestors[root.node] <- dummy.root.node
  descendants <- as.integer(vapply(preorder, descendant, double(2))) # children of each node, a pair at a time, right-left, right-left

  return(CallMorphy(n.char, nTip, descendants, ancestors, matrixAsString))
}
