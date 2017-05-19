CallMorphy <- function (n_char, n_taxa, desc, ancs, rawmatrix) {
  .Call('RMorphy', as.integer(n_char), as.integer(n_taxa), as.integer(desc - 1L),
                  as.integer(ancs - 1L), as.character(rawmatrix))[[1]]
}

MorphyParsimony <- function (tree, data) {
  tip.names <- tree$tip.label
  if (class(data) == 'phyDat') {
    at <- attributes(data)
    tip.names <- tip.names[tip.names %in% at$names]
    weight <- at$weight
    n.char <- sum(weight) # At present, MorphyLib doesn't allow us to leverage this efficient 
                          # solution, so we'll have to be inefficient.
    contrast <- at$contrast
    dim.contrast <- dim(contrast)
    nrow.contrast <- dim.contrast[1]
    contrast <- matrix(as.logical(contrast), nrow.contrast, dim.contrast[2])
    levels <- at$levels
    matrix_as_string <- paste0(c(vapply(tip.names, function (name) paste0(vapply(
      rep(as.integer(data[[name]]), weight), function (x) {
        if (x == nrow.contrast) return('?')
        tokens <- levels[contrast[x ,]]
        if (length(tokens) > 1) return (paste0(c('{', tokens, '}'), collapse=''))
        return (tokens)
    }, character(1)), collapse=''), character(1)), ';'), collapse='')
  } else if (class(data) == 'matrix') {
    tip.names <- tip.names[tip.names %in% rownames(data)]
    dim.dat <- dim(data)
    n.char <- dim.dat[2]
    n.tip  <- dim.dat[1]
    matrix_as_string <- paste0(c(unlist(data), ';'), collapse='')
    weight <- rep(1, n.char)
  }  else if (class(data) == 'list') {
    tip.names <- tip.names[tip.names %in% names(data)]
    dim.dat <- dim(dat.matrix)
    n.char <- dim.dat[1]
    n.tip  <- dim.dat[2]
    dat.matrix <- vapply(dat, as.vector, dat[[1]])
    matrix_as_string <- paste0(c(t(dat.matrix), ';'), collapse='') ##CHECK: we probably don't need the t()
    weight <- rep(1, n.char)
  } else {
    stop ("Unrecognized data format. Try phyDat format; see ?phyDat.")
  }
  n.tip <- length(tip.names)
  if (n.tip < 3) stop ("Sorry, couldn't find enough tips with available data.")
  
  edge <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]

  max.node <- n.tip * 2 - 1
  root.node <- n.tip + 1
  dummy.root.node <- max.node + 1
  
  if (max.node != max(parent)) stop ("Tree must be binary")
  if (root.node != min(parent)) stop ("Root node miscalculated")

  preorder <- root.node:max.node

  ancestor <- function (x) parent[child==x]
  descendant <- function (x) child[parent==x]
  ancestors <- as.integer(c(vapply(1:n.tip, ancestor, double(1)), 0, vapply((n.tip + 2):(n.tip * 2 - 1), ancestor, double(1))))
  ancestors[root.node] <- dummy.root.node
  descendants <- as.integer(vapply(preorder, descendant, double(2))) # children of each node, a pair at a time, right-left, right-left

  return(CallMorphy(n.char, n.tip, descendants, ancestors, matrix_as_string))
}
