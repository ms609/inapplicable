renumber <- function (tree) {
## Numbers the nodes and tips in a tree to conform with the phylo standards.
  tree <- reorder(tree, 'postorder')
  edge <- tree$edge
  nTip <- length(tree$tip.label)
  nNode = nTip - 1L
  edge1 <- edge[,1L]
  edge2 <- edge[,2L]
  NODES <- edge2 > nTip
  TIPS <- !NODES
  
  tip <- edge2[TIPS]
  name <- vector("character", length(tip))
  name[1:nTip] <- tree$tip.label[tip]
  tree$tip.label <- name
  edge2[TIPS] <- 1:nTip
  
  old.node.number <- unique(edge1)
  new.node.number <- (nTip + nNode):(nTip + 1L)
  edge2[NODES] <- new.node.number[match(edge2[NODES], old.node.number)]
  nodeseq <- (1L:nNode) * 2L
  edge1[c(nodeseq, nodeseq-1L)] <- new.node.number
  
  tree$edge[,1] <- edge1; tree$edge[,2] <- edge2
  reorder(tree)
}

single.taxon.tree <- function (label) {
  res <- list(edge=matrix(c(2L,1L), 1, 2), tip.label=label, Nnode=1L)
  class(res) <- 'phylo'
  res
}

extract.clade.robust <- function (phy, node) {
  phy.tip.label <- phy$tip.label
  phy.edge <- phy$edge
  nTip <- length(phy.tip.label)
  if (node <= nTip) return(single.taxon.tree(phy.tip.label[node]))
  if (node == nTip + 1L) return(phy)
  nodes.to.keep <- descendants(phy, node)
  edges.to.keep <- phy.edge[,2L] %in% nodes.to.keep
  phy.edge <- phy.edge[edges.to.keep, ]

  phy.edge1 <- phy.edge[,1L]
  phy.edge2 <- phy.edge[,2L]
  TIPS <- phy.edge2 <= nTip
  tip <- phy.edge2[TIPS]
  name <- vector("character", length(tip))
  name[order(tip)] <- phy.tip.label[tip]
  phy$tip.label <- name
  new.nTip <- length(name)
  phy.edge2[TIPS] <- order(tip)
  
  ## renumber nodes:
  phy.edge2[!TIPS] <- (phy.edge2[!TIPS] - node) + new.nTip + 1L
  phy.edge1 <- (phy.edge1 - node) + new.nTip + 1L
  phy$Nnode <- dim(phy.edge)[1] - new.nTip + 1L
  phy.edge[,1] <- phy.edge1;  phy.edge[,2] <- phy.edge2
  phy$edge <- phy.edge
  
  phy
}
ecr <- extract.clade.robust

add.tip <- function (tree, where, label) {
  nTip <- length(tree$tip.label)
  nNode <- tree$Nnode
  ROOT <- nTip + 1L
  if (where < 1L) where <- ROOT
  new.tip.number <- nTip + 1L
  tree.edge <- tree$edge
  
  ## find the row of 'where' before renumbering
  if (where == ROOT) case <- 1 else {
      insertion.edge <- which(tree.edge[, 2] == where)
      case <- if (where <= nTip) 2 else 3
  }
  ## case = 1 -> y is bound on the root of x
  ## case = 2 -> y is bound on a tip of x
  ## case = 3 -> y is bound on a node of x

### because in all situations internal nodes need to be
### renumbered, they are changed to negatives first, and
### nodes eventually added will be numbered sequentially
  nodes <- tree.edge > nTip
  tree.edge[nodes] <- -(tree.edge[nodes] - nTip)  # -1, ..., -nTip
  next.node <- -nNode - 1L
  ROOT <- -1L # This may change later
  
  switch(case, { # case = 1 -> y is bound on the root of x
      tree.edge <- rbind(c(next.node, tree.edge[1]), tree.edge, c(next.node, new.tip.number))
      ROOT <- next.node
    }, { # case = 2 -> y is bound on a tip of x
      tree.edge[insertion.edge, 2] <- next.node
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, where), c(next.node, new.tip.number), tree.edge[-(1:insertion.edge), ])
    }, { # case = 3 -> y is bound on a node of x
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, tree.edge[insertion.edge, 2]), tree.edge[-(1:insertion.edge), ])
      tree.edge[insertion.edge, 2] <- next.node
      insertion.edge <- insertion.edge + 1L
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, new.tip.number), tree.edge[-(1:insertion.edge), ])
    }
  )
  tree$tip.label <- c(tree$tip.label, label)
  tree$Nnode <- nNode <- nNode + 1L
  
  ## renumber nodes:
  new.numbering <- integer(nNode)
  new.numbering[-ROOT] <- new.tip.number + 1L
  second.col.nodes <- tree.edge[, 2] < 0
  ## executed from right to left, so newNb is modified before x$edge:
  tree.edge[second.col.nodes, 2] <- new.numbering[-tree.edge[second.col.nodes, 2]] <- new.tip.number + 2:nNode
  tree.edge[, 1] <- new.numbering[-tree.edge[, 1]]

  tree$edge <- tree.edge
  tree
  
}

root.robust <- function (tree, outgroup) {
  if (class(tree) != 'phylo') stop ('"tree" must be of class "phylo"')
  if (class(outgroup) == 'character') outgroup <- which(tree$tip.label %in% outgroup)
  if (length(outgroup) < 1) stop ('"outgroup" not specified')
  if (!is.null(tree$edge.length)) {tree$edge.length <- NULL; warning('Edge lengths are not supported and have been dropped.')}
  nTips <- length(tree$tip.label)
  root <- nTips + 1L
  root.children <- Children(tree, root)
  outgroup <- sort(outgroup)
  left.side <- descendants(tree, root.children[1L], TRUE)
  right.side <- descendants(tree, root.children[2L], TRUE)
  if (!any(left.side)) left.side <- root.children[1L]
  if (!any(right.side)) right.side <- root.children[2L]
  if ((length(outgroup) == length(left.side) && outgroup == left.side)
  ||   length(outgroup) == length(right.side) && outgroup == right.side) return (tree)
  if (any(root.children %in% outgroup)) outgroup <- seq_along(tree$tip)[-outgroup] # outgroup straddles root; root on ingroup instead
  
  ancestry <- Ancestors(tree, outgroup)
  if (length(outgroup) > 1) {
    common.ancestors <- Reduce(intersect, ancestry)
    outgroup.root.node <- max(common.ancestors)
  } else {
    common.ancestors <- c(outgroup, ancestry)
    outgroup.root.node <- outgroup
  }
  build.order <- rev(c(outgroup.root.node, Siblings(tree, common.ancestors[-length(common.ancestors)])))
  clades <- mclapply(build.order, function(n) {extract.clade.robust(tree, n)})
  ret <- clades[[1L]]
  if (length(ret$tip.label) == 1) {
    clade2 <- clades[[2L]]
    if (length(clade2$tip.label) == 1) ret <- two.tip.tree(ret$tip.label, clade2$tip.label)
      else ret <- add.tip(clade2, 0, ret$tip.label)
    first.to.add <- 3L
  } else {
    first.to.add <- 2L
  }
  for (i in first.to.add:length(clades)) {
    oClade <- clades[[i]]
    if (length(oClade$tip.label) > 1L) {
      ret$root.edge <- 1L
      oClade$root.edge <- 1L
      ret <- collapse.singles(bind.tree.fast(oClade, ret, position=1L))
    } else ret <- add.tip(ret, 0, oClade$tip.label)
  }
  ret
}

descendants <- function (tree, node, ...) {
# ARGUMENTS:
#   "tree", a phydat object
#   "node", number of an internal node
#   "just.tips", should return value include all nodes or just tips?
# RETURN:
#   vector containing descendant nodes in numerical order
  nTip <- length(tree$tip.label)
  edge <- tree$edge
  edge1 <- edge[,1]
  edge2 <- edge[,2]
  return (which(do.descendants(edge1, edge2, nTip, node, ...)))
}

do.descendants <- function (edge1, edge2, nTip, node, just.tips = FALSE, include.ancestor = FALSE) {
  # ARGUMENTS:
  #   "edge1", parent nodes: from tree$edge[,1]
  #   "edge2", parent nodes: from tree$edge[,2]
  #   "node", number of an internal node
  #   "just.tips", should return value include all nodes or just tips?
  # RETURN:
  #   vector containing descendant nodes in numerical order
  is.descendant <- blank <- logical((nTip * 2) - 1)
#  if (include.ancestor) is.descendant[node] <- TRUE;
  node.children <- function (node, is.descendant) {
    nc <- edge2[edge1 %in% node]
    is.descendant[nc] <- TRUE
    if (length(nc)) is.descendant <- node.children(nc, is.descendant)
    is.descendant
  }
  is.descendant <- node.children(node, is.descendant)
  if (just.tips) return (is.descendant[1:nTip]) else return (is.descendant)
}

two.tip.tree <- function (tip1, tip2) read.tree(text=paste0('(', tip1, ',', tip2, ');'))

bind.tree.fast <- function(x, y, where = "root", position = 0, interactive = FALSE) {
## Copied from ape:::bind.tree; the only change is that I use (x|y).edge in place of (x|y)$edge.

    nx <- length(x$tip.label)
    mx <- x$Nnode
    ROOTx <- nx + 1L
    ny <- length(y$tip.label)
    my <- y$Nnode

    if (interactive) {
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        if (lastPP$type != "phylogram" || lastPP$direction != "rightwards")
            stop("you must plot tree 'x' as a 'rightward phylogram'")
        cat("Click where you want to graft tree 'y'...\n")
        xy <- locator(1)
        d <- abs(xy$y - lastPP$yy)
        d[lastPP$xx - xy$x < 0] <- Inf
        where <- which.min(d)
        position <- lastPP$xx[where] - xy$x
        if (position < 0) position <- 0
        cat("The following parameters are used:\n")
        cat("  where =", where, " position =", position, "\n")
    } else {
        if (where == 0 || where == "root") where <- ROOTx
        if (position < 0) position <- 0
        if (where > nx + mx)
            stop("argument 'where' out of range for tree 'x'")
    }

    ## check whether both trees have branch lengths:
    switch(is.null(x$edge.length) + is.null(y$edge.length) + 1L,
           wbl <- TRUE, {
               x$edge.length <- y$edge.length <- NULL
               wbl <- FALSE
               warning("one tree has no branch lengths, they have been ignored")
           },
           wbl <- FALSE)

    yHasNoRootEdge <- is.null(y$root.edge)
    xHasNoRootEdge <- is.null(x$root.edge)

    x.edge <- x$edge
    y.edge <- y$edge
    ## find the row of 'where' before renumbering
    if (where == ROOTx) case <- 1 else {
        i <- which(x.edge[, 2] == where)
        case <- if (where <= nx) 2 else 3
    }
    ## case = 1 -> y is bound on the root of x
    ## case = 2 -> y is bound on a tip of x
    ## case = 3 -> y is bound on a node of x

    ## check that 'position' is correct
    if (position && wbl) {
### New in ape 3.0-1: this makes possible binding 'y' below
### a node of 'x' thus creating a new node in 'x'
###        if (!wbl)
###            stop("'position' is non-null but trees have no branch lengths")
        if (case == 1) {
            if (xHasNoRootEdge)
                stop("tree 'x' has no root edge")
            if (position > x$root.edge)
                stop("'position' is larger than x's root edge")
        } else {
            if (x$edge.length[i] < position)
                stop("'position' is larger than the branch length")
        }
    }

    ## the special case of substituting two tips:
    if (case == 2 && ny == 1 && !position) {
        x$tip.label[x.edge[i, 2]] <- y$tip.label
        if (wbl)
            x$edge.length[i] <- x$edge.length[i] + y$edge.length
        return(x)
    }

    x <- reorder(x)
    y <- reorder(y)

### because in all situations internal nodes need to be
### renumbered, they are changed to negatives first, and
### nodes eventually added will be numbered sequentially

    nodes <- x.edge > nx
    x.edge[nodes] <- -(x.edge[nodes] - nx) # -1, ..., -mx
    nodes <- y.edge > ny
    y.edge[nodes] <- -(y.edge[nodes] - ny + mx) # -(mx+1), ..., -(mx+my)
    ROOT <- -1L # may change later
    next.node <- -(mx + my) - 1L

    ## renumber now the tips in y:
    new.nx <- if (where <= nx && !position) nx - 1L else nx
    y.edge[!nodes] <- y.edge[!nodes] + new.nx

    ## if 'y' as a root edge, use it:
    if (!yHasNoRootEdge) {
        y.edge <- rbind(c(0, y.edge[1]), y.edge)
        ##                ^ will be filled later
        next.node <- next.node - 1L
        if (wbl) y$edge.length <- c(y$root.edge, y$edge.length)
    }

    switch(case, { # case = 1
        if (position) {
            x$root.edge <- x$root.edge - position
            x.edge <- rbind(c(next.node, x.edge[1]), x.edge)
            ROOT <- next.node
            if (wbl) x$edge.length <- c(position, x$edge.length)
        }
        if (yHasNoRootEdge) {
            j <- which(y.edge[, 1] == y.edge[1])
            y.edge[j, 1] <- ROOT
        } else y.edge[1] <- ROOT
        x.edge <- rbind(x.edge, y.edge)
        if (wbl)
            x$edge.length <- c(x$edge.length, y$edge.length)
    }, { # case = 2
        if (position) {
            x.edge[i, 2] <- next.node
            x.edge <- rbind(x.edge[1:i, ], c(next.node, where), x.edge[-(1:i), ])
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- x.edge[i, 1]
            } else y.edge[1] <- x.edge[i, 1]
        } else {
            if (yHasNoRootEdge) x.edge[i, 2] <- y.edge[1]
            else {
                ## the root edge of y is fused with the terminal edge of x
                if (wbl) y$edge.length[1] <- y$edge.length[1] + x$edge.length[i]
                y.edge[1] <- x.edge[i, 1]
                ## delete i-th edge in x:
                x.edge <- x.edge[-i, ]
                if (wbl) x$edge.length <- x$edge.length[-i]
                i <- i - 1L
            }
            x$tip.label <- x$tip.label[-where]
            ## renumber the tips that need to:
            ii <- which(x.edge[, 2] > where & x.edge[, 2] <= nx)
            x.edge[ii, 2] <- x.edge[ii, 2] - 1L
        }
        x.edge <- rbind(x.edge[1:i, ], y.edge, x.edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    }, { # case = 3
        if (position) {
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- next.node
            } else y.edge[1] <- next.node
            x.edge <- rbind(x.edge[1:i, ], c(next.node, x.edge[i, 2]), x.edge[-(1:i), ])
            x.edge[i, 2] <- next.node
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
        } else {
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- x.edge[i, 2]
            } else y.edge[1] <- x.edge[i, 2]
        }
        x.edge <- rbind(x.edge[1:i, ], y.edge, x.edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    })

    x$tip.label <- c(x$tip.label, y$tip.label)

    if (is.null(x$node.label)) {
        if (!is.null(y$node.label))
            x$node.label <- c(rep(NA, mx), y$node.label)
    } else {
        x$node.label <-
            if (is.null(y$node.label)) c(x$node.label, rep(NA, my))
            else c(x$node.label, y$node.label)
    }

    n <- length(x$tip.label)
    x$Nnode <- dim(x.edge)[1] + 1L - n

    ## update the node labels before renumbering (this adds NA for
    ## the added nodes, and drops the label for those deleted)
    if (!is.null(x$node.label))
        x$node.label <- x$node.label[sort(-unique(x.edge[, 1]))]

    ## renumber nodes:
    newNb <- integer(x$Nnode)
    newNb[-ROOT] <- n + 1L
    sndcol <- x.edge[, 2] < 0
    ## executed from right to left, so newNb is modified before x.edge:
    x.edge[sndcol, 2] <- newNb[-x.edge[sndcol, 2]] <- n + 2:x$Nnode
    x.edge[, 1] <- newNb[-x.edge[, 1]]
    x$edge <- x.edge
    if (!is.null(x$node.label))
        x$node.label <- x$node.label[order(newNb[newNb > 0])]

    x
}

drop.tip.fast <- function(phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(phy), interactive = FALSE) {
# Copied from ape:::drop.dip; edited to avoid excessive calls to $
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (!length(tip)) return(phy)
    phy.edge <- phy$edge
      
    Ntip <- length(phy$tip.label)
    ## find the tips to drop:
    if (interactive) {
        cat("Left-click close to the tips you want to drop; right-click when finished...\n")
        xy <- locator()
        nToDrop <- length(xy$x)
        tip <- integer(nToDrop)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        for (i in 1:nToDrop) {
            d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
            tip[i] <- which.min(d)
        }
    } else {
        if (is.character(tip))
            tip <- which(phy$tip.label %in% tip)
    }
    if (!length(tip)) return(phy)
    if (any(tip > Ntip))
        warning("some tip numbers were higher than the number of tips")

    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }

    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy.edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "pruningwise")
        N <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                as.integer(Nedge), double(Ntip + Nnode),
                DUP = FALSE, PACKAGE = "ape")[[6]]
    }
    wbl <- !is.null(phy$edge.length)
    edge1 <- phy.edge[, 1] # local copies
    edge2 <- phy.edge[, 2] #
    keep <- !logical(Nedge)

    ## delete the terminal edges given by `tip':
    keep[match(tip, edge2)] <- FALSE

    if (trim.internal) {
        ints <- edge2 > Ntip
        ## delete the internal edges that do not have anymore
        ## descendants (ie, they are in the 2nd col of `edge' but
        ## not in the 1st one)
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) break
            keep[sel] <- FALSE
        }
        if (subtree) {
            ## keep the subtending edge(s):
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0) # will store the indices of the edges below the new root
                repeat {
                    i <- which(edge1 == NEWROOT & keep)
                    j <- c(i, j)
                    NEWROOT <- edge2[i]
                    degree <- tabulate(edge1[keep])
                    if (degree[NEWROOT] > 1) break
                }
                keep[j] <- FALSE
                if (length(j) > root.edge) j <- 1:root.edge
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge))
                    NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }

    if (!root.edge) phy$root.edge <- NULL

    ## drop the edges
    phy.edge <- phy.edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]

    ## find the new terminal edges (works whatever 'subtree' and 'trim.internal'):
    TERMS <- !(phy.edge[, 2] %in% phy.edge[, 1])

    ## get the old No. of the nodes and tips that become tips:
    oldNo.ofNewTips <- phy.edge[TERMS, 2]

    ## in case some tips are dropped but kept because of 'subtree = TRUE':
    if (subtree) {
        i <- which(tip %in% oldNo.ofNewTips)
        if (length(i)) {
            phy$tip.label[tip[i]] <- "[1_tip]"
            tip <- tip[-i]
        }
    }

    n <- length(oldNo.ofNewTips) # the new number of tips in the tree

    ## the tips may not be sorted in increasing order in the
    ## 2nd col of edge, so no need to reorder $tip.label
    phy.edge[TERMS, 2] <- rank(phy.edge[TERMS, 2])
    phy$tip.label <- phy$tip.label[-tip]

    ## make new tip labels if necessary:
    if (subtree || !trim.internal) {
        ## get the numbers of the nodes that become tips:
        node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
        new.tip.label <- if (subtree) {
            paste("[", N[node2tip], "_tips]", sep = "")
        } else {
            if (is.null(phy$node.label)) rep("NA", length(node2tip))
            else phy$node.label[node2tip - Ntip]
        }
#        if (!is.null(phy$node.label))
#            phy$node.label <- phy$node.label[-(node2tip - Ntip)]
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }

    phy$Nnode <- dim(phy.edge)[1] - n + 1L # update phy$Nnode

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format, same as in root()
    newNb <- integer(Ntip + Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy.edge[, 2] > n
    ## executed from right to left, so newNb is modified before phy.edge:
    phy.edge[sndcol, 2] <- newNb[phy.edge[sndcol, 2]] <-
        (n + 2):(n + phy$Nnode)
    phy.edge[, 1] <- newNb[phy.edge[, 1]]
    phy$edge <- phy.edge
    storage.mode(phy$edge) <- "integer"
    if (!is.null(phy$node.label)) # update node.label if needed
        phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
    collapse.singles(phy)
}