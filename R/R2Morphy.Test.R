library('phangorn')
dyn.load('../src/RMorphy.dll')  

tree <- read.tree(text="(((a, b), c), (d, (e, f)));")
matrix_as_grid <- t(matrix(c('0', '0', '1', '1', '0', '0',
                             '0', '0', '-', '-', '-', '1',
                             #'0', '+', '?', '-', '1', '1', 
                             '0', '+', '?', '-', '1', '{01}', 
                             '0', '1', '-', '-', '-', '0'), nrow=4, ncol=6, byrow=T))
rownames(matrix_as_grid) <- tree$tip.label
source('R2Morphy.R')
MorphyParsimony(tree, matrix_as_grid)

# devtools::install_github('ms609/inapplicable')
library(inapplicable)
data(SigSut)
phy <- SigSut.phy
tip.names <- names(phy)[-1]
tree <- rtree(length(tip.names), tip.label=tip.names, br=NULL)
MorphyParsimony(tree, phy) 


dyn.unload('../src/RMorphy.dll')  


dyn.load("../src/RMorphyEx.dll")
source("../src/Morphyex.R")


tip.names <- tree$tip.label
at <- attributes(phy)
tip.names <- tip.names[tip.names %in% at$names]
weight <- at$weight
n.char <- sum(weight) # At present, MorphyLib doesn't allow us to leverage this efficient 
                      # solution, so we'll have to be inefficient.
contrast <- at$contrast
dim.contrast <- dim(contrast)
nrow.contrast <- dim.contrast[1]
contrast <- matrix(as.logical(contrast), nrow.contrast, dim.contrast[2])
levels <- at$levels
matrix_as_string <- paste0(c(vapply(tip.names, function (name) paste0(vapply(as.integer(rep(phy[[name]], weight)), function (x) {
    if (x == nrow.contrast) return('?')
    tokens <- levels[contrast[x ,]]
    #if (length(tokens) > 1) return (paste0(c('{', tokens, '}'), collapse=''))
    if (length(tokens) > 1) return ('?')
    return (tokens)
}, character(1)), collapse=''), character(1)), ';'), collapse='')
  
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


morphy <- mpl_new_Morphy()
if (mpl_init_Morphy(n.tip, n.char, morphy)) stop("Error in mpl_init_Morphy")
if (mpl_get_numtaxa(morphy) != n.tip) warning("Taxon count mismatch with mpl_get_numtaxa")
if (mpl_get_num_charac(morphy) != n.char) warning("Character count mismatch in mpl_get_num_charac")
#if (msg <- mpl_attach_symbols(matrix_as_string, morphy)) stop("mpl_attach_symbols Couldn't attach symbols: error", msg)
if (msg <- mpl_attach_rawdata(matrix_as_string, morphy)) stop("mpl_attach_rawdata Couldn't attach symbols: error", msg)
for (i in (1:n.char) - 1) if (msg <- mpl_set_parsim_t(i, tname=as.character("FITCH"), morphy)) stop("mpl_set_parsim_t error with char", i, msg)
if (mpl_get_symbols(morphy) != '+01') stop ("Symbols mismatch in mpl_get_symbols")

if (mpl_set_num_internal_nodes(n.internal, morphy)) stop ("mpl_set_num_internal_nodes failed")
if (mpl_get_num_internal_nodes(morphy) != n.internal) stop('mpl_get_num_internal_nodes mismatch')
if (err <- mpl_apply_tipdata(morphy)) stop("mpl_apply_tipdata failed", err)

score <- 0
for (node in postorder) {
  score <- score + mpl_first_down_recon(node - 1L, left_ids[node] - 1L, right_ids[node] - 1L, morphy)
  #cat("Node ", node, "left: ", left_ids[node], " right, ", right_ids[node], "\n")
}
score
if (mpl_update_lower_root(dummy.root.node - 1L, root.node - 1L, morphy)) stop ("mpl_update_lower_root Error")
if (mpl_update_lower_root(root.node - 1L, root.node - 1L, morphy)) stop ("mpl_update_lower_root Error")
for (node in preorder[-1]) {
  #cat("Node:", node, "left:", left_ids[node], " right:", right_ids[node], "anc:", ancestors[node], "\n")
  score <- score + mpl_first_up_recon(node - 1L,  left_ids[node] - 1L, right_ids[node] - 1L, ancestors[node] - 1L, morphy) 
}
score
for (tip in tips) {
  mpl_update_tip(tip - 1L, ancestors[tip] - 1L, morphy)
}
for (node in postorder) {
  score <- score + mpl_second_down_recon(node - 1L, left_ids[node] - 1L, right_ids[node] - 1L, morphy)
}
if (mpl_update_lower_root(dummy.root.node - 1L, root.node - 1L, morphy)) stop ("Second mpl_update_lower_root Error")
for (node in preorder[-1]) {
  score <- score +mpl_second_up_recon(node - 1L,  left_ids[node] - 1L, right_ids[node] - 1L, ancestors[node] - 1L, morphy) 
}
score

mpl_delete_Morphy(morphy)

