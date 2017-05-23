dyn.load('../src/fitch.dll')
library('phangorn')
source('prepare_data.R')
source('fitch.R')
tree <- read.tree(text='((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));')
plot(tree); nodelabels(12:22); tiplabels(0:11)
data <- strsplit('
10
10
-0
-0
-0
00
00
-0
-0
-0
00
00', '')[[1]]
data <- matrix(data[data != '\n'], nrow=length(tree$tip), byrow=TRUE)
rownames(data) <- tree$tip
phy <- phyDat(data, levels=c(0:3, '-'), type='USER'); mp <- MorphyDat(phy)
result <- InapplicableFitch(tree, mp); names(result) <- c('pscore', 'pars', 'DAT', 'pvec')
result
nodes <- (length(tree$tip) + 1):(2 * length(tree$tip) - 1)
charToExplore <- 1
attach(result)
#plot(tree); nodelabels(DAT[charToExplore, nodes])
plot(tree); nodelabels(pvec[nodes])
tiplabels(data[, charToExplore])
nodelabels(nodes - 1)
detach(result)

dyn.unload('../src/fitch.dll')








#load('../data/SigSut.RData')

#phy <- SigSut.phy
#tree <- rtree(length(phy), tip.label = names(phy), br=NULL)
