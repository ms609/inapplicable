dyn.load('../src/fitch.dll')
library('phangorn')
source('prepare_data.R')
source('fitch.R')
tree <- read.tree(text='((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));')
plot(tree); nodelabels(12:22); tiplabels(0:11)
data <- strsplit('
11
11
-1
-?
-?
01
00
-0
-0
-?
01
01', '')[[1]]
data <- matrix(data[data != '\n'], nrow=length(tree$tip), byrow=TRUE)
rownames(data) <- tree$tip
phy <- phyDat(data, levels=c(0:3, '-'), type='USER')
result <- InapplicableFitch(tree, phy); names(result) <- c('pscore', 'pars', 'DAT', 'pvec', 'APPL')
#result #[c(1, 2, 4)]
charToExplore <- 1
attach(result)
plot(tree); nodelabels(DAT[charToExplore, (length(tree$tip) + 1):ncol(DAT)])
tiplabels(data[, charToExplore])
detach(result)

dyn.unload('../src/fitch.dll')








#load('../data/SigSut.RData')

#phy <- SigSut.phy
#tree <- rtree(length(phy), tip.label = names(phy), br=NULL)
