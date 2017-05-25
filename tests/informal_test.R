dyn.load('../src/fitch.dll')
library('phangorn')
source('prepare_data.R')
source('fitch.R')
library(devtools)
install_github("TGuillerme/Inapp")
library(Inapp)

tree <- read.tree(text='((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));')
mp <- morphyData <- StringToMorphy((data_string<-'1??--??--100'), tree$tip)
plot.states.matrix(apply.reconstruction(tree, data_string, passes = 4, 'NA'),counts=c(1,2), show.labels = c(1,2))
result <- InapplicableFitch(tree, mp); names(result) <- c('pscore', 'pars', 'DAT')
result
nodes <- (length(tree$tip) + 1):(2 * length(tree$tip) - 1)
charToExplore <- 1
attach(result)
tiplabels(data[, charToExplore, drop=FALSE])
detach(result)

dyn.unload('../src/fitch.dll')





#load('../data/SigSut.RData')

#phy <- SigSut.phy
#tree <- rtree(length(phy), tip.label = names(phy), br=NULL)
