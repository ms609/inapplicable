library('phangorn')
SigSut.data <- read.nexus.data('SigSut-fold.nex')
contrast <- matrix(c(
1, 0, 0,
0, 1, 0,
0, 0, 1,
1, 1, 0,
1, 1, 1), ncol=3, byrow=TRUE, dimnames=list(c(0:1, '-', '+', '?'), c(0, 1, '-')))
SigSut.phy <- phyDat(SigSut.data, type='USER', contrast=contrast)
save('SigSut.data', 'SigSut.phy', file='../../data/SigSut.RData')

SuttonEtAl.data <- read.nexus.data('SuttonEtAl.nex')
contrast <- matrix(c(
1, 0, 0, 0, 0,
0, 1, 0, 0, 0,
0, 0, 1, 0, 0,
0, 0, 0, 1, 0,
0, 0, 0, 0, 1,
1, 1, 0, 0, 0, #r
0, 1, 1, 0, 0, #s
0, 1, 1, 1, 0, #t
0, 0, 1, 1, 0, #u
1, 1, 1, 1, 0, #+
1, 1, 1, 1, 1), ncol=5, byrow=TRUE, dimnames=list(c(0:3, '-', 'r', 's', 't', 'u','+', '?'), c(0:3, '-')))
SuttonEtAl.phy <- phyDat(SuttonEtAl.data, type='USER', contrast=contrast)
save('SuttonEtAl.data', 'SuttonEtAl.phy', file='../../data/SuttonEtAl.RData')

  