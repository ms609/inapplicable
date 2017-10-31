[![Build Status](https://travis-ci.org/ms609/inapplicable.svg?branch=master)](https://travis-ci.org/ms609/inapplicable)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

# inapplicable
inapplicable is an R package, built on Phangorn, to allow the logically consistent handling of inapplicable data in parsimony analysis.
It incorporates modifications to phangorn that increase the rate of phylogenetic analysis, and adds support for TBR rearrangements.


You can install inapplicable into R thus:

```r
# Install the devtools package from CRAN, if necessary
if(!require(devtools)) install.packages("devtools")

# Install a working development version of phangorn, version 2.3 or greater
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
# Before you install the inapplicable package from GitHub:
if (!require(inapplicable)) devtools::install_github('ms609/inapplicable')

# Load the inapplicable package into R
library('inapplicable')
```

More details will be added here when details of the algorithm are published.

Here's an example of using the package to conduct tree search:

```r 
library(inapplicable)
# Here we'll use the pre-loaded Lobo matrix:
data(Lobo)
# You can load your own dataset using ape::read.nexus.data.

# Perform a simple search with a random starting tree
best <- TreeSearch(tree=TreeSearch::RandomTree(Lobo.phy, root='Cricocosmia'), dataset=Lobo.phy, Rearrange=TreeSearch::RootedNNI)

# Running a second search from this tree will probably see further improvements:
best <- TreeSearch(best, Lobo.phy)

# Using NNI might help to explore the region of treespace close to the local optimum:
best <- TreeSearch(best, Lobo.phy, Rearrange=TreeSearch::RootedNNI)

# SPR and TBR arrangements help to escape local optima and find better peaks 
# further away in tree space.  Using more hits (maxHits) and more iterations (maxIter)
# means we'll move closer to an optimal tree
best <- TreeSearch(best, Lobo.phy, maxHits=40, maxIter=100000, Rearrange=TreeSearch::RootedSPR, verbosity=2)
best <- TreeSearch(best, Lobo.phy, maxHits=40, maxIter=100000, Rearrange=TreeSearch::RootedTBR, verbosity=2)

# A more comprehensive search of tree space can be accomplished using the Parsimony Ratchet
# It might take a couple of minutes to run.
best <- Ratchet(best, Lobo.phy, verbosity=1)

# Let's view the tree:
plot(best)
```

# Reference

Details of the algorithm have been posted as a pre-print at 

BRAZEAU, M. D., GUILLERME, T. and SMITH, M. R. 2017. [Morphological phylogenetic analysis with inapplicable data](https://www.biorxiv.org/content/early/2017/10/26/209775). BioRχiv. doi:10.1101/209775

