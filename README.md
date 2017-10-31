[![Build Status](https://travis-ci.org/ms609/inapplicable.svg?branch=master)](https://travis-ci.org/ms609/inapplicable)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

# inapplicable
inapplicable is an R package, built on Phangorn, to allow the logically consistent handling of inapplicable data in parsimony analysis.
It incorporates modifications to phangorn that increase the rate of phylogenetic analysis, and adds support for TBR rearrangements.


You can install inapplicable into R thus:

```r
# Install the devtools package from CRAN, if necessary
if(!require(devtools)) install.packages("devtools")

# Install the inapplicable package from github
devtools::install_github('ms609/inapplicable')

# Load the inapplicable package into R
library('inapplicable')
```

More details will be added here when details of the algorithm are published.

Here's an example of using the package to conduct tree search:

```r 
library(inapplicable)
data(Lobo)
best <- TreeSearch(RandomTree(Lobo.phy), phy <- Lobo.phy)
best <- TreeSearch(best, phy)
best <- TreeSearch(best, phy, method='TBR')
best <- TreeSearch(best, phy, maxhits=40, maxiter=100000, method='SPR', verbosity=2)
best <- TreeSearch(best, phy, maxhits=40, maxiter=100000, method='TBR', verbosity=2)
best <- Ratchet(best, phy, outgroup='Cricocosmia', verbosity=1)
plot(Root(best, 'Cricocosmia'))
```
