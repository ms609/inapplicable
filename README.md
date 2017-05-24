# inapplicable
inapplicable is an R package, built on Phangorn, to allow the logically consistent handling of inapplicable data in parsimony analysis.
It incorporates modifications to phangorn that increase the rate of phylogenetic analysis, and adds support for TBR rearrangements.


You can install inapplicable into R thus:

```r
# Install the devtools package from CRAN
install.packages('devtools')

# Install the inapplicable package from github
devtools::install_github('ms609/inapplicable')

# Load the inapplicable package into R
library('inapplicable')
```

More details will be added here when details of the algorithm are published.

Here's an example of using the package to conduct tree search:

```r 
library(inapplicable)
data(SigSut)
taxa <- names(SigSut.phy)
tree <- rtree(length(taxa), tip.label=taxa, br=NULL)
result <- InapplicableFitch(tree, SigSut.phy)
best <- TreeSearch(tree, SigSut.phy)
best <- TreeSearch(best, SigSut.phy)
best <- TreeSearch(best, SigSut.phy, method='TBR')
best <- TreeSearch(best, SigSut.phy, maxhits=40, maxiter=100000, method='SPR', trace=3)
best <- TreeSearch(best, SigSut.phy, maxhits=240, maxiter=100000, method='TBR', trace=3)
plot(Root(best, 'Lingula'))
```
