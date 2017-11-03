[![Build Status](https://travis-ci.org/ms609/inapplicable.svg?branch=master)](https://travis-ci.org/ms609/inapplicable)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/ms609/inapplicable/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/inapplicable)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/inapplicable)](https://cran.r-project.org/package=inapplicable)

# inapplicable
`inapplicable` is an R package that allows parsimony search on morphological datasets that contain
inapplicable data, following the algorithm proposed by Brazeau, Guillerme and Smith (2017).

You can install inapplicable into R thus:

```r
if(!require(inapplicable)) install.packages("inapplicable")
library(inapplicable)
```

If you're feeling brave, you can install the development version thus:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/inapplicable')

# Load the inapplicable package into R
library('inapplicable')
```

Details on how to use the package are provided in the 'Getting started' vignette.

# Reference

Details of the algorithm have been posted as a pre-print at 

BRAZEAU, M. D., GUILLERME, T. and SMITH, M. R. 2017. [Morphological phylogenetic analysis with inapplicable data](https://www.biorxiv.org/content/early/2017/10/26/209775). BioRχiv. doi:10.1101/209775

