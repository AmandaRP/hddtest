
<!-- README.md is generated from README.Rmd. Please edit that file -->
hddtest
=======

Functions for hypothesis testing of high dimensional discrete data. Currently functions for multinomially distributed data are available. Functions for multivariate binary data will be added in the future.

Installation
------------

You can install hddtest from github with:

``` r
# install.packages("devtools")
devtools::install_github("AmandaRP/hddtest")
```

Example
-------

Generate two multinomial count vectors and test whether they come from the same underlying distribution:

``` r
#data <- genMultinomialData(null_hyp=TRUE)
#multinom.test(x=data[[1]],y=data[[2]])
```
