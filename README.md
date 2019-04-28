
<!-- README.md is generated from README.Rmd. Please edit that file -->
    ## Loading hddtest

    ## 
    ## Attaching package: 'testthat'

    ## The following objects are masked from 'package:magrittr':
    ## 
    ##     equals, is_less_than, not

hddtest
=======

Functions for two sample hypothesis testing of high dimensional discrete data, specifically multinomial and multivariate binary data.

Installation
------------

You can install hddtest from github with:

``` r
install.packages("devtools")
devtools::install_github("AmandaRP/hddtest", build_vignettes = TRUE)
library("hddtest")
```

Example
-------

Generate two multinomial count vectors and test whether they come from the same underlying distribution:

``` r
data <- genMultinomialData(null_hyp=FALSE, sample_size = 1)
multinom.test(x=data[[1]], y=data[[2]])
#> $statistic
#> [1] 7.019512
#> 
#> $pvalue
#> [1] 1.113221e-12
```

The last call can also be done using the `magrittr` forward pipe:

``` r
data %>% multinom.test
#> $statistic
#> [1] 7.019512
#> 
#> $pvalue
#> [1] 1.113221e-12
```

Available functions and datasets
--------------------------------

See help documentation on each of the following via `?functionname`

-   multinom.test
-   multinom.neighborhood.test
-   genMultinomialData
-   mvbinary.test
-   genMVBinaryData
-   twoNewsGroups

Vignette
--------

Read more about the multinomial neighborhood test:

``` r
vignette("multinomial-neighborhood-test-vignette")
```

References
----------

\[1\] Plunkett, A. & Park, J. (2018) Two-sample test for sparse high-dimensional multinomial distributions, TEST, doi.org/10.1007/s11749-018-0600-8

\[2\] Plunkett, A. & Park, J. (2017) Two-sample tests for sparse high-dimensional binary data, Communications in Statistics - Theory and Methods, 46:22, 11181-11193, DOI: 10.1080/03610926.2016.1260743
