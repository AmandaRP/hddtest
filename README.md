
<!-- README.md is generated from README.Rmd. Please edit that file -->
    ## Loading hddtest

hddtest
=======

Functions for hypothesis testing of high dimensional discrete data. Currently functions for multinomial data, as described in \[1\], are available. Functions for multivariate binary data will be added in the future.

Installation
------------

You can install hddtest from github with:

``` r
install.packages("devtools")
devtools::install_github("AmandaRP/hddtest")
library("hddtest")
```

Example
-------

Generate two multinomial count vectors and test whether they come from the same underlying distribution:

``` r
data <- genMultinomialData(null_hyp=FALSE,sample_size = 1)
multinom.test(x=data[[1]],y=data[[2]])
#> $statistic
#> [1] 6.422941
#> 
#> $pvalue
#> [1] 6.683298e-11
```

The last call can also be done in a "tidy" fashion using the forward pipe from the `magrittr` package:

``` r
data %>% multinom.test
#> $statistic
#> [1] 6.422941
#> 
#> $pvalue
#> [1] 6.683298e-11
```

References
----------

\[1\] Plunkett, A. & Park, J. (2018) Two-sample test for sparse high-dimensional multinomial distributions, TEST, doi.org/10.1007/s11749-018-0600-8

\[2\] Plunkett, A. & Park, J. (2017) Two-sample tests for sparse high-dimensional binary data, Communications in Statistics - Theory and Methods, 46:22, 11181-11193, DOI: 10.1080/03610926.2016.1260743
