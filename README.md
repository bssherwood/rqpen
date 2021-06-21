
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rqPen

<!-- badges: start -->

[![R-CMD-check](https://github.com/be-green/rqpen/workflows/R-CMD-check/badge.svg)](https://github.com/be-green/rqpen/actions)
<!-- badges: end -->

The rqPen package performs penalized quantile regression for LASSO, SCAD
and MCP functions including group penalties. Provides a function that
automatically generates lambdas and evaluates different models with
cross validation or BIC, including a large p version of BIC.

## Installation

You can install the released version of rqPen from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rqPen")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bssherwood/rqpen")
```

## Penalized Quantile Regression

``` r
library(rqPen)

data("engel", package = "quantreg")

engel <- data.frame(engel, no_relationship = exp(rnorm(nrow(engel))))

x <- model.matrix(foodexp ~ 0 + income + no_relationship, data = engel)

# get a zero for the coefficient which we made up
rq.lasso.fit(
  x = x,
  y = engel$foodexp,
  lambda = 0.3, intercept = T)
#> 
#> Coefficients:
#>       intercept          income no_relationship 
#>     439.2074796       0.1641568       0.0000000
```
