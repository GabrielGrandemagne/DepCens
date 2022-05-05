
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DepCens

<!-- badges: start -->
<!-- badges: end -->

Dependent censoring regression models for survival multivariate data.
These models are based on extensions of the frailty models, capable to
accommodating the dependence between failure and censoring times, with
Weibull and piecewise exponential marginal distributions.

## Installation

You can install the development version of DepCens from
[GitHub](https://github.com/) with:

``` r
#install.packages("devtools")
#devtools::install_github("GabrielGrandemagne/DepCens")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DepCens)
#KidneyMimic is our simulated data frame
delta_t <- ifelse(KidneyMimic$cens==1,1,0)
delta_c <- ifelse(KidneyMimic$cens==2,1,0)
fit <- dependent.censoring(formula = time ~ x1 + x2 | x3 + x1, data=KidneyMimic, delta_t=delta_t,
                           delta_c=delta_c, ident=KidneyMimic$ident, approach = "weibull")
summary_dc(fit)
#> WEIBULL APPROACH:
#> 
#> Name  Estimate    Std. Error  CI INF      CI SUP  
#> Alpha    1.337653    0.3831564   0.5866666   2.088640    
#> Sigma    0.6760496   0.2672546   0.1522306   1.199869    
#> 
#> Coefficients T:
#> 
#> Name  Estimate    Std. Error  CI INF      CI SUP  
#> x1   0.08186471  0.02224396  0.03826655  0.1254629   
#> x2   -1.412163   0.2378199   -1.878290   -0.9460356  
#> 
#> Coefficients C:
#> 
#> Name  Estimate    Std. Error  CI INF      CI SUP  
#> x3   0.2223204   0.1898567   -0.1497987  0.5944395   
#> x1   0.1932663   0.04070445  0.1134856   0.273047    
```
