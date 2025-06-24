# sgeenpe

The goal of sgeenpe is to fit marginal proportional hazard model for clustered survival data under the right censoring mechanism .

## Installation

You can install the development version of sgeenpe from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::devtools("Yangcc123/sgeenpe")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sgeenpe)
data("teeth")
formula111 = Surv(teeth$time, teeth$cens)~Decayed + sAge +Smoking+filled_tooth_sum res = sgeenpe(formula111,teeth)
```
