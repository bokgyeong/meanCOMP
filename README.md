# Supplemental Code for "Fast Bayesian Inference for Spatial Mean-Parameterized Conway--Maxwell--Poisson Models"
Authors: Bokgyeong Kang, John Hughes, and Murali Haran

We provide instructions for fitting spatial mean-parameterized Conway--Maxwell--Poisson (mean-COMP) regression model and spatial zero-inflated mean-parameterized Conway--Maxwell--Poisson (mean-ZICOMP) regression model. 

## Required Packages:
The code has been tested with R version 4.2.2, "Innocent and Trusting."  The following R packages must be installed before the code will run successfully:

- `Rcpp`
- `RcppArmadillo`
- `sitmo`
- `akima`
- `mpcmp`
- `qrng`
- `ngspatial`
- `foreach`
- `coda`
- `batchmeans`
- `tidyverse`

The `mpcmp` can be installed from [GitHub](https://github.com/thomas-fung/mpcmp):
```s
require(devtools)
devtools::install_github("thomas-fung/mpcmp")
```
