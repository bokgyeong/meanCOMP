rm(list=ls())
library(DiceKriging); library(DiceDesign)
library(MASS)
library(foreach)
library(akima)

# install.packages("devtools")
# devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)


# =============================================================================-
# check if Bicubic spline works ----
# =============================================================================-
# num.point = 500
# num.point = 600
# num.point = 700
num.point = 800

load(paste0('~/work/meanZICOMP/appxGP/real_estGP', num.point,  '.RData'))

x.point = data.frame(X1 = log(runif(1000, mu_rng[1], mu_rng[2])), 
                     X2 = log(runif(1000, nu_rng[1], nu_rng[2])))

bs.loglam = interpp(x = th[,1], y = th[,2], z = loglam, 
                    xo = x.point$X1, yo = x.point$X2, linear = FALSE, extrap = TRUE)$z

nr.lam = sapply(1:nrow(x.point), function(i) comp_lambdas(exp(x.point[i,1]), exp(x.point[i,2]), summax = 10000)[[1]])
nr.loglam = log(nr.lam)
nr.logeta = log(nr.lam) / exp(x.point[,2])

# plot(krig.logeta, nr.logeta)
plot(bs.loglam, nr.loglam); abline(a = 0, b = 1, col = 2)
mean( (bs.loglam - nr.loglam)^2 )

save(x.point, bs.loglam, nr.lam, nr.loglam, file = paste0('~/work/meanZICOMP/appxBS/real_test', num.point, '.RData'))

