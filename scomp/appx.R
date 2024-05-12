rm(list=ls())
require(qrng); require(foreach); require(akima)
if(!require(mpcmp)){
  require(devtools)
  devtools::install_github("thomas-fung/mpcmp")
}
require(mpcmp)



# =============================================================================-
# Generate particles
# =============================================================================-
maxmu = 200
num.point = round(log(maxmu) * 5 * 15 / 10) * 10

mu_rng = c(0.01, maxmu)
nu_rng = c(0.01, 5)

Domain = matrix(
  c(log(mu_rng[1]), log(mu_rng[2]), nu_rng[1], nu_rng[2]), # for log mu and nu
  2, 2) 


# Generalized Halton sequence (randomized by definition)
n.boundary = round((Domain[2,1] - Domain[1,1]) / 0.8)*2 + round((Domain[2,2] - Domain[1,2]) / 0.8)*2

if(n.boundary > num.point * 0.2){
  th = ghalton(num.point * 0.8 , 2)
  th[,1] <- (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] <- (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  th = rbind(
    th, # add points around border line
    cbind(seq(Domain[1,1], Domain[2,1], length.out = num.point * 0.1 * 0.6 + 1), Domain[1,2]),
    cbind(seq(Domain[1,1], Domain[2,1], length.out = num.point * 0.1 * 0.6 + 1), Domain[2,2]),
    cbind(Domain[1,1], seq(Domain[1,2], Domain[2,2], length.out = num.point * 0.1 * 0.4 + 1)[-c(1, num.point * 0.1 * 0.4 +1)]),
    cbind(Domain[2,1], seq(Domain[1,2], Domain[2,2], length.out = num.point * 0.1 * 0.4 + 1)[-c(1, num.point * 0.1 * 0.4 +1)])
  )
  
} else {
  th = ghalton(num.point - n.boundary , 2)
  th[,1] <- (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] <- (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  th = rbind(
    th, # add points around border line
    cbind(seq(Domain[1,1], Domain[2,1], length.out = round((Domain[2,1] - Domain[1,1]) / 0.8) + 1), Domain[1,2]),
    cbind(seq(Domain[1,1], Domain[2,1], length.out = round((Domain[2,1] - Domain[1,1]) / 0.8) + 1), Domain[2,2]),
    cbind(Domain[1,1], seq(Domain[1,2], Domain[2,2], length.out = round((Domain[2,2] - Domain[1,2]) / 0.8) + 1)[-c(1, round((Domain[2,2] - Domain[1,2]) / 0.8) +1)]),
    cbind(Domain[2,1], seq(Domain[1,2], Domain[2,2], length.out = round((Domain[2,2] - Domain[1,2]) / 0.8) + 1)[-c(1, round((Domain[2,2] - Domain[1,2]) / 0.8) +1)])
  )
}

# plot(th)



# =============================================================================-
# Find lambda at each particle using Newton-Raphson ----
# =============================================================================-

ptm = proc.time()[3]
lam.nr.qrng = foreach(i = 1:nrow(th), .combine = 'c') %do% {
  comp_lambdas(exp(th[i,1]), th[i,2], summax = max(c(2000, 10 * exp(th[i,1]))))[[1]]
}
rtime.nr.qrng = proc.time()[3] - ptm

loglam = log(lam.nr.qrng)

save(Domain, num.point, th, lam.nr.qrng, rtime.nr.qrng, loglam,
     file = 'appx/set.RData')


