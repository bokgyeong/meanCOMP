rm(list=ls())
require(ngspatial); require(foreach)
require(Rcpp); require(RcppArmadillo)

set.seed(1)
# ==============================================================================
# simulation from the spatial mean-parameterized COMP regression model ----
# ==============================================================================

filename = "data/sim.RData"

# design matrix
n = 30
A = adjacency.matrix(n)
x = y = seq(0, 1, length.out = n)
coord = cbind(rep(x, times = n), rep(y, each = n))
nMonth = 1
nYear = 1
X = foreach(j = 1:nYear, .combine = 'rbind') %:%
  foreach(k = 1:nMonth, .combine = 'rbind') %do% {
    cbind(coord)
  }

t = nMonth * nYear
p = ncol(X); N = nrow(X)
n = rep(nrow(A), t); nt = c(0, cumsum(n))

ones = rep(1, n[1])
Imat = ones %*% solve(t(ones) %*% ones) %*% t(ones)
Ic = diag(1, n[1]) - Imat
Moran = Ic %*% A %*% Ic
eig = eigen(Moran)
M0 = eig$vectors[,1:300]
Q = diag(rowSums(A, n[1])) - A


# true parameters
q = 25
M = M0[,1:q]

beta = c(2, 2)
nu = 1.7
alpha = log(nu)
tau = 0.2

Qs = t(M) %*% Q %*% M
C = chol( solve( tau * Qs ) )
delta = t(C) %*% rnorm(q)
mu = foreach(j = 1:t, .combine = 'c') %do% { exp( as.vector( X[(nt[j]+1):(nt[j+1]),] %*% beta + M %*% delta ) ) }
range(mu)

sourceCpp('src/RcppFtns.cpp')
load(paste0('appx/set.RData'))

logeta = logetaBS(log(mu), rep(nu, length(mu)), th, log(lam.nr.qrng))
eta = exp(logeta)

trpar = list(beta = beta, alpha = alpha, nu = nu, mu = mu, eta = eta,
             delta = delta, tau = tau, q = q)

y = sapply(1:N, function(i)  rCOMP(1, trpar$eta[i], trpar$nu))

save(coord, nMonth, nYear, t, n, nt, N, p, M0, A, Q, X, trpar, y, file = filename)


