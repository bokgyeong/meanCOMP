rm(list=ls())
require(ngspatial); require(foreach)
require(Rcpp); require(RcppArmadillo)


set.seed(1)
# ==============================================================================
# simulation from the spatial mean-parameterized ZICOMP regression model ----
# ==============================================================================

filename = "data/sim.RData"

# design matrix
n = 30
A = adjacency.matrix(n)
x = y = seq(0, 1, length.out = n)
coord = cbind(rep(x, times = n), rep(y, each = n))
nMonth = 6
nYear = 1
X = foreach(j = 1:nYear, .combine = 'rbind') %:%
  foreach(k = 1:nMonth, .combine = 'rbind') %do% {
    cbind(coord)
  }

t = nMonth * nYear
p1 = p2 = ncol(X); N = nrow(X)
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

beta1 = c(-2, 1)
beta2 = c(2, 2)
nu = 1.7
alpha = log(nu)
kappa = 0.2
tau = 0.2
Qs = t(M) %*% Q %*% M

C1 = chol( solve( kappa * Qs ) )
C2 = chol( solve( tau * Qs ) )
gamma = t(C1) %*% rnorm(q)
delta = t(C2) %*% rnorm(q)

eta1 = foreach(j = 1:t, .combine = 'c') %do% { as.vector( X[(nt[j]+1):(nt[j+1]),] %*% beta1 + M %*% gamma ) }
pii = exp(eta1)/(1 + exp(eta1))
mu = foreach(j = 1:t, .combine = 'c') %do% { exp( as.vector( X[(nt[j]+1):(nt[j+1]),] %*% beta2 + M %*% delta ) ) }
range(mu)

sourceCpp('src/RcppFtns.cpp')
load(paste0('appx/set.RData'))
logeta = logetaBS(log(mu), rep(nu, length(mu)), th, log(lam.nr.qrng))
eta = exp(logeta)

trpar = list(beta1 = beta1, beta2 = beta2, alpha = alpha,
             nu = nu, pii = pii, mu = mu, eta = eta,
             gamma = gamma, delta = delta,
             kappa = kappa, tau = tau, q = q)


### Rejection sampling by Chanialidis et al., 2018 -----------
w = rbinom(N, size = 1, prob = trpar$pii)
y = numeric(N)
y[w == 1] = sapply(which(w == 1), function(i)  rCOMP(1, trpar$eta[i], trpar$nu))

save(coord, nMonth, nYear, t, n, nt, N, p1, p2, M0, A, Q, X, trpar, y, w, file = filename)

