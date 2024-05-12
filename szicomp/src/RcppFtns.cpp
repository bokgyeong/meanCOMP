// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <limits>
#include <sitmo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(sitmo)]]
// [[Rcpp::plugins(openmp)]]

#define minimum(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// =============================================================================
//  Chanialidis's rejection sampling for COMP distribution
// =============================================================================
// Evaluate unnormalised density of the COM-poisson distribution
// if fudge is set to lgammafn(mode+1) then the unnormalised density is one at the mode
// If mode and fudge are set to 0 then the usual unnormalised density is computed
// [[Rcpp::export]]
double unnorm_ldcpois(double x, double mu, double nu, double mode, double fudge) {
  return nu*((x-mode)*log(mu)-lgamma(x+1)+fudge);
}


// Sample from a geometric distribution truncated to {0,1,...,n}
// u is a U[0,1] realisation
// [[Rcpp::export]]
double truncated_geo_sample(double u, double logq, double n) {
  double C;
  if(logq > -DBL_EPSILON){
    return 0;
  } else {
    C = -expm1(logq*(n+1));
    return floor(log(1-C*u)/logq);
  }
}


// Sample from a geometric distribution with range {0,1,2,...}
// u is a U[0,1] realisation
// [[Rcpp::export]]
double untruncated_geo_sample(double u, double logq) {
  if (logq > -DBL_EPSILON){
    return 0;
  } else {
    return floor(log(u)/logq);
  }
}


// Compute the finite geometric series (1+q+...+q^deltax)
// [[Rcpp::export]]
double truncated_lweights(double deltax, double logq) {
  if (logq > -DBL_EPSILON)
    return log(deltax+1)+logq;
  return log1p(-exp((deltax+1)*logq)) - log1p(-exp(logq));
}


// Compute the geometric series (1+q+...)
// [[Rcpp::export]]
double untruncated_lweights(double logq) {
  return -log1p(-exp(logq));
}


// Sample from an COMP distribution
// [[Rcpp::export]]
vec rCOMP(int n, double mu, double nu){ // mu = mode
  double negInf = -std::numeric_limits<float>::infinity();;
  
  double logmu, lmode, rmode, fudge, sd, lsd, rsd, maxlweight, logprob, x, u;
  int i, attempts;
  vec ldens(4), lweights(4), logq(4), sweights(4), result(n);
  logmu = log(mu);
  
  // Figure out mode and standard deviation
  lmode = ceil(mu)-1;
  fudge = lgamma(lmode+1);
  rmode = lmode+1;
  fudge = lgamma(lmode+1);
  sd = ceil(sqrt(mu)/sqrt(nu));
  if (sd < 5) {
    sd = 5;
  }
  
  // Set up two points at mode +/- sd
  lsd = round(lmode-sd);
  if (lsd < 0){
    lsd = -1;
  }
  rsd = round(rmode+sd);
  
  // Left most tail
  if (lsd == -1) {
    lweights[0] = negInf;
    logq[0] = 0;
    ldens[0] = negInf;
  } else {
    ldens[0] = unnorm_ldcpois(lsd, mu, nu, lmode, fudge);
    if (lsd == 0) {
      lweights[0] = ldens[0];
      logq[0] = 0;
    } else {
      logq[0] = nu * (-logmu + log(lsd));
      lweights[0] = ldens[0] + truncated_lweights(lsd, logq[0]);
    }
  }
  
  // within 1sd to the left of the mode
  ldens[1] = 0;
  if (lmode == 0) {
    lweights[1] = 0;
    logq[1] = 1;
  } else {
    logq[1] = nu * (-logmu + log(lmode));
    lweights[1] = truncated_lweights(lmode-lsd-1, logq[1]);
  }
  
  // within 1sd to the right of the mode
  logq[2] = nu * (logmu - log(rmode+1));
  ldens[2] = nu * (logmu - log(rmode));
  lweights[2] = ldens[2] + truncated_lweights(rsd-rmode-1, logq[2]);
  
  // right tail
  logq[3] = nu * (logmu - log(rsd+1));
  ldens[3] = unnorm_ldcpois(rsd, mu, nu, lmode, fudge);
  lweights[3] = ldens[3] + untruncated_lweights(logq[3]);
  
  // Find maximum log-weight
  maxlweight = lweights[0];
  for (i = 1; i < 4; i++){
    if (lweights[i] > maxlweight) { maxlweight = lweights[i]; }
  }
  // Compute the cumulative sum of the weights
  for (i = 0; i < 4; i++) {
    lweights[i] = lweights[i]-maxlweight;
    sweights[0] = exp(lweights[0]);
  }
  for (i = 1; i < 4; i++) {
    sweights[i]=sweights[i-1]+exp(lweights[i]);
  }
  
  // Draw the sample by rejection sampling
  attempts = 0;
  for (i = 0; i < n; i++) {
    while (TRUE) {
      attempts = attempts + 1;
      u = randu() * sweights[3]; // *** randu ***
      if (u < sweights[0]) {
        u = u / sweights[0];
        x = truncated_geo_sample(u, logq[0], lsd);
        logprob = ldens[0]+x*logq[0];
        x = lsd-x;
      } else {
        if (u < sweights[1]) {
          u = (u-sweights[0])/(sweights[1]-sweights[0]);
          x = truncated_geo_sample(u, logq[1], lmode-lsd-1);
          logprob = ldens[1]+x*logq[1];
          x = lmode - x;
        } else {
          if (u<sweights[2]) {
            u = (u-sweights[1])/(sweights[2]-sweights[1]);
            x = truncated_geo_sample(u, logq[2], rsd-rmode-1);
            logprob = ldens[2]+x*logq[2];
            x = rmode + x;
          } else {
            u = (u-sweights[2])/(sweights[3]-sweights[2]);
            x = untruncated_geo_sample(u, logq[3]);
            logprob = ldens[3]+x*logq[3];
            x = rsd + x;
          }
        }
      }
      if (log(randu()) < unnorm_ldcpois(x, mu, nu, lmode, fudge) - logprob) { // *** randu ***
        result[i] = x;
        break;
      }
    }
  }
  return result;
}


// [[Rcpp::export]]
vec rCOMP2(vec vecmu, vec vecnu){ // mu = mode
  int n = vecmu.size();
  vec res = zeros(n);
  for(int i = 0; i < n; i ++){
    res[i] = rCOMP(1, vecmu[i], vecnu[i])[0];
  }
  return res;
}



// Sample from an COMP distribution in parallel
// [[Rcpp::export]]
vec rCOMP_parallel(vec vecmu, vec vecnu, vec seeds){ // mu = mode
  double negInf = -std::numeric_limits<float>::infinity();;
  
  unsigned int Nvec = vecmu.size();
  vec result(Nvec);
  unsigned int ncores = seeds.size();
  
#pragma omp parallel num_threads(ncores)
{
  
  uint32_t coreseed = static_cast<uint32_t>(seeds[omp_get_thread_num()]);
  sitmo::prng eng(coreseed);
  double mx = sitmo::prng::max();
  
  double mu, nu, logmu, lmode, rmode, fudge, sd, lsd, rsd, maxlweight, logprob, x, u;
  int i, attempts;
  vec ldens(4), lweights(4), logq(4), sweights(4);
  
#pragma omp for 
  for(unsigned int ii = 0; ii < Nvec; ii ++){
    mu = vecmu[ii], nu = vecnu[ii];
    x = 0;
    
    logmu = log(mu);
    
    // Figure out mode and standard deviation
    lmode = ceil(mu)-1;
    fudge = lgamma(lmode+1);
    rmode = lmode+1;
    fudge = lgamma(lmode+1);
    sd = ceil(sqrt(mu)/sqrt(nu));
    if (sd < 5) {
      sd = 5;
    }
    
    // Set up two points at mode +/- sd
    lsd = round(lmode-sd);
    if (lsd < 0){
      lsd = -1;
    }
    rsd = round(rmode+sd);
    
    // Left most tail
    if (lsd == -1) {
      lweights[0] = negInf;
      logq[0] = 0;
      ldens[0] = negInf;
    } else {
      ldens[0] = unnorm_ldcpois(lsd, mu, nu, lmode, fudge);
      if (lsd == 0) {
        lweights[0] = ldens[0];
        logq[0] = 0;
      } else {
        logq[0] = nu * (-logmu + log(lsd));
        lweights[0] = ldens[0] + truncated_lweights(lsd, logq[0]);
      }
    }
    
    // within 1sd to the left of the mode
    ldens[1] = 0;
    if (lmode == 0) {
      lweights[1] = 0;
      logq[1] = 1;
    } else {
      logq[1] = nu * (-logmu + log(lmode));
      lweights[1] = truncated_lweights(lmode-lsd-1, logq[1]);
    }
    
    // within 1sd to the right of the mode
    logq[2] = nu * (logmu - log(rmode+1));
    ldens[2] = nu * (logmu - log(rmode));
    lweights[2] = ldens[2] + truncated_lweights(rsd-rmode-1, logq[2]);
    
    // right tail
    logq[3] = nu * (logmu - log(rsd+1));
    ldens[3] = unnorm_ldcpois(rsd, mu, nu, lmode, fudge);
    lweights[3] = ldens[3] + untruncated_lweights(logq[3]);
    
    // Find maximum log-weight
    maxlweight = lweights[0];
    for (i = 1; i < 4; i++){
      if (lweights[i] > maxlweight) { maxlweight = lweights[i]; }
    }
    // Compute the cumulative sum of the weights
    for (i = 0; i < 4; i++) {
      lweights[i] = lweights[i]-maxlweight;
      sweights[0] = exp(lweights[0]);
    }
    for (i = 1; i < 4; i++) {
      sweights[i]=sweights[i-1]+exp(lweights[i]);
    }
    
    // Draw the sample by rejection sampling
    attempts = 0;
    while (TRUE) {
      attempts = attempts + 1;
      u = eng() / mx * sweights[3]; // *** randu ***
      if (u < sweights[0]) {
        u = u / sweights[0];
        x = truncated_geo_sample(u, logq[0], lsd);
        logprob = ldens[0]+x*logq[0];
        x = lsd-x;
      } else {
        if (u < sweights[1]) {
          u = (u-sweights[0])/(sweights[1]-sweights[0]);
          x = truncated_geo_sample(u, logq[1], lmode-lsd-1);
          logprob = ldens[1]+x*logq[1];
          x = lmode - x;
        } else {
          if (u<sweights[2]) {
            u = (u-sweights[1])/(sweights[2]-sweights[1]);
            x = truncated_geo_sample(u, logq[2], rsd-rmode-1);
            logprob = ldens[2]+x*logq[2];
            x = rmode + x;
          } else {
            u = (u-sweights[2])/(sweights[3]-sweights[2]);
            x = untruncated_geo_sample(u, logq[3]);
            logprob = ldens[3]+x*logq[3];
            x = rsd + x;
          }
        }
      }
      if (log( eng() / mx ) < unnorm_ldcpois(x, mu, nu, lmode, fudge) - logprob) { // *** randu ***
        result[ii] = x;
        break;
      }
    }
  }
  
}
return result;
} 


// =============================================================================
// Distribution functions
// =============================================================================

// log of negative binomial density
// [[Rcpp::export]]
double NB_logd(double y, double lambda, double nu){
  return lgamma(y + nu) - lgamma(nu) - lgamma(y + 1) + y * ( log(lambda) - log(lambda + nu) ) + nu * ( log(nu) - log(lambda + nu) );
} 

// [[Rcpp::export]]
double modeCOMP_logh_mar(double y, double mode, double nu){
  return nu * ( y * log(mode) -  lgamma(y+1) );
} 

// Unnormalized log likelihood of count component of ZICOMP
// [[Rcpp::export]]
double modeZICOMP_logh_count(vec y, vec w, vec logmode, vec nu){
  return sum( w % nu % ( y % logmode -  lgamma(y+1) ) );
} 


// Unnormalized log of normal
// [[Rcpp::export]]
double Normal_logh_mar(double y, double mu, double sig2){
  return - 0.5 * pow(y - mu, 2) / sig2;
}  

// Unnormalized log multivariate normal
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
} 




// =============================================================================
// Spline interpolation to rate parameter
// =============================================================================

// [[Rcpp::export]]
vec rcpp_interpp(vec x, vec y, vec z, vec xo, vec yo){
  
  // Obtaining namespace of akima package
  Environment pkg = Environment::namespace_env("akima");
  
  // Picking up interpp() function from akima package
  Function f = pkg["interpp"];
  
  // Executing Matrix( m, sparse = TRIE )
  List result = f(x, y, z, xo, yo, Named("linear", false), Named("extrap", true));
  return result[2];
}



// [[Rcpp::export]]
vec logetaBS(vec logmu, vec nu, mat designMat, vec loglam){
  
  vec result = rcpp_interpp(designMat.col(0), designMat.col(1), loglam, logmu, nu); 
  return( result / nu );
}



// =============================================================================
// run the exchange algorithm with spline approximatoin
// =============================================================================

// [[Rcpp::export]]
List szicompExBSnu(int outer, vec y, mat Z, mat X,
                   mat M, mat Qs, vec nt,
                   vec w, vec beta1, vec beta2, double alpha, 
                   vec gamma, vec delta, vec Igamma, vec Idelta,
                   double kappa, double tau, double phi,
                   double shape_s, double rate_s,
                   vec sigma2, mat COVbeta1, mat COVbeta2,
                   double COValpha, mat COVgamma, mat COVdelta, 
                   bool updateCOV, int adaptInterval, double adaptFactorExponent, int adapIter,
                   int preiter, int updateInterval, int thin, mat designMat, vec loglam, vec seeds){
  
  double positiveInf = std::numeric_limits<float>::infinity();;
  double negativeInf = -std::numeric_limits<float>::infinity();;
  
  int p1 = beta1.size(), p2 = beta2.size(), q = gamma.size(), t = nt.size()-1;
  int N = y.size(), n0 = M.n_rows, iter = 0, npart = designMat.n_rows, ncores = seeds.size();
  mat posterior(outer, p1+p2+1+q+q+q+q+1+1), posterior_thined;
  vec postj;
  vec gammaprop(q), deltaprop(q);
  vec Igammaprop(q), Ideltaprop(q), jprops;
  vec pii(N), piiprop(N);
  vec logmu(N), logmuprop(N), logmode(N), logmodeprop(N), mode(N), modeprop(N), eta(N), etaprop(N);
  vec nu(N), nuprop(N);
  vec dummy, aux(N), Zcomp(N);
  vec beta1prop(p1), beta2prop(p2), wprop(N), logprob_w(N);
  double shape, scale, logprob, u, aux_i, alphaprop;
  vec rhat = zeros(5), gamma1 = zeros(5), gamma2 = zeros(5), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234; 
  mat accprob = zeros(outer, 7);
  mat ZB(N, p1), XB(N, p2), MZ(n0, q), MG(n0, q);
  mat ZBprop(N, p1), XBprop(N, p2), MZprop(n0, q), MGprop(n0, q);
  
  mat cholCOVbeta1(p1, p1), cholCOVbeta2(p2, p2), cholCOVgamma(q, q), cholCOVdelta(q, q);
  double cholCOValpha;
  cholCOVbeta1 = trans( chol( sigma2[0] * ( COVbeta1 + 0.00001 * diagmat(ones(p1)) ) ) );
  cholCOVbeta2 = trans( chol( sigma2[1] * ( COVbeta2 + 0.00001 * diagmat(ones(p2)) ) ) );
  cholCOValpha = sqrt( sigma2[2] * COValpha );
  cholCOVgamma = trans( chol( sigma2[3] * ( COVgamma + 0.00001 * diagmat(ones(q)) ) ) );
  cholCOVdelta = trans( chol( sigma2[4] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
  
  int percentile = 0.0025 * npart; 
  mat Domain(2, 2);                                                         
  for(int i = 0; i < 2; i ++){
    vec dummy = sort( designMat.col(i) );
    Domain(0,i) = dummy(percentile);
    Domain(1,i) = dummy(npart - 1 - percentile);
  }
  
  // initial model parameters and latent variables
  ZB = Z * beta1;
  XB = X * beta2;
  MZ = M * (Igamma % gamma);
  MG = M * (Idelta % delta);
  
  for(int m = 0; m < t; m++){
    eta.rows(nt[m], nt[m+1]-1) = ZB.rows(nt[m], nt[m+1]-1) + MZ;
    logmu.rows(nt[m], nt[m+1]-1) = XB.rows(nt[m], nt[m+1]-1) + MG; 
  }
  
  for(int i = 0; i < N; i++){
    if(exp(eta[i]) == positiveInf){ pii[i] = 1; } else { pii[i] = exp(eta[i]) / ( 1 + exp(eta[i]) ); }
  } 
  nu = exp(alpha) * ones(N);
  logmode = logetaBS(logmu, nu, designMat, loglam);
  mode = exp(logmode);
  
  
  // Start of MCMC Chain ======================================================
  for(int k = 0; k < outer; k++){
    
    
    // sampling latent variable w ---------------------------------
    for(int i = 0; i < N; i ++){
      if(y[i] > 0) {
        w[i] = 1;
      } else {
        // exchange algorithm
        wprop[i] = 1 - w[i];
        
        // Using negative binomial
        if(wprop[i] == 0){
          aux_i = rnbinom_mu(1, nu[i], mode[i])[0];
          logprob = modeCOMP_logh_mar(aux_i, mode[i], nu[i]) - NB_logd(aux_i, mode[i], nu[i]) +
            log(1-pii[i]) - log(pii[i]);
        } else {
          aux_i = rCOMP(1, mode[i], nu[i])[0];
          logprob = NB_logd(aux_i, mode[i], nu[i]) - modeCOMP_logh_mar(aux_i, mode[i], nu[i]) +
            log(pii[i]) - log(1-pii[i]);
        }
        
        u = log( randu() );
        if( u < logprob ){
          w[i] = wprop[i];
        }
      }
    } 
    
    
    
    // sampling beta1 ---------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter, c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta1 = COVbeta1 + gamma1[0] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(0, p1-1) ) ) - COVbeta1 );
        cholCOVbeta1 = trans( chol( sigma2[0] * ( COVbeta1 + 0.00001 * diagmat(ones(p1)) ) ) );
      }
    }
    
    beta1prop = beta1 + cholCOVbeta1 * randn(p1);
    ZBprop = Z * beta1prop;
    for(int m = 0; m < t; m++){
      etaprop.rows(nt[m], nt[m+1]-1) = ZBprop.rows(nt[m], nt[m+1]-1) + MZ;
    }
    for(int i = 0; i < N; i++){
      if(exp(etaprop[i]) == positiveInf){ piiprop[i] = 1; } else { piiprop[i] = exp(etaprop[i]) / ( 1 + exp(etaprop[i]) ); }
    }
    
    logprob = sum( w % ( log(piiprop) - log(pii) ) + (1 - w) % ( log(1-piiprop) - log(1-pii) ) ) +
      MVN_logh(beta1prop, zeros(p1), diagmat(ones(p1))/100) - MVN_logh(beta1, zeros(p1), diagmat(ones(p1))/100);
    
    u = log( randu() );
    if( u < logprob ){
      beta1 = beta1prop;
      ZB = ZBprop;
      eta = etaprop;
      pii = piiprop;
      
      accprob(k, 0) = 1;
    }
    
    for(int l = 0; l < p1; l ++){
      posterior(k,l) = beta1[l];
    }
    
    
    
    // sampling beta2 ---------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1);
        rhat[1] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[1] = 1 / pow(adapIter, c1);
        gamma2[1] = c0 * gamma1[1];
        sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
        
        COVbeta2 = COVbeta2 + gamma1[1] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(p1, p1+p2-1) ) ) - COVbeta2 );
        cholCOVbeta2 = trans( chol( sigma2[1] * ( COVbeta2 + 0.00001 * diagmat(ones(p2)) ) ) );
      }
    }
    
    beta2prop = beta2 + cholCOVbeta2 * randn(p2);
    XBprop = X * beta2prop;
    for(int m = 0; m < t; m++){
      logmuprop.rows(nt[m], nt[m+1]-1) = XBprop.rows(nt[m], nt[m+1]-1) + MG;
    }
    
    if( max( logmuprop ) > Domain(1,0) || min( logmuprop ) < Domain(0,0) ){
      logprob = negativeInf;
    } else {
      logmodeprop = logetaBS(logmuprop, nu, designMat, loglam);
      modeprop = exp(logmodeprop);
      
      if( min(modeprop) == 0 ){
        logprob = negativeInf;
        
      } else {
        if(ncores == 1){
          aux = rCOMP2(modeprop, nu);
        } else {
          aux = rCOMP_parallel(modeprop, nu, seeds * randi() );
        }
        
        logprob = modeZICOMP_logh_count(y, w, logmodeprop, nu) - modeZICOMP_logh_count(y, w, logmode, nu) +
          modeZICOMP_logh_count(aux, w, logmode, nu) - modeZICOMP_logh_count(aux, w, logmodeprop, nu) +
          MVN_logh(beta2prop, zeros(p2), diagmat(ones(p2))/100) - MVN_logh(beta2, zeros(p2), diagmat(ones(p2))/100);
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      beta2 = beta2prop;
      XB = XBprop;
      logmu = logmuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 1) = 1;
    }
    
    for(int l = 0; l < p2; l ++){
      posterior(k, p1+l) = beta2[l];
    }
    
    
    
    // sampling alpha ----------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(2);
        rhat[2] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[2] = 1 / pow(adapIter, c1);
        gamma2[2] = c0 * gamma1[2];
        sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
        
        postj = posterior.col(p1+p2);
        COValpha = COValpha + gamma1[2] * ( var( postj.rows(k+1-adaptInterval, k-1) ) - COValpha );
        cholCOValpha = sqrt( sigma2[2] * COValpha );
      }
    }
    
    alphaprop = alpha + cholCOValpha * randn();
    
    if( exp(alphaprop) > Domain(1,1) || exp(alphaprop) < Domain(0,1) ){
      logprob = negativeInf;
      
    } else {
      nuprop = exp(alphaprop) * ones(N);
      logmodeprop = logetaBS(logmu, nu, designMat, loglam);
      modeprop = exp(logmodeprop);
      
      if( min(modeprop) == 0 || min(nuprop) == 0 ){
        logprob = negativeInf;
        
      } else {
        if(ncores == 1){
          aux = rCOMP2(modeprop, nuprop);
        } else {
          aux = rCOMP_parallel(modeprop, nuprop, seeds * randi());
        }
        
        logprob = modeZICOMP_logh_count(y, w, logmodeprop, nuprop) - modeZICOMP_logh_count(y, w, logmode, nu) +
          modeZICOMP_logh_count(aux, w, logmode, nu) - modeZICOMP_logh_count(aux, w, logmodeprop, nuprop) +
          Normal_logh_mar(alphaprop, 0, 100) - Normal_logh_mar(alpha, 0, 100);
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      alpha = alphaprop;
      nu = nuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 2) = 1;
    }
    
    posterior(k,p1+p2) = alpha;
    
    
    // Sampling gamma ------------------------------------------
    if( updateCOV ){
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(3);
        rhat[3] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[3] = 1 / pow(adapIter, c1);
        gamma2[3] = c0 * gamma1[3];
        sigma2[3] = exp( log(sigma2[3]) + gamma2[3] * (rhat[3] - ropt) );
        
        COVgamma = COVgamma + gamma1[3] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(p1+p2+1, p1+p2+1+q-1) ) ) - COVgamma );
        cholCOVgamma = trans( chol( sigma2[3] * ( COVgamma + 0.00001 * diagmat(ones(q)) ) ) );
      }
    }
    
    gammaprop = gamma + cholCOVgamma * randn(q);
    MZprop = M * (Igamma % gammaprop);
    for(int m = 0; m < t; m++){
      etaprop.rows(nt[m], nt[m+1]-1) = ZB.rows(nt[m], nt[m+1]-1) + MZprop;
    }
    for(int i = 0; i < N; i++){
      if(exp(etaprop[i]) == positiveInf){ piiprop[i] = 1; } else { piiprop[i] = exp(etaprop[i]) / ( 1 + exp(etaprop[i]) ); }
    }
    
    logprob = sum( w % ( log(piiprop) - log(pii) ) + (1 - w) % ( log(1-piiprop) - log(1-pii) ) ) +
      MVN_logh(gammaprop, zeros(q), kappa * Qs) - MVN_logh(gamma, zeros(q), kappa * Qs);
    
    u = log( randu() );
    if( u < logprob ){
      gamma = gammaprop;
      MZ = MZprop;
      eta = etaprop;
      pii = piiprop;
      
      accprob(k, 3) = 1;
    }
    
    for(int i = 0; i < q; i ++){
      posterior(k, p1+p2+1+i) = gamma[i];
    }
    
    
    // sampling Igamma ---------------------------------
    if( (preiter+k+1 >= updateInterval) && (preiter+k+1 - (updateInterval * trunc((preiter+k+1) / updateInterval)) == 0) ){
      
      for(int l = 0; l < q; l ++){
        Igammaprop = Igamma;
        Igammaprop[l] = 1 - Igamma[l];
        MZprop = M * (Igammaprop % gamma);
        for(int m = 0; m < t; m++){
          etaprop.rows(nt[m], nt[m+1]-1) = ZB.rows(nt[m], nt[m+1]-1) + MZprop;
        }
        for(int i = 0; i < N; i++){
          if(exp(etaprop[i]) == positiveInf){ piiprop[i] = 1; } else { piiprop[i] = exp(etaprop[i]) / ( 1 + exp(etaprop[i]) ); }
        }
        
        logprob = sum( w % ( log(piiprop) - log(pii) ) + (1 - w) % ( log(1-piiprop) - log(1-pii) ) ) +
          (Igammaprop[l] - Igamma[l]) * (log(phi) + log(1-phi));
        
        u = log( randu() );
        if( u < logprob ){
          Igamma = Igammaprop;
          MZ = MZprop;
          eta = etaprop;
          pii = piiprop;
        }
        
        posterior(k, p1+p2+1+q+q+l) = Igamma[l];
      }
    } else {
      
      for(int l = 0; l < q; l ++){
        posterior(k, p1+p2+1+q+q+l) = Igamma[l];
      }
    }
    
    
    // Sampling delta ------------------------------------------
    if( updateCOV ){
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(4);
        rhat[4] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[4] = 1 / pow(adapIter, c1);
        gamma2[4] = c0 * gamma1[4];
        sigma2[4] = exp( log(sigma2[4]) + gamma2[4] * (rhat[4] - ropt) );
        
        COVdelta = COVdelta + gamma1[4] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(p1+p2+1+q, p1+p2+1+q+q-1) ) ) - COVdelta );
        cholCOVdelta = trans( chol( sigma2[4] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
        
        adapIter = adapIter + 1;
      }
    }
    
    deltaprop = delta + cholCOVdelta * randn(q);
    MGprop = M* (Idelta % deltaprop);
    for(int m = 0; m < t; m++){
      logmuprop.rows(nt[m], nt[m+1]-1) = XB.rows(nt[m], nt[m+1]-1) + MGprop;
    }
    
    if( max( logmuprop ) > Domain(1,0) || min( logmuprop ) < Domain(0,0) ){
      logprob = negativeInf;
    } else {
      logmodeprop = logetaBS(logmuprop, nu, designMat, loglam);
      modeprop = exp(logmodeprop);
      
      if( min(modeprop) == 0 ){
        logprob = negativeInf;
        
      } else {
        if(ncores == 1){
          aux = rCOMP2(modeprop, nu);
        } else {
          aux = rCOMP_parallel(modeprop, nu, seeds * randi() );
        }

        logprob = modeZICOMP_logh_count(y, w, logmodeprop, nu) - modeZICOMP_logh_count(y, w, logmode, nu) +
          modeZICOMP_logh_count(aux, w, logmode, nu) - modeZICOMP_logh_count(aux, w, logmodeprop, nu) +
          MVN_logh(deltaprop, zeros(q), tau * Qs) - MVN_logh(delta, zeros(q), tau * Qs);
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      delta = deltaprop;
      MG = MGprop;
      logmu = logmuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 4) = 1;
    }
    
    for(int i = 0; i < q; i ++){
      posterior(k, p1+p2+1+q+i) = delta[i];
    }
    
    
    
    // sampling Idelta ---------------------------------
    if( (preiter+k+1 >= updateInterval) && (preiter+k+1 - (updateInterval * trunc((preiter+k+1) / updateInterval)) == 0) ){
      
      for(int l = 0; l < q; l ++){
        Ideltaprop = Idelta;
        Ideltaprop[l] = 1 - Idelta[l];
        MGprop = M * (Ideltaprop % delta);
        for(int m = 0; m < t; m++){
          logmuprop.rows(nt[m], nt[m+1]-1) = XB.rows(nt[m], nt[m+1]-1) + MGprop;
        }
        
        if( max( logmuprop ) > Domain(1,0) || min( logmuprop ) < Domain(0,0) ){
          logprob = negativeInf;
        } else {
          logmodeprop = logetaBS(logmuprop, nu, designMat, loglam);
          modeprop = exp(logmodeprop);
          
          if( min(modeprop) == 0 ){
            logprob = negativeInf;
            
          } else {
            if(ncores == 1){
              aux = rCOMP2(modeprop, nu);
            } else {
              aux = rCOMP_parallel(modeprop, nu, seeds * randi() );
            }
            
            logprob = modeZICOMP_logh_count(y, w, logmodeprop, nu) - modeZICOMP_logh_count(y, w, logmode, nu) +
              modeZICOMP_logh_count(aux, w, logmode, nu) - modeZICOMP_logh_count(aux, w, logmodeprop, nu) +
              (Ideltaprop[l] - Idelta[l]) * (log(phi) + log(1-phi));
          }
        }
        
        u = log( randu() );
        if( u < logprob ){
          Idelta = Ideltaprop;
          MG = MGprop;
          logmu = logmuprop;
          logmode = logmodeprop;
          mode = modeprop;
        }
        
        posterior(k, p1+p2+1+q+q+q+l) = Idelta[l];
      }
    } else {
      
      for(int l = 0; l < q; l ++){
        posterior(k, p1+p2+1+q+q+q+l) = Idelta[l];
      }
    }
    
    
    // sampling kappa ---------------------------------
    shape = shape_s + 0.5 * q;
    scale = ( 1 / ( 1/rate_s + 0.5 * trans(gamma) * Qs * gamma ) )[0];
    
    kappa = randg(1, distr_param(shape, scale))[0];
    posterior(k, p1+p2+1+q+q+q+q) = kappa;
    
    accprob(k, 5) = 1;
    
    
    // sampling tau ---------------------------------
    shape = shape_s + 0.5 * q;
    scale = ( 1 / ( 1/rate_s + 0.5 * trans(delta) * Qs * delta ) )[0];
    
    tau = randg(1, distr_param(shape, scale))[0];
    posterior(k, p1+p2+1+q+q+q+q+1) = tau;
    
    accprob(k, 6) = 1;
    
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }
    
    if ( (k+1) % 100 == 0 ) {
      Rprintf("Generated %d samples...\n", k+1);
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta1") = COVbeta1,
                            Rcpp::Named("COVbeta2") = COVbeta2,
                            Rcpp::Named("COValpha") = COValpha,
                            Rcpp::Named("COVgamma") = COVgamma,
                            Rcpp::Named("COVdelta") = COVdelta,
                            Rcpp::Named("w") = w);
}

