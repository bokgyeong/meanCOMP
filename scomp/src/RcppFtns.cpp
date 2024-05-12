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

// [[Rcpp::export]]
double modeCOMP_logh(vec y, vec logeta, vec nu){
  return sum( nu % ( y % logeta -  lgamma(y+1) ) );
}


// Unnormalized log normal
// [[Rcpp::export]]
double Normal_logh_mar(double y, double mu, double sig2){
  double result = - 0.5 * pow(y - mu, 2) / sig2;
  return result;
} 


// Unnormalized log multivariate normal
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
} 



// =============================================================================
// Bicubic spline interpolation to rate parameter
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
// Fit models
// =============================================================================

// [[Rcpp::export]]
List scompExBS(int outer, vec y, mat X, 
               mat M, mat Qs, vec nt, 
               vec beta, double alpha,
               vec delta, vec Idelta, double tau, double phi,
               double shape_s, double rate_s,
               vec sigma2, mat COVbeta, double COValpha, mat COVdelta,
               bool updateCOV, int adaptInterval, double adaptFactorExponent, int adapIter,
               int preiter, int updateInterval, int thin, mat designMat, vec loglam, vec seeds){
  
  double positiveInf = std::numeric_limits<float>::infinity();;
  double negativeInf = -std::numeric_limits<float>::infinity();;
  int p = beta.size(), q = delta.size(), t = nt.size()-1;
  int N = y.size(), n0 = M.n_rows, iter = 0, npart = designMat.n_rows, ncores = seeds.size();
  mat posterior(outer, p+1+q+q+1), posterior_thined;
  vec postj;
  vec deltaprop(q), Ideltaprop(q), jprops;
  vec logmu(N), logmuprop(N), logmode(N), logmodeprop(N), mode(N), modeprop(N);
  vec nu(N), nuprop(N);
  vec dummy, aux(N), Zcomp(N), betaprop(p);
  double shape, scale, logprob, u, aux_i, alphaprop;
  vec rhat = zeros(3), gamma1 = zeros(3), gamma2 = zeros(3), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  mat accprob = zeros(outer, 4);
  mat XB(N, p), MG(n0, q);
  mat XBprop(N, p), MGprop(n0, q);
  
  mat cholCOVbeta(p, p), cholCOVdelta(q, q);
  double cholCOValpha;
  cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.00001 * diagmat(ones(p)) ) ) );
  cholCOValpha = sqrt( sigma2[1] * COValpha );
  cholCOVdelta = trans( chol( sigma2[2] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
  
  int percentile = 0.0025 * npart;
  mat Domain(2, 2);
  for(int i = 0; i < 2; i ++){
    vec dummy = sort( designMat.col(i) );
    Domain(0,i) = dummy(percentile);
    Domain(1,i) = dummy(npart - 1 - percentile);
  }
  
  // initial model parameters and latent variables
  XB = X * beta;
  MG = M * (Idelta % delta);
  
  for(int m = 0; m < t; m++){
    logmu.rows(nt[m], nt[m+1]-1) = XB.rows(nt[m], nt[m+1]-1) + MG;
  }
  nu = exp(alpha) * ones(N);
  logmode = logetaBS(logmu, nu, designMat, loglam);
  mode = exp(logmode);
  
  
  // Start of MCMC Chain ======================================================
  for(int k = 0; k < outer; k++){
    
    
    // sampling beta ---------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter, c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta = COVbeta + gamma1[0] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(0, p-1) ) ) - COVbeta );
        cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.00001 * diagmat(ones(p)) ) ) );
      }
    }
    
    betaprop = beta + cholCOVbeta * randn(p);
    XBprop = X * betaprop;
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
        
        logprob = modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
          MVN_logh(betaprop, zeros(p), diagmat(ones(p))/100) - MVN_logh(beta, zeros(p), diagmat(ones(p))/100);  
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      beta = betaprop;
      XB = XBprop;
      logmu = logmuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 0) = 1;
    }
    
    for(int l = 0; l < p; l ++){
      posterior(k, l) = beta[l];
    }
    
    
    // sampling alpha ----------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1);
        rhat[1] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[1] = 1 / pow(adapIter, c1);
        gamma2[1] = c0 * gamma1[1];
        sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
        
        postj = posterior.col(p);
        COValpha = COValpha + gamma1[1] * ( var( postj.rows(k+1-adaptInterval, k-1) ) - COValpha );
        cholCOValpha = sqrt( sigma2[1] * ( COValpha + 0.00001 ) );
      }
    }
    
    alphaprop = alpha + cholCOValpha * randn();
    
    if( exp(alphaprop) > Domain(1,1) || exp(alphaprop) < Domain(0,1) ){
      logprob = negativeInf;	
      
    } else {
      nuprop = exp(alphaprop) * ones(N);
      logmodeprop = logetaBS(logmu, nuprop, designMat, loglam);
      modeprop = exp(logmodeprop);
      
      if( min(modeprop) == 0 || min(nuprop) == 0 ){
        logprob = negativeInf;
        
      } else {
        if(ncores == 1){
          aux = rCOMP2(modeprop, nuprop);
        } else {
          aux = rCOMP_parallel(modeprop, nuprop, seeds * randi());  
        }
        
        logprob = modeCOMP_logh(y, logmodeprop, nuprop) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nuprop) +
          Normal_logh_mar(alphaprop, 0, 100) - Normal_logh_mar(alpha, 0, 100);
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      alpha = alphaprop;
      nu = nuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 1) = 1;
    }
    
    posterior(k, p) = alpha;
    
    
    // Sampling delta ------------------------------------------
    if( updateCOV ){
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(2);
        rhat[2] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[2] = 1 / pow(adapIter, c1);
        gamma2[2] = c0 * gamma1[2];
        sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
        
        COVdelta = COVdelta + gamma1[2] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(p+1, p+1+q-1) ) ) - COVdelta );
        cholCOVdelta = trans( chol( sigma2[2] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
        
        adapIter = adapIter + 1;
      }
    }
    
    deltaprop = delta + cholCOVdelta * randn(q);
    MGprop = M * (Idelta % deltaprop);
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
        
        logprob = modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
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
      
      accprob(k, 2) = 1;
    }
    
    for(int i = 0; i < q; i ++){
      posterior(k, p+1+i) = delta[i];
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
            
            logprob =  modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
              modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
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
        
        posterior(k, p+1+q+l) = Idelta[l];
      }
    } else {
      
      for(int l = 0; l < q; l ++){
        posterior(k, p+1+q+l) = Idelta[l];
      }
    }
    
    
    // sampling tau ---------------------------------
    shape = shape_s + 0.5 * q;
    scale = ( 1 / ( 1/rate_s + 0.5 * trans(delta) * Qs * delta ) )[0];
    
    tau = randg(1, distr_param(shape, scale))[0];
    posterior(k, p+1+q+q) = tau;
    
    accprob(k, 3) = 1;
    
    
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
                            Rcpp::Named("COVbeta") = COVbeta,
                            Rcpp::Named("COValpha") = COValpha,
                            Rcpp::Named("COVdelta") = COVdelta);
}



// [[Rcpp::export]]
List scompExBSoffset(int outer, vec y, mat X, vec offset,
                     mat M, mat Qs, vec nt, 
                     vec beta, double alpha,
                     vec delta, vec Idelta, double tau, double phi,
                     double shape_s, double rate_s,
                     vec sigma2, mat COVbeta, double COValpha, mat COVdelta,
                     bool updateCOV, int adaptInterval, double adaptFactorExponent, int adapIter,
                     int preiter, int updateInterval, int thin, mat designMat, vec loglam, vec seeds){
  
  double positiveInf = std::numeric_limits<float>::infinity();;
  double negativeInf = -std::numeric_limits<float>::infinity();;
  int p = beta.size(), q = delta.size(), t = nt.size()-1;
  int N = y.size(), n0 = M.n_rows, iter = 0, npart = designMat.n_rows, ncores = seeds.size();
  mat posterior(outer, p+1+q+q+1), posterior_thined;
  vec postj;
  vec deltaprop(q), Ideltaprop(q), jprops;
  vec logmu(N), logmuprop(N), logmode(N), logmodeprop(N), mode(N), modeprop(N);
  vec nu(N), nuprop(N);
  vec dummy, aux(N), Zcomp(N), betaprop(p);
  double shape, scale, logprob, u, aux_i, alphaprop;
  vec rhat = zeros(3), gamma1 = zeros(3), gamma2 = zeros(3), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  mat accprob = zeros(outer, 4);
  mat XB(N, 1), MG(n0, 1);
  mat XBprop(N, 1), MGprop(n0, 1);
  
  mat cholCOVbeta(p, p), cholCOVdelta(q, q);
  double cholCOValpha;
  cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.00001 * diagmat(ones(p)) ) ) );
  cholCOValpha = sqrt( sigma2[1] * COValpha );
  cholCOVdelta = trans( chol( sigma2[2] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
  
  int percentile = 0.0025 * npart;
  mat Domain(2, 2);
  for(int i = 0; i < 2; i ++){
    vec dummy = sort( designMat.col(i) );
    Domain(0,i) = dummy(percentile);
    Domain(1,i) = dummy(npart - 1 - percentile);
  }
  
  // initial model parameters and latent variables
  XB = X * beta + offset;
  MG = M * (Idelta % delta);
  
  for(int m = 0; m < t; m++){
    logmu.rows(nt[m], nt[m+1]-1) = XB.rows(nt[m], nt[m+1]-1) + MG;
  }
  nu = exp(alpha) * ones(N);
  logmode = logetaBS(logmu, nu, designMat, loglam);
  mode = exp(logmode);
  
  
  // Start of MCMC Chain ======================================================
  for(int k = 0; k < outer; k++){
    
    
    // sampling beta ---------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter, c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta = COVbeta + gamma1[0] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(0, p-1) ) ) - COVbeta );
        cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.00001 * diagmat(ones(p)) ) ) );
      }
    }
    
    betaprop = beta + cholCOVbeta * randn(p);
    XBprop = X * betaprop + offset;
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
        
        logprob = modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
          MVN_logh(betaprop, zeros(p), diagmat(ones(p))/100) - MVN_logh(beta, zeros(p), diagmat(ones(p))/100);  
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      beta = betaprop;
      XB = XBprop;
      logmu = logmuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 0) = 1;
    }
    
    for(int l = 0; l < p; l ++){
      posterior(k, l) = beta[l];
    }
    
    
    // sampling alpha ----------------------------------------
    if( updateCOV ){
      
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1);
        rhat[1] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[1] = 1 / pow(adapIter, c1);
        gamma2[1] = c0 * gamma1[1];
        sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
        
        postj = posterior.col(p);
        COValpha = COValpha + gamma1[1] * ( var( postj.rows(k+1-adaptInterval, k-1) ) - COValpha );
        cholCOValpha = sqrt( sigma2[1] * ( COValpha + 0.00001 ) );
      }
    }
    
    alphaprop = alpha + cholCOValpha * randn();
    
    if( exp(alphaprop) > Domain(1,1) || exp(alphaprop) < Domain(0,1) ){
      logprob = negativeInf;	
    } else {
      nuprop = exp(alphaprop) * ones(N);
      logmodeprop = logetaBS(logmu, nuprop, designMat, loglam);
      modeprop = exp(logmodeprop);
      
      if( min(modeprop) == 0 || min(nuprop) == 0 ){
        logprob = negativeInf;
        
      } else {
        if(ncores == 1){
          aux = rCOMP2(modeprop, nuprop);  
        } else {
          aux = rCOMP_parallel(modeprop, nuprop, seeds * randi());
        }
        
        logprob = modeCOMP_logh(y, logmodeprop, nuprop) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nuprop) +
          Normal_logh_mar(alphaprop, 0, 100) - Normal_logh_mar(alpha, 0, 100);
      }
    }
    
    u = log( randu() );
    if( u < logprob ){
      alpha = alphaprop;
      nu = nuprop;
      logmode = logmodeprop;
      mode = modeprop;
      
      accprob(k, 1) = 1;
    }
    
    posterior(k, p) = alpha;
    
    
    // Sampling delta ------------------------------------------
    if( updateCOV ){
      if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(2);
        rhat[2] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
        gamma1[2] = 1 / pow(adapIter, c1);
        gamma2[2] = c0 * gamma1[2];
        sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
        
        COVdelta = COVdelta + gamma1[2] * ( cov( posterior( span(k+1-adaptInterval, k-1), span(p+1, p+1+q-1) ) ) - COVdelta );
        cholCOVdelta = trans( chol( sigma2[2] * ( COVdelta + 0.00001 * diagmat(ones(q)) ) ) );
        
        adapIter = adapIter + 1;
      }
    }
    
    deltaprop = delta + cholCOVdelta * randn(q);
    MGprop = M * (Idelta % deltaprop);
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
        
        logprob = modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
          modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
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
      
      accprob(k, 2) = 1;
    }
    
    for(int i = 0; i < q; i ++){
      posterior(k, p+1+i) = delta[i];
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
            
            logprob =  modeCOMP_logh(y, logmodeprop, nu) - modeCOMP_logh(y, logmode, nu) +
              modeCOMP_logh(aux, logmode, nu) - modeCOMP_logh(aux, logmodeprop, nu) +
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
        
        posterior(k, p+1+q+l) = Idelta[l];
      }
    } else {
      
      for(int l = 0; l < q; l ++){
        posterior(k, p+1+q+l) = Idelta[l];
      }
    }
    
    
    // sampling tau ---------------------------------
    shape = shape_s + 0.5 * q;
    scale = ( 1 / ( 1/rate_s + 0.5 * trans(delta) * Qs * delta ) )[0];
    
    tau = randg(1, distr_param(shape, scale))[0];
    posterior(k, p+1+q+q) = tau;
    
    accprob(k, 3) = 1;
    
    
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
                            Rcpp::Named("COVbeta") = COVbeta,
                            Rcpp::Named("COValpha") = COValpha,
                            Rcpp::Named("COVdelta") = COVdelta);
}



