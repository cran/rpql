#include <RcppArmadillo.h>
// #include <R.h>
// #include <Rmath.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
SEXP newcoeffn(SEXP xstarin, SEXP Qin, SEXP ystarin) {
     mat xstar = as<mat>(xstarin);
     mat Q = as<mat>(Qin);
     vec ystar = as<vec>(ystarin);
     //mat cp = inv(trans(xstar)*xstar + Q);
     mat cp = trans(xstar)*xstar + Q;
     //vec coef = cp*trans(xstar)*ystar; 
     vec b = trans(xstar)*ystar; 
     vec coef = solve(cp, b); 
     return(wrap(coef));
     }


