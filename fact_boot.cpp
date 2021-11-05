// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppNumerical.h>

#include <math.h>

#define _USE_MATH_DEFINES

#include <cmath>

//using namespace Numer;
using namespace std;


//using namespace std;

// Function declaration with export tag
// [[Rcpp::export]]
Rcpp::List bootstrap_cpp(Rcpp::NumericMatrix amats, 
              Rcpp::NumericMatrix etw, Rcpp::NumericMatrix ayts,  
              int B = 500) {
  
  int n1c = etw.nrow(), p1c = etw.ncol();
  int r1c = amats.ncol();
  
  // Output the final X(t,omega)* bootstrapped matrix, b=1,2,...,500
  Rcpp::NumericMatrix xtout(B*n1c,p1c);
  
  // Perform bootstrap
  for(int i = 0; i < B; i++) {

    // Generate indices to select for iteration i
    Rcpp::NumericVector seedsv = floor(Rcpp::runif(n1c, 0, n1c));
    
    Rcpp::NumericMatrix gen_errs(n1c,p1c);   
    for(int j = 0; j < n1c; j++){
      
      int curseed = static_cast<int>( seedsv[j] );
      
    gen_errs( j , _ ) = etw(curseed, _ );
    
    xtout( i*n1c + j , _ ) = ayts(j, _ ) +  gen_errs( j , _ ); 
      
    }   
    
  }
  
  // Return bootstrap results
  Rcpp::List outlist=List::create(xtout);
  return outlist;
}