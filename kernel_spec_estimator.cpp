//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>


#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <cmath>
#include <complex> 
#include <limits>



using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
Rcpp::List kernspec(arma::mat X, NumericVector frqr , double h){
  //inputs
  //X time series with T rows and R columns (column vectors are components of MV time series)
  //h is the bandwidth
  //frqr is the vector of frequencies at which spectral matrix is required
  
  //initialization
  int R=X.n_cols; //number of components of MV time series
  int Ts=X.n_rows; //total length of time series
  int Fs=frqr.length(); //number of required frequencies for spectral estimator
  
  double kwt;

  arma::cx_mat fftmat(Ts,R);
  arma::cx_cube pgram(R,R,Ts);
  arma::cx_cube specm(R,R,Fs);
  
  //Temporary variables
  arma::cx_mat tmp1(Ts,1);
  arma::cx_mat tmp2(Ts,1);
  
  arma::cx_mat tmp3(R,R);
  arma::cx_mat tmp4(R,R);
  for(int ii; ii<R; ii++){
    for(int jj=0; jj<R; jj++){
      tmp3( ii , jj )=0;
      tmp4(ii,jj) = 0;
    }
  }
  
  // Find the Fourier transforms 
  for(int i=0; i<R; i++){
   fftmat.col(i) = fft(X.col(i))/pow(Ts,0.5); ///pow(2*PI*Ts,0.5); 
  }

  // Find the periodogram
  for(int i=0;i<Ts;i++){
    tmp1 = (fftmat.submat( i , 0 , i , R-1 ));
    tmp2 = (fftmat.submat( i , 0 , i , R-1 ));
    
    pgram.slice(i) =  (tmp1.st())*conj(tmp2);
    }
  
  
  //Smooth the periodogram
  for(int k=0; k<Fs; k++){ //loop over frequencies
    
    specm.slice(k) = tmp4; // initialize
    
    double kwts=0; // for summing kernel weights
    
    double kdist1,kdist2,kwt1,kwt2 = 0;
    
    for(int j=0;j<Ts;j++){
      
      if(j<=(Ts/2-1)){
      
        kdist1 = abs( frqr[k] - (j*1.0/Ts))  ;

        if(kdist1<=h) { kwt1 = 0.75*(1 - kdist1*kdist1/pow(h,2))/h/Ts;
          }
        if(kdist1>h) { kwt1 = 0; 
          }
      
      specm.slice(k) = specm.slice(k) +  kwt1*pgram.slice(j);       
      kwts = kwts + kwt1;
      
      }
      
      if(j>=(Ts/2)){
         kdist2 = abs( frqr[k] - (j*1.0/Ts)) ; 
           //min( abs(frqr[k] + (j/Ts-1))  , abs(frqr[k] - (j/Ts)) )  ;
        
        if(kdist2<=h) { kwt2 = 0.75*(1 - kdist2*kdist2/pow(h,2))/h/Ts; 
          }
        if(kdist2>h) { kwt2 = 0; 
          }
        
        specm.slice(k) = specm.slice(k) +  kwt2*pgram.slice(j);
        kwts = kwts + kwt2;
      }
        
    } // end loop over j=0,1,...,Ts-1 time points 
    
   // Rcout<<kwts<<"---";
    
  } // end loop over k=0,1,...,Fs frequencies
  
  Rcpp::List out=List::create(specm);
  return out;
}
