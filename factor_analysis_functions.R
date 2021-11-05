# Factor Model of Multivariate Time Series: A Frequency Components Approach #

# Function definitions #

######################################################################################################

# R libraries #

library(Matrix)
library(pracma)
library(expm)

# Rcpp libraries # 
library(Rcpp)
library(RcppNumerical)
library(RcppArmadillo)

################################################################################

# C++ routines to call #

sourceCpp("kernel_spec_estimator.cpp") # For kernel spectral estimator

sourceCpp("eigenvalues.cpp") # For getting Eigenvalues using a C++ routine

# Nonparametric residual bootstrapping #
sourceCpp("fact_boot.cpp") # Factor model bootstrapping

################################################################################

# Additional ancillary functions
source("ancillary_functions.R")

# Band-pass filtering function
source("passfilt.R")

###############################################################################

# Eigendecomposition of M_w matrices (Eq. 17 of the paper)
# Returns a list of size length(frq.vec)
# Each element of the output list has eigenvalues and eigenvectors of the M_w matrix
fact_freq_eigenanalysis <- function( x , frq.vec , h ){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  # no. of frequencies
  n.frq = length(frq.vec)
  
  # Get list of spectral matrices at frequencies in frq.vec
  kern.sp <- kernspec( x, frq.vec , h ) # Rcpp function 
  
  sum.spec.mat <- list() 
  
  mwmats <- lapply( 1:n.frq , mw_mats , x10 = x , h10 = h ,frq.vec = frq.vec )
  
  eigen_mwmats <- lapply(1:n.frq , eigen_mw_mats , mw_matrices = mwmats )
  
  # returns a list of size n.frq. Each element has eigenvalues and eigenvectors
  # of the M_w matrices
  return(eigen_mwmats)
  
} # end function fact_freq_eigenanalysis

###################################################################################

# Bootstrap tests for factor dimension r # 

# Function to generate 500 instances of X(t,w)^*, b=1,2,...,B=500
# Inputs: frequency index wk, A_w matrix, 
# filter.w is the length of filter window used when estimating X(t,\omega)
# A_w is (p x r_0) matrix
boot.gen.f = function( wk , amat.inp,  xt.inp , frq.vec , filter.w = 0.01  ){
  
  p = dim(xt.inp)[2]
  n = dim(xt.inp)[1] 
  
  # filter xt.inp to get X(t,\omega)
  xt.filt = xt.inp # initialize
  for(pp in 1:p){
    
    xt.filt[,pp] = pass.filt1( xt.inp[,pp] , 
                              W = c(frq.vec[wk]-filter.w,frq.vec[wk]+filter.w) 
                              , type = "pass" , n=1 )
  }
  
  
  # Obtain the estimated \widehat{Y}_{t,\omega}
  yt.wk.est = xt.filt%*%amat.inp # (n x r_0) matrix
  
  # Obtain residuals 
  et_est = xt.filt -  yt.wk.est%*%t(amat.inp) # (n x p) matrix  
  
  # Nonparametric residual bootstrap
  # outputs a matrix of size 500*n x p
  # xt_boot: A matrix of size 500*n x p
   xt_boot = bootstrap_cpp( amat.inp , et_est , yt.wk.est%*%t(amat.inp) )
  
  # Note: small correction added to residual covariance matrix
  # to attain positive definiteness and better condition number
  
  return(xt_boot[[1]])
  
} # end function boot.gen.f


# Function to print bootstrap p-value for a given candidate value of r_0 
# Returns local and global test p-values
boot.pvalue.rinput  = function( rinp , x , frq.vec , h ){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  nfrq = length(frq.vec)
  
  eig.info = fact_freq_eigenanalysis(x , frq.vec , h )
  
  eig.boot.info = list() 
  for(fr in 1:nfrq){
    
    amat.cur = eig.info[[fr]]$vectors[,1:rinp]
    if(rinp==1) amat.cur = matrix(amat.cur)
    
    # get the bootstrapped X(t,\omega)
    # list of size 500, each is a (n x p) matrix
    xboot.info = boot.gen.f( fr , amat.cur , x  , frq.vec )
    
    
    # list of size 500, each element has p eigenvalues
    eig.boot.info[[fr]] = lapply( 1:500 , return_eigvals_boot , 
                                  x.boot=xboot.info , frq.boot.vec=frq.vec , 
                                  h.boot=h , fr.ind = fr , n10 = n )
    
    
  } # end loop over fr
  
  
  # Find bootstrap p-values #
  
  # test statistics will be evaluated using these frequencies
  selec.frqs = 1:nfrq
  
  pvals = rep(0, length(selec.frqs) + 1)
  local.ts = rep(0,length(selec.frqs)) # local test statistic
  local.t.boot = matrix(0,length(selec.frqs),500)
  
  for(fr in 1:length(selec.frqs)){
    
    eglm = matrix(unlist(eig.boot.info[[ selec.frqs[fr] ]]) , 
                  nrow = 500 , 
                  byrow = TRUE)
    
    eglm = (abs(eglm))
    
    obs.evals = ( abs(eig.info[[ selec.frqs[fr] ]]$values) )
    
    tot.eig.sum = sum(abs(obs.evals[1:p]))
    
    local.ts[fr] =  (sum( abs(obs.evals[1:(rinp)]) )/tot.eig.sum ) - ( 
      sum( abs(obs.evals[1:(rinp-1)]) )/tot.eig.sum) 
    
    # Bootstrapped test statistic #
    local.t.boot[fr,] = apply( eglm , 1 , 
                               function(uu) ( sum(abs(uu[1:rinp]))/sum(abs(uu[1:p]))  )  - 
                                 ( sum(abs(uu[1:(rinp-1)]))/sum(abs(uu[1:p]))  )  ) 

    
    pvals[fr] = sum( local.t.boot[fr,] > local.ts[fr]  )/500
    
  } # end loop over fr
  
  
  global.t.boot = apply(local.t.boot , 2 , mean ) # bootstrapped global test statistics
  # the above is a vector of lenth 500
  
  # global test p-value
  pvals[ length(selec.frqs) + 1] = sum( global.t.boot > mean(local.ts)  )/500
  
  
  out.ls = list(pvals , global.t.boot , local.t.boot , 
                local.ts , mean(local.ts) ,
                eig.boot.info , eig.info)
  
  names(out.ls) = c("pvals" , "globalt" , "localt" , 
                    "obs.localt"  ,"obs.globalt" , 
                    "eigen.boot" , "obs.eigen")  
  
  return( out.ls )
} # end function boot.pvalue.rinput


# Wrapper function #
# Returns factor dimension r
# Returns local factor dimensions
# Returns estimated factor series Y_t, t=1,2,...,n
factor_analysis = function(x , frq.vec , h){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  nfrq = length(frq.vec)
  
  filter.w = 0.01 # filter width
  
  pvals.ca = rep(1,p-2)
  out.l.tmp = list()
  for(r0 in 1:(p-2) ){
    
    # pvalue from global test
    out.l.tmp[[r0]]  = boot.pvalue.rinput(r0+1 , x , frq.vec , h)

    # Global test p-value
    pvals.ca[r0] =  sum(out.l.tmp[[r0]]$globalt  > 
                          out.l.tmp[[r0]]$obs.globalt)/500
    
  }
  
  # find the factor dimension
  destf = min(which(pvals.ca>0.01))

  # find local dimensions (optional include)
  #loc.dim = NULL # local dimension estimate at each frequency
  #for(frr in 1:length(frq.vec)){
    
  #  pvalvec.curr = NULL
  #  for(pp in 1:(p-2)) pvalvec.curr = c(pvalvec.curr , out.l.tmp[[pp]]$pvals[frr] )
  #  loc.dim.curr = min(which(pvalvec.curr>0.01))
  #  if(loc.dim.curr<=0) loc.dim.curr = 1
  #  loc.dim = c(loc.dim , loc.dim.curr)
  #}
  
  # Obtain the latent factor series Y_t, t=1,2,...,n
  yt = matrix(0,n,destf)
  for(j in 1:nfrq){
    
    # filter X_t to get X(t,\omega)
    xt.filt = x # initialize
    for(pp in 1:p){
      
      xt.filt[,pp] = pass.filt1( x[,pp] , 
                                 W = c(frq.vec[j]-filter.w,frq.vec[j]+filter.w) 
                                 , type = "pass" , n=1 )
    }
    
    # Obtain the estimated \widehat{Y}_{t,\omega}
    yt = yt + xt.filt%*%out.l.tmp[[destf]]$obs.eigen[[j]]$vectors[,1:destf] # (n x destf) matrix
  } # end loop over j
  
  foutlis = list(destf,yt)
  names(foutlis) = c("r" , "yt")
  return( foutlis )
  
} # end function factor_analysis
