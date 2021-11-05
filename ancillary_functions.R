# Factor Model MTS using Frequency components #

# Ancillary Function Definitions # 

# 4th ancilary function: returns matrix M_w at freq corresponding to indxs4
# See Eq. (13) of paper 
mw_mats <- function(frq.tg , frq.vec , x10 , h10) {
  
  p = dim(x10)[2]
  n = dim(x10)[1]
  
  local.window = frq.vec[frq.tg]+seq(-0.009,0.009,0.003) 
  
  ntotf = length( local.window )
  
  specmts = kernspec(x10 , local.window , h10 )[[1]]
  
  reavg = matrix(0,p,p)
  for(nind in 1:ntotf) reavg = reavg + Re(specmts[,,nind])*0.003/2/0.009
  
  mwout = Im(specmts[,,ceil(ntotf/2)])%*%t(Im(specmts[,,ceil(ntotf/2)])) + 
    ( Re(specmts[,,ceil(ntotf/2)])- reavg)%*%t(Re(specmts[,,ceil(ntotf/2)])- reavg)
  
  return( mwout + 10^(-5)*diag(p) )
  #return( mwout  )
  # Note: small correction added to M_w matrix

} # end function mw_mats()

# Find eigenvalues and eigenvectors of mw matrices at frequency index indxs5
eigen_mw_mats = function(indxs5 , mw_matrices){
  return(eigen(mw_matrices[[indxs5]]))
  
}


#############################################################################

# Ancillary functions for Bootstrap 

return_error_boot = function(pseudo.ind1 , err.mat){
  selected.indices = sample(1:dim(err.mat)[1] , dim(err.mat)[1] , replace=TRUE  )
  return( err.mat[selected.indices , ]  )
}

return_xt_boot = function( pseudo.ind2 , err.mat.list , amat.boot 
                           , yhat.boot ){
  return( yhat.boot%*%t(amat.boot) + err.mat.list[[pseudo.ind2]]    )
}


#############################################################################


# Eigendecomposition of M_w matrix at a given \omega
# eigenvalues and eigenvectors of the M_w matrix at a given \omega
freq_wise_eigenanalysis <- function( x10 , frq.targ.ind , frq.vec , h ){
  
  p = dim(x10)[2]
  n = dim(x10)[1]
  
  nfrq = length(frq.vec)
  
  
  mwmat.fr <- mw_mats( frq.targ.ind , frq.vec , x10 , h )
  
  return( getEigenValues(mwmat.fr) ) # uses C++ routine
  
} # end function freq_wise_eigenanalysis


return_eigvals_boot = function( pseudo.ind2 , x.boot , frq.boot.vec ,  h.boot ,fr.ind , n10 , B10=500 ){
  tmp = freq_wise_eigenanalysis(x.boot[ (pseudo.ind2*n10-n10+1):(pseudo.ind2*n10-n10+1+n10-1),] , fr.ind , 
                                frq.boot.vec , h.boot )
  return( tmp   )
}