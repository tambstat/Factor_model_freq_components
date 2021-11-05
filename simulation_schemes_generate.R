# Functions for generating the series X_t, t=1,2,...,n #

# Inputs #
# n = sample size (should be at least 200)
# p = input dimension
# r = factor dimension (should be 1,2,3 or 4)
# h = bandwidth of the kernel spectral estimator
# sigma2 = variance of the error term
# inp.seed = seed for simulating the stationary time series

###########################################################################
# Scheme 1 #

scheme1_generate = function(n,p,r,h,sigma2,inp.seed){
  
  # Generate skew symmetric matrix 
  set.seed(49)
  A1 = matrix(runif(p^2,-1,1),p,p)
  A1 = A1 - t(A1)
  
  # Generate orthogonal matrix
  set.seed(59)
  A2 = randortho(p)
  
  yt = matrix(0 , n , 4 )
  
  # Frequencies in (0,0.5) #
  fr1 = seq(0.01,0.47,0.02)
  fr2 = seq(0.01,0.47,0.02) + 0.02
  nff  = length(fr1)
  
  # A matrices changing at every frequency #
  Aml = list()
  for(kk in 1:nff) {
    Aml[[kk]] = expm( (fr1[kk]+fr2[kk])/2*A1  )%*%A2[,1:r]
  }
  
  set.seed(inp.seed)
  
  xtemp1 = arima.sim( model = list(ar=c(0.9), ma = c(0.8,-0.2) , sd  = 1 ) ,  n = n+15)
  xtemp3 = arima.sim( model = list(ar=c(1.25,-0.75,0.3) , sd  = 1 )   , n = n+15) 

  
  yt[,1] = xtemp1[1:n] 
  yt[,2] = 0.8*xtemp1[2:(n+1)]
  yt[,3] = xtemp3[1:n] 
  yt[,4] = 0.75*xtemp3[3:(n+2)]
  
  # Noise term
  et = matrix(rnorm(n*p,mean=0,sd=sqrt(sigma2)),n,p)
  
  xt = matrix(0,n,p)
  
  # array for storing filtered Y_t series
  yt.f.c = array(0,c(nff,n,4))
  
  for(find in 1:nff){  
    for(pp in 1:4){
      
      yt.f.c[find,,pp] = pass.filt1(yt[,pp],
                                    W = c(fr1[find],fr2[find]),type="pass" , n=1)
    }
  }
  
  xt = matrix(0,n,p)
  for(find in 1:nff){
    xt = xt + t( Aml[[find]][,1:r]%*%t(yt.f.c[find,,1:r])  )
  }  
  
  xt = xt + et
  
  return(xt)
  
} # end function scheme1_generate()

############################################################################

# Scheme 2 #

scheme2_generate = function(n,p,r,h,sigma2,inp.seed){

# Generate skew symmetric matrix 
set.seed(49)
A1 = matrix(runif(p^2,-1,1),p,p)
A1 = A1 - t(A1)

# Generate orthogonal matrix
set.seed(59)
A2 = randortho(p)

yt = matrix(0 , n , 4 )

# Frequencies in (0,0.5) #
fr1 = seq(0.01,0.47,0.02)
fr2 = seq(0.01,0.47,0.02) + 0.02
nff  = length(fr1)

# A matrices changing at every frequency #
Aml = list()
for(kk in 1:nff) {
  Aml[[kk]] = expm( (fr1[kk]+fr2[kk])/2*A1  )%*%A2[,1:r]
}

set.seed(inp.seed)

xtemp5 = arima.sim( model = list(ar=c(0.8) , sd  = 1 ) ,   n = n+15)
xtemp7 = arima.sim( model = list(ar=c(1.5,-0.75) , sd  = 1 ) ,   n = n+15)

yt[,1] = xtemp5[1:n] 
yt[,2] = 0.8*xtemp5[2:(n+1)]
yt[,3] = xtemp7[1:n] 
yt[,4] = 0.75*xtemp7[3:(n+2)]

# Noise term
et = matrix(rnorm(n*p,mean=0,sd=sqrt(sigma2)),n,p)

xt = matrix(0,n,p)

# array for storing filtered Y_t series
yt.f.c = array(0,c(nff,n,4))

for(find in 1:nff){  
  for(pp in 1:4){
    
    yt.f.c[find,,pp] = pass.filt1(yt[,pp],
                                  W = c(fr1[find],fr2[find]),type="pass" , n=1)
  }
}

xt = matrix(0,n,p)
for(find in 1:nff){
  xt = xt + t( Aml[[find]][,1:r]%*%t(yt.f.c[find,,1:r])  )
}  

xt = xt + et

return(xt)

} # end function scheme2_generate

###############################################################################