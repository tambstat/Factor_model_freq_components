# Factor Model of Multivariate Time Series: A Frequency Components Approach #

# Simulation Example  #

# Read functions for factor analysis
source("factor_analysis_functions.R")


# Read functions for simulating a time series following Shceme 1 or Scheme 2
# See Section 3 of the paper for descriptions #
source("simulation_schemes_generate.R")

n = 500 # sample size
p = 6 # input dimension
r = 2 # factor dimension
sigma2 = 0.5^2  # variance of the noise term
h = n^(-0.2) # bandwidth of the kernel spectral estimator

# Scheme 2 #
xt2 = scheme2_generate(n = n, p = p , r = r , h = h , sigma2=sigma2 , inp.seed=1509)

# Scheme 1 #
xt1 = scheme1_generate(n = n, p = p , r = r , h = h , sigma2=sigma2 , inp.seed=1509)

# Frequencies chosen in (0,0.5) to run the factor analysis #
frq.vec = seq(0.02,0.46,0.04)

################################################################################

# Run factor analysis #
ffact = factor_analysis( xt2 , frq.vec , h)

# Print estimated factor dimension
print(ffact$r)

# Plot the estimated factor series Y_t, t=1,2,...,n
plot.ts(ffact$yt , main = "Factor series")
