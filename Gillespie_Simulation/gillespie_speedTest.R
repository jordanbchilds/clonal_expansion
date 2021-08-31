library("parallel")
library("Rcpp")

# if zero or more than one argument is passed the number 
# of simulations is default N=1000
arg = commandArgs(trailingOnly = TRUE)
if( length(args)!=1 ){
  N.sim = arg
} else {
  N.sim = 500
}

# Conor doesn't like grey
myBlack = rgb(0,0,0, alpha=15, max=255)

# create directory to save pdf output
# full simulation results not saved
dir.create("Simulations", showWarnings = FALSE)
dir.create("Simulations/PDF", showWarnings = FALSE)
dir.create("Simulations/console", showWarnings = FALSE)

# define simulation length and hwo often populations monitored 
T.sim = 80*52*7*24*3600
dt.sim = 7*24*3600

# function to simulate using Gillespie Algo, from SMSB by D. Wilkinson
gillespied_mtDNA = function(N, T.sim, dt.sim, thresh=0.8, ...){
  tt = 0
  n = T.sim%/%dt.sim
  x = N$M
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = ncol(S)
  xmat = matrix(ncol=u, nrow=n)
  i = 1
  target = 0
  C0 = sum(x)
  repeat{
    error = C0 - sum(x)
    h = N$h(x, error)
    h0 = sum(h)
    if( h0 < 1e-10 ) tt = 1e99
    else tt = tt + rexp(1, h0)
    while( tt >= target ){
      xmat[i,] = x
      i = i + 1
      target = target + dt.sim
      if( i>n ){
        C = xmat[,1]+xmat[,2]
        mload = xmat[,2]/C
        thresh_vec = mload > thresh
        return(list("h"=function() lines(ts(mload, start=1/52, frequency=52), col=myBlack),
                    "C"=function() lines(ts(C, start=1/52, frequency=52), col=myBlack),
                    "Cmax"=max(C), 
                    "above_thresh"=thresh_vec) 
        )
      }
    }
    j = sample(v,1,prob=h/h0)
    x = x + S[,j]
  }
}


# generate intial copy number and mutation load
inits = function(n=1){
  C = rnorm(n,200,50) 
  h = rbeta(n,10,18.57)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

# define simulation parameters
# N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
#          Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )
# 
# N$h = function(x, th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025,
#                        d_w=3.06e-8, d_m=3.06e-8, m=0)){
#   th*rep(x,length.out=5)
# }

N$h = function(x, error, K_c=8.99e-9,
               th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025, d_w=3.06e-8, d_m=3.06e-8, m=0)){
  r = c("r_w","r_m")
  # th[r] = th[r] + error*K_c
  if(error>=0){
    with(as.list(c(x, th)),{
      return( th*rep(x,length.out=5) )
    })
  } else {
    # th[r] = 3.06e-8*(th[r]<0)
    with(as.list(c(x, th)), {
      return( th*rep(x,length.out=5) )
    })
  }
}

# create a list of length N (no. of sims) with the same input
# but still draws initial values from init function
gen_N = function(N, N.sim){
  NN = list()
  for(i in 1:N.sim){
    N.temp = N
    N.temp$M = inits()
    NN[[i]] = N.temp
  }
  return(NN)
}

NN = gen_N(N, N.sim)

# run the gillespie algo N.sim times using parallel::mclapply
Gillespie_Sims = mclapply(NN, gillespied_mtDNA, T.sim, dt.sim)

hprop_calc = function(Gillespie_output){
  sim_length = length(Gillespie_output[[1]][["above_thresh"]])
  hmat = matrix(NA, nrow=N.sim, ncol=sim_length)
  for(i in 1:N.sim){
    hmat[i,] = Gillespie_Sims[[i]][["above_thresh"]]
  }
  hprop = ts(colSums(hmat, na.rm=TRUE)/N.sim, start=1/52, frequency=52)
  return(hprop)
}

cppFunction('int sampleC(NumericVector hnorm, double u) {
  int v = hnorm.size(); 
  for( int i=0; i < v; i++ ){
    if( hnorm[i] > u ){
      return i+1;
    }
  }
}')

t = runif(5)
t = cumsum(t)/sum(t)
t
tt = double(1000)
system.time( for(i in 1:10000) tt[i] = sampleC(hnorm=t, u=runif(1)) )
system.time( for(i in 1:10000) sample(1:5,1,prob=t))

