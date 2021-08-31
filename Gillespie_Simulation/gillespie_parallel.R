library("parallel")

closeAllConnections()

# if anything other than one argument is passed the number 
# of simulations is default N=1000
arg = commandArgs(trailingOnly = TRUE)
if( length(args)!=1 ){
  N.sim = arg
} else {
  N.sim = 500
}


numCore = detectCores()-1
cl  = makeCluster(numCore) 

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
gillespied = function(N, T.sim, dt.sim, thresh=0.8, ...){
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
    h = N$h(x)
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
    j = sample(v,1,prob=h)
    x = x + S[,j]
  }
}


# N$h = function(x, error, K_c=8.99e-9,
#                th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025, d_w=3.06e-8, d_m=3.06e-8, m=0)){
#   repo = c("r_w", "r_m")
#   th.1 = th + error*K_c*(names(th)%in%repo)
#   if(error>=0){
#     return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
#   } else {
#     if(th.1["r_w"]<0) th.1["r_w"] = th["r_w"]
#     if(th.1["r_m"]<0) th.1["r_m"] = th["r_m"]
#     return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
#   }
# }

# generate initial copy number and mutation load
inits = function(n=1){
  C = rnorm(n,200,50) 
  h = rbeta(n,10,18.57)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

N$h = function(x, th=c(3.06e-8*0.975, 3.06e-8*1.025,
                       3.06e-8, 3.06e-8, 0)){
  with(as.list(c(x,th)), {
    return(c(th[1]*x[1], th[2]*x[2], th[3]*x[1], th[4]*x[2], th[5]*x[1]))
  })
  
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

clusterExport(cl, c("gillespied", "NN"))

tt = lapply(NN, gillespied, T.sim, dt.sim)

# run the gillespie algo N.sim times using parallel::mclapply
Gillespie_Sims = parLapply(cl, NN, gillespied, T.sim, dt.sim)

stopCluster(cl)

hprop_calc = function(Gillespie_output){
  sim_length = length(Gillespie_output[[1]][["above_thresh"]])
  hmat = matrix(NA, nrow=N.sim, ncol=sim_length)
  for(i in 1:N.sim){
    hmat[i,] = Gillespie_Sims[[i]][["above_thresh"]]
  }
  hprop = ts(colSums(hmat, na.rm=TRUE)/N.sim, start=1/52, frequency=52)
  return(hprop)
}


# save output as pdf 
pdf("Simulations/PDF/simPDF.pdf", width=14, height=8.5)
op = par(mfrow=c(1,2))
plot(1, type='n', ylim=c(0,1), xlim=c(0,80), main="",
     xlab="Time (years)", ylab="Mutation Load", cex.lab=1.4)
for(i in 1:N.sim){
  Gillespie_Sims[[i]][["h"]]()
}
ymax = double(N.sim)
for(i in 1:N.sim) ymax[i] = Gillespie_Sims[[i]][["Cmax"]]
plot(1, type='n', ylim=c(0,max(ymax)),  xlim=c(0,80), main="", 
     xlab="Time (years)", ylab="Copy Number", cex.lab=1.4)
for(i in 1:N.sim){
  Gillespie_Sims[[i]][["C"]]()
}
h_prop = hprop_calc(Gillespie_Sims)
plot(h_prop, ylim=c(0,1), xlim=c(0,80), main="Mutation Load Above 80%",
     ylab="Proportion", cex.lab=1.4, cex.main=1.4)
par(op)
dev.off()















