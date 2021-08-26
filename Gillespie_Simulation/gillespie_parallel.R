library("parallel")

# if zero or more than one argument is passed the number 
# of simulations is default N=1000
arg = commandArgs(trailingOnly = TRUE)
if( length(args)!=1 ){
  N.sim = arg
} else {
  N.sim = 1000
}

# Conor doesn't like grey
myBlack = rgb(0,0,0, alpha=25, max=255)

# create directory to save pdf output
# full simulation results not saved
dir.create("Simulations", showWarnings = FALSE)
dir.create("Simulations/PDF", showWarnings = FALSE)
dir.create("Simulations/console", showWarnings = FALSE)

# saves any output printed to console
conPath = file.path("Simulations/console/sim_gill.txt")
conOutput <- file(conPath, open = "wt")
sink(conOutput)
sink(conOutput, type = "message")

# define simulation length and hwo often populations monitored 
T.sim = 80*52*7*24*3600
dt.sim = 7*24*3600

# function to simulate using Gillespie Algo, from SMSB by D. Wilkinson
gillespied = function(N, T.sim, dt.sim, ...){
  tt = 0
  n = T.sim%/%dt.sim
  x = N$M
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = ncol(S)
  xmat = matrix(ncol=u, nrow=n)
  i = 1
  target = 0
  repeat{
    h = N$h(x, tt, ...)
    h0 = sum(h)
    if( h0 < 1e-10 ) tt = 1e99
    else tt = tt + rexp(1, h0)
    while( tt >= target ){
      xmat[i,] = x
      i = i + 1
      target = target + dt.sim
      if( i>n ){
        C = xmat[,1]+xmat[,2]
        h = xmat[,2]/C
        return(list("h"=function() lines(ts(h, start=1/52, frequency=52), col=myBlack),
                    "C"=function() lines(ts(C, start=1/52, frequency=52), col=myBlack),
                    "Cmax"=max(C) ) 
        )
      }
    }
    j = sample(v,1,prob=h)
    x = x + S[,j]
  }
}


# generate intial copy number and mutation load
inits = function(n=1){
  C = rnorm(n,200,50) 
  h = rbeta(n,18.57,10)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )
N$h = function(x, tt, th=c(r_w=3.06e-8, r_m=3.06e-8, d_w=3.06e-8, d_m=3.06e-8, m=0)){
  with(as.list(c(x,th)),{
    return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
  })
}

# create a list of length N (no. of sims) with the same input
# but still draws initial values from init function
NN = list()
for(i in 1:N.sim){
  N.temp = N
  N.temp$M = inits()
  NN[[i]] = N.temp
}

# run the gillespie algo N.sim times using parallel::mclapply
Gillespie_Sims = lapply(NN, gillespied, T.sim, dt.sim)

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
plot(1, type='n', ylim=c(0,max(ymax)), xlim=c(0,80), main="", 
     xlab="Time (years)", ylab="Copy Number", cex.lab=1.4)
for(i in 1:N.sim){
  Gillespie_Sims[[i]][["C"]]()
}
par(op)
dev.off()
  
sink(type = "message")
sink()

















