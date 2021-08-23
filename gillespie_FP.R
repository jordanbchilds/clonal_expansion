plotfuns = list()
Nsims = 100
T.life = 80*52*7*24*3600
dt.week = 7*24*3600

myDarkGrey = rgb(169,169,169, alpha=50, max=255)

inits = function(n=1){
  C = rnorm(n,200,50) 
  h= rbeta(n,18.57,10)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

N$h = function(x, tt, th=c(r_w=8e-6, r_m=8e-6, d_w=2e-8, d_m=2e-8, m=0)){
  with(as.list(c(x,th)),{
    return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
  })
}

# Visualising a 1D random walk, for example
gillespied = function(N, T=T.life, dt=dt.week, ...){
  # Do some like real heavy computing man
  N$M = inits(n=1)
  tt = 0
  n = T%/%dt
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
      target = target + dt
      if( i > n ){
        xts = ts(xmat, start=0, deltat=dt)
        return( list("plotter" = function() plotGillespie(xts), 
                    "sim"=xts) )
      }
    }
    j = sample(v,1,prob=h)
    x = x + S[,j]
  }
  # Return a new function DEFINITION without executing plotwalk function
}

###############################################
######## Singular GILLESPIES
###############################################
# 
# plotGillespie = function(sim, T=T.life, dt=dt.week) {
#   t = seq(from=1, to=T/dt, by=1)
#   M = sim[,2]
#   W = sim[,1]
#   Mload = M/(M+W)
#   plot(t, Mload, ylab="Mutation Load", type='s', col="darkgrey", ylim=c(0,1), 
#        main=paste0("Simulation ", sprintf("%04d",i)))
#   plot(t, M+W, ylab="C", type='s', col="darkgrey",
#        main=paste0("Simulation ", sprintf("%04d",i)))
# }
# 
# for(i in 1:Nsims){
#   # Without actually executing the plotwalk function, store its output inside another function in e.g. a list
#   plotfuns[[i]] = gillespied(N)[["plotter"]]
# }
# 
# # Four plots per page
# pdf("GillespieSim.pdf")
# op = par(mfrow=c(2,2), mai=c(0.75,0.75,0.5,0.5))
# for(i in 1:Nsims) walk = plotfuns[[i]]()
# par(op)
# dev.off()


###############################################
######## MANY GILLESPIES
###############################################


simGillespies = function(Nsims, T=T, dt=dt){
  n = T%/%dt
  sim_array = array(NA, dim=c(n,2,Nsims))
  for(i in 1:Nsims) sim_array[,,i] = gillespied(N)[["sim"]]
  return( list("plotter" = function() plotGillespies(sim_array), "sim"=sim_array ) ) 
}

plotGillespies = function(sim_array, T=T, dt=dt){
  t = seq(from=dt, to=T, by=dt)
  Nsims = dim(sim_array)[3]
  x.lim = range(sim_array[,1,])
  y.lim = range(sim_array[,2,])
  
  plot(sim_array[,1,], sim_array[,2,], type='n', 
       xlim=c(0,T), ylim=c(0,1), xlab="Time", ylab="h", main="Mutation Load")
  for(i in 1:Nsims){
    W = sim_array[,1,i]
    M = sim_array[,2,i]
    h = M/(M+W)
    lines(t, h, type='s', col=myDarkGrey)
  }

  plot(sim_array[,1,], sim_array[,2,], type='n', 
       xlim=c(0,T), ylim=y.lim, xlab="Time", ylab="C", main="Copy Number") 
  for(i in 1:Nsims){
    W = sim_array[,1,i]
    M = sim_array[,2,i]
    C = M+W
    lines(t, C, col=myDarkGrey, type='s')
  }
}

plotSims = simGillespies(Nsims)[["plotter"]]

pdf("ManyGillespies.pdf", width=14, height=8.5)
op = par(mfrow=c(1,2))
plotSims()
par(op)
dev.off()






