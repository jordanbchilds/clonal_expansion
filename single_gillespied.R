args = commandArgs(trailingOnly = TRUE)
if( length(args)==0 ){
  filePath = file.path("Gillespie_sim.pdf")
} else {
  filePath = file.path(paste0("Gillespie_sim",sprintf("%04s", args), ".txt") ) 
}

dir.create("Simulations", showWarnings = FALSE)

T.life = 80*52*7*24*3600
dt.week = 7*24*3600

myDarkGrey = rgb(169,169,169, alpha=50, max=255)

# random draws for intial values, distributions taken from Joe's diss
inits = function(n=1){
  C = rnorm(n,200,50) 
  h = rbeta(n,18.57,10)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

N$h = function(x, tt, th=c(r_w=8e-6, r_m=8e-6, d_w=8e-6, d_m=8e-6, m=0)){
  with(as.list(c(x,th)),{
    return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
  })
}

# Visualising a 1D random walk, for example
gillespied = function(N, T=T.life, dt=dt.week, ...){
  # Some like real heavy computing man
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
        xts = ts(xmat)
        return( xts )
      }
    }
    j = sample(v,1,prob=h)
    x = x + S[,j]
  }
  # Return a new function DEFINITION without executing plotwalk function
}

###############################################
######## run GILLESPIES
###############################################

Gillespie = gillespied(N)

plot(Gillespie[,2])
plot(Gillespie[,1])

C = Gillespie[,1] + Gillespie[,2]
h = Gillespie[,2]/C

write.table(Gillespie, file=filePath)
  
  
  
  
  
  
  