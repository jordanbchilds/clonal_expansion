myGrey = rgb(169,169,169, alpha=50, max=255)

gillespie = function(N, n, ...){
  tt = 0 
  x = N$M # M: the initial state
  # Pre: LHS matrix
  # Post: RHS matrix
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = ncol(S)
  tvec = vector("numeric", n)
  xmat = matrix(ncol=u, nrow=n+1)
  xmat[1,] = x
  for(i in 1:n){
    h = N$h(x,tt,...)
    tt = tt + rexp(1,sum(h))
    j = sample(v,1,prob=h)
    x = x + S[,j]
    tvec[i] = tt
    xmat[i+1,] = x
  }
  return(list(t=tvec, x=xmat))
}


N = list(M=c(100,100), Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

N$h = function(x, tt, th=c(r_w=2, r_m=2, d_w=2, d_m=2, m=0.1)){
  with(as.list(c(x,th)),{
    return(c(r_w*x[1], r_m*x[2], d_w*x[1], d_m*x[2], m*x[1]))
  })
}

sim = gillespie(N, n=10000)
t = c(0,sim$t)
M = sim$x[,2]
W = sim$x[,1]
Mprop = M/(M+W)
par(mfrow=c(1,2))
plot(t, Mprop, type='s', col='darkgrey', ylim=c(0,1))
plot(t, M+W, type='s', col='darkgrey')


gillespied = function(N, T=100, dt=1, ...){
  tt=0
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
      if( i > n ) return(ts(xmat, start=0, deltat=dt))
    }
    j = sample(v,1,prob=h)
    x = x + S[,j]
  }
}

sim = gillespied(N)
t = seq(1,100, by=1)
M = sim[,2]
W = sim[,1]
Mprop = M/(M+W)
par(mfrow=c(1,2))
plot(t, Mprop, type='s', col='darkgrey', ylim=c(0,1))
plot(t, M+W, type='s', col='darkgrey')






