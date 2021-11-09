# install.packages("smfsb")
library("smfsb")

Moraes_data = read.csv("../moraes_data.csv", header=TRUE, sep=",")

Moraes_data$Time.days. = Moraes_data$Time.days.*3600*24
names(Moraes_data) = c("Time", "Red", "Blue")

N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

N$h = function(x, error, th){
  k = th[1:5]
  Kc = th[6:7]
  th.1 = k
  if(error>=0){
    th.1[1:2] = k[1:2] + error*Kc[1]
    return( th.1*c(x,x,x[1]) )
  } else {
    th.1[1:2] = 2*k[1:2]/(1+exp(-error/Kc[2]))
    return( th.1*c(x,x,x[1]) )
  }
}

N$M = c(150,50)

# function to simulate using Gillespie Algo, from SMSB by D. Wilkinson
gillespied_mtDNA = function(N, Tdata, Tmax, dt, th, ...){
  tt = 0
  n = Tmax%/%dt
  ndata = length(Tdata)
  simT = seq(0, Tmax, by=dt)
  
  x = N$M
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = as.integer( ncol(S) )
  xmat = matrix(ncol=u, nrow=n)
  ymat = matrix(ncol=u, nrow=ndata+1)
  ymat[1,] = x
  
  i = 1L
  j = 1L
  target = 0
  C0 = sum(x)
  repeat{
    if( tt >= Tdata[j] && j<=ndata){
      # print( ymat )
      ymat[j+1,] = x
      j = j + 1L
    }
    error = C0 - sum(x)
    h = N$h(x, error, th)
    h0 = sum(h)
    if( h0 < 1e-10 ) tt = 1e99
    else tt = tt + rexp(1, h0)
    while( tt >= target){
      xmat[i,] = x
      print(x)
      i = i + 1L
      target = target + dt
      if( i>n ){
        mutation_load_fs = xmat[,2]/rowSums(xmat)
        mutation_load_ds = ymat[,2]/rowSums(ymat)
        return(list("h_full" = mutation_load_fs, "h_data" = mutation_load_ds) )
      }
    }
    k = sample(1:v,1L,prob=h)
    x = x + S[,k]
  }
}

abcRun = function(n, rprior, rdist){
  v = vector("list", 1000)
  # cl  = makeCluster(4)
  # clusterExport(cl, c("rprior", "rdist", "distance","summStats", "rmodel", "ssd",
  #                     "gillespied_mtDNA", "Moraes_data", "v", "N"))
    p = lapply(v, function(x){ rprior() })
    d = lapply(p, rdist)
  # stopCluster(cl)
  
  pm = t(sapply(p, identity))
  if(dim(pm[1]) ==1)
    pm = as.vector(pm)
  dm = t(sapply(d, identity))
  if(dim(dm)[1] == 1)
    dm = as.vecto(dm)
  list(param=pm, dist=dm)
}

rprior = function(){ c(rbeta(5, 2,50), rbeta(1, 2,50), runif(1)*1000) }

rmodel = function(th){ gillespied_mtDNA(N, Tdata=Moraes_data$Time, Tmax=50*24*3600, dt=3600, th)[["h_data"]] }

summStats = identity

ssd = summStats(Moraes_data$Blue)

distance = function(s){
  diff = s - ssd
  sqrt(sum(diff*diff))
}

rdist = function(th){ distance(summStats(rmodel(th))) }

abcRun(100, rprior, rdist)








