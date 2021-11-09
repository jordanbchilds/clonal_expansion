
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
  ymat = matrix(ncol=u, nrow=ndata)
  
  i = 1L
  j = 1L
  target = 0
  C0 = sum(x)
  repeat{
    if( tt>=Tdata[j] && j<=ndata ){
      ymat[j,] = x
      j = j + 1L
    }
    error = C0 - sum(x)
    h = N$h(x, error, th)
    h0 = sum(h)
    if( h0 < 1e-10 ) tt = 1e99
    else tt = tt + rexp(1, h0)
    while( tt >= target){
      xmat[i,] = x
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

rprior = function(){ c(runif(5,0,1e-4), runif(1,0, 1e-4), runif(1,200,800)) }

distance = function(data.star, data){
  diff = data.star - data
  sqrt(sum(diff*diff))
}

abcrej = function(n.iters, N, data, Tdata, epsilon){
  sim.data = vector("list", n.iters)
  theta.out = matrix(ncol=length( rprior()), nrow=n.iters)
  count = 1
  
  while( count<=n.iters ){
    print(count)
    theta.star = rprior()
    data.star = gillespied_mtDNA(N, Tdata, Tmax=50*24*3600, dt=3600, theta.star)
    
    dist = distance(data.star[["h_data"]], data)
    
    if( dist<epsilon ){
      theta.out[count, ] = theta.star
      sim.data[[count]] = data.star[["h_full"]]
      count = count + 1
    } 
  }
  return(list("posterior"=theta.out, "simulation"=sim.data))
}

tt = abcrej(n.iters=10, N, data=Moraes_data$Blue, Tdata=Moraes_data$Time, epsilon=2)

plot(NA, xlim=c(0,50), ylim=c(0,1), ylab="Mutation Load", xlab="Time (days)")
for(sim in tt[[2]]){
  lines(ts(sim, start=0, end=50, frequency=24) )
}


gillespied_mtDNA(N, Tdata=1e99, Tmax=80*365*24*3600, dt=24*3600, 
                 th=c(3.06e-8,3.06e-8,3.06e-8,3.06e-8,0.0,8.99e-9,500))










