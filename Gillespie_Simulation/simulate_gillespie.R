library("parallel")

myBlack = rgb(0,0,0, alpha=15, max=255)

# create directory to save pdf output
# full simulation results not saved
dir.create("Simulations", showWarnings=F)

# define simulation length and how often populations monitored 
T.sim = 80*365*24*3600
dt.sim = 24*3600
N.sim = 10000 # nuo. of simulations

# function to simulate using Gillespie Algo, from SMSB by D. Wilkinson
gillespied_mtDNA = function(N, T.sim, dt.sim, ...){
  tt = 0
  n = T.sim%/%dt.sim
  x = N$M
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = as.integer( ncol(S) )
  xmat = matrix(ncol=u, nrow=n)
  i = 1L
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
      i = i + 1L
      target = target + dt.sim
      if( i>n ){
        return( xmat )
      }
    }
    j = sample(1:v,1L,prob=h)
    x = x + S[,j]
  }
}

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

# define simulation parameters
# N$h = function(x, th=c(3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0)){
#   th*rep(x,length.out=5)
# }

N$h = function(x, error, K_c=c(8.99e-9, 1/500),
               th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025, d_w=3.06e-8, d_m=3.06e-8, m=3.06e-9)){
  r = c("r_w","r_m")
  th.1 = th
  if(error>=0){
    th.1[r] = th[r] + error*K_c[1]
    return( th.1*c(x,x,x[1]) )
  } else {
    th.1[r] = 2*th[r]/(1+exp(-error*K_c[2]))
    return( th.1*c(x,x,x[1]) )
  }
}

# generate intial copy number and mutation load
inits = function(n=1){
  C = rnorm(n,200,50)
  h = rbeta(n,10,18.57)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}
# create a list of length N (no. of sims) with the same input
# but still draws initial values from init function
gen_N = function(N, N.sim){
  NN = list()
  for(i in 1:N.sim){
    N.temp = N
    # const initial populations
    N.temp$M = inits()
    # N.temp$M = inits()
    NN[[i]] = N.temp
  }
  return(NN)
}

NN = gen_N(N, N.sim)

gillr_time = system.time({
  cl  = makeCluster(5) 
  clusterExport(cl, c("gillespied_mtDNA", "NN"))
  Gillespie_sims = parLapply(cl, NN, gillespied_mtDNA, T.sim, dt.sim)
  stopCluster(cl)
})
### 1000 simple (equal rates, no mutations) simulations: 79.516 seconds 


#### convert raw population data to copy number and mutation load
raw_to_hC = function(sim){
  n = nrow(sim)
  copy_num = rowSums(sim)
  mut_load = sim[,2]/copy_num
  mut_load[is.na(mut_load)] = 0.0

  return(cbind(copy_num, mut_load))
}

CH_sim = lapply(Gillespie_sims, raw_to_hC)

slice_dist = function(sims, t, T.sim, dt.sim){
  t_tot = seq(dt.sim,T.sim, by=dt.sim)
  sim_t = list()
  sim_t[[1]] = matrix(NA, nrow=length(sims), ncol=length(t))
  sim_t[[2]] =  sim_t[[1]]
  for(i in 1:length(sims)){
    for(j in 1:length(t)){
      sim_t[[1]][i,j] = sims[[i]][t_tot==t[j],1]
      sim_t[[2]][i,j] = sims[[i]][t_tot==t[j],2]
    }
  }
  return(sim_t)
}

dist_sims = slice_dist(CH_sim, 1:8*10*365*24*3600, T.sim, dt.sim)

quantiles = function(sims, p){
  Nsim = length(sims)
  n = nrow(sims[[1]])
  out = list()
  out$copy_num = matrix(NA, nrow=n,ncol=length(p))
  out$mut_load = matrix(NA, nrow=n,ncol=length(p))
  for(t in 1:n){
    copy_t = vector("numeric", Nsim)
    mut_t  = vector("numeric", Nsim)
    for(i in 1:Nsim){
      copy_t[i] = sims[[i]][t,1]
      mut_t[i] = sims[[i]][t,2]
    }
    out$copy_num[t,] = quantile(copy_t, p)
    out$mut_load[t,] = quantile(mut_t, p)
  }
  return( out )
}
  
quant_sim = quantiles(CH_sim, p=c(0.025,0.25,0.5,0.75,0.975))

write.table(quant_sim$copy_num, file="Simulations/CN_qnt_gill_r.txt",
            row.names=F, col.names=F)
write.table(quant_sim$mut_load, file="Simulations/ML_qnt_gill_r.txt",
            row.names=F, col.names=F)
write.table(dist_sims[[1]], file="Simulations/CN_ts_gill_r.txt",
            row.names=F, col.names=F)
write.table(dist_sims[[2]], file="Simulations/ML_ts_gill_r.txt",
            row.names=F, col.names=F)

# CN_ymax = numeric(N.sim)
# for(i in 1:N.sim){
#   CN_ymax[i] = max(CH_sim[[i]][,"copy_num"])
# }
# pdf("Simulations/PDF/gill_complex_sim.pdf", width=14, height=8.5)
# par(mfrow=c(1,2))
# plot(1, type='n', xlab="Time (years)", ylab="", main="Copy Number", 
#      xlim=c(0,80), ylim=c(0,max(CN_ymax)))
# for(i in 1:N.sim){
#   lines(ts(CH_sim[[i]][,"copy_num"], start=0, end=80, frequency=52), col=myBlack)
# }
# plot(1, type='n', xlab="Time (years)", ylab="", main="Mutation Load", 
#      xlim=c(0,80), ylim=c(0,1))
# for(i in 1:N.sim){
#   lines(ts(CH_sim[[i]][,"mut_load"], start=0, end=80, frequency=52), col=myBlack)
# }
# dev.off()



