library("parallel")
##########
## POSSION LEAP PRCOESS
##########
myBlack = rgb(0,0,0, alpha=15, max=255)

dir.create("Simulations", showWarnings=F)

Tmax = 80*365*24*3600
dt = 24*3600
N.sim = 10000 # no. of simulations

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

# define simulation parameters
N$h = function(x, th=c(3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0)){
  th*rep(x,length.out=5)
}

tauleap_simulation = function(N, Tmax, dt){
  n = Tmax%/%dt
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = as.integer( ncol(S) )
  xmat = matrix(ncol=u, nrow=n)
  xmat[1,] = N$M
  
  i = 1L
  while( i<n ){
    h = N$h(xmat[i,])
    if( sum(h)<1e-10 ){
      xmat[(i+1):n,] = c(0,0)
      return(xmat)
    }
    r = rpois(v,h*dt)
    xmat[i+1, ] = xmat[i,] + S%*%r
    i = i+1L
  }
  xmat
}

gen_N = function(N, N.sim){
  NN = list()
  for(i in 1:N.sim){
    N.temp = N
    # const initial populations
    N.temp$M = c(100, 100)
    # N.temp$M = inits()
    NN[[i]] = N.temp
  }
  return(NN)
}

NN = gen_N(N, N.sim)
###
### Nsim simulations
tauLeap_time = system.time({
  cl  = makeCluster(12)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw = parLapply(cl, NN, tauleap_simulation, Tmax, dt)
  stopCluster(cl)
})
# 1000 simulations takes 48.555 seconds


###
### convert to copy number and mutation load
raw_to_hC = function(sim){
  n = nrow(sim)
  copy_num = rowSums(sim)
  mut_load = sim[,2]/copy_num
  mut_load[is.na(mut_load)] = 0.0
  
  return(cbind(copy_num, mut_load))
}
tauLeap_sim = lapply(tauLeap_raw, raw_to_hC)

###
### calculate quantiles
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
tauLeap_qntl = quantiles(tauLeap_sim, p=c(0.025,0.25,0.5,0.75,0.975))

###
### save output
write.table(tauLeap_qntl[[1]], "./Simulations/CN_qnt_tau_r.txt",row.names=F, col.names=F)
write.table(tauLeap_qntl[[2]], "./Simulations/ML_qnt_tau_r.txt",row.names=F, col.names=F)

###
### distribution at given times 
slice_dist = function(sims, t, Tmax, dt){
  t_tot = seq(dt, Tmax, by=dt)
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

sims_dist = slice_dist(tauLeap_sim, 1:8*10*365*24*3600, Tmax, dt)

###
### save time slice distributions
write.table(sims_dist[[1]], file="Simulations/CN_ts_tau_r.txt",
            row.names=F, col.names=F)
write.table(sims_dist[[2]], file="Simulations/ML_ts_tau_r.txt",
            row.names=F, col.names=F)










