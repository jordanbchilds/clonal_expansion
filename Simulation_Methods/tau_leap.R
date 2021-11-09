library("parallel")

CN_gill = read.table("./Simulations/CN_qnt_gill_r.txt", header=FALSE)
ML_gill = read.table("./Simulations/ML_qnt_gill_r.txt", header=FALSE)

##########
## POSSION LEAP PRCOESS
##########
myBlack = rgb(0,0,0, alpha=125, max=255)
myBlue = rgb(0,0,255, alpha=125, max=255)

qntl_lines = function(qntls, Tend, freq, col){
  lines(ts(qntls[,3], start=0, end=Tend, frequency=freq), lwd=3, col=col) 
  lines(ts(qntls[,2], start=0, end=Tend, frequency=freq), lwd=2, col=col)
  lines(ts(qntls[,4], start=0, end=Tend, frequency=freq), lwd=2, col=col)
  lines(ts(qntls[,1], start=0, end=Tend, frequency=freq), lwd=1, col=col)
  lines(ts(qntls[,5], start=0, end=Tend, frequency=freq), lwd=1, col=col)
}

qntl_plotter = function(qntls, ylim, Tend, freq, title, main, comp_qntls=NULL, cols){
  plot(1, type="n", ylim=ylim, xlim=c(0,80), main=title, 
       xlab="Time (years)", ylab="")
  if( is.null(comp_qntls) ){
    qntl_lines(qntls, Tend, freq=52, cols[1])
  } else {
    qntl_lines(qntls, Tend, freq=52, cols[1])
    qntl_lines(comp_qntls, Tend, freq=52, cols[2])
  }
  title(main=main, outer=T, line=-1)
}

dist_comp = function(sim_list, titles, legend){
  n = length(sim_list) # no. of simulations being compared
  p = length(titles) # no. of time slices per simulation
  # eqiuv. p = ncol(sim_list[[i]]) for i=1:p
  
  for(i in 1:p){ # cycle through time slices
    xrange = matrix(NA, nrow=n, ncol=2) # matrix to store min and max values in x axis
    yrange = double(n) # vector to store max y vals (min always zero)
    for(j in 1:n){ # cycle through simulations
      xrange[j,] = as.numeric(range(sim_list[[j]][,i])) # range of simulation j and time slice i
      yrange[j] = max(density(sim_list[[j]][,i])$y) # max density for simulation j and time slice i
    }
    xlim = range(xrange)
    ylim = c(0,max(yrange))
    plot(1, type='n', xlim=xlim, ylim=ylim, xlab="", ylab="Density", 
         main=titles[i] )
    for(k in 1:n){
      lines(density(sim_list[[k]][,i]), col=k, lty=1)
    }
    legend("topright", legend=legend, col=1:n, lty=1)
  }
  
}

dir.create("Simulations", showWarnings=F)

Tmax = 80*365*24*3600
dt = 24*3600
N.sim = 5000 # no. of simulations

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

# define simulation parameters
N$h = function(x, th=c(3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0)){
  x[x<0] = 0
  th*rep(x,length.out=5)
}

# N$h = function(x, error, K_c=c(8.99e-9, 1/500),
#                th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025, d_w=3.06e-8, d_m=3.06e-8, m=3.06e-9)){
#   x[x<0] = 0
#   r = c("r_w","r_m")
#   th.1 = th
#   if(error>=0){
#     th.1[r] = th[r] + error*K_c[1]
#     return( th.1*c(x,x,x[1]) )
#   } else {
#     th.1[r] = 2*th[r]/(1+exp(-error*K_c[2]))
#     return( th.1*c(x,x,x[1]) )
#   }
# }

tauleap_simulation = function(N, Tmax, dt){
  n = Tmax%/%dt
  S = t(N$Post - N$Pre)
  u = nrow(S)
  v = as.integer( ncol(S) )
  xmat = matrix(ncol=u, nrow=n)
  xmat[1,] = N$M
  C0 = sum(xmat[1,])
  
  i = 1L
  while( i<n ){
    error = C0 - sum(xmat[i,])
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

# taur_time = system.time({
#   lapply(NN, tauleap_simulation, Tmax, dt)
# })
# simulations takes ~230 seconds

###
### Nsim simulations
tauLeap_time = system.time({
  cl  = makeCluster(20)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw = parLapply(cl, NN, tauleap_simulation, Tmax, dt)
  stopCluster(cl)
})
# 1000 simulations takes 48.555 seconds
# 1000 simulations with control hazard function: 348.89

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
slice_dist = function(sims, t, Tmax, delta.t){
  t_tot = seq(dt, Tmax, by=delta.t)
  t_dash = t_tot + delta.t
  
  sim_t = list()
  sim_t[[1]] = matrix(NA, nrow=length(sims), ncol=length(t))
  sim_t[[2]] =  sim_t[[1]]
  for(i in 1:length(sims)){
    for(j in 1:length(t)){
      sim_t[[1]][i,j] = sims[[i]][which(t_tot<t[j] & t[j]<=t_dash) ,1]
      sim_t[[2]][i,j] = sims[[i]][which(t_tot<t[j] & t[j]<=t_dash),2]
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


###
### test different delta t
## simulations
tauLeap_time_2dt = system.time({
  cl  = makeCluster(20)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw_2dt = parLapply(cl, NN, tauleap_simulation, Tmax, dt*2)
  stopCluster(cl)
})
# 1000 simulations, 4 cores: 33 seconds

tauLeap_time_3dt = system.time({
  cl  = makeCluster(20)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw_3dt = parLapply(cl, NN, tauleap_simulation, Tmax, dt*3)
  stopCluster(cl)
})
# 1000 simulations, 4 cores: 23 seconds

tauLeap_time_4dt = system.time({
  cl  = makeCluster(20)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw_4dt = parLapply(cl, NN, tauleap_simulation, Tmax, dt*4)
  stopCluster(cl)
})
# 1000 simulations, 4 cores: 17 seconds

tauLeap_time_5dt = system.time({
  cl  = makeCluster(20)
  clusterExport(cl, c("tauleap_simulation", "NN"))
  tauLeap_raw_5dt = parLapply(cl, NN, tauleap_simulation, Tmax, dt*5)
  stopCluster(cl)
})
# 1000 simulations, 4 cores: 14 seconds

# simulation post processing
tauLeap_sim_2dt = lapply(tauLeap_raw_2dt, raw_to_hC)
tauLeap_qntl_2dt = quantiles(tauLeap_sim_2dt, p=c(0.025,0.25,0.5,0.75,0.975))
sims_dist_2dt = slice_dist(tauLeap_sim_2dt, t=1:8*10*365*24*3600, Tmax, delta.t=dt*2)

tauLeap_sim_3dt = lapply(tauLeap_raw_3dt, raw_to_hC)
tauLeap_qntl_3dt = quantiles(tauLeap_sim_3dt, p=c(0.025,0.25,0.5,0.75,0.975))
sims_dist_3dt = slice_dist(tauLeap_sim_3dt, 1:8*10*365*24*3600, Tmax, dt*3)

tauLeap_sim_4dt = lapply(tauLeap_raw_4dt, raw_to_hC)
tauLeap_qntl_4dt = quantiles(tauLeap_sim_4dt, p=c(0.025,0.25,0.5,0.75,0.975))
sims_dist_4dt = slice_dist(tauLeap_sim_4dt, 1:8*10*365*24*3600, Tmax, dt*4)

tauLeap_sim_5dt = lapply(tauLeap_raw_5dt, raw_to_hC)
tauLeap_qntl_5dt = quantiles(tauLeap_sim_5dt, p=c(0.025,0.25,0.5,0.75,0.975))
sims_dist_5dt = slice_dist(tauLeap_sim_5dt, 1:8*10*365*24*3600, Tmax, dt*5)

pdf("Simulations/PDF/TauLeap_dt_test.pdf", width=14, height=8.5)
par(mfrow=c(1,2))
  qntl_plotter(tauLeap_qntl[[1]], ylim=c(0,600), Tend=80, freq=52, title="Copy Number", main="dt = 1 day", comp_qntls=CN_gill, cols=c(myBlack, myBlue) )
  qntl_plotter(tauLeap_qntl[[2]], ylim=c(0,1), Tend=80, freq=52, title="Mutation Load", main="dt = 1 day", comp_qntls=ML_gill, cols=c(myBlack, myBlue)  )
  qntl_plotter(tauLeap_qntl_2dt[[1]], ylim=c(0,600), Tend=80, freq=52, title="Copy Number", main="dt = 2 days", comp_qntls=CN_gill, cols=c(myBlack, myBlue)   )
  qntl_plotter(tauLeap_qntl_2dt[[2]], ylim=c(0,1), Tend=80, freq=52,title="Mutation Load", main="dt = 2 days", comp_qntls=ML_gill, cols=c(myBlack, myBlue)   )
  qntl_plotter(tauLeap_qntl_3dt[[1]], ylim=c(0,600), Tend=80, freq=52,title="Copy Number", main="dt = 3 days", comp_qntls=CN_gill, cols=c(myBlack, myBlue)   )
  qntl_plotter(tauLeap_qntl_3dt[[2]], ylim=c(0,1), Tend=80, freq=52,title="Mutation Load", main="dt = 3 days", comp_qntls=ML_gill, cols=c(myBlack, myBlue)   )
  qntl_plotter(tauLeap_qntl_4dt[[1]], ylim=c(0,600), Tend=80, freq=52,title="Copy Number", main="dt = 4 days",  comp_qntls=CN_gill, cols=c(myBlack, myBlue)  )
  qntl_plotter(tauLeap_qntl_4dt[[2]], ylim=c(0,1), Tend=80, freq=52,title="Mutation Load", main="dt = 4 days",  comp_qntls=ML_gill, cols=c(myBlack, myBlue)  )
  qntl_plotter(tauLeap_qntl_5dt[[1]], ylim=c(0,600), Tend=80, freq=52,title="Copy Number", main="dt = 5 days",  comp_qntls=CN_gill, cols=c(myBlack, myBlue)  )
  qntl_plotter(tauLeap_qntl_5dt[[2]], ylim=c(0,1), Tend=80, freq=52,title="Mutation Load", main="dt = 5 days",  comp_qntls=ML_gill, cols=c(myBlack, myBlue)  )
dev.off()

pdf("Simulations/PDF/CN_Tau_dt_test_ts.pdf", width=14, height=8.5)
  par(mfrow=c(2,4))
  dist_comp(list(sims_dist[[1]], sims_dist_2dt[[1]], sims_dist_3dt[[1]], sims_dist_4dt[[1]], sims_dist_5dt[[1]]),
            titles=paste(1:8*10, "Years"), 
            legend=c("1 day", "2 days", "3 days", "4 days", "5 days"))
  title(main="Copy Number", outer=T, line=-1)
dev.off()

pdf("Simulations/PDF/ML_Tau_dt_test_ts.pdf", width=14, height=8.5)
  par(mfrow=c(2,4))
  dist_comp(list(sims_dist[[2]], sims_dist_2dt[[2]], sims_dist_3dt[[2]], sims_dist_4dt[[2]], sims_dist_5dt[[2]]),
          titles=paste(1:8*10, "Years"), 
          legend=c("1 day", "2 days", "3 days", "4 days", "5 days"))
  title(main="Muation Load", outer=T, line=-1)
dev.off()






