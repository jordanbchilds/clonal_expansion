myBlack = rgb(0,0,0, alpha=15, max=255)

# create directory to save pdf output
# full simulation results not saved
dir.create("Gillespie_Simulation", showWarnings=F)

# define simulation length and how often populations monitored 
T.sim = 80*365*24*3600
dt.sim = 24*3600
N.sim = 1000 # nuo. of simulations

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
    h = N$h(x)
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

# generate intial copy number and mutation load
inits = function(n=1){
  C = rnorm(n,200,50)
  h = rbeta(n,10,18.57)
  inits = round( c( C*(1-h), C*h ) )
  return( inits )
}

# define simulation parameters
N = list(Pre=matrix(c(1,0,0,1,1,0,0,1,1,0), byrow=TRUE, ncol=2, nrow=5),
         Post=matrix(c(2,0,0,2,0,0,0,0,1,1), byrow=TRUE, ncol=2, nrow=5) )

# define simulation parameters
N$h = function(x, th=c(3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0)){
  th*rep(x,length.out=5)
}

# N$h = function(x, error, K_c=8.99e-9,
#                th=c(r_w=3.06e-8*0.975, r_m=3.06e-8*1.025, d_w=3.06e-8, d_m=3.06e-8, m=0)){
#   r = c("r_w","r_m")
#   th.1 = th
#   if(error>=0){
#     th.1[r] = th[r] + error*K_c
#     return( th.1*c(x,x,x[1]) )
#   } else {
#     th.1[r] = 2*th[r]/(1+exp(-error/500))
#     return( th.1*c(x,x,x[1]) )
#   }
# }

# create a list of length N (no. of sims) with the same input
# but still draws initial values from init function
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


######
###### SINGLE RUN 

N$M = c(100,100) # add initial conditions for test
test = gillespied_mtDNA(N, T.sim, dt.sim)

######
######

NN = gen_N(N, N.sim)
cl  = makeCluster(4) 
clusterExport(cl, c("gillespied_mtDNA", "NN"))
Gillespie_sims = parLapply(cl, NN, gillespied_mtDNA, T.sim, dt.sim)
stopCluster(cl)

#### convert raw population data to copy number and mutation load
raw_to_summ = function(sim){
  n = nrow(sim)
  copy_num = rowSums(sim)
  mut_load = sim[,2]/copy_num
  mut_load[is.na(mut_load)] = 0.0

  return(cbind(copy_num, mut_load))
}

CH_sim = lapply(Gillespie_sims, raw_to_summ)

quantiles = function(sims, p){
  Nsim = length(sims)
  n = nrow(sims[[1]])
  out = list()
  out$copy_num = matrix(NA, nrow=n,ncol=length(p))
  out$mut_load = matrix(NA, nrow=n,ncol=length(p))
  for(t in 1:n){
    copy_t = vector("numeric", Nsim)
    mut_t = vector("numeric", Nsim)
    for(i in 1:Nsim){
      copy_t[i] = sims[[i]][t,1]
      mut_t[i] = sims[[i]][t,2]
    }
    out$copy_num[t,] = quantile(copy_t, p)
    out$mut_load[t,] = quantile(mut_t, p)
  }
  return( out )
}
  
quant_sim = quantiles(CH_sim, p=c(0.025,0.1,0.5,0.9,0.975))

write.table(quant_sim$copy_num, file="Gillespie_Simulation/copy_number_quantiles_r.txt",
            row.names=F, col.names=F)
write.table(quant_sim$mut_load, file="Gillespie_Simulation/mutation_load_quantiles_r.txt",
            row.names=F, col.names=F)












