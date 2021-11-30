CN_qnt_gill_jl = read.delim("./Simulations/CN_qnt_gill_jl.txt")
ML_qnt_gill_jl = read.delim("./Simulations/ML_qnt_gill_jl.txt")

CN_qnt_gill_r = read.table("./Simulations/CN_qnt_gill_r.txt")
ML_qnt_gill_r = read.table("./Simulations/ML_qnt_gill_r.txt")

CN_ts_gill_jl = read.table("./Simulations/CN_ts_gill_jl.txt")
ML_ts_gill_jl = read.table("./Simulations/ML_ts_gill_jl.txt")

CN_ts_gill_r = read.table("./Simulations/CN_ts_gill_r.txt")
ML_ts_gill_r = read.table("./Simulations/ML_ts_gill_r.txt")

CN_qnt_tau_jl = read.table("./Simulations/CN_qnt_tau_jl.txt")
ML_qnt_tau_jl = read.table("./Simulations/ML_qnt_tau_jl.txt")

CN_qnt_tau_r = read.table("./Simulations/CN_qnt_tau_r.txt")
ML_qnt_tau_r = read.table("./Simulations/ML_qnt_tau_r.txt")

CN_ts_tau_r = read.table("./Simulations/CN_ts_tau_r.txt")
ML_ts_tau_r = read.table("./Simulations/ML_ts_tau_r.txt")

CN_ts_tau_jl = read.table("./Simulations/CN_ts_tau_jl.txt")
ML_ts_tau_jl = read.table("./Simulations/ML_ts_tau_jl.txt")

CN_qnt_abm_onehour = read.table("./Simulations/CN_qnt_abm_jl_onehour.txt")
ML_qnt_abm_onehour = read.table("./Simulations/ML_qnt_abm_jl_onehour.txt")

CN_qnt_abm_oneday = read.table("./Simulations/CN_qnt_abm_jl_oneday.txt")
ML_qnt_abm_oneday = read.table("./Simulations/ML_qnt_abm_jl_oneday.txt")



Tmax = 80*365*24*3600
dt = 24*3600
t = seq(0,Tmax, by=dt)

dir.create("Simulations", showWarnings = F)
dir.create("Simulations/PDF", showWarnings = F)

plotter = function(x, ylim, t, title, main){
  plot(1, type="n", ylim=ylim, xlim=c(0,80), main=title, 
       xlab="Time (years)", ylab="")
  lines(ts(x[,3], start=0,end=80, frequency=52), lwd=2, col="black")
  lines(ts(x[,2], start=0,end=80, frequency=52), lwd=1, col="black")
  lines(ts(x[,4], start=0,end=80, frequency=52), lwd=1, col="black")
  lines(ts(x[,1], start=0,end=80, frequency=52), lwd=1, col="darkgrey")
  lines(ts(x[,5], start=0,end=80, frequency=52), lwd=1, col="darkgrey")
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

ks_tester = function(sim_list, times, names){
  n = length(sim_list)
  p = ncol(sim_list[[1]])
  out = vector("list", length=p)
  for(t in 1:p){ # cycle though time slices
    out[[paste(times[t])]] = matrix(NA, ncol=n, nrow=n, dimnames=names)
    for(i in 1:n){ # cycle through simulations
      for(j in 1:n){
        out[[paste(times[t])]][i,j] = ks.test(sim_list[[i]][,t], sim_list[[j]][,t])$p.value
      }
    }
  }
  out
}

pdf("./Simulations/PDF/ABM_qntls.pdf", width=14, height=8.5)
par(mfrow=c(1,2))
plotter(CN_qnt_abm_jl, ylim=c(0,max(CN_qnt_gill_jl)), t, title="Copy Number", main="Agent Based Julia")
plotter(ML_qnt_abm_jl, ylim=c(0,1), t, title="Mutation Load", main="Agent Based Julia")
dev.off()


pdf("./Simulations/PDF/sim_comparison.pdf", width=14,height=8.5)
  op = par(mfrow=c(1,2))
  plotter(CN_qnt_gill_jl, ylim=c(0,max(CN_qnt_gill_jl)), t, title="Copy Number", main="Gillespie Julia")
  plotter(ML_qnt_gill_jl, ylim=c(0,1), t, title="Mutation Load", main="Gillespie Julia")
  plotter(CN_qnt_gill_r, ylim=c(0,max(CN_qnt_gill_r)), t, title="Copy Number", main="Gillespie R")
  plotter(ML_qnt_gill_r, ylim=c(0,1), t, title="Mutation Load", main="Gillespie R")
  plotter(CN_qnt_tau_r, ylim=c(0,max(CN_qnt_tau_r)), t, title="Copy Number", main="Tau Leap R")
  plotter(ML_qnt_tau_r, ylim=c(0,1), t, title="Mutation Load", main="Tau Leap R")
  plotter(CN_qnt_tau_jl, ylim=c(0,max(CN_qnt_tau_jl)), t, title="Copy Number", main="Tau Leap Julia")
  plotter(ML_qnt_tau_jl, ylim=c(0,1), t, title="Mutation Load", main="Tau Leap Julia")
  par(op)
dev.off()

time_titles = paste(1:8*10, " Years")
pdf("./Simulations/PDF/CN_ts_comparison.pdf", width=14, height=7)
  op = par(mfrow=c(2,4))
  dist_comp(list(CN_ts_gill_r, CN_ts_gill_jl, CN_ts_tau_r, CN_ts_tau_jl), time_titles,
            legend=c("Gill R", "Gill Jl", "Tau R", "Tau Jl"))
  title(main="Copy Number", outer=T, line=-1)
  par(op)
dev.off()

pdf("./Simulations/PDF/ML_ts_comparison.pdf", width=14, height=7)
  op = par(mfrow=c(2,4))
  dist_comp(list(ML_ts_gill_r, ML_ts_gill_jl, ML_ts_tau_r, ML_ts_tau_jl), time_titles,
            legend=c("Gill R", "Gill Jl", "Tau R", "Tau Jl"))
  title(main="Mutation Load", outer=T, line=-1)
  par(op)
dev.off()

ML_ks_test = ks_tester(list(ML_ts_gill_r, ML_ts_gill_jl, ML_ts_tau_r, ML_ts_tau_jl),
                    times=paste0(1:8*10, "_years"),
                    names=list(c("gill_r", "gill_jl", "tau_r", "tau_jl"),
                               c("gill_r", "gill_jl", "tau_r", "tau_jl")))

CN_ks_test = ks_tester(list(CN_ts_gill_r, CN_ts_gill_jl, CN_ts_tau_r, CN_ts_tau_jl),
                       times=paste0(1:8*10, "_years"),
                       names=list(c("gill_r", "gill_jl", "tau_r", "tau_jl"),
                                  c("gill_r", "gill_jl", "tau_r", "tau_jl")))


#### plot Gillespie in Julia and R for 80 year time slice Mutation Load
plot(1, type='n', xlab="", ylab="", main="Mutation Load", xlim=c(-0.2,1.2), ylim=c(0,2.5))
lines(density(ML_ts_gill_r[,8]), col="black", lwd=2)
lines(density(ML_ts_gill_jl[,8]), col="red", lwd=2)
legend("topleft", lty=1, lwd=2, col=1:2, legend=c("Gill R", "Gill Jl"))


