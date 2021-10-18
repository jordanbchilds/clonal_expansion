copy_num_jl = read.delim("./Gillespie_Simulation/Julia/copy_number.txt")
mut_load_jl = read.delim("./Gillespie_Simulation/Julia/mutation_load.txt")

copy_num_r = read.table("./Gillespie_Simulation/R/copy_number.txt")
mut_load_r = read.table("./Gillespie_Simulation/R/mutation_load.txt")

Tmax = 80*365*24*3600
dt = 24*3600
t = seq(0,Tmax, by=dt)

plotter = function(x, ylim, t, title, main){
  plot(1, type="n", ylim=ylim, xlim=c(0,80), ylab=title, xlab="Time (years)")
  for(i in 1:ncol(x)){
    lines(ts(x[,i], start=0,end=80, frequency=365), lwd=1, col="black")
  }
  title(main=main, outer=T, line=-1)
}


pdf("./Gillespie_Simulation/julia_v_R.pdf", width=14,height=8.5)
op = par(mfrow=c(1,2))
plotter(copy_num_jl, ylim=c(0,max(copy_num_jl)), t, title="Copy Number", main="Julia")
plotter(mut_load_jl, ylim=c(0,1), t, title="Mutation Load", main="Julia")
plotter(copy_num_r, ylim=c(0,max(copy_num_r)), t, title="Copy Number", main="R")
plotter(mut_load_r, ylim=c(0,1), t, title="Mutation Load", main="R")
par(op)
dev.off()
