myBlack = rgb(0,0,0,0.8)
myBlue = rgb(0,0,1, 0.5)
myRed = rgb(1,0,0, 0.5)

steps = c(1,2,5,10,25,50,100,150,175,185,190,200,300,500,1000,5000)

abm_sims = list()
tau_sims = list()
for( step in steps ){
  file_temp = paste0("./Simulations/CN_qnt_abm_",step,"d.txt")
  if( file.exists(file_temp)){ abm_sims[[paste0("CN_",step,"d")]] = read.table(file_temp) }
  file_temp = paste0("./Simulations/ML_qnt_abm_",step,"d.txt")
  if( file.exists(file_temp)){ abm_sims[[paste0("ML_",step,"d")]] = read.table(file_temp) }
  
  file_temp = paste0("./Simulations/CN_qnt_tau_",step,"d.txt")
  if( file.exists(file_temp)){ tau_sims[[paste0("CN_",step,"d")]] = read.table(file_temp) }
  file_temp = paste0("./Simulations/ML_qnt_tau_",step,"d.txt")
  if( file.exists(file_temp)){ tau_sims[[paste0("ML_",step,"d")]] = read.table(file_temp) }
}

CN_qnt_gill = read.table("./Simulations/CN_qnt_gill_jl.txt")
ML_qnt_gill = read.table("./Simulations/ML_qnt_gill_jl.txt")

gill_times = read.table("./Simulations/gill_times.txt")
abm_185d_times = read.table("./Simulations/abm_185d_times.txt")
tau_200d_times = read.table("./Simulations/tau_200d_times.txt")

comp_plotter = function(x1, x2=NULL, ymax, ts, ts2=NULL, title, main){
  plot(1, type='n', xlim=range(ts), ylim=c(0,ymax),
       xlab="Time (years)", ylab=main, cex.lab=1.2, cex.main=1.4)
  qntl_lines(x1, ts, col=myBlack)
  if( !is.null(x2) & is.null(ts2) ){
    qntl_lines(x2, ts, col=rgb(0,0,1,0.5))
  } else {
    qntl_lines(x2, ts2, col=rgb(0,0,1,0.5))
  }

  title(main=title, outer=TRUE, line=-1)
}

qntl_lines = function(x, ts, col){
  if(!is.null(x)){
    lines(ts, x[,3], type='s', lwd=3, col=col)
    lines(ts, x[,2], type='s', lwd=2, col=col)
    lines(ts, x[,4], type='s', lwd=2, col=col)
    lines(ts, x[,1], type='s', lwd=1, col=col)
    lines(ts, x[,5], type='s', lwd=1, col=col)
  }
}


hour = 3600
day = 24*hour
year = 365*day

Tmax = 80*365*24*3600
dtout = 7*24*3600

pdf("./Simulations/PDF/model_comp.pdf", width=14.5, height=8)
par(mfrow=c(1,2))
for(step in steps){
  comp_plotter(CN_qnt_gill, ymax=800, ts=seq(0,Tmax,by=dtout)/year, title=paste("Time Step:",step,"Day"), 
               main="Copy Number")
  qntl_lines(abm_sims[[paste0("CN_",step,"d")]], ts=seq(0,Tmax,by=step*day)/year, col=myBlue)
  qntl_lines(tau_sims[[paste0("CN_",step,"d")]], ts=seq(0,Tmax,by=step*day)/year, col=myRed)
  legend("topleft", lty=1, col=c(myBlack, myRed, myBlue), legend=c("Gillespie","Tau","ABM"))
  
  comp_plotter(ML_qnt_gill, ymax=1, ts=seq(0,Tmax,by=dtout)/year, title="", main="Mutation Load")
  qntl_lines(abm_sims[[paste0("ML_",step,"d")]], ts=seq(0,Tmax,by=step*day)/year, col=myBlue)
  qntl_lines(tau_sims[[paste0("ML_",step,"d")]], ts=seq(0,Tmax,by=step*day)/year, col=myRed)
  legend("topleft", lty=1, col=c(myBlack, myRed, myBlue), legend=c("Gillespie","Tau","ABM"))
}
dev.off()

pdf("./Simulations/PDF/model_times.pdf", width=14.5, height=8)
par(mfrow=c(1,1))
plot(1, type='n', main="Simulation Times", xlim=c(0,0.08), ylim=c(0,1e2), ylab="",
     xlab="Time (seconds)")
lines(density(gill_times[[1]]), col="black", lwd=2)
lines(density(abm_185d_times[[1]]), col="blue", lwd=2)
lines(density(tau_200d_times[[1]]), col="red", lwd=2)
legend("topright", lty=1, col=c("black", "blue", "red"), legend=c("Gillespie", expression("ABM,"~Delta~"t = 185 days"), expression("Tau Leap,"~Delta~"t = 200 days")))
dev.off()










