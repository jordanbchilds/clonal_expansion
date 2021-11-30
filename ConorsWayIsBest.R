# Simulating two schemes for irregular sampling from regular timeseries
deltasim = 0.1
deltarec = 1/3
tmax = 10.0

simtimes = seq(0,tmax,0.1)
rectargets = seq(0,tmax,deltarec)

recindex = 1
nextrec = rectargets[1]
actualrecs = c()
actualrecs_index = c()

J.rec = c()
J.target = rectargets[1]

for (st in simtimes){
  # What I think I'm doing
  if( st>=J.target ){
    J.rec = c(J.rec, st)
    J.target = J.target + deltarec
  }
  # Jordan method
  if(st>=nextrec){
    actualrecs = c(actualrecs,st)
    nextrec = st + deltarec
  }
  # Conor method (better)
  if(st>=rectargets[recindex]){
    actualrecs_index = c(actualrecs_index,st)
    recindex = recindex+1
  }
}

# These points should all lie on the red line, but...
trim = min(c(length(actualrecs),length(rectargets),length(actualrecs_index)))
plot(rectargets[1:trim],actualrecs[1:trim],xlim=c(0,tmax),ylim=c(0,tmax),cex=1.5,pch=16)
points(rectargets[1:trim],actualrecs_index[1:trim],col="blue",cex=1.5,pch=16)
points(rectargets[1:trim], J.rec[1:trim],col="green",cex=1.5,pch=16)
abline(a=0,b=1,col="red",lwd=3)
legend("topleft",c("Jordan's intentions", "Jordan","Conor","Perfection"),col=c("green", "black","blue","red"),pch=16)


## errors
















