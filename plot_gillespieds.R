
myDarkGrey = rgb(169,169,169, alpha=50, max=255)

args = commandArgs(trailingOnly = TRUE)
if( length(args)!=0 ){
  N = args
} else {
  N = NULL
}

if( length(N)>0 ){
  sims = sprintf( "Gillespie_sim", sprintf("%04d",1:N),".pdf")
  
  # list to store output as functions
  h_list = list()
  C_list = list()
  Cmax = double(lengths(N))
  
  # function to plot h (mutation load)
  h_line = function(h){
    t = seq(0,80, by=1/52)
    return( function(h) lines(t, h, col="myDarkGrey"))
  }
  
  # function to plot C (copy number)
  C_line = function(C){
    t = seq(0,80, by=1/52)
    return( function(C) lines(t, C, col="myDarkGrey"))
  }
  
  # loop through datasets and store h_line and C_line functions
  for(i in 1:N){
    df = read.delim(paste0("Simulations/",sims[i]), header=FALSE, stringsAsFactors=FALSE)
    C = df[,1] + df[,2]
    h = df[,2]/C
    C_list[[i]] = C_line(C)
    Cmax[i] = max(C)
    h_list[[i]] = h_line(h)
  }
  
  plot(NA, ylim=c(0,1), xlim=c(0,80), main="Mutation Load",
       xlab="Time (years)", ylab="h")
  for(i in 1:N){
    C_list[[i]]()
  }
  
  plot(NA, ylim=c(0,max(Cmax)), xlim=c(0,80), main="Copy Number", 
       xlab="Time (years)", ylab="C")
  for(i in 1:N){
    h_list[[i]]()
  }
} else {
  error("Number of simulations required.")
}










