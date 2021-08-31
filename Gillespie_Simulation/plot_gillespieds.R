
myDarkGrey = rgb(169,169,169, alpha=20, max=255)

args = commandArgs(trailingOnly = TRUE)
if( length(args)!=0 ){
  N = args
} else {
  N = NULL
}

dir.create("Simulations", showWarnings = FALSE)
dir.create("Simulations/PDF", showWarnings = FALSE)

N=500
if( length(N)>0 ){
  sims = paste0( "Gillespie_sim",str_pad(as.character(1:N), width=4 ),".txt")
  
  # list to store output as functions
  h_list = list()
  C_list = list()
  Cmax = double(lengths(N))
  
  # function to plot h (mutation load)
  h_line = function(h){
    hh = ts(h, start=1/52, frequency=52)
    return( list("plotter"=function() lines(hh, col=myDarkGrey), "ts"=hh ) )
  }
  
  # function to plot C (copy number)
  C_line = function(C){
    CC = ts(C, star=1/52, frequency=52)
    return( list("plotter"=function() lines(CC, col=myDarkGrey), "ts"=CC ) )
  }
  
  # loop through datasets and store h_line and C_line functions
  for(i in 1:N){
    if( file.exists(paste0("Simulations/OUTPUT_optim/",sims[i]))){ # TEMPORARY SOLUTION  
      df = read.delim(paste0("Simulations/OUTPUT_optim/",sims[i]),
                      sep=" ", header=TRUE, stringsAsFactors=FALSE)
      C = df[,1] + df[,2]
      h = df[,2]/C
      C_list[[i]] = C_line(C)
      Cmax[i] = max(C)
      h_list[[i]] = h_line(h)
    } else { 
      C_list[[i]] = C_line(NA)
      h_list[[i]] = h_line(NA)
    }
  }
  
  pdf("Simulations/PDF/simPDF.pdf", width=14, height=8.5)
  op = par(mfrow=c(1,2))
  plot(NA, ylim=c(0,1), xlim=c(0,80), main="Mutation Load",
       xlab="Time (years)", ylab="h")
  for(i in 1:N){
    lines(h_list[[i]]$ts, col=myDarkGrey)
  }
  
  na.Cmax = na.omit(Cmax)
  ymax = sort(na.Cmax)[length(na.Cmax)-1]
  plot(NA, ylim=c(0,ymax), xlim=c(0,80), main="Copy Number", 
       xlab="Time (years)", ylab="C")
  for(i in 1:N){
    lines(C_list[[i]]$ts, col=myDarkGrey)
  }
  par(op)
  dev.off()
  
} else {
  error("Number of simulations required.")
}











