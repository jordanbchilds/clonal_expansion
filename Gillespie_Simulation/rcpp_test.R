install.packages("Rcpp")
library("Rcpp")

cppFunction('int sampleC(NumericVector hnorm, double u) {
  int v = hnorm.size(); 
  for( int i=0; i < v; i++ ){
    if( hnorm[i] > u ){
      return i+1;
    }
  }
}')

t = runif(5)
t = cumsum(t)/sum(t)
t
tt = double(1000)
system.time( for(i in 1:10000) tt[i] = sampleC(hnorm=t, u=runif(1)) )
system.time( for(i in 1:10000) sample(1:5,1,prob=t))

cppFunction('NumericVector cumsumC(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for(int i=0; i<n; i++) {
    float temp = 0;
    for(int j=0; j<i+1; j++){
      temp += x[j];
    }
    out[i] = temp;
  }
  return out;
}')


cppFunction('float sumC(NumericVector x) {
  int n = x.size();
  float sum = 0;
  for(int i = 0; i < n; i++) {
    sum += x[i];
  }
  return sum;
}')

system.time( sumC( tt ) )

