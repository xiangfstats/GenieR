#include <Rcpp.h>
using namespace Rcpp;

//' Function to calculate the loglikelihood of a expansion coalescent model
//' @param parr A parameter vector for population size, growth rate and population size at when t is infinity.
//' @param t A sampling and coalescent time vector.
//' @param A A lineages vector
//' @author Fei Xiang (\email{xf3087@@gmail.com})
//' @return loglikelihood of a expansion coalescent model.
//' @useDynLib genieR
//' @importFrom Rcpp sourceCpp
//' @examples
//' library(ape)
//' t1=rcoal(20)
//' x=att(t1)
//' negllc_expan(c(1,1,10), x$t, x$A)
//' library(minqa)
//' bobyqa(c(1,2,10),negllc_expan,lower=0,upper=Inf,t=x$t,A=x$A)
//' Geniefit(t1,Model="expan",start=c(10,2,2),upper=Inf,lower=0)
//' @export
// [[Rcpp::export]]
double negllc_expan(NumericVector parr, NumericVector t, NumericVector A){

  // Initialization
  int nint = t.size()-1;
  double ll = 0.0;
  double N0=parr[0];
  double alpha=parr[2];
  double r=parr[1];
  double dt, intval;
  int a, ch;
  // Main loop
  for(int i=0;i<nint;i++){
    dt=log(alpha*exp(r*t[i+1])+1-alpha)-log(alpha*exp(r*t[i])+1-alpha);
    a=A[i];
    ch=a*(a-1)/2;
    intval=dt/N0/r/alpha;
    if(A[i+1]==(A[i]-1)){
      ll = ll + log(ch)-log(N0)-log(alpha+(1-alpha)*exp(-r*t[i+1]))-ch*intval;
    }else{
      ll = ll -ch*intval;
    }
  }
  return(-ll);
}
