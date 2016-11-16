#include <Rcpp.h>
using namespace Rcpp;



//' Function to calculate the loglikelihood of a exponential coalescent model
//' @param parr A parameter vector for population size and growth rate.
//' @param t A sampling and coalescent time vector.
//' @param A A lineages vector
//' @author Fei Xiang (\email{xf3087@@gmail.com})
//' @return loglikelihood of a exponential coalescent model.
//' @useDynLib genieR
//' @importFrom Rcpp sourceCpp
//' @examples
//' library(ape)
//' t1=rcoal(20)
//' x=att(t1)
//' negllc_expo(c(1,1), x$t, x$A)
//' library(minqa)
//' bobyqa(c(1,2),negllc_expo,lower=0,upper=Inf,t=x$t,A=x$A)
//' Geniefit(t1,Model="expo",start=c(10,2),upper=Inf,lower=0)
//' @export
// [[Rcpp::export]]
double negllc_expo(NumericVector parr, NumericVector t, NumericVector A){

  // Initialization
  int nint = t.size()-1;
  double ll = 0.0;
  double N0=parr[0];
  double r=parr[1];
  double dt, intval;
  int a, ch;
  // Main loop
  for(int i=0;i<nint;i++){
    dt=exp( r*t[i+1] )-exp( r*t[i] );
    a=A[i];
    ch=a*(a-1)/2;
    intval=dt/N0/r;
    if(A[i+1]==(A[i]-1)){
      ll = ll + log(ch)-log(N0)+r*t[i+1]-ch*intval;
    }else{
      ll = ll -ch*intval;
    }
  }
  return(-ll);
}



