#include <Rcpp.h>
using namespace Rcpp;


//' Function to calculate the logistic of a logistic coalescent model
//' @param parr A parameter vector for population size, growth rate and logistic shape parameter.
//' @param t A sampling and coalescent time vector.
//' @param A A lineages vector
//' @author Fei Xiang (\email{xf3087@@gmail.com})
//' @return loglikelihood of a logistic coalescent model.
//' @useDynLib genieR
//' @importFrom Rcpp sourceCpp
//' @examples
//' library(ape)
//' t1=rcoal(20)
//' x=att(t1)
//' negllc_log(c(1,1,10), x$t, x$A)
//' library(minqa)
//' bobyqa(c(1,2,10),negllc_log,lower=0,upper=Inf,t=x$t,A=x$A)
//' Geniefit(t1,Model="log",start=c(10,20,20),upper=Inf,lower=0)
//' library(dfoptim)
//' nmkb(c(1,2,10),negllc_log,lower=0,upper=Inf,t=x$t,A=x$A)
//' @export
// [[Rcpp::export]]
double negllc_log(NumericVector parr, NumericVector t, NumericVector A){

  // Initialization
  int nint = t.size()-1;
  double ll = 0.0;
  double N0=parr[0];
  double c=parr[2];
  double r=parr[1];
  double dt, intval;
  int a, ch;
  // Main loop
  for(int i=0;i<nint;i++){
    dt=t[i+1]-t[i]+c/r*(exp(r*t[i+1])-exp(r*t[i]));
    a=A[i];
    ch=a*(a-1)/2;
    intval=dt/N0/(1+c);
    if(A[i+1]==(A[i]-1)){
      ll = ll + log(ch)-log(N0)-log(1+c)+log(1+c*exp(r*t[i+1]))-ch*intval;
    }else{
      ll = ll -ch*intval;
    }
  }
  return(-ll);
}
