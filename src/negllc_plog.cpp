#include <Rcpp.h>
using namespace Rcpp;


//' Function to calculate the loglikelihood of a piecewise logistic coalescent model
//' @param parr A parameter vector for population size, exponential growth rate and change time.
//' @param t A sampling and coalescent time vector.
//' @param A A lineages vector
//' @author Fei Xiang (\email{xf3087@@gmail.com})
//' @return loglikelihood of a piecewise logistic coalescent model.
//' @useDynLib genieR
//' @importFrom Rcpp sourceCpp
//' @examples
//' library(ape)
//' t1=rcoal(20)
//' x=att(t1)
//' negllc_plog(c(1,1,.1), x$t, x$A)
//' library(minqa)
//' bobyqa(c(1,2,1),negllc_plog,lower=0,upper=Inf,t=x$t,A=x$A)
//' Geniefit(t1,Model="plog",start=c(10,2,1),upper=Inf,lower=0)
//' library(dfoptim)
//' nmkb(c(1,2,10),negllc_plog,lower=0,upper=Inf,t=x$t,A=x$A)
//' @export
// [[Rcpp::export]]
double negllc_plog(NumericVector parr, NumericVector t, NumericVector A){

  // Initialization
  int nint = t.size()-1;
  double ll = 0.0;
  double N0=parr[0];
  double r=parr[1];
  double x=parr[2];
  double intval, integrand;
  int a, ch;

  // Main loop
  for(int i=0;i<nint;i++){
    a=A[i];
    ch=a*(a-1)/2;
    integrand=int(t[i+1]<x)*N0+int(t[i+1]>x)*exp(-r*(t[i+1]-x))*N0;
    intval=int(t[i+1]<x)*(t[i+1]-t[i])/N0+int(t[i]>x)*(-exp(r*t[i]-r*x)+exp(r*t[i+1]-r*x))/N0/r
      +int(t[i]<x && t[i+1]>x)*((t[i+1]-x)/N0+(-exp(r*t[i]-r*x)+1)/N0/r);



    if(A[i+1]==(A[i]-1)){
      ll = ll + log(ch)-log(integrand)-ch*intval;
    }else{
      ll = ll -ch*intval;
    }
  }
  return(-ll);
}


