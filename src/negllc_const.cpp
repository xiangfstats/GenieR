#include <Rcpp.h>
using namespace Rcpp;


//' Function to calculate the loglikelihood of a constant coalescent model
//' @param N0 A population size.
//' @param t A sampling and coalescent time vector.
//' @param A A lineages vector
//' @author Simon Frost (\email{sdwfrost@@gmail.com})
//' @return Sorted sampling times, coalescent times and sampling lineages.
//' @useDynLib genieR
//' @importFrom Rcpp sourceCpp
//' @examples
//' library(ape)
//' t1=rcoal(20)
//' x=att(t1)
//' negllc_const(1, x$t, x$A)
//' @export
// [[Rcpp::export]]
double negllc_const(double N0, NumericVector t, NumericVector A){

  // Initialization
  int nint = t.size()-1;
  double ll = 0.0;
  double dt, intval;
  int a, ch;
  // Main loop
  for(int i=0;i<nint;i++){
    dt=t[i+1]-t[i];
    a=A[i];
    ch=a*(a-1)/2;
    intval=dt/N0;
    if(A[i+1]==(A[i]-1)){
      ll = ll + log(ch)-log(N0)-ch*intval;
    }else{
      ll = ll -ch*intval;
    }
  }
  return(-ll);
}









