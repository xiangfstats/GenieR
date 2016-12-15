#include <Rcpp.h>
using namespace Rcpp;

//' Cumulative sum.
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @export
//' @return bounded numeric vector
// [[Rcpp::export]]
NumericVector cumsum_bounded(NumericVector x, double low, double high) {
  NumericVector res(x.size());
  double acc = 0;
  for (int i=0; i < x.size(); ++i) {
    acc += x[i];
    if (acc < low)  acc = low;
    else if (acc > high)  acc = high;
    res[i] = acc;
  }
  return res;
}
