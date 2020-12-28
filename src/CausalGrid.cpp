#include <Rcpp.h>
using namespace Rcpp;

//' Determine if vector is constant (or zero length)
//'
//' @param var A numeric vector
//' @keywords internal
// [[Rcpp::export]]
bool const_vect(NumericVector var){
  for (int i = 0, size = var.size(); i < size; ++i) {
    if (var[i] - var[0] > 0 || var[0] - var[i] > 0)
      return false;
  }

  return true;
}
