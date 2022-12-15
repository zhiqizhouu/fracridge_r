#include <Rcpp.h>
using namespace Rcpp;

//' A fast implementation of finding the mean of each vector
//' @param x A size n vector
//' @return mean of vector x
double meanVec(NumericVector x) {
  int num = x.size();
  double total = 0;
  for(int i = 0; i < num; ++i) {
      total = total + x[i];
    }
  return total/num;
}

//' A fast implementation of finding the sample standard deviation of each vector
//' @param x A size n vector
//' @return sample standard deviation of vector x
double stdVec(NumericVector x) {
  int num = x.size();
  double mean_vec = meanVec(x);
  double total = 0;
  for(int i = 0; i < num; ++i) {
    total = total + (x[i] - mean_vec) * (x[i] - mean_vec);
  }
  return total/(num - 1);
}

//' A fast implementation of scaling the matrix using meanVec and stdVec defined in rcpp
//' @param x A size n*p matrix
//' @return scaled matrix x
//' @export
// [[Rcpp::export]]
NumericMatrix scaleMat(NumericMatrix x) {
  int n_col = x.ncol();
  int n_row = x.nrow();
  NumericMatrix out_mat(n_row, n_col);
  for(int i = 0; i < n_col; ++i) {
    double mean_vec = meanVec(x(_ , i));
    double std_vec = stdVec(x(_, i));
    out_mat(_,i) = (x(_ , i) - mean_vec)/(std_vec);
  }
  return out_mat;
}


