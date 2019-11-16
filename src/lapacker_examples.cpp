#include <Rcpp.h>
#include "lapacke.h"

//'@title Example: using lapacke_dgels
//'@description Solve linear systems using lapacke_dgels
//'@examples
//'## Expected output:
//'example_lapacke_dgels()
//'@export
// [[Rcpp::export]]
void example_lapacke_dgels() {
  double a[5][3] = {{1, 1, 1}, {2, 3, 4}, {3, 5, 2}, {4, 2, 5}, {5, 4, 3}};
  double b[5][2] = {{-10, -3}, {12, 14}, {14, 12}, {16, 16}, {18, 16}};
  lapack_int info, m, n, lda, ldb, nrhs;
  int i, j;

  m = 5;
  n = 3;
  nrhs = 2;
  lda = 3;
  ldb = 2;

  info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, *a, lda, *b, ldb);

  for (i = 0; i < n; i++) {
    Rcpp::Rcout << "[ ";
    for (j = 0; j < nrhs; j++) {
      Rcpp::Rcout << b[i][j] << " ";
    }
    Rcpp::Rcout << " ]" << std::endl;
  }
}
