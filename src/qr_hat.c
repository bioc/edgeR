#include <Rmath.h>
#include "R_ext/Applic.h"

/* hat values by QR decomposition */

/*
  inputs:
  x input matrix with n > p
  n number of rows of x
  p number of columns of x
  
  output:
  hat hat values for the matrix x
*/

void QR_hat (double* x, int n, int p, double* hat)
{
  int k, jpvt[p];  // k will be the rank of matrix x
  for(int i=0; i<p; ++i) jpvt[i]=i+1; // pivot vector
  double tol=1e-7, qraux[p], work[2*p];
  F77_CALL(dqrdc2)(x, &n, &n, &p, &tol, &k, qraux, jpvt, work); // QR decomposition
  double y[n*k]; for(int i=0; i<(n*k);++i) y[i]=0; // initialization
  for(int i=0; i<k; ++i) y[i*n+i]=1; // diagonalization
  F77_CALL(dqrqy)(x, &n, &k, qraux, y, &k, y); // QR_econ
  for(int i=0;i<n;++i){
    for(int j=0;j<k;++j)
      *(hat+i) += y[i+j*n]*y[i+j*n]; // compute hat values
  } 
  return;
}
