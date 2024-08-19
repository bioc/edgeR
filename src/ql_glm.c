#include "edgeR.h"

/* compute the adjusted deviance and degree of freedom for QL method */

/*
  inputs:
  y       raw count matrix
  mu      fitted value matrix
  design  design matrix

  disp    dispersion compressMatrix
  weights weight compressMatrix

  prior   average quasi-dispersion

  vector outputs:
  df  adjusted degree of freedom 
  dev adjusted deviance or total.dev
  s2  quasi dispersion estimate by adjusted dev and df

  matrix outputs:
  lvmat hatvalue matrix
  dvmat unit deviance matrix
  dfmat individual df matrix
*/

const double thresholdzero=1e-4;

/* return s2 adjusted deviance and df vectors */
void compute_adjust_vec (cmx *y, cmx *mu, cmx *design, cmx *disp, double prior, cmx *weights, double *df, double *dev, double *s2)
{
  /* number of tags samples, groups */
  int ntag=(y->nrow), nlib=(y->ncol), ncoef=(design->ncol); 

  /* row vectors for y mu disp weights */
  double *yptr = R_Calloc(nlib,double);
  double *uptr = R_Calloc(nlib,double);
  double *wptr = R_Calloc(nlib,double);
  double *dptr = R_Calloc(nlib,double);
  
  // sqrt(W)X, weighted design matrix
  double *xdpt = R_Calloc(nlib*ncoef,double);
  // hatvalues
  double *hptr = R_Calloc(nlib,double);
  // working weights
  double *zwpt = R_Calloc(nlib,double);
  
  /* store the weights of adjusted deviance, complementary hat values */
  double wpt[2], hdp, udp;    

  /* repeat by tags */                              
  for(int tag=0; tag<ntag; ++tag){
    get_row4(y,mu,disp,weights,tag,yptr,uptr,dptr,wptr);
    /*
      The preparation of working weight matrix is separated in two steps:
      I.  compute the working weight zwpt for each sample
      II. compute the sqrt(W)X by multiply the weight zwpt looped by columns (trick [i%nlib])
    */
    for(int lib=0; lib<nlib; ++lib) zwpt[lib]=sqrt(uptr[lib]/(1+(uptr[lib]*dptr[lib]*wptr[lib]/prior))), hptr[lib]=0;
    for(int i=0; i<nlib*ncoef; ++i) xdpt[i]=(design->dmat)[i] * zwpt[i % nlib];
    
    /* compute hat values */ 
    qr_hat(xdpt, nlib, ncoef, hptr);
    
    /* compute s2 deviance and df */ 
    dev[tag]=0, df[tag]=0;
    for(int lib=0; lib<nlib; ++lib){
      compute_weight(uptr[lib], dptr[lib], prior/wptr[lib], wpt);
      udp = compute_unit_nb_deviance(yptr[lib], uptr[lib], dptr[lib]*wptr[lib]/prior);
      hdp = 1.0 - hptr[lib];
      //correction for extremely small leverage
      if(hdp < thresholdzero) udp=0.0, hdp=0.0;
      dev[tag] += (udp*wpt[0])*wptr[lib], df[tag] += hdp*wpt[1];
    }
    /* correct s2 for extremely small df */
    s2[tag] = (df[tag] < thresholdzero) ? 0 : dev[tag]/df[tag];
  }

  R_Free(yptr);
  R_Free(uptr);
  R_Free(wptr);
  R_Free(dptr);
  R_Free(xdpt);
  R_Free(hptr);
  R_Free(zwpt);

  return;
}

/* return the unit matrices: unit deviance, unit df, and leverages */
void compute_adjust_mat (cmx *y, cmx *mu, cmx *design, cmx *disp, double prior, cmx *weights, double *dfmat, double *dvmat, double *lvmat)
{
  int ntag=(y->nrow), nlib=(y->ncol), ncoef=(design->ncol);
  double *yptr = R_Calloc(nlib,double);
  double *uptr = R_Calloc(nlib,double);
  double *wptr = R_Calloc(nlib,double);
  double *dptr = R_Calloc(nlib,double);
  double *xdpt = R_Calloc(nlib*ncoef,double);
  double *hptr = R_Calloc(nlib,double);
  double *zwpt = R_Calloc(nlib,double);
  double wpt[2], hdp, udp;
  double *fptr, *vptr, *lptr;
  for(int tag=0; tag<ntag; ++tag){
    get_row4(y,mu,disp,weights,tag,yptr,uptr,dptr,wptr);
    for(int lib=0; lib<nlib; ++lib) zwpt[lib]=sqrt(uptr[lib]/(1+(uptr[lib]*dptr[lib]*wptr[lib]/prior))), hptr[lib]=0; 
    for(int i=0; i<nlib*ncoef; ++i) xdpt[i]=(design->dmat)[i] * zwpt[i % nlib];
    qr_hat(xdpt, nlib, ncoef, hptr);
    
    /*
    for(int lib=0; lib<nlib; ++lib){
      compute_weight(uptr[lib], dptr[lib], prior/wptr[lib], wpt);
      udp = compute_unit_nb_deviance(yptr[lib], uptr[lib], dptr[lib]*wptr[lib]/prior);
      hdp = 1.0 - hatvalues[lib];
      if(hdp < thresholdzero) udp=0.0, hdp=0.0;
      R_xlen_t ii = (R_xlen_t)(ntag) * lib + tag; 
      dvmat[ii] = udp*wpt[0];     // dev matrix
      dfmat[ii] = hdp*wpt[1];     // df  matrix
      lvmat[ii] = hatvalues[lib]; // hat matrix
    } 
    */
    
    fptr=dfmat+tag;
    vptr=dvmat+tag;
    lptr=lvmat+tag;
    for(int lib=0; lib<nlib; ++lib, fptr+=ntag, vptr+=ntag, lptr+=ntag){
      compute_weight(uptr[lib], dptr[lib], prior/wptr[lib], wpt);
      udp = compute_unit_nb_deviance(yptr[lib], uptr[lib], dptr[lib]*wptr[lib]/prior);
      hdp = 1.0 - hptr[lib];
      if(hdp < thresholdzero) udp=0.0, hdp=0.0;
      (*vptr) = udp*wpt[0];     // dev matrix
      (*fptr) = hdp*wpt[1];     // df  matrix
      (*lptr) = hptr[lib];      // hat matrix
    }
  }

  R_Free(yptr);
  R_Free(uptr);
  R_Free(wptr);
  R_Free(dptr);
  R_Free(xdpt);
  R_Free(hptr);
  R_Free(zwpt);

  return;
}

/* compute the average quasi dispersion by iteration using two updates */
double update_prior(cmx *y, cmx *mu, cmx *design, cmx *disp, cmx *weights, double *avg)
{
  int ntag=(y->nrow);
  double *df  = R_Calloc(ntag, double);
  double *dev = R_Calloc(ntag, double);
  double *s2  = R_Calloc(ntag, double);

  // initialize prior = 1 
  double prior=1.0;
  
  // first update
  compute_adjust_vec(y,mu,design,disp,prior,weights,df,dev,s2);
  prior = compute_prior(avg,s2,df, ntag);
  // second update
  compute_adjust_vec(y,mu,design,disp,prior,weights,df,dev,s2);
  prior = compute_prior(avg,s2,df, ntag);

  R_Free(df);
  R_Free(dev);
  R_Free(s2);

  return prior;
}


/* The following code is used to compute average quasi-dispersion in edgeR
 * a c-version for .computePrior()
 * Created by Lizhong Chen
 * Last revised on 3 May 2024 
 */

/*
 * a wrapper for the input data, preparing for clowess
 * sort x and y using rsort_with_index
 * calculate delta = diff(range(x))
 */

void clowess2 (double *x, double *y, int n, int iter, double f, double *ans)
{
	int *ind = R_Calloc(n, int);
	for(int i=0;i<n;++i) ind[i]=i;
	rsort_with_index(x,ind,n);

	double delta=0.01*(x[n-1]-x[0]);
  double *yy  = R_Calloc(n, double);
  double *rw  = R_Calloc(n, double);
  double *res = R_Calloc(n, double);
	for(int i=0;i<n;++i) yy[i]=y[ind[i]];

	clowess(x, yy, n, f, iter, delta, ans, rw, res);

	R_Free(ind);
	R_Free(yy);
	R_Free(rw);
	R_Free(res);
	return;

}

// compute average quasi-dispersion using lowess trend
double compute_prior (double *ag, double *s2, double *df, int ntag)
{
	// default settings: f=0.5, iter=3L, t=1e-8
	double t=1e-8, f=0.5, out;
	int k=0, iter=3;

  double *xx  = R_Calloc(ntag, double);
  double *yy  = R_Calloc(ntag, double);

	// remove obs with small df and fit on x^1/4 level
	for(int i=0;i<ntag;++i){
		if(df[i]>t) xx[k]=ag[i], yy[k]=sqrt(sqrt(s2[i])), ++k;
	}
  double *ans = R_Calloc(k, double);

	// fit lowess trend
	clowess2(xx,yy,k,iter,f,ans);

	// compute 90% quantile of the trend, converting from quantile.R with type = 7
	double m = (k-1) * 0.9;
	int lo   = (int)(m);
	rPsort(ans,k,lo);
	rPsort(ans,k,lo+1);
	double h = m-lo;
	double p = (1-h)*ans[lo] + h*ans[lo+1];

	// low bound is 1 
	if(p < 1.0) p=1.0;

	out=p*p*p*p;

	R_Free(xx);
	R_Free(yy);
	R_Free(ans);

	return out;
}

/* hat values by QR decomposition
 * this function is extracted from hat() in stats package
 * hat <- function(x, intercept = TRUE)
 * {
 *   if(is.qr(x)) n <- nrow(x$qr)
 *   else {
 *	    if(intercept) x <- cbind(1, x)
 *	    n <- nrow(x)
 *	    x <- qr(x)
 *   }
 *   rowSums(qr.qy(x, diag(1, nrow = n, ncol = x$rank))^2)
 * }
 * 
 * it contains two parts:
 * 1. QR decomposition of x
 * 2. QR_econ for computation of hat values
 */


/*
  inputs:
  x input matrix with n > p
  n number of rows of x
  p number of columns of x
  
  output:
  hat hat values for the matrix x
*/

void qr_hat (double* x, int n, int p, double* hat)
{
  /* k: the rank of matrix x */
  int k; 

  /* jpvt: pivot vector */
  int *jpvt = R_Calloc(p,int);
  for(int i=0; i<p; ++i) jpvt[i]=i+1; 
  
  /* preparation for QR decomposition */
  double tol=1e-7;
  double *qraux = R_Calloc(p,double);
  double *work  = R_Calloc(2*p,double);

  /* call dqrdc2: https://svn.r-project.org/R/branches/R-4-4-branch/src/appl/dqrdc2.f */
  F77_CALL(dqrdc2)(x, &n, &n, &p, &tol, &k, qraux, jpvt, work); 

  /* preparation for QR_econ */

  /*
  double *y = R_Calloc((R_xlen_t)(n)*k,double); 
  for(R_xlen_t i=0; i < (R_xlen_t)(n)*k;++i) y[i]=0;         // initialization
  for(int i=0; i<k; ++i){
    R_xlen_t ii = (R_xlen_t)(n)*i+i;
    y[ii]=1;                                                 // diagonalization
  }
  */

  int nk=n*k;
  double *y = R_Calloc(nk,double); 
  for(int i=0; i<nk; ++i) y[i]=0;                            // initialization

  double *yptr; yptr=y;
  for(int i=0; i<k; ++i, yptr+=n){
    yptr[i]=1;                                               // diagonalization
  }

  /* call dqrqy: https://svn.r-project.org/R/branches/R-4-4-branch/src/appl/dqrutl.f */
  F77_CALL(dqrqy)(x, &n, &k, qraux, y, &k, y); 

  /* compute hat values */

  /*
  for(int i=0;i<n;++i){
    for(int j=0;j<k;++j){
      R_xlen_t ii = (R_xlen_t)(n)*j + i;
      hat[i] += y[ii]*y[ii];
    }
  } 
  */

  for(int i=0;i<n;++i){
    yptr=y+i;
    for(int j=0;j<k;++j,yptr+=n){
      hat[i] += fsquare(*yptr);
    }
  } 

  R_Free(jpvt);
  R_Free(qraux);
  R_Free(work);
  R_Free(y);

  return;
}
