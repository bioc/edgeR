#include "edgeR.h"

/* function to compute the unit deviance and the sum
 *
 * inputs:
 * y        row vector or matrix of counts
 * mu       row vector or matrix of means
 * disp     compressedMatrix of dispersion
 * weights  compressedMatrix of weights
 * design   design matrix
 * do_adust whether to include adjust term
 * 
 * output:
 * output   row vector of adjusted profile likelihood
 * 
 * comment:
 * this is converted from adj_coxreid.cpp and R_compute_apl.cpp written by Aaron
 */

void compute_adj_profile_ll(cmx *y, cmx *mu, cmx *disp, cmx *weights, cmx *design, int do_adjust, double *output) 
{
    const char uplo = 'U';
    const double low_value = 1e-10, log_low_value = log(1e-10);

    int ntag=(y->nrow), nlib=(y->ncol), nvar=(design->ncol);
    double *dm; dm = (design->dmat);

    double *xtwx = R_Calloc(nvar*nvar, double);
    int *pivots  = R_Calloc(nvar, int);
    
    /* We also want to identify the optimal size of the 'work' array 
     * using the ability of the dystrf function to call ILAENV. We then
     * reallocate the work pointer to this value.
     */
    int info=0, lwork=-1;
	double temp_work;
    F77_CALL(dsytrf)(&uplo, &nvar, xtwx, &nvar, pivots, &temp_work, &lwork, &info FCONE);
	if (info) { error("failed to identify optimal size of workspace through ILAENV"); }
    lwork = (int)(temp_work +0.5);
    if (lwork < 1) { lwork = 1; }

    // working space for dsytrf
    double *work = R_Calloc(lwork, double);

    // zwpt: working weights in XtWX, or Wz
    double *zwpt = R_Calloc(nlib, double);

    // row vectors for y mu weights and disp
    double *yptr = R_Calloc(nlib, double);
    double *uptr = R_Calloc(nlib, double);
    double *wptr = R_Calloc(nlib, double);
    double *dptr = R_Calloc(nlib, double);

    for (int tag=0; tag<ntag; ++tag) {
        get_row4(y,mu,disp,weights,tag,yptr,uptr,dptr,wptr);

        output[tag] = 0;
        /* First computing the log-likelihood. */
        for (int lib=0; lib<nlib; ++lib) {
            if (uptr[lib]==0) {
                if (do_adjust) {
                    zwpt[lib] = 0;
                }
                continue; // Mean should only be zero if count is zero, where the log-likelihood would then be 0.
            }

            // Each y is assumed to be the average of 'weights' counts, so we convert
            // from averages to the "original sums" in order to compute NB probabilities.
            double curw = wptr[lib];
            double curu = uptr[lib] * curw;
            double cury = yptr[lib] * curw;
            double curd = dptr[lib] / curw;

            double loglik=0;
            if (curd > 0) {
                // same as loglik <- rowSums(weights*dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
                double r=1/curd;
                double logmur=log(curu+r);
                loglik = cury*log(curu) - cury*logmur + r*log(r) - r*logmur + lgamma(cury+r) - lgamma(cury+1) - lgamma(r);
            } else {
                // same as loglik <- rowSums(weights*dpois(y,lambda=mu,log = TRUE))
                loglik = cury*log(curu) - curu - lgamma(cury+1);
            }
            output[tag] += loglik;

            // Adding the Jacobian, to account for the fact that we actually want the log-likelihood
            // of the _scaled_ NB distribution (after dividing the original sum by the weight).
            // output[tag] += log(curw);

            if (do_adjust) {
                /* Computing 'W', the matrix of negative binomial working weights.
                 * The class computes 'XtWX' and performs an LDL decomposition
                 * to compute the Cox-Reid adjustment factor.
                 */
                zwpt[lib] = curu /(1 + curd * curu);
            }
        }

        if (do_adjust) {
            double adj=0;
            if (nvar==1) {
                for(int lib=0;lib<nlib;++lib) { adj += zwpt[lib]; }
                adj=log(fabs(adj))/2;
            } 
            else {
                compute_xtwx(nlib, nvar, dm, zwpt, xtwx);

                /* DSYTRF computes the factorization of a real symmetric matrix A using
                 * the Bunch-Kaufman diagonal pivoting method.  The form of the factorization is
                 * A = U*D*U**T  or  A = L*D*L**T
                 * where U (or L) is a product of permutation and unit upper (lower)
                 * triangular matrices, and D is symmetric and block diagonal with
                 * 1-by-1 and 2-by-2 diagonal blocks.
                 * 
                 * https://netlib.org/lapack/explore-3.2-html/dsytrf.f.html
                 * 
                 */

                F77_CALL(dsytrf)(&uplo, &nvar, xtwx, &nvar, pivots, work, &lwork, &info FCONE);
                
                if (info<0) { error("LDL factorization failed for XtWX."); }
                    
                // the sum of half log diagonal entries
                for (int i=0; i<nvar; ++i) {
                    double cur_val = xtwx[i*nvar + i];
                    adj = (cur_val < low_value)? adj+log_low_value : adj+log(cur_val)*0.5;
                }
            }
            output[tag] -= adj;
        }
    }

    R_Free(xtwx);
    R_Free(pivots);
    R_Free(zwpt);
    R_Free(work);

    R_Free(yptr);
    R_Free(uptr);
    R_Free(dptr);
    R_Free(wptr);

    return;
}


/* EXPLANATION:
   
   XtWX represents the expected Fisher information. The overall strategy is to compute the 
   log-determinant of this matrix, to compute the adjustment factor for the likelihood (in 
   order to account for uncertainty in the nuisance parameters i.e. the fitted values).

   We want to apply the Cholesky decomposition to the XtWX matrix. However, to be safe,
   we call the routine to do a symmetric indefinite factorisation i.e. A = LDLt. This 
   guarantees factorization for singular matrices when the actual Cholesky decomposition 
   would fail because it would start whining about non-positive eigenvectors. 
  
   We then try to compute the determinant of XtWX. Here we use two facts:

   - For triangular matrices, the determinant is the product of the diagonals.
   - det(LDL*)=det(L)*det(D)*det(L*)
   - All diagonal elements of 'L' are unity.

   Thus, all we need to do is to we sum over all log'd diagonal elements in 'D', which - 
   happily enough - are stored as the diagonal elements of 'xtwx'. (And then 
   divide by two, because that's just how the Cox-Reid adjustment works.)

   'info > 0' indicates that one of the diagonals is zero. We handle this by replacing the
   it with an appropriately small non-zero value, if the diagonal element is zero or NA. This
   is valid because the zero elements correspond to all-zero columns in "WX", which in turn
   only arise when there are fitted values of zero, which will be constant at all dispersions. 
   Thus, any replacement value will eventually cancel out during interpolation to obtain the CRAPLE.

   Note that the LAPACK routine will also do some pivoting, essentially solving PAP* = LDL* for 
   some permutation matrix P. This shouldn't affect anything; the determinant of the permutation 
   is either 1 or -1, but this cancels out, so det(A) = det(PAP*).

   Further note that the routine can theoretically give block diagonals, but this should 
   not occur for positive (semi)definite matrices, which is what XtWX should always be.
*/

// Computes upper-triangular matrix.
void compute_xtwx (int nlib, int nvar, double* X, double* W, double* out) 
{
    const double* xptr1=X;
    for (int coef1=0; coef1<nvar; ++coef1, xptr1+=nlib) {
        const double* xptr2=X;
        for (int coef2=0; coef2<=coef1; ++coef2, xptr2+=nlib) {
            out[coef2]=0;
            for (int lib=0; lib<nlib; ++lib) {
                out[coef2] += xptr1[lib]*xptr2[lib]*W[lib];
            }
        }
        out += nvar;
    }
    return;
}

/* comment for compute_xtwx
 * it is a standard way to loop matrix using pointers
 * compared with index method,
 * it saves time for multiplication
 * only addition is required for pointers
 */