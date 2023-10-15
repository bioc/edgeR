/* a wrapper function for compute_adjusted_deviance */

#include <R.h>
#include <Rinternals.h>

void c_compute_adjust_s2 (double *, double *, int *, int *, double *, int *, double *, double *, double *, double *, double *, double *, double *,double *, double *);

SEXP compute_adjust_s2 (SEXP y, SEXP mu, SEXP de, SEXP di, SEXP wo, SEXP we)
{
    SEXP ans;
    SEXP df, dev, s2, lever, unitdev, unitdf;
    int nt, nl, nv;

    /* extract dimensions of count matrix and design matrix */
    ans = getAttrib(y, R_DimSymbol);
    int *dims = INTEGER(ans);
    nl = dims[0], nt = dims[1];
    ans = getAttrib(de, R_DimSymbol);
    int *dim = INTEGER(ans);
    nv = dim[1];

    /* ensure double input for computation */
	PROTECT(y = coerceVector(y, REALSXP));
	PROTECT(mu = coerceVector(mu, REALSXP));
	PROTECT(de = coerceVector(de, REALSXP));
	PROTECT(di = coerceVector(di, REALSXP));
	PROTECT(wo = coerceVector(wo, REALSXP));
	PROTECT(we = coerceVector(we, REALSXP));

    /* prepare output */
    PROTECT(df  = allocVector(REALSXP, nt));
    PROTECT(dev = allocVector(REALSXP, nt));
    PROTECT(s2  = allocVector(REALSXP, nt));
    PROTECT(lever   = allocMatrix(REALSXP, nl, nt));
    PROTECT(unitdev = allocMatrix(REALSXP, nl, nt));
    PROTECT(unitdf  = allocMatrix(REALSXP, nl, nt));

    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, df);
    SET_VECTOR_ELT(ans, 1, dev);
    SET_VECTOR_ELT(ans, 2, s2);
    SET_VECTOR_ELT(ans, 3, lever);
    SET_VECTOR_ELT(ans, 4, unitdev);
    SET_VECTOR_ELT(ans, 5, unitdf);

    /* call compute function */
    c_compute_adjust_s2(REAL(y), REAL(mu), &nt, &nl, REAL(de), &nv, REAL(di), REAL(wo), REAL(we), REAL(df), REAL(dev), REAL(s2), REAL(lever), REAL(unitdev), REAL(unitdf));       

    /* check the number of protect */
    UNPROTECT(13);
    return ans;
}