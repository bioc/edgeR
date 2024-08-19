#ifndef EDGER_H
#define EDGER_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/* Alternatively, do as R does as from version 3.6.2 and pass the character length(s) from C to Fortran. A portable way to do this is
 *
 * // before any R headers, or define in PKG_CPPFLAGS
 * #ifndef  USE_FC_LEN_T
 * # define USE_FC_LEN_T
 * #endif
 * #include <Rconfig.h>
 * #include <R_ext/BLAS.h>
 * #ifndef FCONE
 * # define FCONE
 * #endif
    
   F77_CALL(dgemm)("N", "T", &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z, &nrx FCONE FCONE);

 * (Note there is no comma before or between the FCONE invocations.) 
 *
 * It is strongly recommended that packages which call from C/C++ BLAS/LAPACK routines
 * with character arguments adopt this approach: packages not using will fail to install as from R 4.3.0.
 * 
 * https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings
 */

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif

#include <Rconfig.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>	/* sort.c */

#ifndef FCONE
#define FCONE
#endif

/* Notes:
 * 1. the memory management is better to to choose
 *    R_Calloc and R_Free in a safe R-style way
 *    #define R_Calloc(n,t)  (t *) R_chk_calloc( (size_t) (n), sizeof(t) )
 * 2. only use wrapper function for R-C interface
 *    refer to lmfit.c in stats package  
 * 3. void function is used for the preparation to wrap
 * 4. use to duplicate() to inherit the data structure
 * 5. use getAttrib() with install() to get attributes from R object
 * 6. use Rf_mkNamed() to create a list with names
 *    char *names[]={"a", ..., ""} last should be the empty
 *    or use setAttrib(), refer to arima.c in stats package
 * 
 * References:
 * 1. Writing R extensions Chapter 5 & 6
 *    https://cran.r-project.org/doc/manuals/R-exts.html
 * 2. Advance R Chapter R’s C interface
 *    http://adv-r.had.co.nz/C-interface.html
 * 3. R internals (Rinternal.h for R-C interface)
 *    https://cran.r-project.org/doc/manuals/R-ints.html
 * 4. Learning from other packages: base and stats
 */

/* --------------------------------------------------------- 
 * C struct for matrix and CompressedMatrix 
 * it can handle different input, integer or double
 * 
 * data structure:
 * 
 * dmat: point to real vector for REALSXP object
 * imat: point to integer vector for INTSXP object
 * nrow: number of nows
 * ncol: number of columns
 * isint: indicate whether it is integer
 * type: 4 cases 
 *   0: common matrix
 *   1: repeated by row
 *   2: repeated by column
 *   3: repeated by row and column
 * 
 * help functions:
 * 
 * SEXPtocmx1: convert a common atrix to cmx with type 0
 * SEXPtocmx2: convert a compressMatrix to cmx with isint = 0
 * 
 * get_row: get a row vector from cmx
 * get_row3: get the same row from 3 cmx, a wrapper function
 * get_row4: get the same row from 4 cmx, a wrapper function
 * 
 * max_cmx: get the maximal entry from cmx
 * 
 * check_row_scalar: check whether a row of cmx is equal to ?
 */

typedef struct 
{
    double *dmat;
    int *imat;
    int nrow, ncol;
    int type, isint;
} cmx;

/* ---------------------------------------------------------------
 * object.c
 */

/* convert SEXP to struct cmx */
cmx SEXPtocmx1 (SEXP);
cmx SEXPtocmx2 (SEXP);

/* helper function for cmx object */

/* extract 1, 3, 4 rows from cmx objects */
void get_row (cmx *, int, double *);
void get_row3(cmx *, cmx *, cmx *, int, double *, double *, double *);
void get_row4 (cmx *, cmx *, cmx *, cmx *, int, double *, double *, double *, double *);

/* max element in cmx object */
double max_cmx (cmx *);

/* check whether a row of cmx equal to 0 or 1 */
int check_row_scalar (cmx *, int, double);

/* -------------------------------------------------------------
 * add_prior_count.c
 */

/* computeadjusted offset and prior */
void compute_offsets (cmx *, cmx *, int, int, int, double *, double *);

/* add prior count to the count matrix */
void add_prior_count_vec(cmx *, cmx *, cmx *, double *, double *);
void add_prior_count_mat(cmx *, cmx *, cmx *, double *, double *);

/* ------------------------------------------------------------
 * clowess.c
 * fsquare fcube and clowess function, copy from stat package   
 */
double fsquare(double x);
double fcube(double x);
void clowess(double *, double *, int , double, int, double, double *, double *, double *);

/* ------------------------------------------------------------
 * compute_apl.c 
 */

/* compute XtWX matrix, only upper triangle */
void compute_xtwx (int, int, double*, double*, double*); 

/* compute adjusted profile likelihood */
void compute_adj_profile_ll(cmx *, cmx *, cmx *, cmx *, cmx *, int, double *);

/* -------------------------------------------------------------
 * calculate_cpm.c
 */

/* compute cpm or logCPM */
void calc_cpm_log(cmx *, cmx *, cmx *, double *);
void calc_cpm_raw(cmx *, cmx *, double *);

/* compute average logCPM */
void average_log_cpm(cmx *, cmx *, cmx *, cmx *, cmx *, int, double, double *); 

/* ------------------------------------------------------------
 * compute_nbdev.c
 */

/* compute unit deviance for poisson or negative binomial distribution */
double compute_unit_nb_deviance (double, double, double);

/* compute negative binomial deviance for matrix, sum by row or not */
void compute_nbdev_sum(cmx *, cmx *, cmx *, cmx *, double *);
void compute_nbdev_unit(cmx *, cmx *, cmx *, double *);

/* -------------------------------------------------------------
 * exact_test_by_dev.c
 * perform exact test by deviance 
 */
void exact_test_by_dev(int *, int *, int, int, int, double *, double *); 

/* -------------------------------------------------------------
 * fmm_spline.c
 * fmm_spline function, copy from stat package
 */
void fmm_spline(int, double *, double *, double *, double *, double *);

/* -------------------------------------------------------------
 * glm.c
 */

/* one group fitting for row vector */
void glm_one_group_vec(int, double *,double *,double *, double *, int, double, double, double *,int *); 

/* compute the the fitted value from estimated coef for levenberg method */
void fit_leven_autofill(int, double *, int, double *, double *, double *); 

/* levenberg fitting for row vector */
void fit_leven_vec(int, double *, double *, double *, double *, int, double *, int, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *);

/* one group fitting for matrix input */
void fit_one_group_mat (cmx *, cmx *, cmx *, cmx *, int, double, double *, double *, int *); 

/* compute fitted value for one group fitting */
void get_one_way_fit(cmx *, cmx *, int *, double *);

/* initialize levenberg fitting, compute starting coefficients */
void get_leven_start (cmx *, cmx *, cmx *, cmx *, cmx *, int, double *);

/* levenberg fitting for matrix input */
void fit_leven(cmx *, cmx *, cmx *, cmx *, cmx *, cmx *, double, int, double *, double *, double *, int *, int *);

/* check poisson bound */
void check_poi_bound(cmx *, cmx *, cmx *, int *);

/* -------------------------------------------------------------
 * good_turing.c
 * simple good turing function
 */
void good_turing (int *, int *, int, double, double *, double *); 

/* -------------------------------------------------------------
 * interpolator.c
 * maximal interpolant function
 */
void max_interpolant(double *, cmx *, double *); 

/* -------------------------------------------------------------
 * loess_by_col.c
 * fit loess curve by columns 
 */
void loess_by_column(double *, cmx *, int, double *, double *);

/* -------------------------------------------------------------
 * ql_glm.c
 */

/* qr decomposition for hat values */
void qr_hat (double*, int, int, double*);

/* wrapper function for clowess */
void clowess2 (double *, double *, int, int, double, double *);

/* compute prior quasi-dispersion*/
double compute_prior (double *, double *, double *, int);

/* update prior quasi-dispersion */
double update_prior(cmx *, cmx *, cmx *, cmx *, cmx *, double *);

/* compute adjusted deviance df quasi-dispersion */
void compute_adjust_vec (cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *);

/* compute unit deviance df leverage matrix */
void compute_adjust_mat (cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *);

/* -------------------------------------------------------------
 * ql_weights.c
 * compute weight function
 */
void compute_weight(double, double, double, double*);

#endif
