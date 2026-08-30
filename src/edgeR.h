#ifndef EDGER_H
#define EDGER_H

/* omp.h before R headers: R remaps match->Rf_match, clashing with
 * clang omp.h declare-variant match() clause (LLVM). */
#ifdef _OPENMP
#include <omp.h>
#endif

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

/* ---------------------------------------------------------------
 * OpenMP support
 *
 * The genewise kernels loop over independent tags and can be run in
 * parallel. The code must still compile and run correctly when the
 * compiler has no OpenMP, so every use of <omp.h> is guarded by
 * _OPENMP and a serial fallback is provided.
 *
 * Per-thread scratch is allocated once before each parallel region
 * (R_Calloc/R_Free are not thread-safe and must never be called from
 * within a parallel region); every kernel below takes an explicit
 * nthreads argument supplied from R.
 */

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
 * 2. Advance R Chapter R's C interface
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
 * nrow: number of (logical) rows -- the subset size when this is a row-subset view
 * ncol: number of columns
 * isint: indicate whether it is integer
 * type: 4 cases
 *   0: common matrix
 *   1: repeated by column
 *   2: repeated by row
 *   3: repeated by row and column
 * pnrow: physical row count = column stride; equals nrow unless a row subset
 * rowsel: NULL, or a length-nrow array of physical row indices (a subset view)
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
    int pnrow;           /* physical row count (column stride); == nrow when not a subset */
    const int *rowsel;   /* NULL, or length-nrow physical row indices for a subset view   */
} cmx;

/* ---------------------------------------------------------------
 * object.c
 */

/* convert SEXP to struct cmx */
cmx SEXPtocmx1 (SEXP);
cmx SEXPtocmx2 (SEXP);

/* build a cmx view over a caller-owned buffer */
cmx make_cmx (double *, int *, int, int, int, int);

/* zero-copy row-subset view of a FULL cmx (src.rowsel must be NULL); rows are
   length-nrows physical row indices (0-based), caller-owned */
cmx subset_cmx (cmx, const int *, int);

/* logical indicator -> row indices, into caller-owned idx (length >= n); returns count */
int which_index (const int *, int, int *);

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

/* computeadjusted offset and prior
 * the trailing two pointers are caller-owned nlib scratch buffers so this
 * helper performs no allocation of its own (safe to call inside a parallel
 * region)
 */
void compute_offsets (cmx *, cmx *, int, int, int, double *, double *, double *, double *);

/* add prior count to the count matrix */
void add_prior_count_vec(cmx *, cmx *, cmx *, double *, double *, int);
void add_prior_count_mat(cmx *, cmx *, cmx *, double *, double *, int);

/* ------------------------------------------------------------
 * utils.c : small general-purpose helpers
 */
int clamp_threads(int);
double fsquare(double);
double fcube(double);
double brent_fmin(double, double, double (*)(double, void *), void *, double);
/* lowess.f (Fortran) : Cleveland LOWESS smoother; replaces the former clowess.c */
void F77_NAME(lowess)(double *, double *, int *, double *, int *, double *, double *, double *, double *);

/* ------------------------------------------------------------
 * compute_apl.c 
 */

/* compute XtWX matrix, only upper triangle */
void compute_xtwx (int, int, double*, double*, double*); 

/* compute adjusted profile likelihood */
void compute_adj_profile_ll(cmx *, cmx *, cmx *, cmx *, cmx *, int, double *, int);

/* -------------------------------------------------------------
 * calculate_cpm.c
 */

/* compute cpm or logCPM */
void calc_cpm_log(cmx *, cmx *, cmx *, double *, int);
void calc_cpm_raw(cmx *, cmx *, double *, int);

/* compute average logCPM */
void average_log_cpm(cmx *, cmx *, cmx *, cmx *, cmx *, int, double, double *, int);

/* ------------------------------------------------------------
 * compute_nbdev.c
 */

/* compute unit deviance for poisson or negative binomial distribution */
double compute_unit_nb_deviance (double, double, double);

/* weighted sum of unit deviances over n libraries: sum_i w[i]*unit_dev(y[i],mu[i],disp[i]) */
double nb_deviance_sum (int, const double *, const double *, const double *, const double *);

/* compute negative binomial deviance for matrix, sum by row or not */
void compute_nbdev_sum(cmx *, cmx *, cmx *, cmx *, double *);
void compute_nbdev_unit(cmx *, cmx *, cmx *, double *);

/* -------------------------------------------------------------
 * exact_test_by_dev.c
 * perform exact test by deviance
 */
void exact_test_by_dev(int *, int *, int, int, int, double *, double *);

/* -------------------------------------------------------------
 * spline.f (Fortran)
 * cubic spline coefficients; replaces the former fmm_spline.c
 */
void F77_NAME(spline)(int *, double *, double *, double *, double *, double *);

/* -------------------------------------------------------------
 * glm.c
 */

/* one group fitting for row vector */
void glm_one_group_vec(int, double *,double *,double *, double *, int, double, double, double *,int *); 

/* levenberg fitting for row vector
 * trailing int* is a LAPACK-failure flag set in place of error() so the
 * caller can raise the error outside any parallel region
 */
void fit_leven_vec(int, double *, double *, double *, double *, int, double *, int, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);

/* one group fitting for matrix input */
void fit_one_group_mat (cmx *, cmx *, cmx *, cmx *, int, double, double *, double *, int *, int);

/* fit one way layout for all groups */
void fit_one_way_mat(cmx *, cmx *, cmx *, cmx *, int *, int, int, double, double *, double *, int *, int);

/* compute fitted value for one group fitting */
void get_one_way_fit(cmx *, cmx *, int *, double *, int);

/* initialize levenberg fitting, compute starting coefficients */
void get_leven_start (cmx *, cmx *, cmx *, cmx *, cmx *, int, double *, int);

/* levenberg fitting for matrix input */
void fit_leven(cmx *, cmx *, cmx *, cmx *, cmx *, cmx *, double, int, double *, double *, double *, int *, int *, int);

/* unified fitting for negative binomial generalized linear models */
void fit_glm_mat(cmx *, cmx *, cmx *, cmx *, cmx *, int, double, double *, double *, double *, double *, int *, int *, int *, int);

/* single-call predictive fold-change fit: add prior counts then fit, return coef */
void pred_fc_mat(cmx *, cmx *, cmx *, cmx *, cmx *, cmx *, int, double, double *, int);

/* check poisson bound */
void check_poi_bound(cmx *, cmx *, cmx *, int *, int);

/* get groups from design matrix */
int get_groups_from_design(cmx *, int *);

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

/* qr decomposition for hat values
 * jpvt (p), qraux (p), work (2*p) and qy (n*p) are caller-owned scratch so
 * this helper allocates nothing and is safe to call inside a parallel region
 */
void qr_hat (double*, int, int, double*, int*, double*, double*, double*);

/* update prior quasi-dispersion */
double update_prior(cmx *, cmx *, cmx *, cmx *, cmx *, double *, int);

/* shared QL adjustment drivers (defined in ql_glm.c, called from ql_bin.c too):
 * one vec + one mat kernel dispatched by distribution.  aux holds the NB dispersion
 * or the binomial coverage; prior is used only for QL_NB. */
enum ql_family { QL_NB, QL_BIN };
void ql_adjust_vec (int, cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *, int);
void ql_adjust_mat (int, cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *, int);

/* compute adjusted deviance df quasi-dispersion for negative binomial distribution */
void compute_adjust_vec (cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *, int);

/* compute unit deviance df leverage matrix negative binomial distribution*/
void compute_adjust_mat (cmx *, cmx *, cmx *, cmx *, double, cmx *, double *, double *, double *, int);

/* -------------------------------------------------------------
 * ql_bin.c
 */

/* compute adjusted deviance df quasi-dispersion for binomial distribution */
void compute_adjust_vec_bin (cmx *, cmx *, cmx *, cmx *, cmx *, double *, double *, double *, int);

/* compute unit deviance df leverage matrix for binomial distribution */
void compute_adjust_mat_bin (cmx *, cmx *, cmx *, cmx *, cmx *, double *, double *, double *, int);

/* -------------------------------------------------------------
 * sample_weights.c
 */

/* empirical per-sample weights from adjusted unit deviances (NB and binomial) */
void sample_weights_nb (cmx *, cmx *, cmx *, cmx *, double, cmx *, const double *, int, int, double *);
void sample_weights_binom (cmx *, cmx *, cmx *, cmx *, cmx *, int, double *);

/* single-call QL binomial fit: fit + leverage-adjusted deviance/df/dispersion */
void bin_ql_fit_mat (cmx *, cmx *, cmx *, int, cmx *, cmx *, int, double, int, double, double *, double *, double *, int *, int *, int *, double *, double *, double *, int);

/* -------------------------------------------------------------
 * ql_weights.c
 * compute weight function for negative binomial and binomial distributions
 */
void compute_weight_negbin (double, double, double, double*);
void compute_weight_binomial (double , int , double *);

/* -------------------------------------------------------------
 * binomial.c
 */

/* helper function for binomial distributions */
double y_log_y(double , double);

/* this function fits one group design of binomial models for matrix input */
void bin_one_group_mat (cmx *, cmx *, cmx *, int , double , double *, double *, double *, int *, int);

/* this function fits binomial models for matrix input using IWLS */
void bin_iwls_mat (cmx *, cmx *, cmx *, cmx *, int , double , double *, double *, double *, int *, int *, int);

/* unified single-call binomial fitter (oneway + general IWLS), mirrors fit_glm */
void bin_fit_mat (cmx *, cmx *, cmx *, cmx *, int, double, int, double, double *, double *, double *, int *, int *, int *, int);

/* -------------------------------------------------------------
 * ftest.c
 */

/* single-call binomial QL F-test with TREAT/UPSHOT p-values */
void bin_ftest_mat (cmx *, cmx *, cmx *, int, cmx *, int *, int, double *, double *, double *, double *, cmx *, int, double *, double *, int);

/* single-call negative-binomial QL F-test with TREAT/UPSHOT p-values */
void glm_ftest_mat (cmx *, cmx *, cmx *, cmx *, cmx *, int *, int, double *, double *, double *, double *, cmx *, int, double *, double *, int);

/* -------------------------------------------------------------
 * coxreid.c
 */
double coxreid_disp_opt(cmx *, cmx *, cmx *, cmx *, double, double, double, int, int);
double coxreid_disp_top_opt(cmx *, cmx *, cmx *, cmx *, const double *, double, double, double, int, int);

/* -------------------------------------------------------------
 * diffsplice.c
 */
void fit_diff_splice_mat(cmx *, cmx *, cmx *, cmx *, cmx *, int, int, int, int, int *, int *, int, int, double, double *, double *, double *, double *, double *, int);
#endif
