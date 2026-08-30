#include "edgeR.h"

/* Single-call predictive fold-change fit: augment the counts and offsets with
 * library-size-scaled prior counts (add_prior_count_mat), then fit the genewise
 * NB GLM on the augmented data (fit_glm_mat), returning only the coefficients.
 * This folds the former R predFC() path (.cxx_add_prior_count then .cxx_fit_glm,
 * run single-threaded) used by glmFit.default into one nthreads-aware C call.
 * predFC divides the coefficients by log(2) and glmFit.default multiplies them
 * back, so the natural-log coefficients from fit_glm_mat are returned as-is.
 *
 * inputs:
 *   y        count matrix, ntag x nlib (dense cmx)
 *   offsets  log offsets (log library sizes), compressed cmx
 *   priors   prior counts, compressed cmx
 *   disp     dispersions, compressed cmx
 *   weights  observation weights, compressed cmx
 *   design   design matrix, nlib x ncoef (dense cmx, coerced double)
 *   maxit / tol   IWLS convergence settings forwarded to the fitter
 *   nthreads thread count forwarded to both workers
 *
 * output (caller-allocated):
 *   coef     coefficients, ntag x ncoef (natural log scale)
 *
 * Memory: yy/adjoff and the throwaway mu/dev/iter/failed are R_Calloc, freed at
 * the single exit.  Caveat: fit_glm_mat can error() (longjmp) and skip the frees,
 * a bounded one-time leak on a fatal path, matching the other combined workers.
 */
void pred_fc_mat (cmx *y, cmx *offsets, cmx *priors, cmx *disp, cmx *weights, cmx *design, int maxit, double tol, double *coef, int nthreads)
{
    int ntag = (y->nrow);
    int nlib = (y->ncol);

    /* scratch (all R_Calloc, freed at the single exit) */
    double *yy = R_Calloc((size_t) ntag * nlib, double);
    double *adjoff = R_Calloc((size_t) ntag * nlib, double);
    double *mu = R_Calloc((size_t) ntag * nlib, double);
    double *dev = R_Calloc(ntag, double);
    int *iter = R_Calloc(ntag, int);
    int *failed = R_Calloc(ntag, int);

    /* augment counts and offsets with library-size-scaled prior counts */
    add_prior_count_mat(y, offsets, priors, yy, adjoff, nthreads);

    /* dense cmx views over the augmented counts and adjusted offsets */
    cmx yymx = make_cmx(yy, NULL, ntag, nlib, 0, 0);
    cmx adjoffmx = make_cmx(adjoff, NULL, ntag, nlib, 0, 0);

    /* fit the genewise NB GLMs on the augmented data (oneway-vs-Levenberg decided
       in C); only coef is kept, mu / dev / iter / failed / method are scratch */
    int method = 0;
    fit_glm_mat(&yymx, &adjoffmx, disp, weights, design, maxit, tol, NULL, coef, mu, dev, iter, failed, &method, nthreads);

    R_Free(yy);
    R_Free(adjoff);
    R_Free(mu);
    R_Free(dev);
    R_Free(iter);
    R_Free(failed);

    return;
}
