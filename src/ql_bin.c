#include "edgeR.h"

/* compute the adjusted deviance and degree of freedom for bionomial distribution */

/*
  inputs:
  y       raw proportion matrix
  mu      fitted value matrix
  cover   coverage matrix
  design  design matrix
  weights weight matrix

  vector outputs:
  df  adjusted degree of freedom 
  dev adjusted deviance or total.dev
  s2  quasi dispersion estimate by adjusted dev and df

  matrix outputs:
  lvmat hatvalue matrix
  dvmat unit deviance matrix
  dfmat individual df matrix
*/

/* The QL adjustment scratch (ql_ws) and the shared vec/mat drivers live in ql_glm.c;
 * the binomial entry points below are thin wrappers over ql_adjust_vec / ql_adjust_mat
 * with fam=QL_BIN (aux holds the coverage matrix; the prior argument is unused). */

/* Compute the QL binomial adjusted total deviance, adjusted residual df, and
 * quasi-dispersion for every gene, parallelised over tags. For each library it
 * forms the working weight sqrt(cover*mu*(1-mu)), builds sqrt(W)X, and takes the
 * hat values from qr_hat, then accumulates the leverage-adjusted deviance and df.
 *
 * inputs:
 *   y        raw proportion matrix (ntag x nlib)
 *   mu       fitted value matrix
 *   design   design matrix (nlib x ncoef)
 *   cover    coverage matrix
 *   weights  weight matrix
 *   nthreads requested OpenMP thread count (clamped)
 *
 * outputs (caller-allocated, length ntag):
 *   df   adjusted degree of freedom per gene
 *   dev  adjusted deviance (total.dev) per gene
 *   s2   quasi dispersion estimate dev/df per gene
 *
 * No return value; writes only into the caller-owned df, dev, s2 arrays.
 * Thin wrapper over the shared ql_adjust_vec / ql_adjust_mat drivers (ql_glm.c), fam=QL_BIN.
 */
/* return s2 adjusted deviance and df vectors */
void compute_adjust_vec_bin (cmx *y, cmx *mu, cmx *design, cmx *cover, cmx *weights, double *df, double *dev, double *s2, int nthreads)
{
  ql_adjust_vec(QL_BIN, y, mu, design, cover, 1.0, weights, df, dev, s2, nthreads);
}

/* Compute the per-observation QL binomial unit matrices for every gene,
 * parallelised over tags. Mirrors compute_adjust_vec_bin's per-library math but
 * stores the unit deviance, unit df, and hat value for each (gene, library)
 * instead of summing them.
 *
 * inputs:
 *   y        raw proportion matrix (ntag x nlib)
 *   mu       fitted value matrix
 *   design   design matrix (nlib x ncoef)
 *   cover    coverage matrix
 *   weights  weight matrix
 *   nthreads requested OpenMP thread count (clamped)
 *
 * outputs (caller-allocated, each ntag x nlib, one library per column block):
 *   dfmat  individual df matrix
 *   dvmat  unit deviance matrix
 *   lvmat  hatvalue (leverage) matrix
 *
 * No return value; writes only into the caller-owned dfmat, dvmat, lvmat arrays.
 * Thin wrapper over the shared ql_adjust_vec / ql_adjust_mat drivers (ql_glm.c), fam=QL_BIN.
 */
/* return the unit matrices: unit deviance, unit df, and leverages */
void compute_adjust_mat_bin (cmx *y, cmx *mu, cmx *design, cmx *cover, cmx *weights, double *dfmat, double *dvmat, double *lvmat, int nthreads)
{
  ql_adjust_mat(QL_BIN, y, mu, design, cover, 1.0, weights, dfmat, dvmat, lvmat, nthreads);
}

/* Single-call QL binomial fit: fit the genewise binomial GLM, then form the
 * leverage-adjusted deviance / df / quasi-dispersion, in one call. Folds the
 * former two-step R path (.cxx_bin_fit then .cxx_compute_adj_vec_bin) used by
 * binQLFit.default. Coverage / proportion / coverage-scaled weights are derived
 * from the raw counts here (mirroring the R lines coverage=y+z,
 * weights0=weights*coverage, prop=y/pmax(coverage,1), and ftest.c); then
 * bin_fit_mat produces the fit and compute_adjust_vec_bin the QL adjustment, so
 * the outputs are identical to the former two calls.
 *
 * inputs:
 *   ycnt     successes (counts) matrix, ntag x nlib (dense cmx, int or double)
 *   zcnt     failures  (counts2) matrix, ntag x nlib
 *   wts      observation weights, ntag x nlib, or unused when has_w == 0
 *   has_w    1 if wts holds weights, 0 if weights were NULL on the R side
 *   offsets  log offset, compressed cmx (passed straight to bin_fit_mat)
 *   design   design matrix, nlib x ncoef (dense cmx, coerced double)
 *   maxit_ow / tol_ow   oneway convergence settings
 *   maxit_iw / tol_iw   IWLS convergence settings
 *   nthreads thread count forwarded to the fitter and the adjustment
 *
 * outputs (caller-allocated):
 *   coef     coefficients, ntag x ncoef
 *   mu       fitted values, ntag x nlib (also read back by the adjustment)
 *   dev      residual deviance, length ntag
 *   iter     IWLS iteration count per gene (from bin_fit_mat), length ntag
 *   failed   convergence-failure flag per gene (from bin_fit_mat), length ntag
 *   method   0 = oneway shortcut, 1 = general IWLS (from bin_fit_mat), scalar
 *   df_adj   adjusted residual df, length ntag
 *   dev_adj  adjusted (total) deviance, length ntag
 *   s2       quasi-dispersion dev_adj/df_adj, length ntag
 *
 * iter / failed / method are passed straight through from bin_fit_mat so the
 * binQLFit object carries the same fit diagnostics as binFit.
 *
 * Memory: prop/cover/w0 and the per-gene row buffers are R_Calloc, freed at the
 * single exit. Caveat: bin_fit_mat can error() (longjmp) and skip the frees, a
 * bounded one-time leak on a fatal path, matching ftest.c.
 */
void bin_ql_fit_mat (cmx *ycnt, cmx *zcnt, cmx *wts, int has_w, cmx *offsets, cmx *design, int maxit_ow, double tol_ow, int maxit_iw, double tol_iw, double *coef, double *mu, double *dev, int *iter, int *failed, int *method, double *df_adj, double *dev_adj, double *s2, int nthreads)
{
  int ntag = ycnt->nrow;
  int nlib = ycnt->ncol;

  /* scratch (all R_Calloc, freed at the single exit) */
  double *prop = R_Calloc((size_t) ntag * nlib, double);
  double *cover = R_Calloc((size_t) ntag * nlib, double);
  double *w0 = R_Calloc((size_t) ntag * nlib, double);
  double *ybuf = R_Calloc(nlib, double);
  double *zbuf = R_Calloc(nlib, double);
  double *wbuf = R_Calloc(nlib, double);

  /* coverage, coverage-scaled weights and proportion per gene */
  for(int g = 0; g < ntag; ++g)
  {
    get_row(ycnt, g, ybuf);
    get_row(zcnt, g, zbuf);
    if(has_w)
    {
      get_row(wts, g, wbuf);
    }
    /* cover/w0/prop walked by ntag over libraries (no g + l*ntag product) */
    /*
    for(int l = 0; l < nlib; ++l)
    {
      double cov = ybuf[l] + zbuf[l];
      cover[g + l * ntag] = cov;
      w0[g + l * ntag] = has_w ? (wbuf[l] * cov) : cov;
      prop[g + l * ntag] = ybuf[l] / fmax2(cov, 1.0);
    }
    */
    double *cvp = cover + g, *wp = w0 + g, *pp = prop + g;
    for(int l = 0; l < nlib; ++l, cvp += ntag, wp += ntag, pp += ntag)
    {
      double cov = ybuf[l] + zbuf[l];
      *cvp = cov;
      *wp = has_w ? (wbuf[l] * cov) : cov;
      *pp = ybuf[l] / fmax2(cov, 1.0);
    }
  }

  /* dense cmx views over prop, coverage-scaled weights, coverage and the fitted
     values (mu is filled by bin_fit_mat then read by the adjustment) */
  cmx propmx = make_cmx(prop, NULL, ntag, nlib, 0, 0);
  cmx w0mx = make_cmx(w0, NULL, ntag, nlib, 0, 0);
  cmx covermx = make_cmx(cover, NULL, ntag, nlib, 0, 0);
  cmx mumx = make_cmx(mu, NULL, ntag, nlib, 0, 0);

  /* fit the tagwise binomial GLMs (oneway-vs-IWLS decided in C), filling
     coef / mu / dev and the fit diagnostics iter / failed / method */
  bin_fit_mat(&propmx, offsets, &w0mx, design, maxit_ow, tol_ow, maxit_iw, tol_iw, coef, mu, dev, iter, failed, method, nthreads);

  /* leverage-adjusted deviance, residual df and quasi-dispersion */
  compute_adjust_vec_bin(&propmx, &mumx, design, &covermx, &w0mx, df_adj, dev_adj, s2, nthreads);

  R_Free(prop);
  R_Free(cover);
  R_Free(w0);
  R_Free(ybuf);
  R_Free(zbuf);
  R_Free(wbuf);

  return;
}
