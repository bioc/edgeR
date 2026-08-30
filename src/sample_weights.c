#include "edgeR.h"
#include <math.h>

/* ============================================================================
 * Empirical sample (quality) weights from the adjusted unit deviances.
 *
 * Folds the two per-class pipelines of R/sampleWeights.R into single C calls:
 * each computes the genes-by-samples adjusted unit deviance / unit df matrices
 * (via the existing compute_adjust_mat / compute_adjust_mat_bin kernels) and
 * reduces them to one weight per sample, so the full matrices are never returned
 * to R. The reduction reproduces the former R code:
 *   w <- colSums(unit.deviance * weights / s2) / colSums(unit.df)
 *   w <- exp(mean(log(w)) - log(w))
 * Sums accumulate in plain double; over the (positive) per-gene terms this is
 * accurate to ~1e-12 relative, so results match the former R path to that level
 * (not bit-identical, since R's colSums/mean use long-double accumulation).
 * ========================================================================== */

/* Reduce the per-observation dev/df matrices to one weight per sample.
 * dvmat, dfmat : genes-by-samples, column-major with ntag rows (dvmat[g+j*ntag]).
 * wmx          : the original observation weights (read by row, rowsel-aware).
 * s2           : NULL for the s2=1 (binomial) case; otherwise the prior variance,
 *                length s2len (1 => scalar, else per-gene), floored at 1 (pmax).
 * Writes w[0..nlib-1]. */
static void reduce_sample_weights(const double *dvmat, const double *dfmat, cmx *wmx, const double *s2, int s2len, int ntag, int nlib, double *w)
{
    double *num = R_Calloc(nlib, double);
    double *den = R_Calloc(nlib, double);
    double *wrow = R_Calloc(nlib, double);

    for(int g=0; g<ntag; ++g)
    {
        get_row(wmx, g, wrow);
        double s2g = 1.0;
        if(s2)
        {
            s2g = (s2len==1) ? s2[0] : s2[g];
            if(s2g < 1.0) s2g = 1.0;
        }
        const double *dvptr = dvmat + g;
        const double *dfptr = dfmat + g;
        for(int j=0; j<nlib; ++j, dvptr+=ntag, dfptr+=ntag)
        {
            double contrib = (*dvptr) * wrow[j];
            if(s2)
            {
                contrib /= s2g;
            }
            num[j] += contrib;
            den[j] += (*dfptr);
        }
    }

    /* w = log(colSums(dev*weights/s2) / colSums(df)), then exp(mean(w) - w) */
    double m = 0;
    for(int j=0; j<nlib; ++j)
    {
        w[j] = log(num[j] / den[j]);
        m += w[j];
    }
    m /= nlib;

    for(int j=0; j<nlib; ++j)
    {
        w[j] = exp(m - w[j]);
    }

    R_Free(num);
    R_Free(den);
    R_Free(wrow);
}

/* Negative-binomial sample weights (DGEGLM). Computes the adjusted unit dev/df
 * with compute_adjust_mat, then reduces with the (pmax'd) prior variance s2.
 *   ymx, umx : counts, fitted values
 *   gmx      : design
 *   dmx      : compressed dispersion; prior = average quasi-dispersion
 *   wmx      : compressed observation weights
 *   s2       : prior variance (s2.prior), length s2len; floored at 1 in the reducer
 * Writes w (length nlib). Allocates and frees its own scratch. */
void sample_weights_nb(cmx *ymx, cmx *umx, cmx *gmx, cmx *dmx, double prior, cmx *wmx, const double *s2, int s2len, int nthreads, double *w)
{
    int ntag = ymx->nrow;
    int nlib = ymx->ncol;
    double *dfmat = R_Calloc((size_t) ntag*nlib, double);
    double *dvmat = R_Calloc((size_t) ntag*nlib, double);
    double *lvmat = R_Calloc((size_t) ntag*nlib, double);

    compute_adjust_mat(ymx, umx, gmx, dmx, prior, wmx, dfmat, dvmat, lvmat, nthreads);
    reduce_sample_weights(dvmat, dfmat, wmx, s2, s2len, ntag, nlib, w);

    R_Free(dfmat);
    R_Free(dvmat);
    R_Free(lvmat);
}

/* Binomial sample weights (DGEBIN). Derives coverage / proportion / scaled
 * weights from the raw successes and failures, computes the adjusted unit dev/df
 * with compute_adjust_mat_bin, then reduces with s2 = 1.
 *   ymx, y2mx : successes (counts), failures (counts2)
 *   umx       : fitted values (proportions)
 *   gmx       : design
 *   wmx       : compressed observation weights (used both to scale and to reduce)
 * Writes w (length nlib). Allocates and frees its own scratch. */
void sample_weights_binom(cmx *ymx, cmx *y2mx, cmx *umx, cmx *gmx, cmx *wmx, int nthreads, double *w)
{
    int ntag = ymx->nrow;
    int nlib = ymx->ncol;

    double *cover = R_Calloc((size_t) ntag*nlib, double);
    double *prop  = R_Calloc((size_t) ntag*nlib, double);
    double *wght0 = R_Calloc((size_t) ntag*nlib, double);
    double *yr  = R_Calloc(nlib, double);
    double *y2r = R_Calloc(nlib, double);
    double *wr  = R_Calloc(nlib, double);

    for(int g=0; g<ntag; ++g)
    {
        get_row(ymx, g, yr);
        get_row(y2mx, g, y2r);
        get_row(wmx, g, wr);
        double *covptr = cover + g;
        double *prpptr = prop + g;
        double *w0ptr  = wght0 + g;
        for(int j=0; j<nlib; ++j, covptr+=ntag, prpptr+=ntag, w0ptr+=ntag)
        {
            double cov = yr[j] + y2r[j];
            *covptr = cov;
            *prpptr = yr[j] / fmax(cov, 1.0);
            *w0ptr  = wr[j] * cov;
        }
    }
    R_Free(yr);
    R_Free(y2r);
    R_Free(wr);

    cmx propmx = make_cmx(prop, NULL, ntag, nlib, 0, 0);
    cmx covmx  = make_cmx(cover, NULL, ntag, nlib, 0, 0);
    cmx w0mx   = make_cmx(wght0, NULL, ntag, nlib, 0, 0);

    double *dfmat = R_Calloc((size_t) ntag*nlib, double);
    double *dvmat = R_Calloc((size_t) ntag*nlib, double);
    double *lvmat = R_Calloc((size_t) ntag*nlib, double);

    compute_adjust_mat_bin(&propmx, umx, gmx, &covmx, &w0mx, dfmat, dvmat, lvmat, nthreads);
    reduce_sample_weights(dvmat, dfmat, wmx, NULL, 0, ntag, nlib, w);

    R_Free(dfmat);
    R_Free(dvmat);
    R_Free(lvmat);
    R_Free(cover);
    R_Free(prop);
    R_Free(wght0);
}
