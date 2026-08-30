#include "edgeR.h"

/* ============================================================================
 * ftest.c
 *
 * Single-call quasi-likelihood F-test with TREAT / UPSHOT p-values, shared by
 * the NB GLM (glm_ftest_mat, glmQLFTest) and the binomial (bin_ftest_mat,
 * binQLFTest) paths.  The TREAT/UPSHOT engine -- base null refit, the QL
 * F-statistic and its p-value, the UPSHOT 8-node Gauss-quadrature loop, the
 * 2^ncoef boundary-corner enumeration with the running-min deviance, and the
 * final qf() back-transform -- is written once in ftest_engine().  Each
 * distribution supplies only its input prep + base fit (in the *_mat entry), its
 * per-corner null fitter, its gather, and its boundary-offset base, through a
 * small ftest_dist vtable.  Merged from the former glm_ftest.c + bin_ftest.c;
 * results are bit-identical to those files.
 *
 * Memory: all scratch is R_Calloc / R_Free (repo convention), allocated up front
 * in each *_mat entry and freed at its single exit.  Caveat: the fitters can
 * error() (longjmp) on a catastrophic numerical failure (singular design, QR /
 * Cholesky / dgesv), which would skip the R_Free below; that one-time bounded
 * leak on a fatal error path matches the fitters' own behaviour and aborts the
 * call anyway.
 * ==========================================================================*/

/* NB fitter convergence settings, identical to the R .Call(.cxx_fit_glm, ...) in
 * glmFit.default (250 iterations, tolerance 1e-6). */
static const int    MAXIT = 250;
static const double TOL   = 1e-6;

/* binomial fitter convergence settings, identical to the R .Call(.cxx_bin_fit,
 * ...) in binFit.default (oneway 50/1e-8, IWLS 100/1e-6). */
static const int    MAXIT_OW = 50;
static const double TOL_OW   = 1e-8;
static const int    MAXIT_IW = 100;
static const double TOL_IW   = 1e-6;

/* UPSHOT Gauss-quadrature rule (8 folded nodes + centre weight), copied verbatim
 * from the former R constants in glmQLFTest / binQLFTest; used only by
 * ftest_engine() below. */
static const double gq_nodes[8]   = {0.1784842, 0.3512318, 0.5126905, 0.6576712, 0.7815140, 0.8802392, 0.9506755, 0.9905755};
static const double gq_weights[8] = {0.1765627, 0.1680041, 0.1540458, 0.1351364, 0.1118839, 0.0850362, 0.0554595, 0.0241483};
static const double gq_centre     = 0.08972310;

/* ----------------------------------------------------------------------------
 * Distribution vtable.  The genuine differences between the NB and binomial
 * F-tests are narrow: (1) input prep + base null fit (done in each *_mat entry
 * before the engine runs), (2) the per-corner null refit (fit), (3) gathering
 * the selected genes' inputs (gather), and (4) the boundary-offset base
 * (offbase: the NB library-size offset gathered per node, or NULL for the
 * binomial which has none).  Everything else lives in ftest_engine().
 * ------------------------------------------------------------------------- */
typedef struct ftest_dist ftest_dist;
struct ftest_dist
{
    /* dimensions */
    int ntag, nlib, ncoef, nthreads;
    /* inputs (shared, from R) */
    double *deviance, *s2post, *df_total, *lfc;
    cmx *logFCt;
    double *design1;
    cmx d0mx;
    /* outputs (length ntag) */
    double *Fout, *Pout;
    /* shared scratch */
    double *pfbase, *nodep, *mindev, *dev_scr, *offset_sub, *sign, *coefval;
    int *sel;
    /* fitter scratch (shared shape); method is the in/out oneway-vs-IWLS flag */
    double *coef_scr, *mu_scr;
    int *iter_scr, *failed_scr, method;
    /* seam (4): per-gene base offset for the eta seed (gathered), or NULL */
    const double *offbase;
    /* seams (2)/(3): refit the null for the gathered subset; gather the subset */
    void (*gather)(ftest_dist *d, int nsub, const int *sel);
    void (*fit)(ftest_dist *d, cmx *offsubmx, double *dev_out);
    int upshot;
};

/* design1 = design[,coef] (coef order); design0 = design[,-coef] (original
 * order).  Byte-identical block shared by both entries. */
static void split_design(cmx *design, const int *coef, int ncoef, int nlib, int nbeta, double *design1, double *design0, int *is_tested)
{
    for (int c = 0; c < nbeta; ++c)
    {
        is_tested[c] = 0;
    }
    for (int c = 0; c < ncoef; ++c)
    {
        int dc = coef[c] - 1;
        is_tested[dc] = 1;
        for (int l = 0; l < nlib; ++l)
        {
            design1[l + c * nlib] = design->dmat[l + dc * nlib];
        }
    }
    int cc = 0;
    for (int c = 0; c < nbeta; ++c)
    {
        if (is_tested[c])
        {
            continue;
        }
        for (int l = 0; l < nlib; ++l)
        {
            design0[l + cc * nlib] = design->dmat[l + c * nlib];
        }
        ++cc;
    }
}

/* ----------------------------------------------------------------------------
 * The shared TREAT / UPSHOT engine.  On entry d->dev_scr holds the base null
 * deviance (from the entry's base fit).  Computes the base QL F-statistic and
 * p-value, then (if any lfc > 0) the UPSHOT quadrature TREAT loop and the qf()
 * back-transform, delegating gather / null-refit to the vtable.
 * ------------------------------------------------------------------------- */
static void ftest_engine(ftest_dist *d)
{
    int ntag = d->ntag, nlib = d->nlib, ncoef = d->ncoef;

    /* Base QL F-statistic and its p-value.  LR0 = null_dev - full_dev is the
       deviance drop from adding the tested coefs, clamped at 0.  F = LR0 / ncoef
       / s2post divides the drop by the number of tested coefs and the
       EB-moderated quasi-dispersion s2post; the p-value is the upper tail of
       F(df1=ncoef, df2=df_total).  The two FALSE args to pf() are
       lower.tail=FALSE (upper tail) and log.p=FALSE. */
    for (int g = 0; g < ntag; ++g)
    {
        double LR0 = fmax2(d->dev_scr[g] - d->deviance[g], 0.0);
        double F0 = LR0 / ncoef / d->s2post[g];
        d->pfbase[g] = pf(F0, (double) ncoef, d->df_total[g], FALSE, FALSE);
        d->Fout[g] = F0;
        d->Pout[g] = d->pfbase[g];
    }

    /* Any positive fold-change threshold?  If every lfc == 0 the plain QL F-test
       above is the final answer and the TREAT/UPSHOT block is skipped. */
    int any_lfc = 0;
    for (int c = 0; c < ncoef; ++c)
    {
        if (d->lfc[c] > 0)
        {
            any_lfc = 1;
        }
    }

    if (any_lfc)
    {
        /* UPSHOT integrates the TREAT p-value over a Uniform[-lfc,lfc] prior on
           the true logFC by Gauss quadrature: 8 folded nodes here plus the centre
           node seeded below.  Plain TREAT (upshot=0) uses the single node s=1 (the
           full threshold).  nboundary = 2^ncoef sign patterns of the tested coefs
           (the corners of the [-lfc,lfc] box). */
        int nnode = d->upshot ? 8 : 1;
        int nboundary = 1 << ncoef;
        if (d->upshot)
        {
            /* centre node (true logFC = 0): weight gq_centre times the base p */
            for (int g = 0; g < ntag; ++g)
            {
                d->Pout[g] = gq_centre * d->pfbase[g];
            }
        }

        for (int i = 0; i < nnode; ++i)
        {
            double s = d->upshot ? gq_nodes[i] : 1.0;

            /* Select genes whose |logFC| exceeds the node threshold lfc*s in at
               least one tested coef; the rest keep node p-value 1 (a gene inside
               the interval cannot reject H0: |logFC| <= lfc).  The *M_LN2 factor
               (log2 -> natural log) appears on both sides and cancels; it is kept
               to mirror the R helpers .glmtreat / .bintreat. */
            int nsub = 0;
            for (int g = 0; g < ntag; ++g)
            {
                d->nodep[g] = 1.0;
                int hit = 0;
                const double *lp = d->logFCt->dmat + g;
                for (int c = 0; c < ncoef; ++c, lp += ntag)
                {
                    if (fabs(*lp) * M_LN2 > d->lfc[c] * s * M_LN2)
                    {
                        hit = 1;
                        break;
                    }
                }
                if (hit)
                {
                    d->sel[nsub] = g;
                    ++nsub;
                }
            }

            if (nsub > 0)
            {
                /* gather the selected genes' inputs (distribution-specific) and
                   view the offset scratch as an nsub x nlib matrix */
                d->gather(d, nsub, d->sel);

                cmx offsubmx = make_cmx(d->offset_sub, NULL, nsub, nlib, 0, 0);

                /* mindev tracks the running minimum deviance over the 2^ncoef
                   corners: TREAT takes the least-conservative admissible boundary
                   point.  Seed with DBL_MAX (a finite value larger than any
                   achievable deviance) so the first corner always replaces it. */
                for (int j = 0; j < nsub; ++j)
                {
                    d->mindev[j] = DBL_MAX;
                }

                /* Enumerate the 2^ncoef corners of the [-lfc,lfc] box.  Bit c of b
                   gives the sign of the c-th tested coef: sign = 2*bit - 1 = +/-1. */
                for (int b = 0; b < nboundary; ++b)
                {
                    for (int c = 0; c < ncoef; ++c)
                    {
                        d->sign[c] = 2.0 * ((b >> c) & 1) - 1.0;
                    }
                    /* per-gene offset = base offset + design1 %*% (pmin(lfc,|logFC|)*sign) */
                    for (int j = 0; j < nsub; ++j)
                    {
                        int g = d->sel[j];
                        /* Project the |logFC| MLE onto [0, lfc*s] (natural-log
                           scale) then apply the corner sign: the tested coef is
                           pinned at min(lfc, |logFC|)*sign, the closest admissible
                           boundary point (pmin(lfc, |logFC|)*boundary in R). */
                        const double *lp = d->logFCt->dmat + g;
                        for (int c = 0; c < ncoef; ++c, lp += ntag)
                        {
                            double thr = fmin2(d->lfc[c] * s * M_LN2, fabs(*lp) * M_LN2);
                            d->coefval[c] = thr * d->sign[c];
                        }
                        /* eta = design1 %*% coefval pins the tested coefs; the NB
                           path seeds the accumulator with the library-size offset
                           (d->offbase), the binomial with 0 (d->offbase == NULL).
                           The base offset MUST stay the accumulator seed -- moving
                           it out of the sum would change the FP rounding.  The
                           pointer-walks (stride nsub) avoid the j + l*nsub int
                           product overflowing. */
                        const double *obp = d->offbase ? d->offbase + j : NULL;
                        double *osp = d->offset_sub + j;
                        for (int l = 0; l < nlib; ++l, osp += nsub)
                        {
                            double e = obp ? *obp : 0.0;
                            for (int c = 0; c < ncoef; ++c)
                            {
                                e += d->design1[l + c * nlib] * d->coefval[c];
                            }
                            *osp = e;
                            if (obp)
                            {
                                obp += nsub;
                            }
                        }
                    }

                    /* refit the null (design0) for the selected genes with this
                       corner's offset; dev_scr is the boundary deviance */
                    d->fit(d, &offsubmx, d->dev_scr);

                    /* keep the smallest boundary deviance seen so far per gene */
                    for (int j = 0; j < nsub; ++j)
                    {
                        d->mindev[j] = fmin2(d->mindev[j], d->dev_scr[j]);
                    }
                }

                /* Node p-value: the same QL F-test as the base, but against the
                   best (minimum-deviance) boundary null instead of the point null. */
                for (int j = 0; j < nsub; ++j)
                {
                    int g = d->sel[j];
                    double LR = fmax2(d->mindev[j] - d->deviance[g], 0.0);
                    double Fs = LR / ncoef / d->s2post[g];
                    d->nodep[g] = pf(Fs, (double) ncoef, d->df_total[g], FALSE, FALSE);
                }
            }

            /* Accumulate this node into the p-value.  UPSHOT: add the quadrature
               term gq_weights[i]*nodep (the sum already holds the centre term).
               Plain TREAT: average the base and worst-case (s=1, single node)
               p-values. */
            if (d->upshot)
            {
                for (int g = 0; g < ntag; ++g)
                {
                    d->Pout[g] += gq_weights[i] * d->nodep[g];
                }
            }
            else
            {
                for (int g = 0; g < ntag; ++g)
                {
                    d->Pout[g] = (d->pfbase[g] + d->nodep[g]) / 2.0;
                }
            }
        }

        /* The UPSHOT/TREAT p-value has no natural F statistic, so recover one by
           inverting the upper-tail F cdf: Fout = qf(p, ncoef, df_total). */
        for (int g = 0; g < ntag; ++g)
        {
            d->Fout[g] = qf(d->Pout[g], (double) ncoef, d->df_total[g], FALSE, FALSE);
        }
    }
}

/* ------------------------------ NB GLM path ------------------------------- */
typedef struct
{
    ftest_dist d;
    cmx *ycnt, *offset, *disp, *weights;
    double *counts_sub, *disp_sub, *weights_sub, *offbase_sub;
    double *ybuf, *obuf, *dbuf, *wbuf;
    cmx csubmx, dsubmx, wsubmx;
} glm_dist;

/* gather selected genes into dense type-0 rows before the subset refit. The counts
   are densified (NOT a subset_cmx view) because the NB deviance kernel
   (compute_nbdev.c) reads the count matrix directly with an nrow-strided pointer
   walk, which a subset view (physical stride != nrow) would break; dispersion,
   weights and the library-size offset base are densified here too. */
static void glm_gather(ftest_dist *dd, int nsub, const int *sel)
{
    glm_dist *g = (glm_dist *) dd;
    int nlib = dd->nlib;
    for (int j = 0; j < nsub; ++j)
    {
        int gg = sel[j];
        get_row4(g->ycnt, g->offset, g->disp, g->weights, gg, g->ybuf, g->obuf, g->dbuf, g->wbuf);
        double *cp = g->counts_sub + j, *ob = g->offbase_sub + j, *dp = g->disp_sub + j, *wp = g->weights_sub + j;
        for (int l = 0; l < nlib; ++l, cp += nsub, ob += nsub, dp += nsub, wp += nsub)
        {
            *cp = g->ybuf[l];
            *ob = g->obuf[l];
            *dp = g->dbuf[l];
            *wp = g->wbuf[l];
        }
    }
    g->csubmx = make_cmx(g->counts_sub, NULL, nsub, nlib, 0, 0);
    g->dsubmx = make_cmx(g->disp_sub, NULL, nsub, nlib, 0, 0);
    g->wsubmx = make_cmx(g->weights_sub, NULL, nsub, nlib, 0, 0);
}

static void glm_fit(ftest_dist *dd, cmx *offsubmx, double *dev_out)
{
    glm_dist *g = (glm_dist *) dd;
    fit_glm_mat(&g->csubmx, offsubmx, &g->dsubmx, &g->wsubmx, &dd->d0mx, MAXIT, TOL, NULL, dd->coef_scr, dd->mu_scr, dev_out, dd->iter_scr, dd->failed_scr, &dd->method, dd->nthreads);
}

/* Single-call negative-binomial QL F-test with TREAT/UPSHOT p-values.
 *
 * inputs:
 *   ycnt     count matrix, ntag x nlib (dense cmx, int or double)
 *   offset   log offset (lib size), compressed cmx (row-repeat form)
 *   disp     NB dispersion, compressed cmx (per-gene column-repeat, or scalar)
 *   weights  observation weights, compressed cmx (all-1 when R weights were NULL)
 *   design   design matrix, nlib x nbeta (dense cmx, coerced double)
 *   coef     1-based indices of the tested design columns, length ncoef
 *   ncoef    number of tested coefficients
 *   deviance full-model residual deviance per gene, length ntag
 *   s2post   EB-moderated quasi-dispersion per gene, length ntag
 *   df_total total degrees of freedom per gene, length ntag
 *   lfc      log2 fold-change thresholds aligned to coef, length ncoef
 *   logFCt   log2 fold-changes, ntag x ncoef (dense cmx), aligned to coef
 *   upshot   1 for the UPSHOT quadrature, 0 for the plain averaged TREAT
 *   nthreads thread count forwarded to fit_glm_mat
 *
 * outputs (caller-allocated, length ntag):
 *   Fout     QL F-statistic
 *   Pout     QL / TREAT p-value
 */
void glm_ftest_mat (cmx *ycnt, cmx *offset, cmx *disp, cmx *weights, cmx *design, int *coef, int ncoef, double *deviance, double *s2post, double *df_total, double *lfc, cmx *logFCt, int upshot, double *Fout, double *Pout, int nthreads)
{
    int ntag = ycnt->nrow;
    int nlib = ycnt->ncol;
    int nbeta = design->ncol;
    int ncoef0 = nbeta - ncoef;

    /* scratch (all R_Calloc, freed at the single exit) */
    double *design1 = R_Calloc(nlib * ncoef, double);
    double *design0 = R_Calloc(nlib * ncoef0, double);
    int *is_tested = R_Calloc(nbeta, int);
    double *counts_sub = R_Calloc((size_t) ntag * nlib, double);
    double *disp_sub = R_Calloc((size_t) ntag * nlib, double);
    double *weights_sub = R_Calloc((size_t) ntag * nlib, double);
    double *offbase_sub = R_Calloc((size_t) ntag * nlib, double);
    double *offset_sub = R_Calloc((size_t) ntag * nlib, double);
    double *pfbase = R_Calloc(ntag, double);
    double *nodep = R_Calloc(ntag, double);
    double *mindev = R_Calloc(ntag, double);
    double *dev_scr = R_Calloc(ntag, double);
    double *mu_scr = R_Calloc((size_t) ntag * nlib, double);
    double *coef_scr = R_Calloc((size_t) ntag * ncoef0, double);
    int *iter_scr = R_Calloc(ntag, int);
    int *failed_scr = R_Calloc(ntag, int);
    int *sel = R_Calloc(ntag, int);
    double *ybuf = R_Calloc(nlib, double);
    double *obuf = R_Calloc(nlib, double);
    double *dbuf = R_Calloc(nlib, double);
    double *wbuf = R_Calloc(nlib, double);
    double *sign = R_Calloc(ncoef, double);
    double *coefval = R_Calloc(ncoef, double);

    split_design(design, coef, ncoef, nlib, nbeta, design1, design0, is_tested);

    cmx d0mx = make_cmx(design0, NULL, nlib, ncoef0, 0, 0);

    /* base null refit over all genes (folds glmQLFTest.R:249-258); the compressed
       offset / dispersion / weights cmx are passed straight through, exactly as
       glmFit's base call feeds fit_glm_mat.  method (in/out) chooses oneway
       closed form vs Levenberg IWLS and carries into the corner refits. */
    int method = 0;
    fit_glm_mat(ycnt, offset, disp, weights, &d0mx, MAXIT, TOL, NULL, coef_scr, mu_scr, dev_scr, iter_scr, failed_scr, &method, nthreads);

    glm_dist g;
    g.d.ntag = ntag;
    g.d.nlib = nlib;
    g.d.ncoef = ncoef;
    g.d.nthreads = nthreads;
    g.d.deviance = deviance;
    g.d.s2post = s2post;
    g.d.df_total = df_total;
    g.d.lfc = lfc;
    g.d.logFCt = logFCt;
    g.d.design1 = design1;
    g.d.d0mx = d0mx;
    g.d.Fout = Fout;
    g.d.Pout = Pout;
    g.d.pfbase = pfbase;
    g.d.nodep = nodep;
    g.d.mindev = mindev;
    g.d.dev_scr = dev_scr;
    g.d.offset_sub = offset_sub;
    g.d.sign = sign;
    g.d.coefval = coefval;
    g.d.sel = sel;
    g.d.coef_scr = coef_scr;
    g.d.mu_scr = mu_scr;
    g.d.iter_scr = iter_scr;
    g.d.failed_scr = failed_scr;
    g.d.method = method;
    g.d.offbase = offbase_sub;
    g.d.gather = glm_gather;
    g.d.fit = glm_fit;
    g.d.upshot = upshot;
    g.ycnt = ycnt;
    g.offset = offset;
    g.disp = disp;
    g.weights = weights;
    g.counts_sub = counts_sub;
    g.disp_sub = disp_sub;
    g.weights_sub = weights_sub;
    g.offbase_sub = offbase_sub;
    g.ybuf = ybuf;
    g.obuf = obuf;
    g.dbuf = dbuf;
    g.wbuf = wbuf;

    ftest_engine(&g.d);

    R_Free(design1);
    R_Free(design0);
    R_Free(is_tested);
    R_Free(counts_sub);
    R_Free(disp_sub);
    R_Free(weights_sub);
    R_Free(offbase_sub);
    R_Free(offset_sub);
    R_Free(pfbase);
    R_Free(nodep);
    R_Free(mindev);
    R_Free(dev_scr);
    R_Free(mu_scr);
    R_Free(coef_scr);
    R_Free(iter_scr);
    R_Free(failed_scr);
    R_Free(sel);
    R_Free(ybuf);
    R_Free(obuf);
    R_Free(dbuf);
    R_Free(wbuf);
    R_Free(sign);
    R_Free(coefval);

    return;
}

/* ----------------------------- binomial path ------------------------------ */
typedef struct
{
    ftest_dist d;
    double *prop_full, *w0_full;
    cmx propsubmx, w0submx;   /* zero-copy row-subset views of prop_full / w0_full */
} bin_dist;

/* selected-gene proportions and coverage-weights as zero-copy row-subset views
   (subset_cmx) over the full dense prop / weight matrices */
static void bin_gather(ftest_dist *dd, int nsub, const int *sel)
{
    bin_dist *b = (bin_dist *) dd;
    int nlib = dd->nlib, ntag = dd->ntag;
    b->propsubmx = subset_cmx(make_cmx(b->prop_full, NULL, ntag, nlib, 0, 0), sel, nsub);
    b->w0submx   = subset_cmx(make_cmx(b->w0_full,   NULL, ntag, nlib, 0, 0), sel, nsub);
}

static void bin_fit(ftest_dist *dd, cmx *offsubmx, double *dev_out)
{
    bin_dist *b = (bin_dist *) dd;
    bin_fit_mat(&b->propsubmx, offsubmx, &b->w0submx, &dd->d0mx, MAXIT_OW, TOL_OW, MAXIT_IW, TOL_IW, dd->coef_scr, dd->mu_scr, dev_out, dd->iter_scr, dd->failed_scr, &dd->method, dd->nthreads);
}

/* Single-call binomial QL F-test with TREAT/UPSHOT p-values.
 *
 * inputs:
 *   ycnt     successes (counts) matrix, ntag x nlib (dense cmx, int or double)
 *   zcnt     failures  (counts2) matrix, ntag x nlib
 *   wts      observation weights, ntag x nlib, or unused when has_w == 0
 *   has_w    1 if wts holds weights, 0 if weights were NULL on the R side
 *   design   design matrix, nlib x nbeta (dense cmx, coerced double)
 *   coef     1-based indices of the tested design columns, length ncoef
 *   ncoef    number of tested coefficients
 *   deviance full-model residual deviance per gene, length ntag
 *   s2post   EB-moderated quasi-dispersion per gene, length ntag
 *   df_total total degrees of freedom per gene, length ntag
 *   lfc      log2 fold-change thresholds aligned to coef, length ncoef
 *   logFCt   log2 fold-changes, ntag x ncoef (dense cmx), aligned to coef
 *   upshot   1 for the UPSHOT quadrature, 0 for the plain averaged TREAT
 *   nthreads thread count forwarded to bin_fit_mat
 *
 * outputs (caller-allocated, length ntag):
 *   Fout     QL F-statistic
 *   Pout     QL / TREAT p-value
 */
void bin_ftest_mat (cmx *ycnt, cmx *zcnt, cmx *wts, int has_w, cmx *design, int *coef, int ncoef, double *deviance, double *s2post, double *df_total, double *lfc, cmx *logFCt, int upshot, double *Fout, double *Pout, int nthreads)
{
    int ntag = ycnt->nrow;
    int nlib = ycnt->ncol;
    int nbeta = design->ncol;
    int ncoef0 = nbeta - ncoef;

    /* scratch (all R_Calloc, freed at the single exit) */
    double *design1 = R_Calloc(nlib * ncoef, double);
    double *design0 = R_Calloc(nlib * ncoef0, double);
    int *is_tested = R_Calloc(nbeta, int);
    double *prop_full = R_Calloc((size_t) ntag * nlib, double);
    double *w0_full = R_Calloc((size_t) ntag * nlib, double);
    double *off_sub = R_Calloc((size_t) ntag * nlib, double);
    double *pfbase = R_Calloc(ntag, double);
    double *nodep = R_Calloc(ntag, double);
    double *mindev = R_Calloc(ntag, double);
    double *dev_scr = R_Calloc(ntag, double);
    double *mu_scr = R_Calloc((size_t) ntag * nlib, double);
    double *coef_scr = R_Calloc((size_t) ntag * ncoef0, double);
    int *iter_scr = R_Calloc(ntag, int);
    int *failed_scr = R_Calloc(ntag, int);
    int *sel = R_Calloc(ntag, int);
    double *ybuf = R_Calloc(nlib, double);
    double *zbuf = R_Calloc(nlib, double);
    double *wbuf = R_Calloc(nlib, double);
    double *sign = R_Calloc(ncoef, double);
    double *coefval = R_Calloc(ncoef, double);

    split_design(design, coef, ncoef, nlib, nbeta, design1, design0, is_tested);

    /* proportions and coverage-scaled weights, per gene (mirrors binFit.default).
       coverage = successes + failures; the binomial weight is the prior weight
       times coverage (coverage alone when R weights were NULL); the response is
       the proportion p = y / coverage.  fmax2(cov, 1) guards a zero-coverage
       observation (its p is 0 and weight 0).  pointer-walk strides by ntag so the
       g + l*ntag int product cannot overflow. */
    for (int g = 0; g < ntag; ++g)
    {
        get_row(ycnt, g, ybuf);
        get_row(zcnt, g, zbuf);
        if (has_w)
        {
            get_row(wts, g, wbuf);
        }
        double *pp = prop_full + g, *wp = w0_full + g;
        for (int l = 0; l < nlib; ++l, pp += ntag, wp += ntag)
        {
            double cov = ybuf[l] + zbuf[l];
            double w = has_w ? (wbuf[l] * cov) : cov;
            double p = ybuf[l] / fmax2(cov, 1.0);
            *pp = p;
            *wp = w;
        }
    }

    /* dense cmx views over the full-gene proportions / weights and design0 */
    cmx propmx = make_cmx(prop_full, NULL, ntag, nlib, 0, 0);
    cmx w0mx = make_cmx(w0_full, NULL, ntag, nlib, 0, 0);
    cmx d0mx = make_cmx(design0, NULL, nlib, ncoef0, 0, 0);

    /* zero offset as a type-3 (repeat row and column) cmx over a single 0.0 */
    double dzero = 0.0;
    cmx zoffmx = make_cmx(&dzero, NULL, ntag, nlib, 3, 0);

    /* Base null model: fit design0 (the untested columns) to every gene with a
       zero offset and coverage weights; dev_scr receives the null deviance.
       method (in/out) carries into the corner refits. */
    int method = 0;
    bin_fit_mat(&propmx, &zoffmx, &w0mx, &d0mx, MAXIT_OW, TOL_OW, MAXIT_IW, TOL_IW, coef_scr, mu_scr, dev_scr, iter_scr, failed_scr, &method, nthreads);

    bin_dist b;
    b.d.ntag = ntag;
    b.d.nlib = nlib;
    b.d.ncoef = ncoef;
    b.d.nthreads = nthreads;
    b.d.deviance = deviance;
    b.d.s2post = s2post;
    b.d.df_total = df_total;
    b.d.lfc = lfc;
    b.d.logFCt = logFCt;
    b.d.design1 = design1;
    b.d.d0mx = d0mx;
    b.d.Fout = Fout;
    b.d.Pout = Pout;
    b.d.pfbase = pfbase;
    b.d.nodep = nodep;
    b.d.mindev = mindev;
    b.d.dev_scr = dev_scr;
    b.d.offset_sub = off_sub;
    b.d.sign = sign;
    b.d.coefval = coefval;
    b.d.sel = sel;
    b.d.coef_scr = coef_scr;
    b.d.mu_scr = mu_scr;
    b.d.iter_scr = iter_scr;
    b.d.failed_scr = failed_scr;
    b.d.method = method;
    b.d.offbase = NULL;
    b.d.gather = bin_gather;
    b.d.fit = bin_fit;
    b.d.upshot = upshot;
    b.prop_full = prop_full;
    b.w0_full = w0_full;

    ftest_engine(&b.d);

    R_Free(design1);
    R_Free(design0);
    R_Free(is_tested);
    R_Free(prop_full);
    R_Free(w0_full);
    R_Free(off_sub);
    R_Free(pfbase);
    R_Free(nodep);
    R_Free(mindev);
    R_Free(dev_scr);
    R_Free(mu_scr);
    R_Free(coef_scr);
    R_Free(iter_scr);
    R_Free(failed_scr);
    R_Free(sel);
    R_Free(ybuf);
    R_Free(zbuf);
    R_Free(wbuf);
    R_Free(sign);
    R_Free(coefval);

    return;
}
