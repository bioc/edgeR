#include "edgeR.h"
#include <math.h>
#include <stdlib.h>

/* ============================================================================
 * Cox-Reid common-dispersion estimation in a single C call.
 *
 * Maximizes sum_tag adjustedProfileLik(disp) over the dispersion using a 1-D
 * Brent search (in the par = disp^0.25 parametrization, matching dispCoxReid).
 * Each objective evaluation fits all genes (fit_glm_mat) and sums the adjusted
 * profile likelihood (compute_adj_profile_ll). Fits are cold-started: an earlier
 * version warm-started from the previous dispersion's coefficients, but that
 * corrupted fits on degenerate (zero-count) genes and produced wrong
 * dispersions, so warm-starting was removed. This replaces the R optimize() loop
 * in dispCoxReid() (Brent here matches R's Brent_fmin) and is faster, especially
 * multi-threaded. See R/dispCoxReid.R.
 * ========================================================================== */

/* GLM-fit controls, matching adjustedProfileLik()/glmFit(); distinct from the
 * Brent (optimize) tolerance, which is supplied by the caller. */
#define COXREID_FIT_MAXIT 250
#define COXREID_FIT_TOL   1e-6

/* Optimization context. The dispersion is a scalar shared across all tags,
 * represented as a cmx of type 3 (read as dmat[0]); changing disp_val between
 * evaluations needs no re-allocation. */
typedef struct {
    cmx *ymx, *omx, *wmx, *gmx;   /* inputs (counts/offset/weights/design)   */
    cmx dmx;                      /* scalar-dispersion cmx, dmat = &disp_val */
    double disp_val;              /* current dispersion (= par^4)            */

    int ntag, nlib, ncoef;        /* rows, cols, coefficients                */
    int do_adjust, nthreads;      /* CR adjustment flag; fit threads         */

    /* work buffers */
    double *coef;                 /* coefficient output buffer               */
    double *mu, *apl;             /* fitted means; adjusted profile lik.     */
    double *dev;                  /* throwaway fit_glm_mat outputs  */
    int *iter;
    int *failed;
} cr_ctx;

/* One objective evaluation: sum of adjusted profile likelihood at dispersion
 * par^4. Each fit is cold-started (coef_start = NULL): warm-starting from the
 * previous dispersion's coefficients corrupts fits on degenerate (zero-count)
 * genes and gave wrong dispersions, so it is intentionally not used. */
static double eval_sum_apl(cr_ctx *c, double par)
{
    c->disp_val = par*par*par*par;   /* dispersion = par^4 */
    int method_val = 0;

    fit_glm_mat(c->ymx, c->omx, &c->dmx, c->wmx, c->gmx, COXREID_FIT_MAXIT, COXREID_FIT_TOL, NULL, c->coef, c->mu, c->dev, c->iter, c->failed, &method_val, c->nthreads);

    /* wrap the freshly fitted means as a plain (type 0) cmx for the APL */
    cmx umx = make_cmx(c->mu, NULL, c->ntag, c->nlib, 0, 0);

    compute_adj_profile_ll(c->ymx, &umx, &c->dmx, c->wmx, c->gmx, c->do_adjust, c->apl, c->nthreads);

    double s = 0;
    for (int tag=0; tag<c->ntag; ++tag)
    {
        s += c->apl[tag];
    }

    return s;
}

/* Brent minimizes, so negate (we maximize the summed APL). */
static double neg_eval(double par, void *info)
{
    return -eval_sum_apl((cr_ctx*)info, par);
}

/* Optimize the Cox-Reid common dispersion by a 1-D Brent search over
 * par = disp^0.25. Plain-C worker behind the .Call shim coxreid_disp (R_exports.c).
 *   ymx, omx, wmx, gmx : counts, compressed offsets/weights, design (as cmx)
 *   lower, upper       : Brent bounds in par = disp^0.25 space
 *   tol                : Brent (optimize) tolerance
 *   do_adjust          : include the Cox-Reid adjustment
 *   nthreads           : threads for the per-evaluation GLM fits
 * Allocates its own scratch, runs Brent, frees, and returns the optimal
 * dispersion (par_opt^4). */
double coxreid_disp_opt(cmx *ymx, cmx *omx, cmx *wmx, cmx *gmx, double lower, double upper, double tol, int do_adjust, int nthreads)
{
    cr_ctx c;
    c.ymx = ymx;
    c.gmx = gmx;
    c.omx = omx;
    c.wmx = wmx;

    c.ntag = ymx->nrow;
    c.nlib = ymx->ncol;
    c.ncoef = gmx->ncol;
    c.do_adjust = do_adjust;
    c.nthreads = nthreads;

    /* scalar-dispersion cmx (single value, repeated by row and column) */
    c.dmx = make_cmx(&c.disp_val, NULL, c.ntag, c.nlib, 3, 0);

    /* Scratch freed on the normal path below. A callee error() in fit_glm_mat or
     * compute_adj_profile_ll would longjmp past those R_Free (a one-time bounded
     * leak), but none can fire for a full-rank design: the rank-dependent solves
     * (dgesv, dtrtrs) require a singular design, and the remaining LAPACK calls are
     * illegal-argument guards -- a singular per-gene XtWX is factored (dsytrf
     * linfo>0) and floored, not errored. Matches the sibling workers' convention. */
    c.coef   = R_Calloc((size_t) c.ntag * c.ncoef, double);
    c.mu     = R_Calloc((size_t) c.ntag * c.nlib, double);
    c.apl    = R_Calloc(c.ntag, double);
    c.dev    = R_Calloc(c.ntag, double);
    c.iter   = R_Calloc(c.ntag, int);
    c.failed = R_Calloc(c.ntag, int);

    double par_opt = brent_fmin(lower, upper, neg_eval, &c, tol);

    R_Free(c.coef);
    R_Free(c.mu);
    R_Free(c.apl);
    R_Free(c.dev);
    R_Free(c.iter);
    R_Free(c.failed);

    return par_opt*par_opt*par_opt*par_opt;   /* dispersion = par^4 */
}

/* {abundance, original index} pair used to rank genes for top-gene selection. */
typedef struct
{
    double score;
    int idx;
} scored_row;

/* qsort comparator: descending score, ascending original index on ties. This
 * reproduces R's stable order(AveLogCPM, decreasing=TRUE), where equal scores
 * keep their ascending original-index order. */
static int cmp_scored_row(const void *a, const void *b)
{
    const scored_row *ra = (const scored_row*) a;
    const scored_row *rb = (const scored_row*) b;
    if(ra->score > rb->score) return -1;
    if(ra->score < rb->score) return 1;
    if(ra->idx < rb->idx) return -1;
    if(ra->idx > rb->idx) return 1;
    return 0;
}

/* Cox-Reid common dispersion over the most-abundant genes, in a single call.
 *
 * Folds the R block in glmQLFit.default (formerly lines 111-116 of
 * glmQLFTest.R) into C: it chooses the top fraction of genes by abundance and
 * runs the Cox-Reid Brent search (coxreid_disp_opt) over just those genes.
 * Reproduces
 *   df.residual <- ncol(y) - ncol(design)
 *   top.prop    <- chooseLowessSpan(ngenes*sqrt(df.residual), 20, 0.02)  (power 1/3)
 *   ntop        <- ceiling(top.prop * ngenes)
 *   i           <- order(AveLogCPM, decreasing=TRUE)[1:ntop]
 *   dispCoxReid(y[i,], design, offset[i,], weights[i,])
 * but deliberately omits dispCoxReid's min.row.sum filtering, systematicSubset
 * subsetting and nonEstimable design check (unnecessary on this path: the
 * selected genes are the highest-abundance genes and the design is validated
 * upstream by glmFit). chooseLowessSpan's closed form (limma) is
 * pmin(min.span + (1-min.span)*(small.n/n)^power, 1).
 *
 * inputs:
 *   ymx, omx, wmx : counts, compressed offsets, compressed weights (full cmx)
 *   gmx           : design (library-indexed; not subset by gene)
 *   avelogcpm     : length-nrow(ymx) abundance used to rank genes
 *   lower, upper  : Brent bounds in par = disp^0.25 space
 *   tol           : Brent (optimize) tolerance
 *   do_adjust     : include the Cox-Reid adjustment
 *   nthreads      : threads for the per-evaluation GLM fits
 *
 * output: the optimal common dispersion (par_opt^4).
 * side effects: allocates and frees its own selection scratch and a contiguous
 *               double copy of the selected counts; inputs read only.
 */
double coxreid_disp_top_opt(cmx *ymx, cmx *omx, cmx *wmx, cmx *gmx, const double *avelogcpm, double lower, double upper, double tol, int do_adjust, int nthreads)
{
    int ngenes = ymx->nrow;
    int nlib = ymx->ncol;
    int ncoef = gmx->ncol;
    int df_residual = nlib - ncoef;

    double top_prop = fmin(0.02 + 0.98*pow(20.0/((double)ngenes*sqrt((double)df_residual)), 1.0/3.0), 1.0);
    int ntop = (int) ceil(top_prop*(double)ngenes);
    if(ntop < 1) ntop = 1;
    if(ntop > ngenes) ntop = ngenes;

    /* rank genes by abundance and keep the top ntop physical row indices */
    scored_row *sr = R_Calloc(ngenes, scored_row);
    for(int g=0; g<ngenes; ++g)
    {
        sr[g].score = avelogcpm[g];
        sr[g].idx = g;
    }
    qsort(sr, (size_t) ngenes, sizeof(scored_row), cmp_scored_row);

    int *rowsel = R_Calloc(ntop, int);
    for(int k=0; k<ntop; ++k)
    {
        rowsel[k] = sr[k].idx;
    }
    R_Free(sr);

    /* Materialize the selected counts into a contiguous double matrix. Offsets
     * and weights are only ever read through get_row (rowsel-aware), so they can
     * stay zero-copy subset views; but compute_nbdev_sum() reads y->dmat
     * directly (bypassing get_row), so the counts must be a plain, contiguous,
     * double type-0 matrix -- a row-subset view is read there at the wrong
     * stride, and an integer count matrix has no dmat at all. get_row converts
     * integer counts to double, matching the former R path (which passed
     * y[i,,drop=FALSE] as a fresh matrix to dispCoxReid). */
    double *ycopy = R_Calloc((size_t) ntop*nlib, double);
    double *rowbuf = R_Calloc(nlib, double);
    for(int k=0; k<ntop; ++k)
    {
        get_row(ymx, rowsel[k], rowbuf);
        for(int lib=0; lib<nlib; ++lib)
        {
            ycopy[k + (size_t) lib*ntop] = rowbuf[lib];
        }
    }
    R_Free(rowbuf);

    cmx ymat = make_cmx(ycopy, NULL, ntop, nlib, 0, 0);
    cmx osub = subset_cmx(*omx, rowsel, ntop);  /* design is not gene-indexed */
    cmx wsub = subset_cmx(*wmx, rowsel, ntop);

    double disp = coxreid_disp_opt(&ymat, &osub, &wsub, gmx, lower, upper, tol, do_adjust, nthreads);

    R_Free(ycopy);
    R_Free(rowsel);

    return disp;
}
