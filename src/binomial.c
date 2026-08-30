
#include "edgeR.h"

static const double THRESH = 30.0;
static const double MTHRESH = -30.0;
static const double DOUEPS = 1e-16;
static const double INVEPS = 1e16;

/** 
 * Helper functions for Binomial distribution
 * convert from family.c in package stat 
 */

/**
 * Evaluate x/(1 - x). 
 * An inline function is used so that x is evaluated once only.
 */
static double domx(double x)
{
    return x /(1 - x);
}

/**
 * Evaluate x/(1 + x). 
 */
static double dopx(double x) 
{
    return x /(1 + x);
}

/**
 * Evaluate logit(x) = log(x/(1 - x)). 
 */
static double logit(double mu) 
{
    return log(domx(mu));
}

/**
 * Evaluate bounded logit(x): log(x/(1 - x)), clamped to MTHRESH or THRESH
 * when mu is within DOUEPS of 0 or 1.
 */
static double logit2(double mu)
{
    return (fabs(mu) < DOUEPS) ? MTHRESH : ( (fabs(1-mu) < DOUEPS) ? THRESH : log(domx(mu)));
}

/**
 * Evaluate sigmod(x) = exp(x)/(1 + exp(x)), which is the inverse of logit. 
 */
static double sigmod(double eta)
{
    double tmp;
    tmp = (eta < MTHRESH) ? DOUEPS : ((eta > THRESH) ? INVEPS : exp(eta));
    return dopx(tmp);
}

/**
 * Evaluate derivate of sigmod(x) = exp(x)/(1 + exp(x))
 */
static double dsigmod(double eta)
{
    double tmp, opexp;

    if(eta < MTHRESH || eta > THRESH)
    {
        tmp = DOUEPS;
    }
    else
    {
        tmp   = exp(eta);
        opexp = 1 + tmp;
        tmp  = tmp / (opexp * opexp);
    }

    return tmp;
}

/**
 * Evaluate y * log(y)
 */
double y_log_y(double y, double mu)
{
    return (fabs(y) > DOUEPS && fabs(mu) > DOUEPS) ? (y * log(y/mu)) : 0;
}

/* this function fits one group design of binomial models for matrix input
 * 
 * common input:
 * y       matrix of proportions
 * weights compressed matrix of total coverage
 * offsets compressed matrix of offsets, only used in treat
 * 
 * maxit   max iteration
 * tol     tolerance
 * 
 * outputs:
 * prob  fitted probabilities
 * coef  fitted coefficients
 * dev   deviance for one group
 * conv  index of convergence
 *
 * No return value; parallelises over tags and writes only into the
 * caller-owned output arrays (prob, coef, dev, conv).
 */
/* Fit one binomial group/gene from the length-n gathered proportions y, offsets o and
 * coverage weights w, plus the group's offset-is-zero flag. Writes the fitted proportions
 * mu_out[n], the group logit *coef_out, the convergence flag *conv_out (0 on hitting maxit),
 * and the group's unclamped residual deviance *dev_out. Shared by bin_one_group_mat and
 * bin_one_way_mat; each caller scatters mu_out/coef/conv/dev into its own output layout and
 * applies its own fmax2(dev,0) clamp. offzero || degenerate-p takes the closed-form weighted
 * mean; otherwise Newton-Raphson on z=logit(p) with offsets. */
static void bin_group_fit(const double *y, const double *o, const double *w, int n, int offzero, int maxit, double tol, double *mu_out, double *coef_out, int *conv_out, double *dev_out)
{
    double sum_counts = 0;
    double sum_weights = 0;
    double dev = 0;
    int conv = 1;
    for (int lib = 0; lib < n; ++lib)
    {
        sum_counts += y[lib] * w[lib];
        sum_weights += w[lib];
    }
    double p = (fabs(sum_weights) < DOUEPS) ? 0.5 : sum_counts / sum_weights;

    if (offzero || fabs(p) < DOUEPS || fabs(1 - p) < DOUEPS || fabs(sum_weights) < DOUEPS)
    {
        for (int lib = 0; lib < n; ++lib)
        {
            mu_out[lib] = p;
            dev += 2 * w[lib] * (y_log_y(y[lib], p) + y_log_y(1 - y[lib], 1 - p));
        }
        *coef_out = logit2(p);
    }
    else
    {
        int iter = 1;
        double eta = 0;
        double step = 0;
        double z = logit(p);
        while (++iter)
        {
            double dl = 0;
            double info = 0;
            for (int lib = 0; lib < n; ++lib)
            {
                eta = z + o[lib];
                dl += w[lib] * (y[lib] - sigmod(eta));
                info += w[lib] * dsigmod(eta);
            }
            step = dl / info;
            z += step;
            if (fabs(step / z) < tol)
            {
                break;
            }
            if (iter == maxit)
            {
                conv = 0;
                break;
            }
        }
        *coef_out = z;
        for (int lib = 0; lib < n; ++lib)
        {
            mu_out[lib] = sigmod(z + o[lib]);
            dev += 2 * w[lib] * (y_log_y(y[lib], mu_out[lib]) + y_log_y(1 - y[lib], 1 - mu_out[lib]));
        }
    }
    *conv_out = conv;
    *dev_out = dev;
}

/* three nlib row buffers per thread for the one-group binomial fitter */
typedef struct {
    double *yptr, *optr, *wptr, *uptr;   /* nlib rows: proportions, offsets, weights, fitted mu */
} bone_ws;

void bin_one_group_mat (cmx *y, cmx *offsets, cmx *weights, int maxit, double tol, double *prob, double *coef, double *dev, int *conv, int nthreads)
{
    int ntag = (y->nrow), nlib = (y->ncol);

    int nth = clamp_threads(nthreads);

    /* row vectors for y offsets weights, one set per thread */
    bone_ws *ws = R_Calloc(nth, bone_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].yptr = R_Calloc(nlib,double);
        ws[t].optr = R_Calloc(nlib,double);
        ws[t].wptr = R_Calloc(nlib,double);
        ws[t].uptr = R_Calloc(nlib,double);
    }

    // Check whether offset is zero.
    int offset_is_zero=1, tag_start=0;
    double zero=0;
    offset_is_zero = check_row_scalar(offsets,tag_start,zero) && (offsets->type >= 2);

    // Iterating through tags and fitting.
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        bone_ws *w = &ws[tid];
        double *yptr = w->yptr, *optr = w->optr, *wptr = w->wptr;
        get_row(y,tag,yptr);
        get_row(offsets,tag,optr);
        get_row(weights,tag,wptr);

        // fit the whole row as a single group (closed-form or Newton-with-offset)
        bin_group_fit(yptr, optr, wptr, nlib, offset_is_zero, maxit, tol, w->uptr, &coef[tag], &conv[tag], &dev[tag]);

        double *pptr = prob+tag;
        for(int lib=0;lib<nlib;++lib,pptr+=ntag)
        {
            (*pptr) = w->uptr[lib];
        }
        dev[tag]=fmax2(dev[tag],0.0);
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].uptr);
    }
    R_Free(ws);

    return;
}

/* this function fits binomial regression using IWLS method
 *
 * inputs:
 * y        matrix of proportions
 * offsets  compressed matrix of offsets
 * weights  compressed matrix of total coverage
 * design   design matrix
 * maxit    max iteration
 * tol      tolerance
 * 
 * outputs:
 * mu      fitted values
 * beta    coefficients
 * dev     deviance 
 * iter    number of iteration
 * fail    index of convergence if failed
 *
 * No return value; parallelises over tags, giving each thread its own
 * biwls_ws scratch, and writes only into the caller-owned output arrays
 * (mu, beta, dev, iter, failed).
 */

/* per-thread scratch for the IWLS binomial fitter (mirrors the arrays
 * passed to bin_iwls_vec) */
typedef struct {
    double *yptr, *optr, *vptr;   /* nlib row buffers              */
    double *bptr;                 /* ncoef coefficients            */
    double *pptr;                 /* nlib fitted values            */
    double *wptr, *lptr;          /* nlib working weights / z       */
    double *xdm;                  /* nlib*ncoef design copy        */
    double *resid, *effects;      /* nlib                          */
    double *qraux;                /* ncoef                         */
    double *work;                 /* 2*ncoef                       */
    int    *pivot;                /* ncoef                         */
} biwls_ws;

/* file-local: fits one gene's binomial GLM by IWLS; defined below, used by bin_iwls_mat */
static void bin_iwls_vec(int, int, int, double *, double *, double *, int, int, double *, int, double, double *, double *, double *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *);

void bin_iwls_mat (cmx *y, cmx *offsets, cmx *weights, cmx *design, int maxit, double tol, double *mu, double *beta, double *dev, int *iter, int *failed, int nthreads)
{
    int ntag = (y->nrow), nlib = (y->ncol), ncoef = (design->ncol);
    int nsize = nlib * ncoef;
    double *dm;
    dm = (design->dmat);

    int nth = clamp_threads(nthreads);

    /* one full scratch workspace per thread (pivot is reset to 1..ncoef as
     * in the serial code; dqrdc2 reinitialises it internally so there is no
     * dependency across tags) */
    biwls_ws *ws = R_Calloc(nth, biwls_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].yptr    = R_Calloc(nlib,double);
        ws[t].optr    = R_Calloc(nlib,double);
        ws[t].vptr    = R_Calloc(nlib,double);
        ws[t].bptr    = R_Calloc(ncoef,double);
        ws[t].pptr    = R_Calloc(nlib,double);
        ws[t].wptr    = R_Calloc(nlib,double);
        ws[t].lptr    = R_Calloc(nlib,double);
        ws[t].xdm     = R_Calloc(nsize,double);
        ws[t].resid   = R_Calloc(nlib,double);
        ws[t].effects = R_Calloc(nlib,double);
        ws[t].qraux   = R_Calloc(ncoef,double);
        ws[t].work    = R_Calloc(ncoef*2,double);
        ws[t].pivot   = R_Calloc(ncoef,int);
        for(int i=0; i<ncoef; ++i)
        {
            ws[t].pivot[i]=i+1;
        }
    }

    // prepared for dqrls
    int ny=1, rank=1;

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        biwls_ws *w = &ws[tid];
        int oiter, ofail;
        double odev;

        get_row(y,tag,w->yptr);
        get_row(offsets,tag,w->optr);
        get_row(weights,tag,w->vptr);

        bin_iwls_vec(nlib,ny,rank,w->yptr,w->optr,w->vptr,ncoef,nsize,dm,maxit,tol,w->bptr,w->pptr,&odev,&oiter,&ofail,w->wptr,w->lptr,w->xdm,w->resid,w->effects,w->qraux,w->work,w->pivot);

        double *uupt=mu+tag;
        for(int lib=0;lib<nlib;++lib,uupt+=ntag)
        {
            (*uupt) = w->pptr[lib];
        }
        double *bbpt=beta+tag;
        for(int coef=0;coef<ncoef;++coef,bbpt+=ntag)
        {
            (*bbpt) = w->bptr[coef];
        }

        dev[tag]    = fmax2(odev,0.0);
        iter[tag]   = oiter;
        failed[tag] = ofail;
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].vptr);
        R_Free(ws[t].bptr);
        R_Free(ws[t].pptr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].lptr);
        R_Free(ws[t].xdm);
        R_Free(ws[t].resid);
        R_Free(ws[t].effects);
        R_Free(ws[t].qraux);
        R_Free(ws[t].work);
        R_Free(ws[t].pivot);
    }
    R_Free(ws);

    return;
}

/* this function fits binomial regression using IWLS method for one row
 *
 * inputs:
 * nlib     number of libraries
 * ny       = 1
 * rank     = 1
 * y        raw prpbabilities
 * offsets  offsets
 * weights  weights
 * ncoef    number of covariates
 * nsize    = nlib * ncoef
 * dm       design matrix
 * maxit    maximal iteration
 * tol      tolerance
 *  
 * outputs:
 * bptr     updated coefficents
 * pptr     fitted values
 * odevi    deviance 
 * oiter    number of iteration
 * ofail    index of convergence if failed
 * 
 * others:
 * wptr     working weights
 * lptr     eta - offsets
 * 
 * used for dqrls
 * xdm      copy of design matrix
 * resid    residuals
 * effects  effects
 * qraux    qraux
 * work     work space
 * pivot    pivot
 */
static void bin_iwls_vec(int nlib, int ny, int rank, double *y, double *offsets, double *weights, int ncoef, int nsize, double *dm, int maxit, double tol, 
                  double *bptr, double *pptr, double *odevi, int *oiter, int *ofail,
                  double *wptr, double *lptr, double *xdm, double *resid, double *effects, double *qraux, double *work, int *pivot)
{
    // starting values, empirical probabilities (pptr), eta-offsets (lptr)
    for(int lib=0;lib<nlib;++lib)
    {
        pptr[lib] = (y[lib] * weights[lib] + 0.5) / (weights[lib] + 1);
        lptr[lib] = logit(pptr[lib]) - offsets[lib];
    }

    // make a copy for design matrix
    for(int i=0; i<nsize; ++i)
    {
        xdm[i]=dm[i];
    }

    // fit X %*% beta = eta - offsets
    /* dqrls (LINPACK) performs a weighted least-squares solve by Householder
     * QR factorisation: it decomposes the design copy xdm and solves the linear
     * system X * b = y for the coefficients b that minimise the residual sum of
     * squares.
     * edgeR calls it here for the initial IWLS step of the binomial GLM, fitting
     * the working response (eta - offsets) against the design so the subsequent
     * Newton / IWLS updates start from an ordinary least-squares fit.
     * Argument mapping (Fortran name = local):
     *   x      = xdm      design copy (n x p), overwritten by its QR factors
     *   n/ldx  = nlib     number of observations (rows and leading dimension)
     *   p      = ncoef    number of coefficients (columns)
     *   y      = lptr     working response (eta - offsets), n x ny
     *   ny     = ny       number of response columns (= 1)
     *   tol    = tol      tolerance for rank determination
     *   b      = bptr     output coefficients (p x ny)
     *   rsd    = resid    output residuals (n x ny)
     *   qty    = effects  output Q' * y (n x ny)
     *   k      = rank     output numerical rank
     *   jpvt   = pivot    column pivots (length p)
     *   qraux  = qraux    auxiliary QR information (length p)
     *   work   = work     workspace (length 2 * p)
     * Result: bptr holds the least-squares coefficients and resid the residuals,
     * used just below to form eta and the starting deviance.
     * Equivalent R operation: lm.fit(x = design, y = eta - offsets)$coefficients
     * (equivalently qr.solve(design, eta - offsets)).
     * Netlib references: dqrls is a LINPACK-based routine shipped with R in
     * src/appl/dqrls.f (calls LINPACK dqrdc2 / dqrsl); see
     * https://netlib.org/linpack/ . */
    F77_CALL(dqrls)(xdm, &nlib, &ncoef, lptr, &ny, &tol, bptr, resid, effects, &rank, pivot, qraux, work);

    int iter=1, failed=0;
    double dev=0, varp=0, lmtol=1e-12;

    // initialize pptr, lptr, deviance
    for(int lib=0;lib<nlib;++lib)
    {
        lptr[lib]=lptr[lib]-resid[lib];
        pptr[lib]=sigmod(lptr[lib]+offsets[lib]);
        dev += 2 * weights[lib] * (y_log_y(y[lib],pptr[lib]) + y_log_y(1-y[lib],1-pptr[lib]));
    }

    // IWLS method
    while(++iter)
    {
        // compute w and sqrt(w) * z
        for(int lib=0;lib<nlib;++lib)
        {
            varp= ((fabs(pptr[lib]) < DOUEPS) || (fabs(1-pptr[lib]) < DOUEPS)) ? DOUEPS : pptr[lib]*(1-pptr[lib]);
            wptr[lib]=sqrt(weights[lib]*varp);
            lptr[lib]=(fabs(weights[lib]) < DOUEPS) ? DOUEPS : wptr[lib] * ((y[lib]-pptr[lib])/varp + lptr[lib]);
        }

        // compute sqrt(w) * design
        for(int i=0; i<nsize; ++i)
        {
            xdm[i]= dm[i] * wptr[i % nlib];
        }

        // solve sqrt(w) * design * beta = sqrt(w) * z, and update coefficients (bptr)
        /* dqrls (LINPACK) again solves a weighted least-squares problem by
         * Householder QR: here xdm already holds sqrt(w) * design and lptr holds
         * sqrt(w) * z, so the solve delivers the reweighted least-squares update
         * of the coefficients.
         * edgeR calls it on every IWLS iteration to update the binomial GLM
         * coefficients from the current working weights and working response.
         * Argument mapping (Fortran name = local):
         *   x      = xdm      sqrt(w)-scaled design (n x p), overwritten by QR factors
         *   n/ldx  = nlib     number of observations (rows and leading dimension)
         *   p      = ncoef    number of coefficients (columns)
         *   y      = lptr     sqrt(w)-scaled working response, n x ny
         *   ny     = ny       number of response columns (= 1)
         *   tol    = lmtol    tolerance for rank determination (1e-12)
         *   b      = bptr     output updated coefficients (p x ny)
         *   rsd    = resid    output residuals (n x ny)
         *   qty    = effects  output Q' * y (n x ny)
         *   k      = rank     output numerical rank
         *   jpvt   = pivot    column pivots (length p)
         *   qraux  = qraux    auxiliary QR information (length p)
         *   work   = work     workspace (length 2 * p)
         * Result: bptr holds the updated coefficients; resid is used below to
         * recover eta = design * beta without an explicit matrix product.
         * Equivalent R operation: lm.fit(x = sqrt(w) * design, y = sqrt(w) * z)
         * (one weighted step of glm.fit / IRLS).
         * Netlib references: dqrls is a LINPACK-based routine shipped with R in
         * src/appl/dqrls.f (calls LINPACK dqrdc2 / dqrsl); see
         * https://netlib.org/linpack/ . */
        F77_CALL(dqrls)(xdm, &nlib, &ncoef, lptr, &ny, &lmtol, bptr, resid, effects, &rank, pivot, qraux, work);

        // update pptr, lptr, deviance
        double ndev=0;
        for(int lib=0;lib<nlib;++lib)
        {
            // trick to compute design * beta, avoiding matrix product
            lptr[lib]=(fabs(weights[lib]) < DOUEPS) ? DOUEPS : (lptr[lib]-resid[lib])/wptr[lib];
            pptr[lib]=sigmod(lptr[lib]+offsets[lib]);
            ndev += 2 * weights[lib] * (y_log_y(y[lib],pptr[lib]) + y_log_y(1-y[lib],1-pptr[lib]));
        }

        if(fabs(dev-ndev)/(fabs(dev)+1) < tol)
        {
            break;
        }
        dev=ndev;

        if(iter > maxit)
        {
            failed=1;
            break;
        }
    }

    // save results
    (*oiter) = iter;
    (*odevi) = dev;
    (*ofail) = failed;

    return;
}

/* ----------------------------------------------------------------------------
 * Unified single-call binomial fitter (mirrors fit_glm / fit_glm_mat in glm.c).
 * Detects a oneway layout from the design and either fits groupwise and solves
 * back with LAPACK dgesv, or falls back to the general IWLS solver.
 * ------------------------------------------------------------------------- */

/* per-thread scratch for the oneway binomial fitter */
typedef struct
{
    double *ysptr;   /* nlib subset of the proportions        */
    double *osptr;   /* nlib subset of the offsets            */
    double *wsptr;   /* nlib subset of the coverage weights   */
    double *usptr;   /* nlib subset fitted mu                 */
    double *yptr;    /* nlib full row of the proportions      */
    double *optr;    /* nlib full row of the offsets          */
    double *wptr;    /* nlib full row of the coverage weights */
} bonew_ws;

/* Fit a oneway-layout binomial GLM for every gene, parallelised over tags.
 * Reuses the per-gene math of bin_one_group_mat for each group of libraries.
 *
 * inputs:
 *   y        proportions (dense matrix)
 *   offsets  compressed matrix of offsets
 *   weights  coverage-scaled weights (dense matrix)
 *   group    library-to-group index, length nlib
 *   ngroups  number of groups
 *   maxit    maximum Newton iterations
 *   tol      convergence tolerance
 *
 * outputs (caller-allocated):
 *   mu    fitted probabilities (ntag x nlib)
 *   beta  group logits (ntag x ngroups)
 *   dev   residual deviance per gene, summed over groups (ntag)
 *   conv  convergence indicator per gene and group (ntag x ngroups)
 *
 * No return value; writes only into the caller-owned output arrays.
 */
static void bin_one_way_mat (cmx *y, cmx *offsets, cmx *weights, int *group, int ngroups, int maxit, double tol, double *mu, double *beta, double *dev, int *conv, int nthreads)
{
    int ntag = (y->nrow);
    int nlibs = (y->ncol);
    int nth = clamp_threads(nthreads);

    /* one set of subset and full-row buffers per thread */
    bonew_ws *ws = R_Calloc(nth, bonew_ws);
    for (int t = 0; t < nth; ++t)
    {
        ws[t].ysptr = R_Calloc(nlibs, double);
        ws[t].osptr = R_Calloc(nlibs, double);
        ws[t].wsptr = R_Calloc(nlibs, double);
        ws[t].usptr = R_Calloc(nlibs, double);
        ws[t].yptr = R_Calloc(nlibs, double);
        ws[t].optr = R_Calloc(nlibs, double);
        ws[t].wptr = R_Calloc(nlibs, double);
    }

    /* group partition (computed once, read-only inside the parallel loop) */
    int *gnlib = R_Calloc(ngroups, int);
    int *g_idx = R_Calloc(nlibs, int);
    int *ogptr = R_Calloc(ngroups + 1, int);
    for (int lib = 0; lib < nlibs; ++lib)
    {
        int g = group[lib];
        if (g >= 0 && g < ngroups)
        {
            gnlib[g]++;
        }
    }
    ogptr[0] = 0;
    for (int g = 0; g < ngroups; ++g)
    {
        ogptr[g + 1] = ogptr[g] + gnlib[g];
    }
    int *cpos = R_Calloc(ngroups, int);
    for (int g = 0; g < ngroups; ++g)
    {
        cpos[g] = ogptr[g];
    }
    for (int lib = 0; lib < nlibs; ++lib)
    {
        int g = group[lib];
        if (g >= 0 && g < ngroups)
        {
            g_idx[cpos[g]++] = lib;
        }
    }
    R_Free(cpos);

    /* per-group offset-is-zero flag: row 0 is representative only when offsets
     * is row-repeated (type >= 2), matching bin_one_group_mat's guard */
    int *offzero = R_Calloc(ngroups, int);
    double *o0 = R_Calloc(nlibs, double);
    get_row(offsets, 0, o0);
    for (int g = 0; g < ngroups; ++g)
    {
        int oz = 1;
        for (int k = 0; k < gnlib[g]; ++k)
        {
            if (o0[g_idx[ogptr[g] + k]] != 0.0)
            {
                oz = 0;
                break;
            }
        }
        offzero[g] = oz && (offsets->type >= 2);
    }
    R_Free(o0);

    /* column base pointers for mu (ntag x nlibs): the hot loop writes
       mucol[lib][tag] instead of mu[lib*ntag+tag].  The column bases are built by
       pointer increment (mcol += ntag), so no lib*ntag product is formed at all. */
    double **mucol = R_Calloc(nlibs, double*);
    double *mcol = mu;
    for (int l = 0; l < nlibs; ++l, mcol += ntag)
    {
        mucol[l] = mcol;
    }

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag = 0; tag < ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        bonew_ws *w = &ws[tid];
        get_row(y, tag, w->yptr);
        get_row(offsets, tag, w->optr);
        get_row(weights, tag, w->wptr);

        dev[tag] = 0;

        /* conv[b_idx]/beta[b_idx] with b_idx = g*ntag+tag walked by ntag over
           groups (cvp/btp); mu[lib*ntag+tag] -> mucol[lib][tag] (table above). */
        int *cvp = conv + tag;
        double *btp = beta + tag;
        for (int g = 0; g < ngroups; ++g, cvp += ntag, btp += ntag)
        {
            int nlibs_g = gnlib[g];
            int s_idx = ogptr[g];
            for (int k = 0; k < nlibs_g; ++k)
            {
                int lib = g_idx[s_idx + k];
                w->ysptr[k] = w->yptr[lib];
                w->osptr[k] = w->optr[lib];
                w->wsptr[k] = w->wptr[lib];
            }

            double gdev = 0;
            bin_group_fit(w->ysptr, w->osptr, w->wsptr, nlibs_g, offzero[g], maxit, tol, w->usptr, btp, cvp, &gdev);
            for (int k = 0; k < nlibs_g; ++k)
            {
                mucol[g_idx[s_idx + k]][tag] = w->usptr[k];
            }
            dev[tag] += fmax2(gdev, 0.0);
        }
    }

    R_Free(mucol);
    R_Free(offzero);
    R_Free(gnlib);
    R_Free(g_idx);
    R_Free(ogptr);
    for (int t = 0; t < nth; ++t)
    {
        R_Free(ws[t].ysptr);
        R_Free(ws[t].osptr);
        R_Free(ws[t].wsptr);
        R_Free(ws[t].usptr);
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
    }
    R_Free(ws);

    return;
}

/* Unified single-call binomial fitter, parallel to fit_glm_mat in glm.c.
 * Detects a oneway layout from the design.  For a oneway layout it fits the
 * groups with bin_one_way_mat and maps the group logits back to design
 * coefficients with LAPACK dgesv; otherwise it runs the general IWLS solver
 * bin_iwls_mat.  Two (maxit, tol) pairs preserve the distinct convergence
 * settings of the former mbinOneWay (50, 1e-8) and mbinIWLS (100, 1e-6).
 *
 * inputs:
 *   y        proportions (dense matrix)
 *   offsets  compressed matrix of offsets
 *   weights  coverage-scaled weights (dense matrix)
 *   design   design matrix (dense)
 *   maxit_ow, tol_ow      oneway Newton settings
 *   maxit_iwls, tol_iwls  IWLS settings
 *
 * outputs (caller-allocated):
 *   coef    coefficients on the logit scale (ntag x ncoef)
 *   mu      fitted probabilities (ntag x nlib)
 *   dev     residual deviance per gene (ntag)
 *   iter    IWLS iterations per gene, 0 for the oneway path (ntag)
 *   failed  non-convergence indicator per gene (ntag)
 *   method  set to 0 for the oneway path, 1 for IWLS
 *
 * No return value; writes only into the caller-owned output arrays.
 */
void bin_fit_mat (cmx *y, cmx *offsets, cmx *weights, cmx *design, int maxit_ow, double tol_ow, int maxit_iwls, double tol_iwls, double *coef, double *mu, double *dev, int *iter, int *failed, int *method, int nthreads)
{
    int ntag = (y->nrow);
    int nlibs = (y->ncol);
    int ncoef = (design->ncol);

    int *group = R_Calloc(nlibs, int);
    int ngroups = get_groups_from_design(design, group);

    if (ngroups == ncoef)
    {
        /* oneway layout: fit groups, then solve back to design coefficients.
         * designunique holds the design row of the first library of each group */
        double *designunique = R_Calloc(ncoef * ncoef, double);
        int *first_lib_of_group = R_Calloc(ncoef, int);
        for (int g = 0; g < ncoef; ++g)
        {
            first_lib_of_group[g] = -1;
        }
        for (int lib = 0; lib < nlibs; ++lib)
        {
            int g = group[lib];
            if (g >= 0 && g < ncoef && first_lib_of_group[g] == -1)
            {
                first_lib_of_group[g] = lib;
            }
        }
        for (int c = 0; c < ncoef; ++c)
        {
            for (int g = 0; g < ncoef; ++g)
            {
                designunique[c * ncoef + g] = design->dmat[c * nlibs + first_lib_of_group[g]];
            }
        }
        R_Free(first_lib_of_group);

        /* group-level binomial fit: mu, coef (group logits), dev, conv */
        int *conv = R_Calloc((size_t) ntag * ngroups, int);
        bin_one_way_mat(y, offsets, weights, group, ngroups, maxit_ow, tol_ow, mu, coef, dev, conv, nthreads);
        /* conv[g*ntag+tag] walked by ntag over groups */
        /*
        for (int tag = 0; tag < ntag; ++tag)
        {
            failed[tag] = 0;
            for (int g = 0; g < ngroups; ++g)
            {
                if (!conv[g * ntag + tag])
                {
                    failed[tag] = 1;
                    break;
                }
            }
            iter[tag] = 0;
        }
        */
        for (int tag = 0; tag < ntag; ++tag)
        {
            failed[tag] = 0;
            const int *cvp = conv + tag;
            for (int g = 0; g < ngroups; ++g, cvp += ntag)
            {
                if (!*cvp)
                {
                    failed[tag] = 1;
                    break;
                }
            }
            iter[tag] = 0;
        }
        R_Free(conv);

        /* copy the group logits into the right-hand side, one gene per column */
        double *rhs = R_Calloc((size_t) ncoef * ntag, double);
        /* rhs[tag*ncoef+g] = coef[g*ntag+tag]: rhs base (size_t)tag*ncoef then +1;
           coef walked by ntag over coefs. */
        /*
        for (int tag = 0; tag < ntag; ++tag)
        {
            for (int g = 0; g < ncoef; ++g)
            {
                rhs[tag * ncoef + g] = coef[g * ntag + tag];
            }
        }
        */
        for (int tag = 0; tag < ntag; ++tag)
        {
            double *rp = rhs + (size_t) tag * ncoef;
            const double *cp = coef + tag;
            for (int g = 0; g < ncoef; ++g, ++rp, cp += ntag)
            {
                *rp = *cp;
            }
        }

        int info = 0;
        int *ipiv = R_Calloc(ncoef, int);

        /* dgesv is the LAPACK double general linear solver.  It LU-factorises the
         * square matrix designunique and overwrites rhs with the solution X of
         * designunique %*% X = rhs.  edgeR calls it here to map each gene's group
         * logits back to design-scale coefficients for a oneway layout.
         * Argument mapping: N = ncoef (system order), NRHS = ntag (one gene per
         * column of rhs), A = designunique (N x N, overwritten by its LU factors),
         * LDA = ncoef, IPIV = ipiv (row pivots), B = rhs (N x NRHS, overwritten by
         * the solution), LDB = ncoef, INFO = info (0 on success, non-zero singular).
         * The solution gives the genewise coefficients on the natural-log-odds scale.
         * Equivalent R operation: beta <- t(solve(designunique, t(beta))).
         * Netlib references: https://netlib.org/lapack/explore-html/ (dgesv). */
        F77_CALL(dgesv)(&ncoef, &ntag, designunique, &ncoef, ipiv, rhs, &ncoef, &info);
        if (info != 0)
        {
            R_Free(ipiv);
            R_Free(designunique);
            R_Free(rhs);
            R_Free(group);
            error("LAPACK dgesv failed to solve the design matrix transformation");
        }
        R_Free(ipiv);
        R_Free(designunique);
        /* coef[g*ntag+tag] = rhs[tag*ncoef+g] (inverse of the pack above) */
        /*
        for (int tag = 0; tag < ntag; ++tag)
        {
            for (int g = 0; g < ncoef; ++g)
            {
                coef[g * ntag + tag] = rhs[tag * ncoef + g];
            }
        }
        */
        for (int tag = 0; tag < ntag; ++tag)
        {
            double *cp = coef + tag;
            const double *rp = rhs + (size_t) tag * ncoef;
            for (int g = 0; g < ncoef; ++g, cp += ntag, ++rp)
            {
                *cp = *rp;
            }
        }
        R_Free(rhs);

        *method = 0;
    }
    else
    {
        /* general design: iteratively reweighted least squares */
        bin_iwls_mat(y, offsets, weights, design, maxit_iwls, tol_iwls, mu, coef, dev, iter, failed, nthreads);
        *method = 1;
    }

    R_Free(group);
    return;
}