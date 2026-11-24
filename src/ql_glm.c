#include "edgeR.h"

/* compute the adjusted deviance and degree of freedom for QL method */

/*
  inputs:
  y       raw count matrix
  mu      fitted value matrix
  design  design matrix

  disp    dispersion compressMatrix
  weights weight compressMatrix

  prior   average quasi-dispersion

  vector outputs:
  df  adjusted degree of freedom
  dev adjusted deviance or total.dev
  s2  quasi dispersion estimate by adjusted dev and df

  matrix outputs:
  lvmat hatvalue matrix
  dvmat unit deviance matrix
  dfmat individual df matrix
*/

static const double thresholdzero=1e-4;

/* Per-thread scratch for the QL adjustment kernels (NB and binomial). The `aux`
 * buffer holds the NB dispersion or the binomial coverage; the trailing four fields
 * are the scratch qr_hat needs, hoisted here so qr_hat performs no allocation inside
 * the parallel region (qraux also holds the group sums of the one-way hat values). */
typedef struct {
    double *yptr, *uptr, *wptr, *aux;    /* nlib row buffers: y, mu, weights, disp/cover */
    double *xdpt;                        /* nlib*ncoef weighted design */
    double *hptr;                        /* nlib hat values           */
    double *zwpt;                        /* nlib working weights       */
    int    *jpvt;                        /* ncoef  qr_hat pivots       */
    double *qraux;                       /* ncoef  qr_hat              */
    double *qrwork;                      /* 2*ncoef qr_hat work        */
    double *vrow;                        /* ncoef  qr_hat solve vector */
} ql_ws;

/* Allocate one ql_ws scratch workspace per thread; release with ql_ws_free. */
static ql_ws *ql_ws_alloc(int nth, int nlib, int ncoef)
{
    ql_ws *ws = R_Calloc(nth, ql_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].yptr   = R_Calloc(nlib,double);
        ws[t].uptr   = R_Calloc(nlib,double);
        ws[t].wptr   = R_Calloc(nlib,double);
        ws[t].aux    = R_Calloc(nlib,double);
        ws[t].xdpt   = R_Calloc(nlib*ncoef,double);
        ws[t].hptr   = R_Calloc(nlib,double);
        ws[t].zwpt   = R_Calloc(nlib,double);
        ws[t].jpvt   = R_Calloc(ncoef,int);
        ws[t].qraux  = R_Calloc(ncoef,double);
        ws[t].qrwork = R_Calloc(2*ncoef,double);
        ws[t].vrow   = R_Calloc(ncoef,double);
    }
    return ws;
}

/* Free a ql_ws array allocated by ql_ws_alloc. */
static void ql_ws_free(ql_ws *ws, int nth)
{
    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].uptr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].aux);
        R_Free(ws[t].xdpt);
        R_Free(ws[t].hptr);
        R_Free(ws[t].zwpt);
        R_Free(ws[t].jpvt);
        R_Free(ws[t].qraux);
        R_Free(ws[t].qrwork);
        R_Free(ws[t].vrow);
    }
    R_Free(ws);
}

/* One-way layout check for the QL hat values, done once per call before the gene loop.
 *
 * inputs:
 *   design  nlib x ncoef design matrix
 *
 * outputs:
 *   group   length nlib: index (0..ncoef-1) of the distinct design row of each library;
 *           meaningful only when the return value is 1
 *
 * returns:
 *   1 when the design has exactly ncoef distinct rows (equal within 1e-12, as in
 *   get_groups_from_design) and the ncoef x ncoef matrix of those rows has full rank: the
 *   column space of the design is then spanned by the group indicators, and the hat value of
 *   library i in a gene is w_i / (sum of w over the libraries of its group). 0 otherwise; the
 *   scan stops as soon as a row ncoef+1 is found, so a design with a continuous covariate costs
 *   O(nlib ncoef^2), not O(nlib^2).
 */
static int ql_oneway_groups(const cmx *design, int *group)
{
    int nlib = design->nrow;
    int ncoef = design->ncol;
    const double *dmat = design->dmat;
    int *rep = R_Calloc(ncoef, int);
    int ngroups = 0;
    int ok = 1;
    for(int i=0; i<nlib && ok; ++i)
    {
        int g;
        for(g=0; g<ngroups; ++g)
        {
            int match = 1;
            for(int c=0; c<ncoef; ++c)
            {
                if(fabs(dmat[c*nlib+i] - dmat[c*nlib+rep[g]]) > 1e-12)
                {
                    match = 0;
                    break;
                }
            }
            if(match)
            {
                break;
            }
        }
        if(g==ngroups)
        {
            if(ngroups==ncoef)
            {
                ok = 0;
                break;
            }
            rep[ngroups] = i;
            ngroups++;
        }
        group[i] = g;
    }
    if(ok && ngroups==ncoef)
    {
        double *u = R_Calloc(ncoef*ncoef, double);
        double *qraux = R_Calloc(ncoef, double);
        double *work = R_Calloc(2*ncoef, double);
        int *jpvt = R_Calloc(ncoef, int);
        double tol = 1e-7;
        int k;
        for(int g=0; g<ncoef; ++g)
        {
            jpvt[g] = g+1;
            for(int c=0; c<ncoef; ++c)
            {
                u[c*ncoef+g] = dmat[c*nlib+rep[g]];
            }
        }
        /* dqrdc2 is R's modified LINPACK QR routine (src/appl/dqrdc2.f): Householder QR of the
         * ncoef x ncoef matrix u of distinct design rows with tolerance-based limited column
         * pivoting, returning its numerical rank k.
         * edgeR calls it here only to confirm that the distinct rows have full rank, so that the
         * design and the group indicators span the same column space.
         * Argument mapping:
         *   x     = u      (ncoef x ncoef distinct rows; overwritten by the QR factors)
         *   ldx   = &ncoef (leading dimension of u)
         *   n     = &ncoef (number of rows)
         *   p     = &ncoef (number of columns)
         *   tol   = &tol   (rank-detection tolerance, 1e-7, as in qr_hat)
         *   k     = &k     (output: numerical rank)
         *   qraux = qraux  (output, length ncoef: Householder scalars, unused here)
         *   jpvt  = jpvt   (in/out, length ncoef: pivots, preset to 1..ncoef)
         *   work  = work   (workspace, length 2*ncoef)
         * The design is treated as a one-way layout only when k == ncoef.
         * Equivalent R operation: qr(u, tol = 1e-7)$rank == ncol(u).
         * Netlib references: https://netlib.org/linpack/dqrdc.f */
        F77_CALL(dqrdc2)(u, &ncoef, &ncoef, &ncoef, &tol, &k, qraux, jpvt, work);
        ok = (k==ncoef);
        R_Free(u);
        R_Free(qraux);
        R_Free(work);
        R_Free(jpvt);
    }
    else
    {
        ok = 0;
    }
    R_Free(rep);
    return ok;
}

/* Hat values of the sqrt(W)-weighted design for one gene, written to w->hptr. The working
 * weight zwpt is the only distribution-specific part: NB uses sqrt(mu*w/(1+mu*disp/
 * prior)), binomial uses sqrt(cover*mu*(1-mu)). When group is not NULL (one-way layout, see
 * ql_oneway_groups) the hat value is zw_i^2 / (sum of zw^2 over the group), 0 for a group whose
 * weights are all 0; otherwise sqrt(W)X is formed in w->xdpt and passed to qr_hat. */
static void ql_tag_hat(ql_ws *w, const double *dmat, int nlib, int ncoef, int fam, double prior, const int *group)
{
    for(int lib=0; lib<nlib; ++lib)
    {
        if(fam==QL_NB)
        {
            w->zwpt[lib]=sqrt(w->uptr[lib]*w->wptr[lib]/(1+(w->uptr[lib]*w->aux[lib]/prior)));
        }
        else
        {
            w->zwpt[lib]=sqrt(w->aux[lib]*w->uptr[lib]*(1-w->uptr[lib]));
        }
        w->hptr[lib]=0;
    }
    if(group)
    {
        double *gsum = w->qraux;
        for(int g=0; g<ncoef; ++g)
        {
            gsum[g] = 0;
        }
        for(int lib=0; lib<nlib; ++lib)
        {
            gsum[group[lib]] += w->zwpt[lib] * w->zwpt[lib];
        }
        for(int lib=0; lib<nlib; ++lib)
        {
            double sg = gsum[group[lib]];
            w->hptr[lib] = (sg > 0) ? w->zwpt[lib] * w->zwpt[lib] / sg : 0.0;
        }
        return;
    }
    for(int i=0; i<nlib*ncoef; ++i)
    {
        w->xdpt[i]=dmat[i] * w->zwpt[i % nlib];
    }
    qr_hat(w->xdpt, nlib, ncoef, dmat, w->zwpt, w->hptr, w->jpvt, w->qraux, w->qrwork, w->vrow);
}

/* Per-library leverage-adjusted unit contributions for one gene, distribution-specific.
 * Writes the QL weight pair wpt[2], the unit deviance *udp, the leverage complement *hdp
 * (1-hat, floored to 0 below thresholdzero), and the binomial coverage mask *wdp (1.0 for
 * NB, so multiplying by it leaves the NB path unchanged). */
static void ql_unit(const ql_ws *w, int lib, int fam, double prior, double *udp, double *wpt, double *hdp, double *wdp)
{
    if(fam==QL_NB)
    {
        compute_weight_negbin(w->uptr[lib], w->aux[lib], prior/w->wptr[lib], wpt);
        *udp = compute_unit_nb_deviance(w->yptr[lib], w->uptr[lib], w->aux[lib]*w->wptr[lib]/prior);
        *wdp = 1.0;
    }
    else
    {
        compute_weight_binomial(w->uptr[lib], (int) (w->aux[lib]), wpt);
        *udp = 2 * (y_log_y(w->yptr[lib],w->uptr[lib]) + y_log_y(1-w->yptr[lib],1-w->uptr[lib]));
        *wdp = (w->aux[lib] > 0.5) ? 1.0 : 0.0;
    }
    *hdp = 1.0 - w->hptr[lib];
    if(*hdp < thresholdzero)
    {
        *udp=0.0;
        *hdp=0.0;
    }
}

/* Shared QL adjusted total-deviance / df / quasi-dispersion driver (NB or binomial).
 * aux holds disp (NB) or cover (binomial); prior is used only for QL_NB. */
void ql_adjust_vec (int fam, cmx *y, cmx *mu, cmx *design, cmx *aux, double prior, cmx *weights, double *df, double *dev, double *s2, int nthreads)
{
    int ntag=(y->nrow), nlib=(y->ncol), ncoef=(design->ncol);

    int nth = clamp_threads(nthreads);
    ql_ws *ws = ql_ws_alloc(nth, nlib, ncoef);
    int *group = R_Calloc(nlib, int);
    if(!ql_oneway_groups(design, group))
    {
        R_Free(group);
        group = NULL;
    }

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        ql_ws *w = &ws[tid];
        double wpt[2], hdp, udp, wdp;
        get_row4(y,mu,aux,weights,tag,w->yptr,w->uptr,w->aux,w->wptr);
        ql_tag_hat(w, design->dmat, nlib, ncoef, fam, prior, group);

        dev[tag]=0;
        df[tag]=0;
        for(int lib=0; lib<nlib; ++lib)
        {
            ql_unit(w, lib, fam, prior, &udp, wpt, &hdp, &wdp);
            dev[tag] += (udp*wpt[0])*w->wptr[lib];
            df[tag] += hdp*wpt[1]*wdp;
        }
        if(fam==QL_BIN)
        {
            dev[tag] = fmax2(dev[tag], 0.0);
        }
        s2[tag] = (df[tag] < thresholdzero) ? 0.0 : dev[tag]/df[tag];
    }

    ql_ws_free(ws, nth);
    if(group)
    {
        R_Free(group);
    }

    return;
}

/* Shared QL per-observation unit-matrix driver (NB or binomial): unit deviance, unit df,
 * and leverage for every (gene, library). */
void ql_adjust_mat (int fam, cmx *y, cmx *mu, cmx *design, cmx *aux, double prior, cmx *weights, double *dfmat, double *dvmat, double *lvmat, int nthreads)
{
    int ntag=(y->nrow), nlib=(y->ncol), ncoef=(design->ncol);

    int nth = clamp_threads(nthreads);
    ql_ws *ws = ql_ws_alloc(nth, nlib, ncoef);
    int *group = R_Calloc(nlib, int);
    if(!ql_oneway_groups(design, group))
    {
        R_Free(group);
        group = NULL;
    }

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        ql_ws *w = &ws[tid];
        double wpt[2], hdp, udp, wdp;
        get_row4(y,mu,aux,weights,tag,w->yptr,w->uptr,w->aux,w->wptr);
        ql_tag_hat(w, design->dmat, nlib, ncoef, fam, prior, group);

        double *fptr=dfmat+tag;
        double *vptr=dvmat+tag;
        double *lptr=lvmat+tag;
        for(int lib=0; lib<nlib; ++lib, fptr+=ntag, vptr+=ntag, lptr+=ntag)
        {
            ql_unit(w, lib, fam, prior, &udp, wpt, &hdp, &wdp);
            (*vptr) = udp*wpt[0];         // dev matrix
            (*fptr) = hdp*wpt[1]*wdp;     // df  matrix
            (*lptr) = w->hptr[lib]*wdp;   // hat matrix
        }
    }

    ql_ws_free(ws, nth);
    if(group)
    {
        R_Free(group);
    }

    return;
}

/* return s2 adjusted deviance and df vectors */
void compute_adjust_vec (cmx *y, cmx *mu, cmx *design, cmx *disp, double prior, cmx *weights, double *df, double *dev, double *s2, int nthreads)
{
    ql_adjust_vec(QL_NB, y, mu, design, disp, prior, weights, df, dev, s2, nthreads);
}

/* return the unit matrices: unit deviance, unit df, and leverages */
void compute_adjust_mat (cmx *y, cmx *mu, cmx *design, cmx *disp, double prior, cmx *weights, double *dfmat, double *dvmat, double *lvmat, int nthreads)
{
    ql_adjust_mat(QL_NB, y, mu, design, disp, prior, weights, dfmat, dvmat, lvmat, nthreads);
}

/* compute the average quasi dispersion by iteration using two updates */
/* file-local: 90% quantile of the lowess quasi-dispersion trend; defined below */
static double compute_prior (double *, double *, double *, int);

double update_prior(cmx *y, cmx *mu, cmx *design, cmx *disp, cmx *weights, double *avg, int nthreads)
{
    int ntag=(y->nrow);
    double *df  = R_Calloc(ntag, double);
    double *dev = R_Calloc(ntag, double);
    double *s2  = R_Calloc(ntag, double);

    // initialize prior = 1
    double prior=1.0;

    // first update
    compute_adjust_vec(y,mu,design,disp,prior,weights,df,dev,s2,nthreads);
    prior = compute_prior(avg,s2,df, ntag);
    // second update
    compute_adjust_vec(y,mu,design,disp,prior,weights,df,dev,s2,nthreads);
    prior = compute_prior(avg,s2,df, ntag);

    R_Free(df);
    R_Free(dev);
    R_Free(s2);

    return prior;
}


/* The following code is used to compute average quasi-dispersion in edgeR
 * a c-version for .computePrior()
 * Created by Lizhong Chen
 * Last revised on 3 May 2024
 */

/*
 * a wrapper for the input data, preparing for clowess
 * sort x and y using R_qsort_I
 * calculate delta = diff(range(x))
 */

static void clowess2 (double *x, double *y, int n, int iter, double f, double *ans)
{
    int *ind = R_Calloc(n, int);
    for(int i=0;i<n;++i)
    {
        ind[i]=i;
    }
    R_qsort_I(x,ind,1,n);

    double delta=0.01*(x[n-1]-x[0]);
    double *yy  = R_Calloc(n, double);
    double *rw  = R_Calloc(n, double);
    double *res = R_Calloc(n, double);
    for(int i=0;i<n;++i)
    {
        yy[i]=y[ind[i]];
    }

    /* lowess: Cleveland LOWESS smoother (lowess.f, the Fortran ancestor of R stats'
     * lowess.c; this copy carries R's double-precision +1e-7 neighbourhood and the
     * cmad<1e-7*sc guard). Fits a robust locally-weighted trend ys through sorted (x, yy).
     * Called locally to build the quasi-dispersion prior trend for glmQLFit.
     * Args -> Fortran: x,yy=sorted data (in); &n=n; &f=span; &iter=nsteps;
     *   &delta=interpolation threshold; ans=ys smoothed fit (out); rw,res=workspace.
     * Equivalent R operation: stats::lowess(x, y, f, iter).
     * Netlib references: Cleveland (1979) JASA 74:829-836; netlib lowess.f.
     */
    F77_CALL(lowess)(x, yy, &n, &f, &iter, &delta, ans, rw, res);

    R_Free(ind);
    R_Free(yy);
    R_Free(rw);
    R_Free(res);
    return;

}

// compute average quasi-dispersion using lowess trend
static double compute_prior (double *ag, double *s2, double *df, int ntag)
{
    // default settings: f=0.5, iter=3L, t=1e-8
    double t=1e-8, f=0.5, out;
    int k=0, iter=3;

    double *xx  = R_Calloc(ntag, double);
    double *yy  = R_Calloc(ntag, double);

    // remove obs with small df and fit on x^1/4 level
    for(int i=0;i<ntag;++i)
    {
        if(df[i]>t)
        {
            xx[k]=ag[i];
            yy[k]=sqrt(sqrt(s2[i]));
            ++k;
        }
    }
    double p;
    if(k > 1)
    {
        double *ans = R_Calloc(k, double);

        // fit lowess trend
        clowess2(xx,yy,k,iter,f,ans);

        // compute 90% quantile of the trend, converting from quantile.R with type = 7
        double m = (k-1) * 0.9;
        int lo   = (int)(m);
        rPsort(ans,k,lo);
        rPsort(ans,k,lo+1);
        double h = m-lo;
        p = (1-h)*ans[lo] + h*ans[lo+1];

        R_Free(ans);
    }
    else
    {
        // k==1: one informative gene, no trend to fit; the 90% quantile is its own
        // value (h==0 in the k==1 quantile), avoiding the ans[lo+1] OOB read.
        p = yy[0];
    }

    // low bound is 1
    if(p < 1.0)
    {
        p=1.0;
    }

    out=p*p*p*p;

    R_Free(xx);
    R_Free(yy);

    return out;
}

/* Hat values of an n x p matrix x (here the sqrt(W)-weighted design of one gene), as in hat() of
 * the stats package:
 *   rowSums(qr.qy(qr(x), diag(1, nrow = n, ncol = qr(x)$rank))^2)
 * The pivoted QR gives X_k = Q_k R_k for the leading k (rank) pivoted columns, so row i of Q_k is
 * z = R_k^{-T} x_i and hat_i = sum(z^2). z is found by one forward solve per row (k^2/2
 * operations) instead of forming Q_k with dqrqy (about 2nk^2), with the same pivoting and rank
 * decision as before.
 *
 * inputs:
 *   x      n x p weighted design, overwritten by dqrdc2
 *   n, p   dimensions of x (n >= p)
 *   dmat   n x p unweighted design and
 *   zw     length-n weights: row i of x is zw[i] * dmat[i, ]; used to rebuild the rows of x
 *          after dqrdc2 has overwritten it
 *
 * outputs:
 *   hat    length n: hat values (each row's sum of squares of Q_k)
 *
 * scratch (caller-owned): jpvt (p), qraux (p), work (2p), vrow (p). No return value.
 */
void qr_hat (double* x, int n, int p, const double *dmat, const double *zw, double* hat, int *jpvt, double *qraux, double *work, double *vrow)
{
    /* k: the rank of matrix x */
    int k;
    double tol=1e-7;

    for(int i=0; i<p; ++i)
    {
        jpvt[i]=i+1;
    }

    /* dqrdc2 is R's modified LINPACK QR routine (src/appl/dqrdc2.f): it forms the
     * QR factorization of the n-by-p matrix x by Householder transformations, using
     * a tolerance-based limited column pivoting to detect and handle rank deficiency.
     * edgeR calls it to factor the sqrt(W)-weighted design so the QL GLM leverages
     * (hat values) used in the adjusted deviance can be recovered from R.
     * Argument mapping:
     *   x     = x     (n x p weighted design; overwritten in place by R and the
     *                  Householder vectors of the QR factorization)
     *   ldx   = &n    (leading dimension of x)
     *   n     = &n    (number of rows / libraries)
     *   p     = &p    (number of columns / coefficients)
     *   tol   = &tol  (rank-detection tolerance, 1e-7)
     *   k     = &k    (output: computed numerical rank of x)
     *   qraux = qraux (output, length p: auxiliary Householder scalars, unused here)
     *   jpvt  = jpvt  (in/out, length p: column pivot indices, preset to 1..p)
     *   work  = work  (workspace, length 2*p)
     * The leading k columns of the upper triangle of x then hold R_k, and jpvt[0..k-1]
     * name the design columns they belong to.
     * Equivalent R operation: qr(x, tol = 1e-7) (the decomposition inside hat()).
     * Netlib references: https://netlib.org/linpack/dqrdc.f */
    /* call dqrdc2: https://svn.r-project.org/R/branches/R-4-4-branch/src/appl/dqrdc2.f */
    F77_CALL(dqrdc2)(x, &n, &n, &p, &tol, &k, qraux, jpvt, work);

    /* forward solve R_k^T z = x_i for each row i; R_k[l, j] = x[l + n*j] for l <= j */
    for(int i=0; i<n; ++i)
    {
        double h = 0;
        for(int j=0; j<k; ++j)
        {
            double s = dmat[i + (size_t) n*(jpvt[j]-1)] * zw[i];
            const double *rptr = x + (size_t) n*j;
            for(int l=0; l<j; ++l)
            {
                s -= rptr[l] * vrow[l];
            }
            vrow[j] = s / rptr[j];
            h += vrow[j] * vrow[j];
        }
        hat[i] = h;
    }
}
