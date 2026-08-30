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
 * the parallel region. */
typedef struct {
    double *yptr, *uptr, *wptr, *aux;    /* nlib row buffers: y, mu, weights, disp/cover */
    double *xdpt;                        /* nlib*ncoef weighted design */
    double *hptr;                        /* nlib hat values           */
    double *zwpt;                        /* nlib working weights       */
    int    *jpvt;                        /* ncoef  qr_hat pivots       */
    double *qraux;                       /* ncoef  qr_hat              */
    double *qrwork;                      /* 2*ncoef qr_hat work        */
    double *qy;                          /* nlib*ncoef qr_hat Q*y      */
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
        ws[t].qy     = R_Calloc(nlib*ncoef,double);
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
        R_Free(ws[t].qy);
    }
    R_Free(ws);
}

/* Build the sqrt(W)-weighted design and take its hat values for one gene. The working
 * weight zwpt is the only distribution-specific part: NB uses sqrt(mu*w/(1+mu*disp/
 * prior)), binomial uses sqrt(cover*mu*(1-mu)). */
static void ql_tag_hat(ql_ws *w, const double *dmat, int nlib, int ncoef, int fam, double prior)
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
    for(int i=0; i<nlib*ncoef; ++i)
    {
        w->xdpt[i]=dmat[i] * w->zwpt[i % nlib];
    }
    qr_hat(w->xdpt, nlib, ncoef, w->hptr, w->jpvt, w->qraux, w->qrwork, w->qy);
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
        ql_tag_hat(w, design->dmat, nlib, ncoef, fam, prior);

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

    return;
}

/* Shared QL per-observation unit-matrix driver (NB or binomial): unit deviance, unit df,
 * and leverage for every (gene, library). */
void ql_adjust_mat (int fam, cmx *y, cmx *mu, cmx *design, cmx *aux, double prior, cmx *weights, double *dfmat, double *dvmat, double *lvmat, int nthreads)
{
    int ntag=(y->nrow), nlib=(y->ncol), ncoef=(design->ncol);

    int nth = clamp_threads(nthreads);
    ql_ws *ws = ql_ws_alloc(nth, nlib, ncoef);

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
        ql_tag_hat(w, design->dmat, nlib, ncoef, fam, prior);

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
 * sort x and y using rsort_with_index
 * calculate delta = diff(range(x))
 */

static void clowess2 (double *x, double *y, int n, int iter, double f, double *ans)
{
    int *ind = R_Calloc(n, int);
    for(int i=0;i<n;++i)
    {
        ind[i]=i;
    }
    rsort_with_index(x,ind,n);

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

/* hat values by QR decomposition
 * this function is extracted from hat() in stats package
 * hat <- function(x, intercept = TRUE)
 * {
 *   if(is.qr(x)) n <- nrow(x$qr)
 *   else {
 *	    if(intercept) x <- cbind(1, x)
 *	    n <- nrow(x)
 *	    x <- qr(x)
 *   }
 *   rowSums(qr.qy(x, diag(1, nrow = n, ncol = x$rank))^2)
 * }
 *
 * it contains two parts:
 * 1. QR decomposition of x
 * 2. QR_econ for computation of hat values
 */


/*
  inputs:
  x input matrix with n > p
  n number of rows of x
  p number of columns of x

  output:
  hat hat values for the matrix x
*/

void qr_hat (double* x, int n, int p, double* hat, int *jpvt, double *qraux, double *work, double *qy)
{
    /* k: the rank of matrix x */
    int k;

    /* jpvt: pivot vector (caller-owned scratch, size p) */
    for(int i=0; i<p; ++i)
    {
        jpvt[i]=i+1;
    }

    /* preparation for QR decomposition
     * qraux: size p, work: size 2*p (both caller-owned) */
    double tol=1e-7;

    /* dqrdc2 is R's modified LINPACK QR routine (src/appl/dqrdc2.f): it forms the
     * QR factorization of the n-by-p matrix x by Householder transformations, using
     * a tolerance-based limited column pivoting to detect and handle rank deficiency.
     * edgeR calls it to factor the sqrt(W)-weighted design so the QL GLM leverages
     * (hat values) used in the adjusted deviance can be recovered from Q.
     * Argument mapping:
     *   x     = x     (n x p weighted design; overwritten in place by R and the
     *                  Householder vectors of the QR factorization)
     *   ldx   = &n    (leading dimension of x)
     *   n     = &n    (number of rows / libraries)
     *   p     = &p    (number of columns / coefficients)
     *   tol   = &tol  (rank-detection tolerance, 1e-7)
     *   k     = &k    (output: computed numerical rank of x)
     *   qraux = qraux (output, length p: auxiliary Householder scalars defining Q)
     *   jpvt  = jpvt  (in/out, length p: column pivot indices, preset to 1..p)
     *   work  = work  (workspace, length 2*p)
     * The factors and rank k feed the following dqrqy multiply that yields the hats.
     * Equivalent R operation: qr(x, tol = 1e-7) (the decomposition inside hat()).
     * Netlib references: https://netlib.org/linpack/dqrdc.f */
    /* call dqrdc2: https://svn.r-project.org/R/branches/R-4-4-branch/src/appl/dqrdc2.f */
    F77_CALL(dqrdc2)(x, &n, &n, &p, &tol, &k, qraux, jpvt, work);

    /* preparation for QR_econ: qy is caller-owned scratch of size n*p (>= n*k) */
    int nk=n*k;
    double *y = qy;
    for(int i=0; i<nk; ++i)
    {
        y[i]=0;                            // initialization
    }

    double *yptr;
    yptr=y;
    for(int i=0; i<k; ++i, yptr+=n)
    {
        yptr[i]=1;                                               // diagonalization
    }

    /* dqrqy is R's LINPACK-derived routine (src/appl/dqrutl.f): it applies the
     * orthogonal factor Q from dqrdc2 to a matrix, forming qy = Q %*% y without
     * ever assembling Q explicitly.
     * edgeR calls it with the first k unit columns so that Q's leading k columns
     * are produced; their squared row sums are the QL GLM leverages (hat values).
     * Argument mapping:
     *   x     = x     (QR factors from dqrdc2, n x k)
     *   n     = &n    (number of rows / libraries)
     *   k     = &k    (rank: number of Householder steps / columns of Q applied)
     *   qraux = qraux (Householder scalars from dqrdc2)
     *   y     = y     (input, n x k: the first k columns of the identity)
     *   ny    = &k    (number of columns of y)
     *   qy    = y     (output: Q %*% y, written in place over the same buffer)
     * The result holds Q(:,1:k); hat[i] = sum_j Q[i,j]^2 is formed in the loop below.
     * Equivalent R operation: qr.qy(qr(x), diag(1, nrow = n, ncol = k)).
     * Netlib references: https://netlib.org/linpack/dqrsl.f */
    /* call dqrqy: https://svn.r-project.org/R/branches/R-4-4-branch/src/appl/dqrutl.f */
    F77_CALL(dqrqy)(x, &n, &k, qraux, y, &k, y);

    /* compute hat values */

    /*
    for(int i=0;i<n;++i){
      for(int j=0;j<k;++j){
        R_xlen_t ii = (R_xlen_t)(n)*j + i;
        hat[i] += y[ii]*y[ii];
      }
    }
    */

    for(int i=0;i<n;++i)
    {
        yptr=y+i;
        for(int j=0;j<k;++j,yptr+=n)
        {
            hat[i] += fsquare(*yptr);
        }
    }

    return;
}
