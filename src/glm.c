#include "edgeR.h"

/* common inputs:
 * y       count matrix
 * offsets compressedMatrix
 * disp    compressedMatrix
 * weigts  compressedMatrix
 * design  design matrix
 * beta    start values
 *
 * maxit   max iteration
 * tol     tolerance
 *
 * nthreads number of OpenMP threads for the genewise loop
 */

/* Per-thread scratch packed into workspace structs. One workspace per
 * thread is allocated before the parallel region and freed afterwards so
 * that no R_Calloc/R_Free happens inside the region. Each thread writes
 * outputs only at its own tag index, so the shared output pointers need no
 * synchronisation.
 */

/* four nlib row buffers shared by the one-group fitters */
typedef struct {
    double *yptr, *optr, *wptr, *dptr;   /* nlib rows: counts, offsets, weights, disp */
} row4_ws;

/* scratch for fit_leven (per tag) */
typedef struct {
    double *yptr, *optr, *wptr, *dptr;   /* nlib row buffers              */
    double *bptr, *nbeta, *dl, *dbeta;   /* ncoef                         */
    double *xtwx, *xtwxc;                /* ncoef*ncoef Fisher info + copy */
    double *uptr, *nmu;                  /* nlib fitted values            */
    double *zwpt, *drvt;                 /* nlib working weight + deriv    */
} leven_ws;

/* scratch for get_leven_start (per tag) */
typedef struct {
    double *xdpt;                        /* nlib*ncoef design / QR factor  */
    double *tau;                         /* ncoef                          */
    double *effects;                     /* nlib                           */
    double *work_geqrf, *work_ormqr;     /* LAPACK workspaces              */
    double *yptr, *optr, *wptr, *dptr;   /* nlib row buffers               */
} lstart_ws;

/* this function fits one group mean for one row
 * this is converted from glm_one_group.cpp written by Aaron
 */
void glm_one_group_vec(int nlib, double* counts, double* offset, double* disp, double* weights, int maxit, double tolerance, double cur_beta, double *beta, int *conv)
{
    /* Setting up initial values for beta as the log of the mean of the ratio of counts to offsets.
     * This is the exact solution for the gamma distribution (which is the limit of the NB as
     * the dispersion goes to infinity. However, if cur_beta is not NA, then we assume it's good.
     */

    const double low_value = 1e-10;
    int allzero = 1;

    if (ISNA(cur_beta))
    {
        cur_beta=0;
        double totweight=0;
        for (int lib=0; lib<nlib; ++lib)
        {
            double cur_val = counts[lib];
            if (cur_val > low_value)
            {
                cur_beta += cur_val / exp(offset[lib]) * weights[lib];
                allzero = 0;
            }
            totweight += weights[lib];
        }
        cur_beta = log(cur_beta/totweight);
    }
    else
    {
        for (int lib=0; lib<nlib; ++lib)
        {
            if (counts[lib] > low_value)
            {
                allzero = 0;
                break;
            }
        }
    }

    // Skipping to a result for all-zero rows.
    // Use a large finite value (-1e8) rather than -Inf so the oneway dgesv
    // back-transform in fit_glm_mat stays finite (avoids NaN design coefficients);
    // exp(-1e8+offset)=0, so fitted values and deviance are unchanged.
    if (allzero)
    {
        (*beta) = -1e8;
        (*conv) = 1;
        return;
    }

    // Newton-Raphson iterations to converge to mean.
    (*conv)=0;
    for (int i=0; i<maxit; ++i)
    {
        double dl=0;
        double info=0;
        for (int lib=0; lib<nlib; ++lib)
        {
            double mu=exp(cur_beta+offset[lib]), denominator=1+mu*disp[lib];
            dl+=(counts[lib]-mu)/denominator * weights[lib];
            info+=mu/denominator * weights[lib];
        }
        double step=dl/info;
        cur_beta+=step;
        if (fabs(step)<tolerance)
        {
            (*beta)=cur_beta;
            (*conv)=1;
            break;
        }
    }

    return;
}

/* this function fits one group design for matrix input
 *
 * input:
 * beta  starting coefficient
 *
 * outputs:
 * coef  fitted coefficients
 * conv  index of convergence
 *
 * comment:
 * this is converted from R_fit_one_group.cpp written by Aaron
 */
void fit_one_group_mat (cmx *y, cmx *offsets, cmx *disp, cmx *weights, int maxit, double tol, double *beta, double *coef, int *conv, int nthreads)
{
    int ntag = (y->nrow), nlib = (y->ncol);

    int nth = clamp_threads(nthreads);

    /* one set of row vectors (y offsets disp weights) per thread */
    row4_ws *ws = R_Calloc(nth, row4_ws);
    for (int t=0; t<nth; ++t)
    {
        ws[t].yptr = R_Calloc(nlib,double);
        ws[t].optr = R_Calloc(nlib,double);
        ws[t].wptr = R_Calloc(nlib,double);
        ws[t].dptr = R_Calloc(nlib,double);
    }

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
        row4_ws *w = &ws[tid];
        get_row4(y,offsets,disp,weights,tag,w->yptr,w->optr,w->dptr,w->wptr);

        double ocoef;
        int oconv;

        glm_one_group_vec(nlib, w->yptr, w->optr, w->dptr, w->wptr, maxit, tol, beta[tag], &ocoef, &oconv);

        coef[tag]=ocoef;
        conv[tag]=oconv;
    }

    for (int t=0; t<nth; ++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].dptr);
    }
    R_Free(ws);

    return;
}

/* this function computes fitted values for one group fitting
 *
 * inputs:
 * group   index vector to show group info
 *
 * outputs:
 * mu      fitted value matrix
 *
 * comment:
 * This is converted from R_get_one_way_fitted.cpp written by Aaron
*/

void get_one_way_fit(cmx *beta, cmx *offsets, int *group, double *mu, int nthreads)
{
    int ntag=(offsets->nrow), nlib=(offsets->ncol), nbeta=(beta->ncol);

    int nth = clamp_threads(nthreads);

    /* one optr (nlib) and bptr (nbeta) row buffer per thread */
    double **optr_t = R_Calloc(nth, double*);
    double **bptr_t = R_Calloc(nth, double*);
    for(int t=0;t<nth;++t)
    {
        optr_t[t] = R_Calloc(nlib, double);
        bptr_t[t] = R_Calloc(nbeta, double);
    }

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0;tag<ntag;++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *optr = optr_t[tid];
        double *bptr = bptr_t[tid];
        get_row(offsets,tag,optr);
        get_row(beta,tag,bptr);

        /*
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            mu[ii]      = exp(optr[lib]+bptr[group[lib]]);
        }
        */
        double *uptr=mu+tag;
        for(int lib=0;lib<nlib;++lib,uptr+=ntag)
        {
            (*uptr) = exp(optr[lib]+bptr[group[lib]]);
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(optr_t[t]);
        R_Free(bptr_t[t]);
    }
    R_Free(optr_t);
    R_Free(bptr_t);

    return;
}

/* free the per-thread lstart_ws array plus the shared efail buffer; one teardown
 * shared by get_leven_start's three exit paths (the two QR-failure branches and the
 * normal exit) */
static void free_lstart_ws(lstart_ws *ws, int nth, int *efail)
{
    for (int t=0; t<nth; ++t)
    {
        R_Free(ws[t].xdpt);
        R_Free(ws[t].tau);
        R_Free(ws[t].effects);
        R_Free(ws[t].work_geqrf);
        R_Free(ws[t].work_ormqr);
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].dptr);
    }
    R_Free(ws);
    R_Free(efail);
}

/* Factor the design the caller has just written into ws[0].xdpt (plain design for the
 * null fit, sqrt(W)*design for repeated weights) with dgeqrf, then replicate the R
 * factor + reflectors into every thread workspace.
 * dgeqrf (LAPACK) = QR factorization A = Q*R via Householder reflectors; the (weighted)
 * design is identical for every gene, so edgeR factorizes it once and reuses the factor
 * across tags to project each gene's working response onto the design's column space.
 * Argument mapping (Fortran = local): M = nlib rows; N = ncoef cols; A = ws[0].xdpt in:
 * (sqrt(W)*)design out: R (upper) + reflectors (below); LDA = nlib; TAU = ws[0].tau
 * reflector scalars; WORK/LWORK = ws[0].work_geqrf/lwork_geqrf; INFO nonzero ->
 * "QR decomposition failed". Equivalent R operation: qr(design) / qr(sqrt(w)*design).
 * Netlib references: https://netlib.org/lapack/explore-html/ (dgeqrf) */
static void factor_and_replicate(lstart_ws *ws, int nth, int nlib, int ncoef, int lwork_geqrf, int *efail)
{
    int info = 0;
    F77_CALL(dgeqrf)(&nlib, &ncoef, ws[0].xdpt, &nlib, ws[0].tau, ws[0].work_geqrf, &lwork_geqrf, &info);
    if (info)
    {
        free_lstart_ws(ws, nth, efail);
        error("QR decomposition failed");
    }
    for (int t=1; t<nth; ++t)
    {
        for (int i=0; i<nlib*ncoef; ++i)
        {
            ws[t].xdpt[i] = ws[0].xdpt[i];
        }
        for (int i=0; i<ncoef; ++i)
        {
            ws[t].tau[i] = ws[0].tau[i];
        }
    }
}

/* this function computes starting coefficient for levenberg fitting
 *
 * inputs:
 * use_null index of using null method
 *
 * outputs:
 * beta     starting coefficient matrix
 *
 * comment:
 * This is converted from R_initialize_levenberg.cpp written by Aaron
*/
void get_leven_start (cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design, int use_null, double *beta, int nthreads)
{
    // char for fortran call
    const char side='L';
    const char trans_o='T';
    const char uplo='U';
    const char trans_t='N';
    const char diag='N';
    const int unity=1;

    int ntag = (y->nrow), nlib = (y->ncol), ncoef = (design->ncol);
    int lwork_geqrf = -1, lwork_ormqr = -1, info;

    int nth = clamp_threads(nthreads);

    // Setting up the workspace for dgeqrf and dormqr with optimal WORK (the
    // query is independent of the matrix contents, so a scratch buffer is fine)
    double tmpwork;
    double *qbuf = R_Calloc(nlib*ncoef, double);
    double *ebuf = R_Calloc(nlib, double);
    double *tbuf = R_Calloc(ncoef, double);
    /* dgeqrf (LAPACK) = QR factorization A = Q*R of a real M-by-N matrix via
     * Householder reflectors. This call is a workspace-size query only: with
     * LWORK = -1 the routine returns the optimal WORK length in WORK(1) and
     * does not touch the factorization.
     * edgeR runs the query once so the per-thread dgeqrf calls below can size
     * their WORK buffers optimally; the query is independent of the contents.
     * Argument mapping (Fortran = local):
     *   M     = nlib          -> rows of A (libraries)
     *   N     = ncoef         -> columns of A (coefficients)
     *   A     = qbuf          -> scratch M-by-N buffer (untouched by a query)
     *   LDA   = nlib          -> leading dimension of A
     *   TAU   = tbuf          -> reflector scalars (untouched by a query)
     *   WORK  = &tmpwork      -> on exit WORK(1) holds the optimal LWORK
     *   LWORK = lwork_geqrf   -> -1 requests the workspace-size query
     *   INFO  = info          -> 0 on success
     * Result: tmpwork is read into lwork_geqrf to size the real factorizations.
     * A workspace-size query has no algorithmic R equivalent.
     * Netlib references: https://netlib.org/lapack/explore-html/ (dgeqrf)
     */
    F77_CALL(dgeqrf)(&nlib, &ncoef, qbuf, &nlib, tbuf, &tmpwork, &lwork_geqrf, &info);
    lwork_geqrf=(int)(tmpwork+0.5);
    if (lwork_geqrf < 1)
    {
        lwork_geqrf = 1;
    }

    /* dormqr (LAPACK) overwrites C with Q**T * C using the reflectors from
     * dgeqrf. This call is a workspace-size query only: with LWORK = -1 it
     * returns the optimal WORK length in WORK(1) and does not touch C.
     * edgeR runs the query once so the per-thread dormqr calls below can size
     * their WORK buffers optimally; the query is independent of the contents.
     * Argument mapping (Fortran = local):
     *   SIDE  = side ('L')    -> apply Q**T from the left
     *   TRANS = trans_o ('T') -> use Q**T (transpose)
     *   M     = nlib          -> rows of C (libraries)
     *   N     = unity (1)     -> columns of C (single right-hand side)
     *   K     = ncoef         -> number of reflectors
     *   A     = qbuf          -> reflectors from dgeqrf (untouched by a query)
     *   LDA   = nlib          -> leading dimension of A
     *   TAU   = tbuf          -> reflector scalars (untouched by a query)
     *   C     = ebuf          -> scratch M-by-N buffer (untouched by a query)
     *   LDC   = nlib          -> leading dimension of C
     *   WORK  = &tmpwork      -> on exit WORK(1) holds the optimal LWORK
     *   LWORK = lwork_ormqr   -> -1 requests the workspace-size query
     *   INFO  = info          -> 0 on success
     *   trailing FCONE FCONE  -> hidden Fortran lengths of SIDE and TRANS
     * Result: tmpwork is read into lwork_ormqr to size the real Q**T products.
     * A workspace-size query has no algorithmic R equivalent.
     * Netlib references: https://netlib.org/lapack/explore-html/ (dormqr)
     */
    F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, qbuf, &nlib, tbuf, ebuf, &nlib, &tmpwork, &lwork_ormqr, &info FCONE FCONE);
    lwork_ormqr=(int)(tmpwork+0.5);
    if (lwork_ormqr < 1)
    {
        lwork_ormqr = 1;
    }
    R_Free(qbuf);
    R_Free(ebuf);
    R_Free(tbuf);

    /* one workspace per thread */
    lstart_ws *ws = R_Calloc(nth, lstart_ws);
    for (int t=0; t<nth; ++t)
    {
        ws[t].xdpt       = R_Calloc(nlib*ncoef, double);
        ws[t].tau        = R_Calloc(ncoef, double);
        ws[t].effects    = R_Calloc(nlib, double);
        ws[t].work_geqrf = R_Calloc(lwork_geqrf, double);
        ws[t].work_ormqr = R_Calloc(lwork_ormqr, double);
        ws[t].yptr       = R_Calloc(nlib,double);
        ws[t].optr       = R_Calloc(nlib,double);
        ws[t].wptr       = R_Calloc(nlib,double);
        ws[t].dptr       = R_Calloc(nlib,double);
    }

    /* per-thread LAPACK failure stage: 1=QR, 2=Q**T multiply, 3=triangular solve */
    int *efail = R_Calloc(nth, int);

    if(use_null)
    {
        // make a copy of design matrix and factor it once (shared, read-only in
        // dormqr/dtrtrs); replicate the factor into every thread workspace
        for(int i=0;i<nlib*ncoef;++i)
        {
            ws[0].xdpt[i] = (design->dmat)[i];
        }

        factor_and_replicate(ws, nth, nlib, ncoef, lwork_geqrf, efail);

        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nth) schedule(static)
        #endif
        for(int tag=0;tag<ntag;++tag)
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            lstart_ws *w = &ws[tid];
            int linfo=0;
            get_row4(y,offsets,disp,weights,tag,w->yptr,w->optr,w->dptr,w->wptr);

            // Computing weighted average of the count:library size ratios.
            double sum_weight=0, sum_exprs=0;
            for (int lib=0; lib<nlib; ++lib)
            {
                double curN=exp(w->optr[lib]);
                double curweight=w->wptr[lib]*curN/(1 + w->dptr[lib] * curN);
                sum_exprs  += w->yptr[lib] * curweight / curN;
                sum_weight += curweight;
            }

            for (int lib=0; lib<nlib; ++lib)
            {
                w->effects[lib]=log(sum_exprs/sum_weight);
            }
            /*  DORMQR overwrites the general real M-by-N matrix C with Q**T * C
             *  using the reflectors returned by DGEQRF.
             *  https://netlib.org/lapack/explore-3.1.1-html/dormqr.f.html
             *
             *  dormqr (LAPACK) overwrites C with Q**T * C using the reflectors from dgeqrf.
             *  Here it projects this gene's working response 'effects' onto the QR basis of
             *  the design, i.e. forms Q**T * effects ready for the triangular solve.
             *  Argument mapping (Fortran = local):
             *    SIDE  = side ('L')    -> apply Q**T from the left
             *    TRANS = trans_o ('T') -> use Q**T (transpose)
             *    M     = nlib          -> rows of C (libraries)
             *    N     = unity (1)     -> single right-hand side
             *    K     = ncoef         -> number of reflectors
             *    A     = w->xdpt       -> QR factor from dgeqrf
             *    LDA   = nlib          -> leading dimension of A
             *    TAU   = w->tau        -> reflector scalars
             *    C     = w->effects    -> in: working response; out: Q**T * effects (overwritten)
             *    LDC   = nlib          -> leading dimension of C
             *    WORK  = w->work_ormqr -> workspace
             *    LWORK = lwork_ormqr   -> workspace length from the query above
             *    INFO  = linfo         -> 0 on success; nonzero sets efail[tid]=2
             *    trailing FCONE FCONE  -> hidden Fortran lengths of SIDE and TRANS
             *  Result: effects holds Q**T * (working response) for the back-substitution.
             *  Equivalent R operation: qr.qty(qr(design), effects).
             *  Netlib references: https://netlib.org/lapack/explore-html/ (dormqr)
             */
            F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, w->xdpt, &nlib, w->tau, w->effects, &nlib, w->work_ormqr, &lwork_ormqr, &linfo FCONE FCONE);
            if (linfo)
            {
                efail[tid]=2;
                continue;
            }

            /*  DTRTRS solves a triangular system A * X = B.
             *  https://netlib.org/lapack/explore-3.1.1-html/dtrtrs.f.html
             *
             *  dtrtrs (LAPACK) solves the triangular system A * X = B in place.
             *  Here A is the upper-triangular factor R from dgeqrf, so this back-solves
             *  R * beta = Q**T * effects to recover the starting coefficients.
             *  Argument mapping (Fortran = local):
             *    UPLO  = uplo ('U')    -> A is upper triangular (R)
             *    TRANS = trans_t ('N') -> solve A * X = B (no transpose)
             *    DIAG  = diag ('N')    -> A is non-unit triangular
             *    N     = ncoef         -> order of A
             *    NRHS  = unity (1)     -> single right-hand side
             *    A     = w->xdpt       -> R stored in the QR factor
             *    LDA   = nlib          -> leading dimension of A
             *    B     = w->effects    -> in: Q**T * response; out: solution beta (overwritten)
             *    LDB   = nlib          -> leading dimension of B
             *    INFO  = linfo         -> 0 on success; nonzero sets efail[tid]=3
             *    trailing FCONE x3     -> hidden Fortran lengths of UPLO, TRANS, DIAG
             *  Result: effects holds the starting coefficients for this gene.
             *  Equivalent R operation: backsolve(R, qty).
             *  Netlib references: https://netlib.org/lapack/explore-html/ (dtrtrs)
             */
            F77_CALL(dtrtrs)(&uplo, &trans_t, &diag, &ncoef, &unity, w->xdpt, &nlib, w->effects, &nlib, &linfo FCONE FCONE FCONE);
            if (linfo)
            {
                efail[tid]=3;
                continue;
            }

            double *bptr = beta+tag;
            for(int var=0;var<ncoef;++var,bptr+=ntag)
            {
                (*bptr) = w->effects[var];
            }
        }
    }
    else
    {
        // Finding the delta
        double delta = max_cmx(y);
        delta = fmin(delta, 1.0/6);

        int weights_repeated = ((weights->type) >= 2);

        // check whether weights is row repeated: factor once and replicate
        if(weights_repeated)
        {
            int tag_start = 0;
            get_row(weights,tag_start,ws[0].wptr);

            for(int i=0;i<nlib*ncoef;++i)
            {
                ws[0].xdpt[i]=(design->dmat)[i]*sqrt(ws[0].wptr[i % nlib]);
            }
            factor_and_replicate(ws, nth, nlib, ncoef, lwork_geqrf, efail);
        }

        #ifdef _OPENMP
        #pragma omp parallel for num_threads(nth) schedule(static)
        #endif
        for(int tag=0;tag<ntag;++tag)
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            lstart_ws *w = &ws[tid];
            int linfo=0;
            get_row4(y,offsets,disp,weights,tag,w->yptr,w->optr,w->dptr,w->wptr);

            if(!weights_repeated)
            {
                for(int i=0;i<nlib*ncoef;++i)
                {
                    w->xdpt[i]=(design->dmat)[i]*sqrt(w->wptr[i % nlib]);
                }
                /* dgeqrf (LAPACK) = QR factorization A = Q*R via Householder reflectors.
                 * With per-gene weights the sqrt(W)-weighted design differs by tag, so each
                 * thread factorizes its own weighted design copy before the Q**T product and
                 * triangular solve below.
                 * Argument mapping (Fortran = local):
                 *   M     = nlib            -> rows of A (libraries)
                 *   N     = ncoef           -> columns of A (coefficients)
                 *   A     = w->xdpt         -> in: sqrt(W)*design; out: R (upper) + reflectors (below)
                 *   LDA   = nlib            -> leading dimension of A
                 *   TAU   = w->tau          -> out: Householder reflector scalars
                 *   WORK  = w->work_geqrf   -> workspace
                 *   LWORK = lwork_geqrf     -> workspace length from the query above
                 *   INFO  = linfo           -> 0 on success; nonzero sets efail[tid]=1
                 * Result: w->xdpt/tau hold this gene's weighted-design QR factor.
                 * Equivalent R operation: qr(sqrt(w) * design).
                 * Netlib references: https://netlib.org/lapack/explore-html/ (dgeqrf)
                 */
                F77_CALL(dgeqrf)(&nlib, &ncoef, w->xdpt, &nlib, w->tau, w->work_geqrf, &lwork_geqrf, &linfo);
                if (linfo)
                {
                    efail[tid]=1;
                    continue;
                }
            }

            // Computing normalized log-expression values.
            for (int lib=0; lib<nlib; ++lib)
            {
                w->yptr[lib]=log(fmax(delta, w->yptr[lib])) - w->optr[lib];
            }

            for (int lib=0; lib<nlib; ++lib)
            {
                w->effects[lib]=w->yptr[lib]*sqrt(w->wptr[lib]);
            }

            /* dormqr (LAPACK) overwrites C with Q**T * C using the reflectors from dgeqrf.
             * Here it forms Q**T * effects, projecting this gene's sqrt(W)-weighted working
             * response onto the QR basis of the weighted design ready for the triangular solve.
             * Argument mapping (Fortran = local):
             *   SIDE  = side ('L')    -> apply Q**T from the left
             *   TRANS = trans_o ('T') -> use Q**T (transpose)
             *   M     = nlib          -> rows of C (libraries)
             *   N     = unity (1)     -> single right-hand side
             *   K     = ncoef         -> number of reflectors
             *   A     = w->xdpt       -> QR factor from dgeqrf
             *   LDA   = nlib          -> leading dimension of A
             *   TAU   = w->tau        -> reflector scalars
             *   C     = w->effects    -> in: sqrt(W)*response; out: Q**T * effects (overwritten)
             *   LDC   = nlib          -> leading dimension of C
             *   WORK  = w->work_ormqr -> workspace
             *   LWORK = lwork_ormqr   -> workspace length from the query above
             *   INFO  = linfo         -> 0 on success; nonzero sets efail[tid]=2
             *   trailing FCONE FCONE  -> hidden Fortran lengths of SIDE and TRANS
             * Result: effects holds Q**T * (weighted response) for the back-substitution.
             * Equivalent R operation: qr.qty(qr(sqrt(w) * design), effects).
             * Netlib references: https://netlib.org/lapack/explore-html/ (dormqr)
             */
            F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, w->xdpt, &nlib, w->tau, w->effects, &nlib, w->work_ormqr, &lwork_ormqr, &linfo FCONE FCONE);
            if (linfo)
            {
                efail[tid]=2;
                continue;
            }

            /* dtrtrs (LAPACK) solves a triangular system A * X = B in place.
             * Here A is the upper-triangular factor R from dgeqrf, so this back-solves
             * R * beta = Q**T * (weighted response) to recover the starting coefficients.
             * Argument mapping (Fortran = local):
             *   UPLO  = uplo ('U')    -> A is upper triangular (R)
             *   TRANS = trans_t ('N') -> solve A * X = B (no transpose)
             *   DIAG  = diag ('N')    -> A is non-unit triangular
             *   N     = ncoef         -> order of A
             *   NRHS  = unity (1)     -> single right-hand side
             *   A     = w->xdpt       -> R stored in the QR factor
             *   LDA   = nlib          -> leading dimension of A
             *   B     = w->effects    -> in: Q**T * response; out: solution beta (overwritten)
             *   LDB   = nlib          -> leading dimension of B
             *   INFO  = linfo         -> 0 on success; nonzero sets efail[tid]=3
             *   trailing FCONE x3     -> hidden Fortran lengths of UPLO, TRANS, DIAG
             * Result: effects holds the starting coefficients for this gene.
             * Equivalent R operation: backsolve(R, qty).
             * Netlib references: https://netlib.org/lapack/explore-html/ (dtrtrs)
             */
            F77_CALL(dtrtrs)(&uplo, &trans_t, &diag, &ncoef, &unity, w->xdpt, &nlib, w->effects, &nlib, &linfo FCONE FCONE FCONE);
            if (linfo)
            {
                efail[tid]=3;
                continue;
            }

            double *bptr=beta+tag;
            for(int coef=0;coef<ncoef;++coef,bptr+=ntag)
            {
                (*bptr) = w->effects[coef];
            }
        }
    }

    int stage=0;
    for(int t=0;t<nth;++t)
    {
        if(efail[t])
        {
            stage=efail[t];
        }
    }

    free_lstart_ws(ws, nth, efail);

    if (stage==1)
    {
        error("QR decomposition failed");
    }
    if (stage==2)
    {
        error("Q**T multiplication failed");
    }
    if (stage==3)
    {
        error("failed to solve the triangular system");
    }

    return;
}

/* this function fits levenberg method
 *
 * outputs:
 * mbeta   updated coefficents
 * mu      fitted values
 * dev     deviance
 * iter    number of iteration
 * fail    index of convergence if failed

 * comment:
 * This is converted from R_fit_levenberg.cpp written by Aaron
 */

void fit_leven (cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design, cmx *beta, double tol, int maxit,
                double *mu, double *mbeta, double *dev, int *iter, int *failed, int nthreads)
{
    int ntag = (y->nrow), nlib = (y->ncol), ncoef = (design->ncol);
    double *dm;
    dm = (design->dmat);

    int nth = clamp_threads(nthreads);

    /* one full scratch workspace per thread */
    leven_ws *ws = R_Calloc(nth, leven_ws);
    for (int t=0; t<nth; ++t)
    {
        ws[t].yptr  = R_Calloc(nlib,double);
        ws[t].optr  = R_Calloc(nlib,double);
        ws[t].wptr  = R_Calloc(nlib,double);
        ws[t].dptr  = R_Calloc(nlib,double);
        ws[t].bptr  = R_Calloc(ncoef,double);
        ws[t].nbeta = R_Calloc(ncoef,double);
        ws[t].dl    = R_Calloc(ncoef,double);
        ws[t].dbeta = R_Calloc(ncoef,double);
        ws[t].xtwx  = R_Calloc(ncoef*ncoef,double);
        ws[t].xtwxc = R_Calloc(ncoef*ncoef,double);
        ws[t].uptr  = R_Calloc(nlib,double);
        ws[t].nmu   = R_Calloc(nlib,double);
        ws[t].zwpt  = R_Calloc(nlib,double);
        ws[t].drvt  = R_Calloc(nlib,double);
    }

    /* per-thread Cholesky-solve failure flag (raised after the region) */
    int *cfail = R_Calloc(nth, int);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        leven_ws *w = &ws[tid];
        int oiter, ofail, oerr=0;
        double odev;

        get_row4(y,offsets,disp,weights,tag,w->yptr,w->optr,w->dptr,w->wptr);
        get_row(beta,tag,w->bptr);

        fit_leven_vec(nlib,w->yptr,w->optr,w->dptr,w->wptr,ncoef,dm,maxit,tol,w->zwpt,w->drvt,w->dl,w->dbeta,w->xtwx,w->xtwxc,w->nbeta,w->nmu,w->bptr,w->uptr,&odev,&oiter,&ofail,&oerr);
        if (oerr)
        {
            cfail[tid]=1;
        }

        double *uupt=mu+tag;
        for(int lib=0;lib<nlib;++lib,uupt+=ntag)
        {
            (*uupt) = w->uptr[lib];
        }
        double *bbpt=mbeta+tag;
        for(int coef=0;coef<ncoef;++coef,bbpt+=ntag)
        {
            (*bbpt) = w->bptr[coef];
        }

        dev[tag]    = odev;
        iter[tag]   = oiter;
        failed[tag] = ofail;
    }

    int solve_failed=0;
    for(int t=0;t<nth;++t)
    {
        solve_failed |= cfail[t];
    }

    for (int t=0; t<nth; ++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].dptr);
        R_Free(ws[t].bptr);
        R_Free(ws[t].nbeta);
        R_Free(ws[t].dl);
        R_Free(ws[t].dbeta);
        R_Free(ws[t].xtwx);
        R_Free(ws[t].xtwxc);
        R_Free(ws[t].uptr);
        R_Free(ws[t].nmu);
        R_Free(ws[t].zwpt);
        R_Free(ws[t].drvt);
    }
    R_Free(ws);
    R_Free(cfail);

    if (solve_failed)
    {
        error("solution using the Cholesky decomposition failed");
    }

    return;
}

/* this function fits levenberg for one row
 *
 * inputs:
 * nlib    number of samples
 * y       count vector
 * offset  offset vector
 * disp    dispersion vector
 * w       weight vector
 *
 * ncoef    number of coefficients
 * dm      design matrix
 *
 * maxit   number of maximal iteration
 * tol     tolerance
 *
 * working space:
 *
 * zwpt    working weight vector
 * drvt    derivate vector
 * dl      difference of likelihood
 * db      difference of beta
 *
 * xtwx    fisher info matrix
 * xtwc    copy of xtwx
 *
 * nbt     updated beta
 * nmu     updated fitted values
 *
 * outputs:
 * obt     coefficents
 * omu     fitted values
 * odev    deviance
 * oiter   number of iteration
 * ofail   failure of convergence
 *
 * comment:
 * This is converted from glm_levenberg.cpp written by Aaron
 */

/* Working weights and the derivative contribution for one Newton-Raphson step.
 * For each library: denom = 1 + mu*disp is the NB variance-to-mean adjustment;
 * 'zwpt' is the working weight mu/denom (times the user weight) that fills the
 * diagonal of W in XtWX, and 'drvt' is (y-mu)/denom (times the user weight),
 * the per-observation contribution to the log-likelihood derivative 'dl'.
 * Factored out of the outer loop for readability; arithmetic is unchanged. */
static void compute_working_weights(int nlib, const double *omu, const double *y,
                                    const double *disp, const double *w,
                                    double *zwpt, double *drvt)
{
    for (int lib=0; lib<nlib; ++lib)
    {
        double cur_mu=omu[lib];
        double denom=(1+cur_mu*disp[lib]);
        zwpt[lib]=cur_mu/denom*w[lib];
        drvt[lib]=(y[lib]-cur_mu)/denom*w[lib];
    }
}

/* file-local: recompute fitted means from coefficients; defined below */
static void fit_leven_autofill(int, double *, int, double *, double *, double *);

void fit_leven_vec(int nlib, double *y, double *offset, double *disp, double *w, int ncoef, double *dm, int maxit, double tol,
                   double *zwpt, double *drvt, double *dl, double *db, double *xtwx, double *xtwc, double *nbt, double *nmu,
                   double *obt, double *omu, double *odev, int *oiter, int *ofail, int *oerr)
{
    // const values
    const double low_value = 1e-10;
    const double one_millionth = 1e-6;
    const double supremely_low_value = 1e-13;
    const double ridiculously_low_value = 1e-100;

    const char uplo='U';
    const int nrhs=1;

    /* LAPACK failure flag: replaces an in-loop error() so the caller can
     * raise it outside any parallel region */
    (*oerr) = 0;

    // We expect 'beta' to be supplied. We then check the maximum value of the counts.
    double ymax=0;
    for (int lib=0; lib<nlib; ++lib)
    {
        ymax = (y[lib]>ymax)? y[lib] : ymax;
    }

    // If we start off with all entries at zero, there's really no point continuing.
    if (ymax<low_value)
    {
        for(int coef=0;coef<ncoef; ++coef)
        {
            obt[coef] = NA_REAL;
        }
        for(int lib=0;lib<nlib; ++lib)
        {
            omu[lib] = 0;
        }
        (*odev)  = 0;
        (*oiter) = 0;
        (*ofail) = 0;
        return;
    }

    // Otherwise, we compute 'mu' based on 'beta'. Returning if there are no coefficients!
    fit_leven_autofill(ncoef,obt,nlib,offset,omu,dm);

    double dev=nb_deviance_sum(nlib,y,omu,disp,w);

    // Iterating using reweighted least squares; setting up assorted temporary objects.
    double max_info=-1, lambda=0;
    int info=0, iter=0, failed=0;
    while ((++iter) <= maxit)
    {

        /* Here we set up the matrix XtWX i.e. the Fisher information matrix. X is the design matrix and W is a diagonal matrix
         * with the working weights for each observation (i.e. library). The working weights are part of the first derivative of
         * the log-likelihood for a given coefficient, multiplied by any user-specified weights. When multiplied by two covariates
         * in the design matrix, you get the Fisher information (i.e. variance of the log-likelihood) for that pair. This takes
         * the role of the second derivative of the log-likelihood. The working weights are formed by taking the reciprocal of the
         * product of the variance (in terms of the mean) and the square of the derivative of the link function.
         *
         * We also set up the actual derivative of the log likelihoods in 'dl'. This is done by multiplying each covariate by the
         * difference between the mu and observation and dividing by the variance and derivative of the link function. This is
         * then summed across all observations for each coefficient. The aim is to solve (XtWX)(dbeta)=dl for 'dbeta'. As XtWX
         * is the second derivative, and dl is the first, you can see that we are effectively performing a multivariate
         * Newton-Raphson procedure with 'dbeta' as the step.
         */
        compute_working_weights(nlib, omu, y, disp, w, zwpt, drvt);

        compute_xtwx(nlib, ncoef, dm, zwpt, xtwx);

        double *dmc, *xtwxIt, *xtwcIt;
        dmc=dm;
        xtwxIt=xtwx;
        for (int coef=0; coef<ncoef; ++coef, dmc+=nlib, xtwxIt+=ncoef)
        {
            dl[coef]=0;
            for(int lib=0;lib<nlib;++lib)
            {
                dl[coef] += drvt[lib]*dmc[lib];
            }
            if (xtwxIt[coef]>max_info)
            {
                max_info=xtwxIt[coef];
            }
        }
        if (iter==1)
        {
            lambda=max_info*one_millionth;
            if (lambda < supremely_low_value)
            {
                lambda=supremely_low_value;
            }
        }

        /* Levenberg/Marquardt damping reduces step size until the deviance increases or no
         * step can be found that increases the deviance. In short, increases in the deviance
         * are enforced to avoid problems with convergence.
         */
        int lev=0, low_dev=0;
        while (++lev)
        {
            do
            {
                /* We need to set up copies as the decomposition routine overwrites the originals, and
                 * we want the originals in case we don't like the latest step. For efficiency, we only
                 * refer to the upper triangular for the XtWX copy (as it should be symmetrical). We also add
                 * 'lambda' to the diagonals. This reduces the step size as the second derivative is increased.
                 */
                xtwxIt=xtwx;
                xtwcIt=xtwc;
                for (int coef1=0; coef1<ncoef; ++coef1, xtwxIt+=ncoef, xtwcIt+=ncoef)
                {
                    for(int coef2=0;coef2<=coef1;++coef2)
                    {
                        xtwcIt[coef2]=xtwxIt[coef2];
                    }
                    xtwcIt[coef1] += lambda;
                }

                // Cholesky decomposition, and then use of the decomposition to solve for dbeta in (XtWX)dbeta = dl.
                /*  DPOTRF computes the Cholesky factorization of a real symmetric
                 *  positive definite matrix A.
                 *
                 *  The factorization has the form
                 *     A = U**T * U,  if UPLO = 'U', or
                 *     A = L  * L**T,  if UPLO = 'L',
                 *  where U is an upper triangular matrix and L is lower triangular.
                 *
                 *  https://www.netlib.org/lapack/lapack-3.1.1/html/dpotrf.f.html
                 *
                 *  Role here: factorize the damped information matrix so the next
                 *  Newton-Raphson step can be solved cheaply by back-substitution.
                 *  Arguments in this context:
                 *    uplo  = 'U'   -> reference/produce the upper triangular factor U
                 *    ncoef         -> order of the matrix (number of coefficients)
                 *    xtwc          -> in: damped XtWX (upper triangle); out: factor U (overwritten)
                 *    ncoef         -> leading dimension of xtwc
                 *    info          -> 0 on success; >0 means the matrix was not positive definite
                 *  The trailing FCONE passes the hidden Fortran length of the single
                 *  character argument 'uplo' (required under R's USE_FC_LEN_T; one
                 *  FCONE per char* arg). See the FCONE setup in edgeR.h.
                 *  Equivalent R operation: chol(XtWX_damped) (the upper factor U).
                 *  Netlib references: https://netlib.org/lapack/explore-html/ (dpotrf)
                 */
                F77_CALL(dpotrf)(&uplo, &ncoef, xtwc, &ncoef, &info FCONE);
                if (info!=0)
                {
                    /* If it fails, it MUST mean that the matrix is singular due to numerical imprecision
                     * as all the diagonal entries of the XtWX matrix must be positive. This occurs because of
                     * fitted values being exactly zero; thus, the coefficients attempt to converge to negative
                     * infinity. This generally forces the step size to be larger (i.e. lambda lower) in order to
                     * get to infinity faster (which is impossible). Low lambda leads to numerical instability
                     * and effective singularity. To solve this, we actually increase lambda; this avoids code breakage
                     * to give the other coefficients a chance to converge. Failure of convergence for the zero-
                     * fitted values isn't a problem as the change in deviance from small --> smaller coefficients isn't
                     * that great when the true value is negative inifinity.
                     */
                    lambda*=10;
                    if (lambda <= 0)
                    {
                        lambda=ridiculously_low_value; // Just to make sure it actually increases.
                    }
                }
                else
                {
                    break;
                }
            } while (1);

            for(int coef=0;coef<ncoef;++coef)
            {
                db[coef]=dl[coef];
            }

            /*  DPOTRS solves a system of linear equations A*X = B with a symmetric
             *  positive definite matrix A using the Cholesky factorization
             *  A = U**T*U or A = L*L**T computed by DPOTRF.
             *
             *  https://www.netlib.org/lapack/lapack-3.1.1/html/dpotrs.f.html
             *
             *  Role here: solve (XtWX_damped) * dbeta = dl for the NR update 'dbeta',
             *  using the Cholesky factor produced by dpotrf above.
             *  Arguments in this context:
             *    uplo  = 'U'   -> the factor stored in xtwc is the upper triangle U
             *    ncoef         -> order of the system (number of coefficients)
             *    nrhs  = 1     -> a single right-hand side
             *    xtwc          -> Cholesky factor U from dpotrf
             *    ncoef         -> leading dimension of xtwc
             *    db            -> in: RHS 'dl'; out: solution 'dbeta' (overwritten)
             *    ncoef         -> leading dimension of db
             *    info          -> 0 on success
             *  The trailing FCONE passes the hidden Fortran length of the single
             *  character argument 'uplo' (required under R's USE_FC_LEN_T; one
             *  FCONE per char* arg). See the FCONE setup in edgeR.h.
             *  Equivalent R operation: backsolve()/solve() using the Cholesky factor.
             *  Netlib references: https://netlib.org/lapack/explore-html/ (dpotrs)
             */
            F77_CALL(dpotrs)(&uplo, &ncoef, &nrhs, xtwc, &ncoef, db, &ncoef, &info FCONE);
            if (info!=0)
            {
                /* signal the failure and bail out; the caller raises the
                 * error outside the parallel region */
                (*oerr)  = 1;
                (*odev)  = dev;
                (*oiter) = iter;
                (*ofail) = 1;
                return;
            }

            // Updating beta and the means. 'dbeta' stores 'Y' from the solution of (X*VX)Y=dl, corresponding to a NR step.
            for (int coef=0; coef<ncoef; ++coef)
            {
                nbt[coef]=obt[coef]+db[coef];
            }

            fit_leven_autofill(ncoef,nbt,nlib,offset,nmu,dm);

            /* Checking if the deviance has decreased or if it's too small to care about. Either case is good
             * and means that we'll be using the updated fitted values and coefficients. Otherwise, if we have
             * to repeat the inner loop, then we want to do so from the original values (as we'll be scaling
             * lambda up so we want to retake the step from where we were before). This is why we don't modify the values
             * in-place until we're sure we want to take the step.
             */
            double ndev=nb_deviance_sum(nlib,y,nmu,disp,w);

            if (ndev/ymax < supremely_low_value)
            {
                low_dev=1;
            }
            if (ndev <= dev || low_dev)
            {
                for(int coef=0;coef<ncoef;++coef)
                {
                    obt[coef]=nbt[coef];
                }
                for(int lib=0;lib<nlib;++lib)
                {
                    omu[lib]=nmu[lib];
                }
                dev=ndev;
                break;
            }

            // Increasing lambda, to increase damping. Again, we have to make sure it's not zero.
            lambda*=2;
            if (lambda <= 0)
            {
                lambda=ridiculously_low_value;
            }

            // Excessive damping; steps get so small that it's pointless to continue.
            if (lambda/max_info > 1/supremely_low_value)
            {
                failed=1;
                break;
            }
        }

        /* Terminating if we failed, if divergence from the exact solution is acceptably low
         * (cross-product of dbeta with the log-likelihood derivative) or if the actual deviance
         * of the fit is acceptably low.
         */
        double divergence=0;
        for(int coef=0;coef<ncoef;++coef)
        {
            divergence += dl[coef]*db[coef];
        }

        if ( failed || low_dev || (divergence<tol))
        {
            break;
        }

        /* If we quit the inner levenberg loop immediately and survived all the break conditions above, that means that deviance is decreasing
         * substantially. Thus, we need larger steps to get there faster. To do so, we decrease the damping factor. Note that this only applies
         * if we didn't decrease the damping factor in the inner levenberg loop, as that would indicate that we need to slow down.
         */
        if (lev==1)
        {
            lambda/=10;
        }
    }

    /* Both exit routes (the convergence/failure break above and natural
     * exhaustion of maxit) land here, so report the outputs once. On a maxit
     * exit, dev/iter/failed hold the last accepted values; failed stays 0 as it
     * is reserved for the Levenberg damping failure, not for hitting maxit. */
    (*odev)  = dev;
    (*oiter) = iter;
    (*ofail) = failed;
    return;
}

/* Compute the fitted mean for one gene from its coefficients (log link).
 * Forms mu = exp(dm * beta + offset) via a single DGEMV plus an exp() loop.
 *
 * inputs:
 * ncoef   number of coefficients
 * beta    coefficient vector (length ncoef)
 * nlib    number of libraries
 * offset  offset vector (length nlib)
 * dm      design matrix (nlib x ncoef, column-major)
 *
 * output:
 * mu      fitted values (length nlib), overwritten
 *
 * No return value; writes only into the caller-owned 'mu' buffer.
 */
// Computing updated mean = beta %*% design + offset.
static void fit_leven_autofill(int ncoef, double *beta, int nlib, double *offset, double *mu, double *dm)
{
    for(int lib=0; lib<nlib; ++lib)
    {
        mu[lib]=offset[lib];
    }

    const char trans = 'N';
    const int incx=1, incy=1;
    const double first_scaling=1, second_scaling=1;

    /*  DGEMV  performs one of the matrix-vector operations
     *
     *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are vectors and A is an
     *  m by n matrix.
     *
     *  https://www.netlib.org/lapack/lapack-3.1.1/html/dgemv.f.html
     *
     *  Role here: form the linear predictor mu = dm * beta + offset. 'mu' is
     *  pre-seeded with 'offset' above, so we use beta(=1)*mu to keep it; the
     *  exp() loop below then applies the log link.
     *  Arguments in this context:
     *    trans = 'N'        -> no transpose; compute A*x
     *    nlib               -> rows of A (m): number of libraries
     *    ncoef              -> cols of A (n): number of coefficients
     *    first_scaling = 1  -> alpha
     *    dm                 -> design matrix A (nlib x ncoef, column-major)
     *    nlib               -> leading dimension of dm
     *    beta               -> coefficient vector x
     *    incx = 1           -> stride of x
     *    second_scaling = 1 -> beta scalar: keep the offset already in mu
     *    mu                 -> in: offset (y); out: offset + dm*beta
     *    incy = 1           -> stride of y
     *  The trailing FCONE passes the hidden Fortran length of the single
     *  character argument 'trans' (required under R's USE_FC_LEN_T; one FCONE
     *  per char* arg). See the FCONE setup in edgeR.h.
     *  Equivalent R operation: mu <- design %*% beta + offset.
     *  Netlib references: https://netlib.org/lapack/explore-html/ (dgemv)
     */
    F77_CALL(dgemv)(&trans, &nlib, &ncoef, &first_scaling, dm, &nlib, beta, &incx, &second_scaling, mu, &incy FCONE);

    for (int lib=0; lib<nlib; ++lib)
    {
        mu[lib]=exp(mu[lib]);
    }

    return;
}

/* this function checks poisson bound for legacy QL method
 *
 * inputs:
 * mu   matrix of fitted values
 * s2   compressedMatrix of quasi dispersion
 *
 * output:
 * out  index vector
 *
 * comment:
 * this is converted from R_check_poisson_bound.cpp written by Aaron
 */

void check_poi_bound(cmx *mu, cmx *disp, cmx *s2, int *out, int nthreads)
{
    int ntag=(mu->nrow), nlib=(mu->ncol);

    int nth = clamp_threads(nthreads);

    /* one dptr/sptr (nlib) row buffer per thread */
    double **dptr_t = R_Calloc(nth, double*);
    double **sptr_t = R_Calloc(nth, double*);
    for(int t=0;t<nth;++t)
    {
        dptr_t[t] = R_Calloc(nlib, double);
        sptr_t[t] = R_Calloc(nlib, double);
    }

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0;tag<ntag;++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *dptr = dptr_t[tid];
        double *sptr = sptr_t[tid];
        get_row(disp,tag,dptr);
        get_row(s2,tag,sptr);

        /* a gene is below poisson bound if existing i,
         * such that s2*(1+u_i*d_i) < 1, that is
         * s2 < (1+u_i*d_i)^{-1}, which means
         * the quasi dispersion is too small
         */
        out[tag]=0;
        double *uptr=(mu->dmat)+tag;
        for(int lib=0;lib<nlib;++lib,uptr+=ntag)
        {
            if(sptr[lib]*(1+dptr[lib]*(*uptr)) < 1)
            {
                out[tag]=1;
                break;
            }
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(dptr_t[t]);
        R_Free(sptr_t[t]);
    }
    R_Free(dptr_t);
    R_Free(sptr_t);

    return;
}

/* per-thread scratch for fit_one_way_mat: per-group working vectors and the
 * four row buffers */
typedef struct {
    double *ysptr, *osptr, *dsptr, *wsptr;   /* nlibs per-group working */
    double *yptr, *optr, *dptr, *wptr;       /* nlibs row buffers       */
} oneway_ws;

/* Fit a one-way-layout NB GLM group mean for every gene, parallelised over tags.
 * Each group's libraries are fitted with glm_one_group_vec.
 *
 * inputs:
 * y, offsets, disp, weights  genewise data (compressed matrices)
 * group       per-library group index (length nlibs)
 * ngroups     number of groups
 * maxit, tol  Newton-Raphson controls
 * coef_start  optional starting coefficients (ntag x ngroups) or NULL
 *
 * outputs (caller-allocated):
 * coef  fitted group coefficients (ntag x ngroups)
 * conv  per-gene, per-group convergence indicator (ntag x ngroups)
 *
 * No return value; writes only into the caller-owned output arrays.
 */
/* fit one-way layout for all groups, genewise (parallelised over tags) */
void fit_one_way_mat(cmx *y, cmx *offsets, cmx *disp, cmx *weights, int *group,
                     int ngroups, int maxit, double tol, double *coef_start,
                     double *coef, int *conv, int nthreads)
{
    int ntag = (y->nrow), nlibs = (y->ncol);

    int nth = clamp_threads(nthreads);

    /* one set of working/row buffers per thread */
    oneway_ws *ws = R_Calloc(nth, oneway_ws);
    for (int t = 0; t < nth; ++t)
    {
        ws[t].ysptr = R_Calloc(nlibs, double);
        ws[t].osptr = R_Calloc(nlibs, double);
        ws[t].dsptr = R_Calloc(nlibs, double);
        ws[t].wsptr = R_Calloc(nlibs, double);
        ws[t].yptr = R_Calloc(nlibs, double);
        ws[t].optr = R_Calloc(nlibs, double);
        ws[t].dptr = R_Calloc(nlibs, double);
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

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag = 0; tag < ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        oneway_ws *w = &ws[tid];
        get_row4(y, offsets, disp, weights, tag, w->yptr, w->optr, w->dptr, w->wptr);

        for (int g = 0; g < ngroups; ++g)
        {
            int nlibs_g = gnlib[g];
            int s_idx = ogptr[g];

            for (int k = 0; k < nlibs_g; ++k)
            {
                int lib = g_idx[s_idx + k];
                w->ysptr[k] = w->yptr[lib];
                w->osptr[k] = w->optr[lib];
                w->dsptr[k] = w->dptr[lib];
                w->wsptr[k] = w->wptr[lib];
            }

            int b_idx = g * ntag + tag;

            double ocoef;
            int oconv;
            double cur_beta = coef_start ? coef_start[b_idx] : NA_REAL;

            glm_one_group_vec(nlibs_g, w->ysptr, w->osptr, w->dsptr, w->wsptr, maxit, tol,
                              cur_beta, &ocoef, &oconv);

            coef[b_idx] = ocoef;
            conv[b_idx] = oconv;
        }
    }

    R_Free(gnlib);
    R_Free(g_idx);
    R_Free(ogptr);

    for (int t = 0; t < nth; ++t)
    {
        R_Free(ws[t].ysptr);
        R_Free(ws[t].osptr);
        R_Free(ws[t].dsptr);
        R_Free(ws[t].wsptr);
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].dptr);
        R_Free(ws[t].wptr);
    }
    R_Free(ws);
}

/* Detect the unique groups (distinct design rows) of a design matrix.
 * Two libraries share a group when all their design-matrix entries match to
 * within 1e-12. Groups are numbered 0,1,... in order of first appearance.
 *
 * inputs:
 * design  design matrix (nlib x ncoef)
 *
 * output:
 * group   per-library group index (length nlib), filled in place
 *
 * returns: the number of distinct groups found.
 */
int get_groups_from_design(cmx *design, int *group)
{
    int nlibs = design->nrow;
    int ncoef = design->ncol;
    double *dmat = design->dmat;
    int ngroups = 0;

    for (int i = 0; i < nlibs; ++i)
    {
        group[i] = -1;
    }

    for (int i = 0; i < nlibs; ++i)
    {
        if (group[i] != -1)
        {
            continue;
        }

        int g = ngroups++;
        group[i] = g;

        for (int j = i + 1; j < nlibs; ++j)
        {
            if (group[j] != -1)
            {
                continue;
            }

            int match = 1;
            for (int c = 0; c < ncoef; ++c)
            {
                if (fabs(dmat[c * nlibs + i] - dmat[c * nlibs + j]) > 1e-12)
                {
                    match = 0;
                    break;
                }
            }
            if (match)
            {
                group[j] = g;
            }
        }
    }

    return ngroups;
}

/* one-way-layout branch of fit_glm_mat: fit group means with fit_one_way_mat, then map
 * the group coefficients back to design-scale coefficients via a dgesv solve. The caller
 * owns `group`; this returns 1 on a singular-design dgesv failure (having freed only its
 * own scratch) so the caller frees `group` and raises the error, and 0 on success. */
static int fit_glm_oneway(cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design, int *group, int ngroups, int ntag, int nlibs, int ncoef, int maxit, double tol, double *coef_start, double *coef, double *mu, double *dev, int *iter, int *failed, int nthreads)
{
    // 1. Construct designunique (size: ncoef x ncoef)
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
            designunique[c * ncoef + g] =
                design->dmat[c * nlibs + first_lib_of_group[g]];
        }
    }
    R_Free(first_lib_of_group);

    // 2. If starting coefficients are provided, transform them: coef.start %*%
    // t(designunique)
    double *coef_start_mat = NULL;
    if (coef_start != NULL)
    {
        coef_start_mat = R_Calloc((size_t) ntag * ncoef, double);
        /* coef_start[c*ntag+tag] / coef_start_mat[g*ntag+tag] walked by ntag;
           designunique[c*ncoef+g] walked by ncoef -- no int product overflows. */
        for (int tag = 0; tag < ntag; ++tag)
        {
            double *csmp = coef_start_mat + tag;
            for (int g = 0; g < ncoef; ++g, csmp += ntag)
            {
                double val = 0.0;
                const double *csp = coef_start + tag;
                const double *dup = designunique + g;
                for (int c = 0; c < ncoef; ++c, csp += ntag, dup += ncoef)
                {
                    val += (*csp) * (*dup);
                }
                *csmp = val;
            }
        }
    }

    // 3. Fit group-level coefficients
    int *conv = R_Calloc((size_t) ntag * ngroups, int);
    fit_one_way_mat(y, offsets, disp, weights, group, ngroups, maxit, tol,
                    coef_start_mat, coef, conv, nthreads);
    /* conv[g*ntag+tag] walked by ntag over groups */
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
    }
    R_Free(conv);
    if (coef_start_mat)
    {
        R_Free(coef_start_mat);
    }

    // 4. Compute fitted values (mu)
    cmx coef_cmx = make_cmx(coef, NULL, ntag, ncoef, 0, 0);
    get_one_way_fit(&coef_cmx, offsets, group, mu, nthreads);

    // 5. Compute deviances
    cmx mu_cmx = make_cmx(mu, NULL, ntag, nlibs, 0, 0);
    compute_nbdev_sum(y, &mu_cmx, disp, weights, dev);

    // 6. Map group coefficients to design coefficients: beta =
    // t(solve(designunique, t(beta))) Solve: designunique * X = t(beta)
    double *rhs = R_Calloc((size_t) ncoef * ntag, double);
    /* rhs[tag*ncoef+g] = coef[g*ntag+tag]: rhs base (size_t)tag*ncoef then +1;
       coef walked by ntag over coefs. */
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
    /* dgesv (LAPACK) is the double general linear solver: it LU-factorizes the
     * square matrix designunique and overwrites rhs with the solution X of
     * designunique * X = rhs. edgeR calls it here to map each gene's group-level
     * fitted coefficients back to design-scale coefficients for a one-way layout.
     * Argument mapping (Fortran = local):
     *   N    = ncoef        -> order of the system (number of coefficients/groups)
     *   NRHS = ntag         -> one gene per column of rhs
     *   A    = designunique -> N x N matrix, overwritten by its LU factors
     *   LDA  = ncoef        -> leading dimension of A
     *   IPIV = ipiv         -> row pivots from the LU factorization
     *   B    = rhs          -> in: group coefficients; out: solved design coefficients
     *   LDB  = ncoef        -> leading dimension of B
     *   INFO = info         -> 0 on success; nonzero means a singular design
     * Result: rhs holds the genewise coefficients on the design (log-mean) scale.
     * Equivalent R operation: beta <- t(solve(designunique, t(beta))).
     * Netlib references: https://netlib.org/lapack/explore-html/ (dgesv)
     */
    F77_CALL(dgesv)(&ncoef, &ntag, designunique, &ncoef, ipiv, rhs, &ncoef, &info);
    if (info != 0)
    {
        R_Free(ipiv);
        R_Free(designunique);
        R_Free(rhs);
        return 1;
    }
    R_Free(ipiv);
    R_Free(designunique);

    // Write solved coefficients back to coef
    /* coef[g*ntag+tag] = rhs[tag*ncoef+g] (inverse of the pack above) */
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

    // 7. Fill iteration indicators
    for (int i = 0; i < ntag; ++i)
    {
        iter[i] = 0;
    }
    return 0;
}

/* unified fitting function for negative binomial generalized linear models.
 * Automatically detects if the design matrix corresponds to a one-way layout.
 * If yes, fits using fast group-wise scoring and solves back using LAPACK dgesv.
 * Otherwise, falls back to the general Levenberg-Marquardt solver.
 * nthreads is forwarded to the genewise sub-kernels.
 */
void fit_glm_mat(cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design,
                 int maxit, double tol, double *coef_start, double *coef,
                 double *mu, double *dev, int *iter, int *failed, int *method,
                 int nthreads)
{
    int ntag = y->nrow;
    int nlibs = y->ncol;
    int ncoef = design->ncol;

    // Detect unique groups from design matrix
    int *group = R_Calloc(nlibs, int);
    int ngroups = get_groups_from_design(design, group);

    if (ngroups == ncoef)
    {
        int oneway_fail = fit_glm_oneway(y, offsets, disp, weights, design, group, ngroups, ntag, nlibs, ncoef, maxit, tol, coef_start, coef, mu, dev, iter, failed, nthreads);
        if (oneway_fail)
        {
            R_Free(group);
            error("LAPACK dgesv failed to solve the design matrix transformation");
        }
        *method = 0; // "oneway"
    }
    else
    {
        // --- General Design Layout Case ---
        // 1. Compute start values if null
        if (coef_start == NULL)
        {
            get_leven_start(y, offsets, disp, weights, design, 1, coef, nthreads);
        }
        else
        {
            memcpy(coef, coef_start, sizeof(double) * ntag * ncoef);
        }

        // 2. Run Levenberg-Marquardt fitting
        cmx beta_cmx = make_cmx(coef, NULL, ntag, ncoef, 0, 0);
        fit_leven(y, offsets, disp, weights, design, &beta_cmx, tol, maxit, mu,
                  coef, dev, iter, failed, nthreads);
        *method = 1; // "levenberg"
    }

    R_Free(group);
}
