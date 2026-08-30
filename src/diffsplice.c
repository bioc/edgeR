#include "edgeR.h"

/* -------------------------------------------------------------------------
 * diffSplice negative-binomial GLM solver.
 *
 * For a gene with n exons, the null model constrains all exons to share the
 * tested coefficient ("betabar") while each exon keeps its own nuisance
 * coefficients. The design is therefore highly structured:
 *
 *   - column 0  : the tested covariate, replicated across all exons
 *                 -> the single shared coefficient betabar
 *   - the rest  : block-diagonal, exon e's (ncoef-1) nuisance covariates D
 *                 act only on exon e's rows -> coefficients beta_e
 *
 * Hence XtWX is a bordered block-diagonal ("arrowhead") matrix
 *
 *     [ a     b_1^T ... b_n^T ]   a   : scalar     = sum_{e,s} w * x_coef^2
 *     [ b_1   A_1        0    ]   b_e : (ncoef-1)  = sum_s    w * x_coef * D
 *     [  :          .         ]   A_e : (ncoef-1)^2= sum_s    w * D^T D   (SPD)
 *     [ b_n    0     ...  A_n ]   off-diagonal exon blocks are exactly zero
 *
 * which is solved by a Schur complement on the scalar border in
 * O(n*(ncoef-1)^3) instead of dense Cholesky's O((n*(ncoef-1))^3). The dense
 * n*nsamples by (1+n*(ncoef-1)) design is never formed. Everything else
 * (working weights, deviance, Levenberg damping, convergence) mirrors the
 * general kernel fit_leven_vec in glm.c.
 * ------------------------------------------------------------------------- */

/* Per-thread scratch, all buffers sized for the full model (N = max_nexon *
 * nsamples, P = 1 + max_nexon * q). The reduced (leave-one-out) fits reuse the
 * same buffers with a smaller leading sub-problem, so nothing is allocated
 * inside the parallel region. */
typedef struct {
    /* flattened gene data (full and reduced) */
    double *y_vec, *off_vec, *disp_vec, *w_vec;   /* N           */
    double *fitted_mu;                            /* N           */
    double *start_vec;                            /* P           */
    double *y_red, *off_red, *disp_red, *w_red;   /* N           */
    double *fitted_mu_red;                        /* N           */
    double *start_red;                            /* P           */
    /* fit_splice_leven_vec working space */
    double *zwpt, *drvt, *nmu;                    /* N           */
    double *dl, *db, *nbt;                        /* P           */
    double *b, *u, *v;                            /* n*q         */
    double *A, *Ad;                               /* n*q*q       */
} splice_ws;

/* Fill the fitted means of the splice null model from the current coefficients.
 *
 * mu = exp(eta), eta_{e,s} = betabar*x_coef[s] + D[s,]*beta_e + offset
 *
 * inputs:
 *   n        number of exons in the gene
 *   nsamples number of samples (columns) per exon
 *   q        nuisance dimension per exon (ncoef - 1)
 *   x_coef   tested covariate column, length nsamples
 *   Dcols    q pointers to the nuisance design columns, each length nsamples
 *   beta     coefficient vector, length 1 + n*q; beta[0] is the shared betabar
 *   offset   offsets, length n*nsamples (exon-major)
 *
 * output:
 *   mu       fitted means, length n*nsamples (exon-major), overwritten
 *
 * No return value; writes only into mu.
 */
static void splice_autofill(int n, int nsamples, int q, const double *x_coef, double *const *Dcols, const double *beta, const double *offset, double *mu)
{
    double betabar = beta[0];
    for (int e = 0; e < n; ++e)
    {
        const double *be = beta + 1 + e * q;
        int base = e * nsamples;
        for (int s = 0; s < nsamples; ++s)
        {
            double eta = offset[base + s] + betabar * x_coef[s];
            for (int j = 0; j < q; ++j)
            {
                eta += be[j] * Dcols[j][s];
            }
            mu[base + s] = exp(eta);
        }
    }
}

/* Solve (XtWX + lambda I) db = dl using the arrowhead structure.
 * Returns 1 if the damped system is not positive definite (caller increases
 * lambda and retries, exactly like the dpotrf retry in fit_leven_vec), else 0
 * with the Newton step written to db. */
static int splice_solve(int n, int q, double a, const double *b, const double *A, const double *dl, double lambda, double *Ad, double *u, double *v, double *db)
{
    const char uplo = 'U';
    const int one = 1;
    double schur = a + lambda;   /* Schur complement of the block-diagonal part */
    double rhs0  = dl[0];

    for (int e = 0; e < n; ++e)
    {
        const double *be = b + e * q;
        const double *Ae = A + e * (size_t) q * q;
        const double *ge = dl + 1 + e * q;
        double *Ade = Ad + e * (size_t) q * q;
        double *ue = u + e * q;
        double *ve = v + e * q;

        /* Ade = A_e + lambda I */
        for (int idx = 0; idx < q * q; ++idx)
        {
            Ade[idx] = Ae[idx];
        }
        for (int j = 0; j < q; ++j)
        {
            Ade[j * q + j] += lambda;
        }

        if (q == 1)
        {
            if (Ade[0] <= 0)
            {
                return 1;
            }
            ue[0] = be[0] / Ade[0];
            ve[0] = ge[0] / Ade[0];
        }
        else
        {
            int info = 0, qq = q;
            /* dpotrf is the LAPACK double symmetric-positive-definite Cholesky
             * factorization: it factors the upper triangle (UPLO='U') of the
             * damped nuisance block Ade = A_e + lambda I as A_e = U^T U and
             * overwrites Ade with the factor U. The two dpotrs calls are the
             * matching LAPACK triangular solves; each reuses the Cholesky factor
             * to solve (A_e + lambda I) X = B for a single right-hand side,
             * yielding A_e^{-1} b_e (in ue) and A_e^{-1} g_e (in ve).
             *
             * edgeR calls these to invert each exon's nuisance block A_e of the
             * arrowhead XtWX in the diffSplice GLM/F-test covariance solve; the
             * resulting products form the Schur complement on the shared tested
             * coefficient betabar.
             *
             * Argument mapping:
             *   dpotrf: UPLO = &uplo ('U', upper triangle), N = &qq (order q of
             *     the block), A = Ade (q x q, overwritten by its Cholesky factor
             *     U), LDA = &qq (leading dimension), INFO = &info (0 on success,
             *     >0 means the block is not positive definite).
             *   dpotrs (called twice): UPLO = &uplo ('U'), N = &qq (order q),
             *     NRHS = &one (one column), A = Ade (the factor U from dpotrf),
             *     LDA = &qq, B = ue then ve (right-hand side on entry b_e / g_e,
             *     solution on exit), LDB = &qq, INFO = &info.
             *   FCONE passes the hidden Fortran length of the UPLO character.
             *
             * A non-zero info makes splice_solve return 1 so the caller damps
             * lambda harder and retries, mirroring the dpotrf retry in
             * fit_leven_vec.
             *
             * Equivalent R operation:
             *   R  <- chol(A_e + lambda * diag(q))
             *   ue <- backsolve(R, backsolve(R, b_e, transpose = TRUE))
             *   ve <- backsolve(R, backsolve(R, g_e, transpose = TRUE))
             *   i.e. solve(A_e + lambda * diag(q), cbind(b_e, g_e)).
             *
             * Netlib references: https://netlib.org/lapack/explore-html/
             *   (dpotrf, dpotrs). */
            F77_CALL(dpotrf)(&uplo, &qq, Ade, &qq, &info FCONE);
            if (info != 0)
            {
                return 1;
            }
            for (int j = 0; j < q; ++j)
            {
                ue[j] = be[j];
                ve[j] = ge[j];
            }
            F77_CALL(dpotrs)(&uplo, &qq, &one, Ade, &qq, ue, &qq, &info FCONE);
            if (info != 0)
            {
                return 1;
            }
            F77_CALL(dpotrs)(&uplo, &qq, &one, Ade, &qq, ve, &qq, &info FCONE);
            if (info != 0)
            {
                return 1;
            }
        }

        double bu = 0, bv = 0;
        for (int j = 0; j < q; ++j)
        {
            bu += be[j] * ue[j];
            bv += be[j] * ve[j];
        }
        schur -= bu;    /* a + lambda - sum b_e^T A_e^{-1} b_e */
        rhs0  -= bv;    /* g0        - sum b_e^T A_e^{-1} g_e */
    }

    if (schur <= 0)
    {
        return 1;   /* indefinite from round-off: damp harder */
    }

    double d0 = rhs0 / schur;
    db[0] = d0;
    for (int e = 0; e < n; ++e)
    {
        const double *ue = u + e * q;
        const double *ve = v + e * q;
        double *de = db + 1 + e * q;
        for (int j = 0; j < q; ++j)
        {
            de[j] = ve[j] - d0 * ue[j];  /* A_e^{-1}(g_e - b_e d0) */
        }
    }
    return 0;
}

/* Levenberg-damped IRLS for one gene's splice null model, exploiting the
 * arrowhead structure. Mirrors fit_leven_vec's iteration logic. obt is the
 * starting coefficient vector (length P) and receives the fit; omu (length N)
 * receives the fitted means; *odev the deviance.
 *
 * *oerr is reserved for an unrecoverable failure but ALWAYS stays 0: the damped
 * solve self-recovers -- splice_solve() below is retried with lambda increased
 * 10x until the damped system is positive definite (mirroring the dpotrf retry
 * in fit_leven_vec), so no LAPACK error can escape. The caller's oerr-gated
 * dfail/error() path is therefore dead code (kept as commented reference there).
 *
 * The separate local `failed` (Levenberg damping exhausted = non-convergence,
 * set at the lambda/max_info ceiling below) is intentionally NOT fatal and NOT
 * reported: the loop breaks and the best-so-far coefficients/deviance are
 * returned. This differs from the NB kernel fit_leven, which flags the same
 * non-convergence to R via fit$failed; diffSplice has no per-gene status output,
 * so such a (rare) gene yields an approximate LR silently. */
static void fit_splice_leven_vec(splice_ws *ws, int n, int nsamples, int q, const double *x_coef, double *const *Dcols, const double *y, const double *offset, const double *disp, const double *wt, int maxit, double tol, double *obt, double *omu, double *odev, int *oerr)
{
    const double low_value = 1e-10;
    const double one_millionth = 1e-6;
    const double supremely_low_value = 1e-13;
    const double ridiculously_low_value = 1e-100;

    int N = n * nsamples, P = 1 + n * q;
    *oerr = 0;

    double ymax = 0;
    for (int i = 0; i < N; ++i)
    {
        ymax = (y[i] > ymax) ? y[i] : ymax;
    }
    if (ymax < low_value)
    {
        for (int c = 0; c < P; ++c)
        {
            obt[c] = NA_REAL;
        }
        for (int i = 0; i < N; ++i)
        {
            omu[i] = 0;
        }
        *odev = 0;
        return;
    }

    splice_autofill(n, nsamples, q, x_coef, Dcols, obt, offset, omu);
    double dev = nb_deviance_sum(N, y, omu, disp, wt);

    double *zwpt = ws->zwpt, *drvt = ws->drvt, *nmu = ws->nmu;
    double *dl = ws->dl, *db = ws->db, *nbt = ws->nbt;
    double *b = ws->b, *A = ws->A, *Ad = ws->Ad, *u = ws->u, *v = ws->v;

    double max_info = -1, lambda = 0;
    int iter = 0, failed = 0;

    while ((++iter) <= maxit)
    {
        /* working weights and gradient pieces */
        for (int i = 0; i < N; ++i)
        {
            double cur_mu = omu[i], denom = 1 + cur_mu * disp[i];
            zwpt[i] = cur_mu / denom * wt[i];
            drvt[i] = (y[i] - cur_mu) / denom * wt[i];
        }

        /* structured assembly of a, b_e, A_e and the gradient dl = X^T drvt */
        double a = 0;
        dl[0] = 0;
        for (int e = 0; e < n; ++e)
        {
            double *be = b + e * q, *Ae = A + e * (size_t) q * q, *ge = dl + 1 + e * q;
            for (int j = 0; j < q; ++j)
            {
                be[j] = 0;
                ge[j] = 0;
            }
            for (int idx = 0; idx < q * q; ++idx)
            {
                Ae[idx] = 0;
            }
            int base = e * nsamples;
            for (int s = 0; s < nsamples; ++s)
            {
                double zw = zwpt[base + s], dr = drvt[base + s], xc = x_coef[s];
                a += zw * xc * xc;
                dl[0] += dr * xc;
                for (int j = 0; j < q; ++j)
                {
                    double dj = Dcols[j][s];
                    be[j] += zw * xc * dj;
                    ge[j] += dr * dj;
                    for (int k = 0; k < q; ++k)
                    {
                        Ae[j * q + k] += zw * dj * Dcols[k][s];
                    }
                }
            }
        }
        if (a > max_info)
        {
            max_info = a;
        }
        for (int e = 0; e < n; ++e)
        {
            double *Ae = A + e * (size_t) q * q;
            for (int j = 0; j < q; ++j)
            {
                if (Ae[j * q + j] > max_info)
                {
                    max_info = Ae[j * q + j];
                }
            }
        }

        if (iter == 1)
        {
            lambda = max_info * one_millionth;
            if (lambda < supremely_low_value)
            {
                lambda = supremely_low_value;
            }
        }

        int lev = 0, low_dev = 0;
        while (++lev)
        {
            /* damped solve, increasing lambda until positive definite */
            while (splice_solve(n, q, a, b, A, dl, lambda, Ad, u, v, db))
            {
                lambda *= 10;
                if (lambda <= 0)
                {
                    lambda = ridiculously_low_value;
                }
            }

            for (int c = 0; c < P; ++c)
            {
                nbt[c] = obt[c] + db[c];
            }
            splice_autofill(n, nsamples, q, x_coef, Dcols, nbt, offset, nmu);

            double ndev = nb_deviance_sum(N, y, nmu, disp, wt);

            if (ndev / ymax < supremely_low_value)
            {
                low_dev = 1;
            }
            if (ndev <= dev || low_dev)
            {
                for (int c = 0; c < P; ++c)
                {
                    obt[c] = nbt[c];
                }
                for (int i = 0; i < N; ++i)
                {
                    omu[i] = nmu[i];
                }
                dev = ndev;
                break;
            }

            lambda *= 2;
            if (lambda <= 0)
            {
                lambda = ridiculously_low_value;
            }
            /* damping exhausted: non-convergence. Deliberately not fatal and not
               propagated to *oerr -- break and return the best-so-far fit (see the
               header note); the caller's oerr-gated error() path stays dead. */
            if (lambda / max_info > 1 / supremely_low_value)
            {
                failed = 1;
                break;
            }
        }

        double divergence = 0;
        for (int c = 0; c < P; ++c)
        {
            divergence += dl[c] * db[c];
        }
        if (failed || low_dev || (divergence < tol))
        {
            break;
        }

        if (lev == 1)
        {
            lambda /= 10;
        }
    }

    *odev = dev;
}

/* Core driver: fit the splice null models for every gene.
 *
 * For every gene the shared-coefficient null model is fitted by
 * fit_splice_leven_vec, giving the gene-level likelihood-ratio statistic
 * against the alternative fit. Exon-level statistics are then produced by one
 * of three cases depending on the number of exons n:
 *   Case A (n > nexons_approx): approximate each exon from the full-null means.
 *   Case B (2 < n <= nexons_approx): exact leave-one-out reduced fits.
 *   Case C (n == 2): the closed-form two-exon result.
 * Genes are processed in one flattened OpenMP loop; each writes disjoint
 * output positions.
 *
 * inputs:
 *   y             exon counts, compressed matrix (nexons x nsamples)
 *   offsets       compressed matrix of offsets
 *   disp          compressed matrix of dispersions
 *   weights       compressed matrix of prior weights
 *   design        design matrix, dense (nsamples x ncoef)
 *   ngenes        number of genes
 *   nsamples      number of samples (design rows)
 *   ncoef         number of design coefficients
 *   coef          index of the tested design column
 *   gene_nexons   number of exons per gene, length ngenes
 *   gene_firstexon first exon (row) index of each gene, length ngenes
 *   nexons_approx exon-count threshold above which Case A is used
 *   maxit         maximum IRLS iterations per fit
 *   tol           IRLS convergence tolerance
 *   coef_start    starting coefficients from the alternative fit (ncoef x nexons)
 *   exon_dev      per-exon alternative-model deviance, length nexons
 *   nthreads      requested thread count (clamped by clamp_threads)
 *
 * outputs (caller-allocated):
 *   exon_LR   per-exon likelihood-ratio statistic, length nexons
 *   exon_coef per-exon tested-coefficient contrast, length nexons
 *   gene_LR   per-gene likelihood-ratio statistic, length ngenes
 *
 * No return value; writes only into the caller-owned output arrays. Allocates
 * per-thread scratch that is freed before returning. If any thread's Cholesky
 * solve fails, a flag is set and, after the parallel region and after all
 * scratch is freed, error() is called.
 */
void fit_diff_splice_mat(cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design, int ngenes, int nsamples, int ncoef, int coef, int *gene_nexons, int *gene_firstexon, int nexons_approx, int maxit, double tol, double *coef_start, double *exon_dev, double *exon_LR, double *exon_coef, double *gene_LR, int nthreads)
{
    int nexons = y->nrow;
    int q = ncoef - 1;                 /* nuisance dimension per exon */

    int max_nexon = 0;
    for (int g = 0; g < ngenes; ++g)
    {
        if (gene_nexons[g] > max_nexon)
        {
            max_nexon = gene_nexons[g];
        }
    }
    if (max_nexon < 2)
    {
        return;         /* nothing testable */
    }

    int nth = clamp_threads(nthreads);

    /* shared read-only mapping from nuisance index to original design column,
     * plus pointers to the tested column and the nuisance columns */
    int *col_map = R_Calloc(q, int);
    double **Dcols = R_Calloc(q, double*);
    for (int j = 0; j < q; ++j)
    {
        col_map[j] = (j < coef) ? j : (j + 1);
        Dcols[j] = design->dmat + (size_t) col_map[j] * nsamples;
    }
    double *x_coef = design->dmat + (size_t) coef * nsamples;

    int Nmax = nsamples * max_nexon;
    int Pmax = 1 + max_nexon * q;

    /* one scratch workspace per thread, sized for the largest gene */
    splice_ws *ws = R_Calloc(nth, splice_ws);
    for (int t = 0; t < nth; ++t)
    {
        ws[t].y_vec = R_Calloc(Nmax, double);
        ws[t].off_vec = R_Calloc(Nmax, double);
        ws[t].disp_vec = R_Calloc(Nmax, double);
        ws[t].w_vec = R_Calloc(Nmax, double);
        ws[t].fitted_mu = R_Calloc(Nmax, double);
        ws[t].start_vec = R_Calloc(Pmax, double);
        ws[t].y_red = R_Calloc(Nmax, double);
        ws[t].off_red = R_Calloc(Nmax, double);
        ws[t].disp_red = R_Calloc(Nmax, double);
        ws[t].w_red = R_Calloc(Nmax, double);
        ws[t].fitted_mu_red = R_Calloc(Nmax, double);
        ws[t].start_red = R_Calloc(Pmax, double);
        ws[t].zwpt = R_Calloc(Nmax, double);
        ws[t].drvt = R_Calloc(Nmax, double);
        ws[t].nmu = R_Calloc(Nmax, double);
        ws[t].dl = R_Calloc(Pmax, double);
        ws[t].db = R_Calloc(Pmax, double);
        ws[t].nbt = R_Calloc(Pmax, double);
        ws[t].b = R_Calloc((size_t) max_nexon * q, double);
        ws[t].u = R_Calloc((size_t) max_nexon * q, double);
        ws[t].v = R_Calloc((size_t) max_nexon * q, double);
        ws[t].A = R_Calloc((size_t) max_nexon * q * q, double);
        ws[t].Ad = R_Calloc((size_t) max_nexon * q * q, double);
    }

    /* --- Dead failure-report scaffold (retained as commented reference) ---
       *oerr from fit_splice_leven_vec is never set to 1 (splice_solve self-
       recovers by increasing lambda), so dfail is never raised, splice_failed
       stays 0, and the error() at the end can never fire. Kept commented to
       mirror the live oerr/error path in fit_leven_vec (glm.c).
    int *dfail = R_Calloc(nth, int);
    */

    /* One flattened, load-balanced loop over all genes (each is independent and
     * writes disjoint output positions). */
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(dynamic)
    #endif
    for (int g = 0; g < ngenes; ++g)
    {
        int n = gene_nexons[g];
        if (n < 2)
        {
            continue;
        }

        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        splice_ws *w = &ws[tid];
        int first_exon = gene_firstexon[g];

        double gene_alt_dev = 0.0;
        for (int e = 0; e < n; ++e)
        {
            gene_alt_dev += exon_dev[first_exon + e];
        }

        /* flatten this gene's data (exon-major: index e*nsamples + s) */
        for (int e = 0; e < n; ++e)
        {
            int exon_idx = first_exon + e;
            get_row4(y, offsets, disp, weights, exon_idx, w->y_vec + e * nsamples, w->off_vec + e * nsamples, w->disp_vec + e * nsamples, w->w_vec + e * nsamples);
        }

        /* starting values for the full null model */
        double beta1 = 0.0;
        for (int e = 0; e < n; ++e)
        {
            beta1 += coef_start[(size_t) coef * nexons + (first_exon + e)];
        }
        beta1 /= n;
        w->start_vec[0] = beta1;
        for (int e = 0; e < n; ++e)
        {
            for (int j = 0; j < q; ++j)
            {
                w->start_vec[1 + e * q + j] = coef_start[(size_t) col_map[j] * nexons + (first_exon + e)];
            }
        }

        double null_dev;
        int oerr = 0;
        fit_splice_leven_vec(w, n, nsamples, q, x_coef, Dcols, w->y_vec, w->off_vec, w->disp_vec, w->w_vec, maxit, tol, w->start_vec, w->fitted_mu, &null_dev, &oerr);
        /* oerr is always 0 (see fit_splice_leven_vec); dead report path:
        if (oerr)
        {
            dfail[tid] = 1;
        }
        */

        gene_LR[g] = null_dev - gene_alt_dev;
        double betabar = w->start_vec[0];

        if (n > nexons_approx)
        {
            /* Case A: approximation from the full-null fitted means */
            for (int e = 0; e < n; ++e)
            {
                int exon_idx = first_exon + e;
                double exon_null_dev = 0.0;
                for (int s = 0; s < nsamples; ++s)
                {
                    double mu_val = w->fitted_mu[e * nsamples + s];
                    exon_null_dev += w->w_vec[e * nsamples + s] * compute_unit_nb_deviance(w->y_vec[e * nsamples + s], mu_val, w->disp_vec[e * nsamples + s]);
                }
                exon_LR[exon_idx] = exon_null_dev - exon_dev[exon_idx];
                exon_coef[exon_idx] = coef_start[(size_t) coef * nexons + exon_idx] - betabar;
            }
        }
        else if (n > 2)
        {
            /* Case B: exact leave-one-out, warm-started from the full-null fit */
            int nr = n - 1;
            for (int k = 0; k < n; ++k)
            {
                int pos = 0;
                for (int e = 0; e < n; ++e)
                {
                    if (e == k)
                    {
                        continue;
                    }
                    for (int s = 0; s < nsamples; ++s)
                    {
                        w->y_red[pos * nsamples + s] = w->y_vec[e * nsamples + s];
                        w->off_red[pos * nsamples + s] = w->off_vec[e * nsamples + s];
                        w->disp_red[pos * nsamples + s] = w->disp_vec[e * nsamples + s];
                        w->w_red[pos * nsamples + s] = w->w_vec[e * nsamples + s];
                    }
                    ++pos;
                }

                /* starting values for the reduced null model */
                double beta1_red = 0.0;
                for (int e = 0; e < n; ++e)
                {
                    if (e == k)
                    {
                        continue;
                    }
                    beta1_red += coef_start[(size_t) coef * nexons + (first_exon + e)];
                }
                beta1_red /= (n - 1);
                w->start_red[0] = beta1_red;
                pos = 0;
                for (int e = 0; e < n; ++e)
                {
                    if (e == k)
                    {
                        continue;
                    }
                    for (int j = 0; j < q; ++j)
                    {
                        w->start_red[1 + pos * q + j] = coef_start[(size_t) col_map[j] * nexons + (first_exon + e)];
                    }
                    ++pos;
                }

                double null_dev_red;
                int oerr_red = 0;
                fit_splice_leven_vec(w, nr, nsamples, q, x_coef, Dcols, w->y_red, w->off_red, w->disp_red, w->w_red, maxit, tol, w->start_red, w->fitted_mu_red, &null_dev_red, &oerr_red);
                /* oerr_red is always 0 (see fit_splice_leven_vec); dead report path:
                if (oerr_red)
                {
                    dfail[tid] = 1;
                }
                */

                int exon_idx = first_exon + k;
                exon_LR[exon_idx] = null_dev - (exon_dev[exon_idx] + null_dev_red);
                exon_coef[exon_idx] = coef_start[(size_t) coef * nexons + exon_idx] - w->start_red[0];
            }
        }
        else
        {
            /* Case C: exactly two exons (closed form) */
            exon_LR[first_exon] = gene_LR[g];
            exon_LR[first_exon + 1] = gene_LR[g];
            exon_coef[first_exon] = coef_start[(size_t) coef * nexons + first_exon] - coef_start[(size_t) coef * nexons + first_exon + 1];
            exon_coef[first_exon + 1] = coef_start[(size_t) coef * nexons + first_exon + 1] - coef_start[(size_t) coef * nexons + first_exon];
        }
    }

    /* dead: dfail is never raised (see the scaffold note above)
    int splice_failed = 0;
    for (int t = 0; t < nth; ++t)
    {
        if (dfail[t])
        {
            splice_failed = 1;
        }
    }
    */

    for (int t = 0; t < nth; ++t)
    {
        R_Free(ws[t].y_vec);
        R_Free(ws[t].off_vec);
        R_Free(ws[t].disp_vec);
        R_Free(ws[t].w_vec);
        R_Free(ws[t].fitted_mu);
        R_Free(ws[t].start_vec);
        R_Free(ws[t].y_red);
        R_Free(ws[t].off_red);
        R_Free(ws[t].disp_red);
        R_Free(ws[t].w_red);
        R_Free(ws[t].fitted_mu_red);
        R_Free(ws[t].start_red);
        R_Free(ws[t].zwpt);
        R_Free(ws[t].drvt);
        R_Free(ws[t].nmu);
        R_Free(ws[t].dl);
        R_Free(ws[t].db);
        R_Free(ws[t].nbt);
        R_Free(ws[t].b);
        R_Free(ws[t].u);
        R_Free(ws[t].v);
        R_Free(ws[t].A);
        R_Free(ws[t].Ad);
    }
    R_Free(ws);
    /* R_Free(dfail); -- dfail alloc commented out above */
    R_Free(Dcols);
    R_Free(col_map);

    /* dead: splice_failed can never be set (see the scaffold note above)
    if (splice_failed)
    {
        error("solution using the Cholesky decomposition failed");
    }
    */
}
