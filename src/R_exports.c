#include "edgeR.h"


/* -----------------------------------------------------------------------
 * from add_prior_count.c
 * Adding a prior count to each observation. 
 * wrapper function for adding prior count 
 */

SEXP add_prior_count(SEXP y, SEXP offsets, SEXP priors, SEXP nthreads)
{
    SEXP yy, offset, ans;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx pmx = SEXPtocmx2(priors);

    int nthr = asInteger(nthreads);

    PROTECT(yy = coerceVector(duplicate(y), REALSXP));

    /* prepare the offset output
     * if both offsets and priors are row repeated
     * the adjusted offset will be row repeated
     * and the return a row vector
     * else return a matrix
     */
    int k = (omx.type>=2 && pmx.type>=2)? 1 : ymx.nrow;
    PROTECT(offset  = allocMatrix(REALSXP,k,ymx.ncol));

    if(omx.type>=2 && pmx.type>=2)
    {
        add_prior_count_vec(&ymx,&omx,&pmx,REAL(yy),REAL(offset),nthr);
    }
    else
    {
        add_prior_count_mat(&ymx,&omx,&pmx,REAL(yy),REAL(offset),nthr);
    }

    const char *names[] = {"y","offset",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, yy);
    SET_VECTOR_ELT(ans, 1, offset);

    UNPROTECT(3);
    return ans;
}

/* -----------------------------------------------------------------------
 * from compute_apl.c
 * wrapper function for computing adjusted profile likelihood
 */
SEXP compute_apl(SEXP y, SEXP mu, SEXP disp, SEXP weights, SEXP adjust, SEXP design, SEXP nthreads)
{
    SEXP ans;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx umx = SEXPtocmx1(mu);
    cmx gmx = SEXPtocmx1(design);

    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int do_adjust = asLogical(adjust);
    int nthr = asInteger(nthreads);

    PROTECT(ans=allocVector(REALSXP, ymx.nrow));

    compute_adj_profile_ll(&ymx, &umx, &dmx, &wmx, &gmx, do_adjust, REAL(ans), nthr);

    UNPROTECT(2);

    return ans;
}

/* -----------------------------------------------------------------------
 * combined GLM fit + adjusted profile likelihood
 *
 * Fits the tagwise GLMs (fit_glm_mat) to obtain the fitted values 'mu' and the
 * coefficients, then feeds 'mu' straight into compute_adj_profile_ll without
 * round-tripping through R. Replaces the glmFit() + .cxx_compute_apl pair in
 * adjustedProfileLik(). Returns a named list: apl (per-tag) and coefficients.
 */
SEXP fit_apl (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP max_iterations, SEXP tolerance, SEXP coef_start, SEXP adjust, SEXP nthreads)
{
    SEXP ans, coef, mu, apl;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int do_adjust = asLogical(adjust);
    int nthr = asInteger(nthreads);

    int protect_count = 1; // design protected
    double *c_start = NULL;
    if (coef_start != R_NilValue)
    {
        PROTECT(coef_start = coerceVector(coef_start, REALSXP));
        c_start = REAL(coef_start);
        protect_count++;
    }

    int ntag = ymx.nrow;
    int ncoef = gmx.ncol;
    int nlibs = ymx.ncol;

    // Outputs we keep: coefficients (returned) and the fitted values 'mu'
    // (consumed below by the APL computation).
    PROTECT(coef = allocMatrix(REALSXP, ntag, ncoef));
    PROTECT(mu = allocMatrix(REALSXP, ntag, nlibs));
    protect_count += 2;

    // fit_glm_mat also reports deviance/iter/failed, which the APL does not
    // need; use throwaway C buffers rather than allocating R objects.
    double *dev = R_Calloc(ntag, double);
    int *iter = R_Calloc(ntag, int);
    int *failed = R_Calloc(ntag, int);
    int method_val = 0;

    fit_glm_mat(&ymx, &omx, &dmx, &wmx, &gmx, maxit, tol, c_start, REAL(coef), REAL(mu), dev, iter, failed, &method_val, nthr);

    R_Free(dev);
    R_Free(iter);
    R_Free(failed);

    // Compute the adjusted profile likelihood directly from the fitted 'mu'.
    cmx umx = SEXPtocmx1(mu);
    PROTECT(apl = allocVector(REALSXP, ntag));
    protect_count++;

    compute_adj_profile_ll(&ymx, &umx, &dmx, &wmx, &gmx, do_adjust, REAL(apl), nthr);

    const char *names[] = {"apl", "coefficients", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    protect_count++;

    SET_VECTOR_ELT(ans, 0, apl);
    SET_VECTOR_ELT(ans, 1, coef);

    UNPROTECT(protect_count);
    return ans;
}

/* -----------------------------------------------------------------------
 * from compute_cpm.c
 * wrapper function for computing CPM and logCPM 
 * and average log CPM 
 */

SEXP calculate_cpm_log (SEXP y, SEXP libsizes, SEXP priors, SEXP nthreads)
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx lmx = SEXPtocmx2(libsizes);
    cmx pmx = SEXPtocmx2(priors);

    int nthr = asInteger(nthreads);

    PROTECT(ans = coerceVector(duplicate(y), REALSXP));

    calc_cpm_log(&ymx, &lmx, &pmx, REAL(ans), nthr);

    UNPROTECT(1);

    return ans;
}

/* wrapper function for computing raw CPM
 * inputs: y count matrix, libsizes library sizes (compressMatrix), nthreads
 * returns: matrix of counts per million, same dimensions as y
 */
SEXP calculate_cpm_raw (SEXP y, SEXP libsizes, SEXP nthreads)
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx lmx = SEXPtocmx2(libsizes);

    int nthr = asInteger(nthreads);

    PROTECT(ans = coerceVector(duplicate(y), REALSXP));

    calc_cpm_raw(&ymx, &lmx, REAL(ans), nthr);

    UNPROTECT(1);

    return ans;
}

/* wrapper function for computing the average log2-CPM per tag
 * inputs: y count matrix, offsets/priors/disp/weights (compressMatrix),
 *   max_iterations, tolerance, nthreads
 * returns: numeric vector of average log-CPM, one value per row of y
 */
SEXP ave_log_cpm(SEXP y, SEXP offsets, SEXP priors, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP nthreads)
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx pmx = SEXPtocmx2(priors);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    PROTECT(ans=allocVector(REALSXP,ymx.nrow));

    average_log_cpm(&ymx,&omx,&pmx,&dmx,&wmx,maxit,tol,REAL(ans),nthr);

    UNPROTECT(1);
    return ans;
}

/* --------------------------------------------------------------------
 * from compute_nbdev.c
 * wrapper function for computing negative binomial deviance 
 */
SEXP compute_nbdev (SEXP y, SEXP mu, SEXP disp, SEXP weights, SEXP dosum)
{
    SEXP ans;

    PROTECT(y=coerceVector(y,REALSXP));
    PROTECT(mu=coerceVector(mu,REALSXP));

    cmx ymx = SEXPtocmx1(y);
    cmx umx = SEXPtocmx1(mu);
    cmx dmx = SEXPtocmx2(disp);

    int do_sum = asLogical(dosum);

    if(do_sum)
    {
        cmx wmx = SEXPtocmx2(weights);
        PROTECT(ans=allocVector(REALSXP,ymx.nrow));

        compute_nbdev_sum(&ymx, &umx, &dmx, &wmx, REAL(ans));

        UNPROTECT(3);
        return ans;
    }
    else
    {
        PROTECT(ans=duplicate(y));

        compute_nbdev_unit(&ymx, &umx, &dmx, REAL(ans));

        UNPROTECT(3);
        return ans;
    }
}

/* --------------------------------------------------------------------
 * from exact_test_by_dev.c
 * exact test by deviance 
 */
SEXP exact_test_by_deviance(SEXP sums_1, SEXP sums_2, SEXP n_1, SEXP n_2, SEXP disp)
{
    SEXP ans;

    int ntags=length(sums_1);

    int n1 = asInteger(n_1);
    int n2 = asInteger(n_2);

    PROTECT(disp=coerceVector(disp,REALSXP));
    PROTECT(ans=duplicate(disp));

    exact_test_by_dev(INTEGER(sums_1),INTEGER(sums_2),ntags,n1,n2,REAL(disp),REAL(ans));

    UNPROTECT(2);

    return ans;
}
/* --------------------------------------------------------------------
 * from glm.c
 * */

/* fit one group mean */

SEXP fit_one_group (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP beta, SEXP nthreads)
{
    SEXP ans, coef, conv;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    PROTECT(coef=allocVector(REALSXP,ymx.nrow));
    PROTECT(conv=allocVector(LGLSXP,ymx.nrow));

    fit_one_group_mat(&ymx,&omx,&dmx,&wmx,maxit,tol,REAL(beta),REAL(coef),INTEGER(conv),nthr);

    const char *names[] = {"coef","convergence",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, conv);

    UNPROTECT(3);
    return ans;
}

/* compute the fitted values without a lot of temporary matrices. */

SEXP get_one_way_fitted (SEXP beta, SEXP offset, SEXP groups, SEXP nthreads)
{
    SEXP ans;

    cmx bmx = SEXPtocmx1(beta);
    cmx omx = SEXPtocmx2(offset);

    int nthr = asInteger(nthreads);

    PROTECT(ans=allocMatrix(REALSXP,omx.nrow,omx.ncol));
    PROTECT(groups=coerceVector(groups,INTSXP));

    get_one_way_fit(&bmx, &omx, INTEGER(groups), REAL(ans), nthr);

    UNPROTECT(2);

    return ans;
}

/* intialize the starting coefficents */
SEXP get_levenberg_start(SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP use_null, SEXP nthreads)
{
    SEXP ans;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);

    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int u_null = asInteger(use_null);
    int nthr = asInteger(nthreads);

    PROTECT(ans=allocMatrix(REALSXP,ymx.nrow,gmx.ncol));

    get_leven_start(&ymx,&omx,&dmx,&wmx,&gmx,u_null,REAL(ans),nthr);

    UNPROTECT(2);

    return ans;
}

/* fit levenberg method for complex design */
SEXP fit_levenberg (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP beta, SEXP tolerance, SEXP max_iteration, SEXP nthreads)
{
    SEXP ans, mbeta, mu, dev, iter, failed;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx bmx = SEXPtocmx1(beta);

    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    // Setting up scalars.
    int maxit  = asInteger(max_iteration);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    // prepare outputs
    PROTECT(mbeta=duplicate(beta));
    PROTECT(mu=coerceVector(duplicate(y),REALSXP));
    PROTECT(dev=allocVector(REALSXP,ymx.nrow));
    PROTECT(iter=allocVector(INTSXP,ymx.nrow));
    PROTECT(failed=allocVector(LGLSXP,ymx.nrow));

    const char *names[] = {"coefficients","fitted.values","deviance","iter","failed",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, mbeta);
    SET_VECTOR_ELT(ans, 1, mu);
    SET_VECTOR_ELT(ans, 2, dev);   
    SET_VECTOR_ELT(ans, 3, iter);   
    SET_VECTOR_ELT(ans, 4, failed);   

    fit_leven(&ymx,&omx,&dmx,&wmx,&gmx,&bmx,tol,maxit,REAL(mu),REAL(mbeta),REAL(dev),INTEGER(iter),INTEGER(failed),nthr);

    UNPROTECT(7);
    return ans;
}

/* check whether the variance is below the Poisson bound. */

SEXP check_poisson_bound (SEXP mu, SEXP disp, SEXP s2, SEXP nthreads)
{
    SEXP ans;

    cmx umx = SEXPtocmx1(mu);
    cmx dmx = SEXPtocmx2(disp);
    cmx smx = SEXPtocmx2(s2);

    int nthr = asInteger(nthreads);

    PROTECT(ans=allocVector(LGLSXP, umx.nrow));
    check_poi_bound(&umx, &dmx, &smx, INTEGER(ans), nthr);

    UNPROTECT(1);

    return ans;
}

/* -----------------------------------------------------------
 * from good_turing.c
 * simple good turing function
 */
SEXP simple_good_turing (SEXP obs, SEXP freq, SEXP conf) 
{
    SEXP ans, prop, pzero;

    int nrows = length(obs);
    double conff = asReal(conf);
    
    PROTECT(pzero  = allocVector(REALSXP, 1));
    PROTECT(prop   = allocVector(REALSXP, nrows));

    good_turing(INTEGER(obs),INTEGER(freq),nrows,conff,REAL(pzero),REAL(prop));

    const char *names[] = {"P0","proportion",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, pzero);
    SET_VECTOR_ELT(ans, 1, prop);

    UNPROTECT(3);

    return ans;
}

/* --------------------------------------------------------------------
 * from interpolator.c
 * maximize interpolant function 
 */
SEXP maximize_interpolant(SEXP spts, SEXP ll)
{
    SEXP ans;

    PROTECT(spts = coerceVector(spts, REALSXP));
    PROTECT(ll   = coerceVector(ll, REALSXP));

    cmx lmx = SEXPtocmx1(ll);

    PROTECT(ans  = allocVector(REALSXP, lmx.nrow));

    max_interpolant(REAL(spts), &lmx, REAL(ans));

    UNPROTECT(3);

    return(ans);    
}

/* ------------------------------------------------------------------
 * from loess_by_col.c
 * wrapper function for loess_by_col 
 * SEXP n_cols was removed and n_cols was extracted from y
 */

SEXP loess_by_col(SEXP x, SEXP y, SEXP s) 
{
    SEXP ans, yy, lv;

    PROTECT(x=coerceVector(x,REALSXP));
    PROTECT(y=coerceVector(y,REALSXP));
    
    cmx ymx = SEXPtocmx1(y);
    int span = asInteger(s);

    PROTECT(yy=duplicate(y));
    PROTECT(lv=duplicate(x));

    loess_by_column(REAL(x),&ymx,span,REAL(yy),REAL(lv));

    const char *names[] = {"fitted.values", "leverages",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, yy);
    SET_VECTOR_ELT(ans, 1, lv);

    UNPROTECT(5);

    return ans;
}

/*---------------------------------------------------------------
 * from ql_glm.c
 */

/* input SEXP variables:
 * ym count or pesudo count matrix
 * um fitted value matrix
 * gm design matrix
 * dm dispersion compressMatrix
 * qd average quasi dispersion scalar
 * wm weights compressMatrix 
 * 
 * outputs:
 * s2 quasi-dispersion
 * dv adjusted deviance
 * df adjusted degree of freedom
 * 
 * dvmat unit deviance matrix
 * dfmat unit degree of freedom matrix
 * lvmat leverage matrix
 */

SEXP compute_adj_vec (SEXP ym, SEXP um, SEXP gm, SEXP dm, SEXP qd, SEXP wm, SEXP nthreads)
{
    SEXP ans, s2, dv, df;

    /* ensure double input for design */
    PROTECT(gm = coerceVector(gm, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(ym);
    cmx umx = SEXPtocmx1(um);
    cmx gmx = SEXPtocmx1(gm);

    cmx dmx = SEXPtocmx2(dm);
    cmx wmx = SEXPtocmx2(wm);

    /* average quasi dispersion */
    double prior = asReal(qd);
    int nthr = asInteger(nthreads);

    /* prepare output */
    PROTECT(df  = allocVector(REALSXP, ymx.nrow));
    PROTECT(dv  = allocVector(REALSXP, ymx.nrow));
    PROTECT(s2  = allocVector(REALSXP, ymx.nrow));

    const char *names[] = {"df","deviance","s2",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, df);
    SET_VECTOR_ELT(ans, 1, dv);
    SET_VECTOR_ELT(ans, 2, s2);

    /* call compute function */
    compute_adjust_vec(&ymx, &umx, &gmx, &dmx, prior, &wmx, REAL(df), REAL(dv), REAL(s2), nthr);

    /* check the number of protect */
    UNPROTECT(5);
    return ans;
}

/* NB empirical sample weights: thin .Call shim over sample_weights_nb
 * (sample_weights.c). Computes the adjusted unit dev/df internally and returns
 * the per-sample weight vector directly (length ncol). s2 is the prior variance
 * (s2.prior), floored at 1 inside the worker. */
SEXP sample_weights (SEXP y, SEXP mu, SEXP design, SEXP disp, SEXP qd, SEXP weights, SEXP s2, SEXP nthreads)
{
    PROTECT(design = coerceVector(design, REALSXP));
    PROTECT(s2 = coerceVector(s2, REALSXP));

    cmx ymx = SEXPtocmx1(y);
    cmx umx = SEXPtocmx1(mu);
    cmx gmx = SEXPtocmx1(design);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    double prior = asReal(qd);
    int nthr = asInteger(nthreads);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, ymx.ncol));
    sample_weights_nb(&ymx, &umx, &gmx, &dmx, prior, &wmx, REAL(s2), LENGTH(s2), nthr, REAL(ans));
    UNPROTECT(3);
    return ans;
}

/* compute average quasi-dispersion */
SEXP compute_ave_qd (SEXP ym, SEXP um, SEXP gm, SEXP dm, SEXP ag, SEXP wm, SEXP nthreads)
{
    SEXP ans;

    /* ensure double input for y and design */
    PROTECT(gm = coerceVector(gm, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(ym);
    cmx umx = SEXPtocmx1(um);
    cmx gmx = SEXPtocmx1(gm);

    cmx dmx = SEXPtocmx2(dm);
    cmx wmx = SEXPtocmx2(wm);

    int nthr = asInteger(nthreads);

    /*prepare output*/
    PROTECT(ans = allocVector(REALSXP, 1));

    /* call compute function */
    REAL(ans)[0]=update_prior(&ymx, &umx, &gmx, &dmx, &wmx, REAL(ag), nthr);

    /* check the number of protect */
    UNPROTECT(2);
    return ans;
}

/* --------------------------------------------------------------------
 * from binomial.c
 * */

/* fit one group mean for binomial models */

SEXP bin_one_group (SEXP y, SEXP offsets, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP nthreads)
{
    SEXP ans, prob, coef, devi, conv;

    PROTECT(y       = coerceVector(y, REALSXP));
    PROTECT(offsets = coerceVector(offsets, REALSXP));

    cmx ymx = SEXPtocmx1(y);
    cmx wmx = SEXPtocmx1(weights);
    cmx omx = SEXPtocmx2(offsets);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    PROTECT(prob=duplicate(y));
    PROTECT(coef=allocVector(REALSXP,ymx.nrow));
    PROTECT(devi=allocVector(REALSXP,ymx.nrow));
    PROTECT(conv=allocVector(LGLSXP,ymx.nrow));

    bin_one_group_mat(&ymx,&omx,&wmx,maxit,tol,REAL(prob),REAL(coef),REAL(devi),INTEGER(conv),nthr);

    const char *names[] = {"fitted.values","coef","deviance","convergence",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, prob);
    SET_VECTOR_ELT(ans, 1, coef);
    SET_VECTOR_ELT(ans, 2, devi);
    SET_VECTOR_ELT(ans, 3, conv);

    UNPROTECT(7);
    return ans;
}

/* wrapper function for the binomial IWLS fitter (general design)
 * inputs: y proportion matrix, offsets (compressMatrix), weights matrix,
 *   design matrix, max_iteration, tolerance, nthreads
 * returns a named list: coefficients, fitted.values, deviance, iter, failed
 */
SEXP bin_fit_iwls (SEXP y, SEXP offsets, SEXP weights, SEXP design, SEXP max_iteration, SEXP tolerance, SEXP nthreads)
{
    SEXP ans, beta, mu, dev, iter, failed;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx wmx = SEXPtocmx1(weights);
    cmx omx = SEXPtocmx2(offsets);

    // Setting up scalars.
    int maxit  = asInteger(max_iteration);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    // prepare outputs
    PROTECT(mu=coerceVector(duplicate(y),REALSXP));
    PROTECT(beta=allocMatrix(REALSXP,ymx.nrow,gmx.ncol));
    PROTECT(dev=allocVector(REALSXP,ymx.nrow));
    PROTECT(iter=allocVector(INTSXP,ymx.nrow));
    PROTECT(failed=allocVector(LGLSXP,ymx.nrow));

    const char *names[] = {"coefficients","fitted.values","deviance","iter","failed",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, beta);
    SET_VECTOR_ELT(ans, 1, mu);
    SET_VECTOR_ELT(ans, 2, dev);   
    SET_VECTOR_ELT(ans, 3, iter);   
    SET_VECTOR_ELT(ans, 4, failed);   

    bin_iwls_mat(&ymx,&omx,&wmx,&gmx,maxit,tol,REAL(mu),REAL(beta),REAL(dev),INTEGER(iter),INTEGER(failed),nthr);

    UNPROTECT(7);
    return ans;
}

/* unified single-call binomial fitter: oneway or general IWLS chosen in C.
 * Mirrors fit_glm; two (maxit,tol) pairs keep the former oneway/IWLS settings. */
SEXP bin_fit (SEXP y, SEXP offsets, SEXP weights, SEXP design, SEXP maxit_oneway, SEXP tol_oneway, SEXP maxit_iwls, SEXP tol_iwls, SEXP nthreads)
{
    SEXP ans, coef, mu, dev, iter, failed, method;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx (binomial convention:
     * y/design/weights dense, offsets compressed) */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx wmx = SEXPtocmx1(weights);
    cmx omx = SEXPtocmx2(offsets);

    int maxit_ow = asInteger(maxit_oneway);
    double tol_ow = asReal(tol_oneway);
    int maxit_iw = asInteger(maxit_iwls);
    double tol_iw = asReal(tol_iwls);
    int nthr = asInteger(nthreads);

    int ntag = ymx.nrow;
    int ncoef = gmx.ncol;
    int nlibs = ymx.ncol;

    PROTECT(coef   = allocMatrix(REALSXP, ntag, ncoef));
    PROTECT(mu     = allocMatrix(REALSXP, ntag, nlibs));
    PROTECT(dev    = allocVector(REALSXP, ntag));
    PROTECT(iter   = allocVector(INTSXP,  ntag));
    PROTECT(failed = allocVector(LGLSXP,  ntag));

    int method_val = 0;
    bin_fit_mat(&ymx, &omx, &wmx, &gmx, maxit_ow, tol_ow, maxit_iw, tol_iw, REAL(coef), REAL(mu), REAL(dev), INTEGER(iter), INTEGER(failed), &method_val, nthr);

    PROTECT(method = allocVector(STRSXP, 1));
    if (method_val == 0)
    {
        SET_STRING_ELT(method, 0, mkChar("oneway"));
    }
    else
    {
        SET_STRING_ELT(method, 0, mkChar("IWLS"));
    }

    const char *names[] = {"coefficients", "fitted.values", "deviance", "iter", "failed", "method", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, mu);
    SET_VECTOR_ELT(ans, 2, dev);
    SET_VECTOR_ELT(ans, 3, iter);
    SET_VECTOR_ELT(ans, 4, failed);
    SET_VECTOR_ELT(ans, 5, method);

    UNPROTECT(8);
    return ans;
}

/* single-call QL binomial fit: bin_fit_mat then compute_adjust_vec_bin, so
 * binQLFit needs one .Call instead of .cxx_bin_fit + .cxx_compute_adj_vec_bin. */
SEXP bin_ql_fit (SEXP y, SEXP z, SEXP weights, SEXP offsets, SEXP design, SEXP maxit_oneway, SEXP tol_oneway, SEXP maxit_iwls, SEXP tol_iwls, SEXP nthreads)
{
    SEXP ans, coef, mu, dev, iter, failed, method, df_adj, dev_adj, s2;

    /* ensure double design (binomial convention: counts/design/weights dense,
     * offsets compressed; weights may be NULL) */
    PROTECT(design = coerceVector(design, REALSXP));

    cmx ymx = SEXPtocmx1(y);
    cmx zmx = SEXPtocmx1(z);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);

    int has_w = !isNull(weights);
    cmx wmx;
    if (has_w)
    {
        wmx = SEXPtocmx1(weights);
    }

    int maxit_ow = asInteger(maxit_oneway);
    double tol_ow = asReal(tol_oneway);
    int maxit_iw = asInteger(maxit_iwls);
    double tol_iw = asReal(tol_iwls);
    int nthr = asInteger(nthreads);

    int ntag = ymx.nrow;
    int ncoef = gmx.ncol;
    int nlibs = ymx.ncol;

    PROTECT(coef    = allocMatrix(REALSXP, ntag, ncoef));
    PROTECT(mu      = allocMatrix(REALSXP, ntag, nlibs));
    PROTECT(dev     = allocVector(REALSXP, ntag));
    PROTECT(iter    = allocVector(INTSXP,  ntag));
    PROTECT(failed  = allocVector(LGLSXP,  ntag));
    PROTECT(df_adj  = allocVector(REALSXP, ntag));
    PROTECT(dev_adj = allocVector(REALSXP, ntag));
    PROTECT(s2      = allocVector(REALSXP, ntag));

    int method_val = 0;
    bin_ql_fit_mat(&ymx, &zmx, &wmx, has_w, &omx, &gmx, maxit_ow, tol_ow, maxit_iw, tol_iw, REAL(coef), REAL(mu), REAL(dev), INTEGER(iter), INTEGER(failed), &method_val, REAL(df_adj), REAL(dev_adj), REAL(s2), nthr);

    /* fit diagnostics, so binQLFit carries the same iter/failed/method as binFit */
    PROTECT(method = allocVector(STRSXP, 1));
    if (method_val == 0)
    {
        SET_STRING_ELT(method, 0, mkChar("oneway"));
    }
    else
    {
        SET_STRING_ELT(method, 0, mkChar("IWLS"));
    }

    const char *names[] = {"coefficients", "fitted.values", "deviance", "iter", "failed", "method", "df.residual.adj", "deviance.adj", "s2", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, mu);
    SET_VECTOR_ELT(ans, 2, dev);
    SET_VECTOR_ELT(ans, 3, iter);
    SET_VECTOR_ELT(ans, 4, failed);
    SET_VECTOR_ELT(ans, 5, method);
    SET_VECTOR_ELT(ans, 6, df_adj);
    SET_VECTOR_ELT(ans, 7, dev_adj);
    SET_VECTOR_ELT(ans, 8, s2);

    UNPROTECT(11);
    return ans;
}

/*---------------------------------------------------------------
 * from ql_bin.c
 */

/* input SEXP variables:
 * ym proportion matrix
 * um fitted value matrix
 * gm design matrix
 * cm coverage Matrix
 * wm weights Matrix 
 * 
 * outputs:
 * s2 quasi-dispersion
 * dv adjusted deviance
 * df adjusted degree of freedom
 * 
 * dvmat unit deviance matrix
 * dfmat unit degree of freedom matrix
 * lvmat leverage matrix
 */

/* binomial empirical sample weights: thin .Call shim over sample_weights_bin
 * (sample_weights.c). Takes the raw successes (y) and failures (y2); coverage,
 * proportion and scaled weights are derived inside the worker. Returns the
 * per-sample weight vector (length ncol). */
SEXP sample_weights_bin (SEXP y, SEXP y2, SEXP mu, SEXP design, SEXP weights, SEXP nthreads)
{
    PROTECT(design = coerceVector(design, REALSXP));

    cmx ymx  = SEXPtocmx1(y);
    cmx y2mx = SEXPtocmx1(y2);
    cmx umx  = SEXPtocmx1(mu);
    cmx gmx  = SEXPtocmx1(design);
    cmx wmx  = SEXPtocmx2(weights);

    int nthr = asInteger(nthreads);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, ymx.ncol));
    sample_weights_binom(&ymx, &y2mx, &umx, &gmx, &wmx, nthr, REAL(ans));
    UNPROTECT(2);
    return ans;
}

/*---------------------------------------------------------------
 * from ftest.c
 */

/* input SEXP variables:
 * counts   successes matrix (ntag x nlib)
 * counts2  failures matrix (ntag x nlib)
 * weights  observation weights (ntag x nlib) or NULL
 * design   design matrix (nlib x nbeta)
 * coef     1-based indices of the tested design columns (length ncoef)
 * deviance full-model residual deviance per gene (ntag)
 * s2post   EB-moderated quasi-dispersion per gene (ntag)
 * df_total total degrees of freedom per gene (ntag)
 * lfc      log2 fold-change thresholds aligned to coef (ncoef)
 * logFCt   log2 fold-changes (ntag x ncoef) aligned to coef
 * upshot   TRUE for the UPSHOT quadrature, FALSE for averaged TREAT
 *
 * outputs:
 * F        QL F-statistic (ntag)
 * PValue   QL / TREAT p-value (ntag)
 */

SEXP bin_ftest (SEXP counts, SEXP counts2, SEXP weights, SEXP design, SEXP coef, SEXP deviance, SEXP s2post, SEXP df_total, SEXP lfc, SEXP logFCt, SEXP upshot, SEXP nthreads)
{
    SEXP ans, Fout, Pout;

    /* ensure double input for the design, logFCt and per-gene vectors; integer for coef */
    PROTECT(design = coerceVector(design, REALSXP));
    PROTECT(logFCt = coerceVector(logFCt, REALSXP));
    PROTECT(deviance = coerceVector(deviance, REALSXP));
    PROTECT(s2post = coerceVector(s2post, REALSXP));
    PROTECT(df_total = coerceVector(df_total, REALSXP));
    PROTECT(lfc = coerceVector(lfc, REALSXP));
    PROTECT(coef = coerceVector(coef, INTSXP));

    /* convert matrices to cmx (counts may be integer; weights may be NULL) */
    cmx ymx = SEXPtocmx1(counts);
    cmx zmx = SEXPtocmx1(counts2);
    cmx gmx = SEXPtocmx1(design);
    cmx lmx = SEXPtocmx1(logFCt);

    int has_w = !isNull(weights);
    cmx wmx;
    if (has_w)
    {
        wmx = SEXPtocmx1(weights);
    }

    int ncoef = LENGTH(coef);
    int upsh = asLogical(upshot);
    int nthr = asInteger(nthreads);
    int ntag = ymx.nrow;

    PROTECT(Fout = allocVector(REALSXP, ntag));
    PROTECT(Pout = allocVector(REALSXP, ntag));

    bin_ftest_mat(&ymx, &zmx, &wmx, has_w, &gmx, INTEGER(coef), ncoef, REAL(deviance), REAL(s2post), REAL(df_total), REAL(lfc), &lmx, upsh, REAL(Fout), REAL(Pout), nthr);

    const char *names[] = {"F", "PValue", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, Fout);
    SET_VECTOR_ELT(ans, 1, Pout);

    UNPROTECT(10);
    return ans;
}

/*---------------------------------------------------------------
 * from ftest.c
 */

/* input SEXP variables:
 * counts     count matrix (ntag x nlib)
 * offset     log offset compressedMatrix (row-repeat)
 * dispersion NB dispersion compressedMatrix (per-gene column-repeat, or scalar)
 * weights    observation weights compressedMatrix (all-1 when R weights were NULL)
 * design     design matrix (nlib x nbeta)
 * coef       1-based indices of the tested design columns (length ncoef)
 * deviance   full-model residual deviance per gene (ntag)
 * s2post     EB-moderated quasi-dispersion per gene (ntag)
 * df_total   total degrees of freedom per gene (ntag)
 * lfc        log2 fold-change thresholds aligned to coef (ncoef)
 * logFCt     log2 fold-changes (ntag x ncoef) aligned to coef
 * upshot     TRUE for the UPSHOT quadrature, FALSE for averaged TREAT
 *
 * outputs:
 * F          QL F-statistic (ntag)
 * PValue     QL / TREAT p-value (ntag)
 */

SEXP glm_ftest (SEXP counts, SEXP offset, SEXP dispersion, SEXP weights, SEXP design, SEXP coef, SEXP deviance, SEXP s2post, SEXP df_total, SEXP lfc, SEXP logFCt, SEXP upshot, SEXP nthreads)
{
    SEXP ans, Fout, Pout;

    /* ensure double input for the design, logFCt and per-gene vectors; integer for coef */
    PROTECT(design = coerceVector(design, REALSXP));
    PROTECT(logFCt = coerceVector(logFCt, REALSXP));
    PROTECT(deviance = coerceVector(deviance, REALSXP));
    PROTECT(s2post = coerceVector(s2post, REALSXP));
    PROTECT(df_total = coerceVector(df_total, REALSXP));
    PROTECT(lfc = coerceVector(lfc, REALSXP));
    PROTECT(coef = coerceVector(coef, INTSXP));

    /* convert matrices to cmx (counts may be integer; offset/disp/weights are compressed) */
    cmx ymx = SEXPtocmx1(counts);
    cmx gmx = SEXPtocmx1(design);
    cmx lmx = SEXPtocmx1(logFCt);
    cmx omx = SEXPtocmx2(offset);
    cmx dmx = SEXPtocmx2(dispersion);
    cmx wmx = SEXPtocmx2(weights);

    int ncoef = LENGTH(coef);
    int upsh = asLogical(upshot);
    int nthr = asInteger(nthreads);
    int ntag = ymx.nrow;

    PROTECT(Fout = allocVector(REALSXP, ntag));
    PROTECT(Pout = allocVector(REALSXP, ntag));

    glm_ftest_mat(&ymx, &omx, &dmx, &wmx, &gmx, INTEGER(coef), ncoef, REAL(deviance), REAL(s2post), REAL(df_total), REAL(lfc), &lmx, upsh, REAL(Fout), REAL(Pout), nthr);

    const char *names[] = {"F", "PValue", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, Fout);
    SET_VECTOR_ELT(ans, 1, Pout);

    UNPROTECT(10);
    return ans;
}

/* -----------------------------------------------------------------------
 * from glm.c / diffsplice.c
 * extra fitters retained in the C backend (one-way layout, unified GLM,
 * differential splicing); each takes a trailing nthreads argument.
 */

/* fit one way layout for all groups */
SEXP fit_one_way (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP group, SEXP ngroups, SEXP max_iterations, SEXP tolerance, SEXP coef_start, SEXP nthreads)
{
    SEXP ans, coef, conv;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int n_groups = asInteger(ngroups);
    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    PROTECT(group = coerceVector(group, INTSXP));
    PROTECT(coef_start = coerceVector(coef_start, REALSXP));

    PROTECT(coef = allocMatrix(REALSXP, ymx.nrow, n_groups));
    PROTECT(conv = allocMatrix(LGLSXP, ymx.nrow, n_groups));

    fit_one_way_mat(&ymx, &omx, &dmx, &wmx, INTEGER(group), n_groups, maxit, tol, REAL(coef_start), REAL(coef), INTEGER(conv), nthr);

    const char *names[] = {"coefficients", "convergence", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, conv);

    UNPROTECT(5);
    return ans;
}

/* unified GLM fitting */
SEXP fit_glm (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP max_iterations, SEXP tolerance, SEXP coef_start, SEXP nthreads)
{
    SEXP ans, coef, mu, dev, iter, failed, method;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    int protect_count = 1; // design protected
    double *c_start = NULL;
    if (coef_start != R_NilValue)
    {
        PROTECT(coef_start = coerceVector(coef_start, REALSXP));
        c_start = REAL(coef_start);
        protect_count++;
    }

    int ntag = ymx.nrow;
    int ncoef = gmx.ncol;
    int nlibs = ymx.ncol;

    // Allocate outputs
    PROTECT(coef = allocMatrix(REALSXP, ntag, ncoef));
    PROTECT(mu = allocMatrix(REALSXP, ntag, nlibs));
    PROTECT(dev = allocVector(REALSXP, ntag));
    PROTECT(iter = allocVector(INTSXP, ntag));
    PROTECT(failed = allocVector(LGLSXP, ntag));
    protect_count += 5;

    int method_val = 0;
    fit_glm_mat(&ymx, &omx, &dmx, &wmx, &gmx, maxit, tol, c_start, REAL(coef), REAL(mu), REAL(dev), INTEGER(iter), INTEGER(failed), &method_val, nthr);

    PROTECT(method = allocVector(STRSXP, 1));
    if (method_val == 0)
    {
        SET_STRING_ELT(method, 0, mkChar("oneway"));
    }
    else
    {
        SET_STRING_ELT(method, 0, mkChar("levenberg"));
    }
    protect_count++;

    const char *names[] = {"coefficients", "fitted.values", "deviance", "iter", "failed", "method", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    protect_count++;

    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, mu);
    SET_VECTOR_ELT(ans, 2, dev);
    SET_VECTOR_ELT(ans, 3, iter);
    SET_VECTOR_ELT(ans, 4, failed);
    SET_VECTOR_ELT(ans, 5, method);

    UNPROTECT(protect_count);
    return ans;
}

/* single-call predictive fold-change fit: add library-size-scaled prior counts
 * then fit the NB GLM on the augmented data, returning only the coefficients.
 * Replaces the R predFC() used by glmFit.default and honours nthreads. */
SEXP pred_fc (SEXP y, SEXP offsets, SEXP priors, SEXP disp, SEXP weights, SEXP design, SEXP max_iterations, SEXP tolerance, SEXP nthreads)
{
    SEXP coef;

    /* ensure double input for design */
    PROTECT(design = coerceVector(design, REALSXP));

    /* convert matrix and compressMatrix to cmx (as in fit_glm: y/design dense,
     * offsets/priors/disp/weights compressed) */
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx pmx = SEXPtocmx2(priors);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    int ntag = ymx.nrow;
    int ncoef = gmx.ncol;

    PROTECT(coef = allocMatrix(REALSXP, ntag, ncoef));

    pred_fc_mat(&ymx, &omx, &pmx, &dmx, &wmx, &gmx, maxit, tol, REAL(coef), nthr);

    UNPROTECT(2);
    return coef;
}

/* Cox-Reid common-dispersion optimization: thin .Call shim over coxreid_disp_opt
 * (coxreid.c). Marshals inputs, runs the Brent search, returns the optimal
 * dispersion (par_opt^4) as a length-1 numeric. */
SEXP coxreid_disp(SEXP y, SEXP offsets, SEXP weights, SEXP design, SEXP lower, SEXP upper, SEXP tol, SEXP adjust, SEXP nthreads)
{
    PROTECT(design = coerceVector(design, REALSXP));
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx wmx = SEXPtocmx2(weights);
    double disp = coxreid_disp_opt(&ymx, &omx, &wmx, &gmx, asReal(lower), asReal(upper), asReal(tol), asLogical(adjust), asInteger(nthreads));
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = disp;
    UNPROTECT(2);
    return ans;
}

/* Cox-Reid common-dispersion over the most-abundant genes: thin .Call shim over
 * coxreid_disp_top_opt (coxreid.c), which chooses the top fraction of genes by
 * avelogcpm and runs the Brent search over that subset. Marshals inputs and
 * returns the optimal dispersion as a length-1 numeric. Replaces the R block in
 * glmQLFit.default that selected top genes and called dispCoxReid(). */
SEXP coxreid_disp_top(SEXP y, SEXP offsets, SEXP weights, SEXP design, SEXP avelogcpm, SEXP lower, SEXP upper, SEXP tol, SEXP adjust, SEXP nthreads)
{
    PROTECT(design = coerceVector(design, REALSXP));
    PROTECT(avelogcpm = coerceVector(avelogcpm, REALSXP));
    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx wmx = SEXPtocmx2(weights);
    double disp = coxreid_disp_top_opt(&ymx, &omx, &wmx, &gmx, REAL(avelogcpm), asReal(lower), asReal(upper), asReal(tol), asLogical(adjust), asInteger(nthreads));
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = disp;
    UNPROTECT(3);
    return ans;
}

/* differential splicing fit */
SEXP fit_diff_splice(SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP coef_idx, SEXP gene_nexons, SEXP gene_firstexon, SEXP nexons_approx, SEXP max_iterations, SEXP tolerance, SEXP coef_start, SEXP exon_dev, SEXP nthreads)
{
    SEXP ans, exon_LR, exon_coef, gene_LR;

    /* coerce design */
    PROTECT(design = coerceVector(design, REALSXP));

    cmx ymx = SEXPtocmx1(y);
    cmx gmx = SEXPtocmx1(design);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int coef = asInteger(coef_idx) - 1; // Convert to 0-indexed
    int nex_approx = asInteger(nexons_approx);
    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);
    int nthr = asInteger(nthreads);

    int protect_count = 1; // design protected
    double *c_start = NULL;
    if (coef_start != R_NilValue)
    {
        PROTECT(coef_start = coerceVector(coef_start, REALSXP));
        c_start = REAL(coef_start);
        protect_count++;
    }
    PROTECT(exon_dev = coerceVector(exon_dev, REALSXP));
    protect_count++;

    int nexons = ymx.nrow;
    int ngenes = length(gene_nexons);
    int nlibs = ymx.ncol;
    int ncoef = gmx.ncol;

    // Allocate outputs
    PROTECT(exon_LR = allocVector(REALSXP, nexons));
    PROTECT(exon_coef = allocVector(REALSXP, nexons));
    PROTECT(gene_LR = allocVector(REALSXP, ngenes));
    protect_count += 3;

    // We can also extract pointer to gene_nexons and gene_firstexon
    PROTECT(gene_nexons = coerceVector(gene_nexons, INTSXP));
    PROTECT(gene_firstexon = coerceVector(gene_firstexon, INTSXP));
    protect_count += 2;

    fit_diff_splice_mat(&ymx, &omx, &dmx, &wmx, &gmx, ngenes, nlibs, ncoef, coef, INTEGER(gene_nexons), INTEGER(gene_firstexon), nex_approx, maxit, tol, c_start, REAL(exon_dev), REAL(exon_LR), REAL(exon_coef), REAL(gene_LR), nthr);

    // Build named list
    const char *names[] = {"exon.LR", "exon.coef", "gene.LR", ""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, exon_LR);
    SET_VECTOR_ELT(ans, 1, exon_coef);
    SET_VECTOR_ELT(ans, 2, gene_LR);
    protect_count++;

    UNPROTECT(protect_count);
    return ans;
}
