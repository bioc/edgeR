#include "edgeR.h"


/* -----------------------------------------------------------------------
 * from add_prior_count.c
 * Adding a prior count to each observation. 
 * wrapper function for adding prior count 
 */

SEXP add_prior_count(SEXP y, SEXP offsets, SEXP priors) 
{
    SEXP yy, offset, ans;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx pmx = SEXPtocmx2(priors);

	PROTECT(yy = coerceVector(duplicate(y), REALSXP));

    /* prepare the offset output 
     * if both offsets and priors are row repeated
     * the adjusted offset will be row repeated
     * and the return a row vector
     * else return a matrix
     */
    int k = (omx.type>=2 && pmx.type>=2)? 1 : ymx.nrow;
    PROTECT(offset  = allocMatrix(REALSXP,k,ymx.ncol));

    if(omx.type>=2 && pmx.type>=2){
        add_prior_count_vec(&ymx,&omx,&pmx,REAL(yy),REAL(offset));
    }else{
        add_prior_count_mat(&ymx,&omx,&pmx,REAL(yy),REAL(offset));
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
SEXP compute_apl(SEXP y, SEXP mu, SEXP disp, SEXP weights, SEXP adjust, SEXP design) 
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

    PROTECT(ans=allocVector(REALSXP, ymx.nrow));

    compute_adj_profile_ll(&ymx, &umx, &dmx, &wmx, &gmx, do_adjust, REAL(ans));

    UNPROTECT(2);

    return ans;
}

/* -----------------------------------------------------------------------
 * from compute_cpm.c
 * wrapper function for computing CPM and logCPM 
 * and average log CPM 
 */

SEXP calculate_cpm_log (SEXP y, SEXP libsizes, SEXP priors) 
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx lmx = SEXPtocmx2(libsizes);
    cmx pmx = SEXPtocmx2(priors);

    PROTECT(ans = coerceVector(duplicate(y), REALSXP));

    calc_cpm_log(&ymx, &lmx, &pmx, REAL(ans));

    UNPROTECT(1);

    return ans;
}

SEXP calculate_cpm_raw (SEXP y, SEXP libsizes) 
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx lmx = SEXPtocmx2(libsizes);

    PROTECT(ans = coerceVector(duplicate(y), REALSXP));

    calc_cpm_raw(&ymx, &lmx, REAL(ans));

    UNPROTECT(1);

    return ans;
}

SEXP ave_log_cpm(SEXP y, SEXP offsets, SEXP priors, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance) 
{
    SEXP ans;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx pmx = SEXPtocmx2(priors);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);

    PROTECT(ans=allocVector(REALSXP,ymx.nrow));

    average_log_cpm(&ymx,&omx,&pmx,&dmx,&wmx,maxit,tol,REAL(ans));

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

SEXP fit_one_group (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP beta) 
{
    SEXP ans, coef, conv;

    cmx ymx = SEXPtocmx1(y);
    cmx omx = SEXPtocmx2(offsets);
    cmx dmx = SEXPtocmx2(disp);
    cmx wmx = SEXPtocmx2(weights);

    int maxit = asInteger(max_iterations);
    double tol = asReal(tolerance);

    PROTECT(coef=allocVector(REALSXP,ymx.nrow));
    PROTECT(conv=allocVector(LGLSXP,ymx.nrow));

    fit_one_group_mat(&ymx,&omx,&dmx,&wmx,maxit,tol,REAL(beta),REAL(coef),INTEGER(conv));

    const char *names[] = {"coef","convergence",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, coef);
    SET_VECTOR_ELT(ans, 1, conv);

    UNPROTECT(3);
    return ans;
}

/* compute the fitted values without a lot of temporary matrices. */

SEXP get_one_way_fitted (SEXP beta, SEXP offset, SEXP groups) 
{ 
    SEXP ans;

    cmx bmx = SEXPtocmx1(beta);
    cmx omx = SEXPtocmx2(offset);

    PROTECT(ans=allocMatrix(REALSXP,omx.nrow,omx.ncol));
    PROTECT(groups=coerceVector(groups,INTSXP));

    get_one_way_fit(&bmx, &omx, INTEGER(groups), REAL(ans));

    UNPROTECT(2);

    return ans;
}

/* intialize the starting coefficents */
SEXP get_levenberg_start(SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP use_null) 
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

    PROTECT(ans=allocMatrix(REALSXP,ymx.nrow,gmx.ncol));

    get_leven_start(&ymx,&omx,&dmx,&wmx,&gmx,u_null,REAL(ans));

    UNPROTECT(2);

    return ans;
}

/* fit levenberg method for complex design */
SEXP fit_levenberg (SEXP y, SEXP offsets, SEXP disp, SEXP weights, SEXP design, SEXP beta, SEXP tolerance, SEXP max_iteration) 
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

    fit_leven(&ymx,&omx,&dmx,&wmx,&gmx,&bmx,tol,maxit,REAL(mu),REAL(mbeta),REAL(dev),INTEGER(iter),INTEGER(failed));

    UNPROTECT(7);
    return ans;
}

/* check whether the variance is below the Poisson bound. */

SEXP check_poisson_bound (SEXP mu, SEXP disp, SEXP s2) 
{
    SEXP ans;

    cmx umx = SEXPtocmx1(mu);
    cmx dmx = SEXPtocmx2(disp);
    cmx smx = SEXPtocmx2(s2);

    PROTECT(ans=allocVector(LGLSXP, umx.nrow));
    check_poi_bound(&umx, &dmx, &smx, INTEGER(ans));

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

SEXP compute_adj_vec (SEXP ym, SEXP um, SEXP gm, SEXP dm, SEXP qd, SEXP wm)
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
    compute_adjust_vec(&ymx, &umx, &gmx, &dmx, prior, &wmx, REAL(df), REAL(dv), REAL(s2));

    /* check the number of protect */
    UNPROTECT(5);
    return ans;
}

SEXP compute_adj_mat (SEXP ym, SEXP um, SEXP gm, SEXP dm, SEXP qd, SEXP wm)
{
    SEXP ans, dfmat, dvmat, lvmat;

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

    /* prepare output */
    PROTECT(dfmat  = duplicate(um));
    PROTECT(dvmat  = duplicate(um));
    PROTECT(lvmat  = duplicate(um));

    const char *names[] = {"unit.df","unit.deviance","leverage",""};
    PROTECT(ans = Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(ans, 0, dfmat);
    SET_VECTOR_ELT(ans, 1, dvmat);
    SET_VECTOR_ELT(ans, 2, lvmat);

    /* call compute function */
    compute_adjust_mat(&ymx, &umx, &gmx, &dmx, prior, &wmx, REAL(dfmat), REAL(dvmat), REAL(lvmat));

    /* check the number of protect */
    UNPROTECT(5);
    return ans;
}

/* compute average quasi-dispersion */
SEXP compute_ave_qd (SEXP ym, SEXP um, SEXP gm, SEXP dm, SEXP ag, SEXP wm)
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

    /*prepare output*/
    PROTECT(ans = allocVector(REALSXP, 1));

    /* call compute function */
    REAL(ans)[0]=update_prior(&ymx, &umx, &gmx, &dmx, &wmx, REAL(ag));

    /* check the number of protect */
    UNPROTECT(2);
    return ans;
}


