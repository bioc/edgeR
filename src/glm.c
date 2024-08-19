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
 */

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

	if (ISNA(cur_beta)) {
		cur_beta=0;
 	   	double totweight=0;
		for (int lib=0; lib<nlib; ++lib) {
			double cur_val = counts[lib];
			if (cur_val > low_value) {
				cur_beta += cur_val / exp(offset[lib]) * weights[lib];
				allzero = 0;
			}
			totweight += weights[lib];
		}
		cur_beta = log(cur_beta/totweight);
	} else {
		for (int lib=0; lib<nlib; ++lib) { 
			if (counts[lib] > low_value) { 
                allzero = 0; 
                break; 
            }
		}
	}

	// Skipping to a result for all-zero rows.
	if (allzero) 
    {
        (*beta) = R_NegInf;
        (*conv) = 1; 
        return;
    }

	// Newton-Raphson iterations to converge to mean.
    (*conv)=0;
	for (int i=0; i<maxit; ++i) {
		double dl=0;
 	    double info=0;
		for (int lib=0; lib<nlib; ++lib) {
			double mu=exp(cur_beta+offset[lib]), denominator=1+mu*disp[lib];
			dl+=(counts[lib]-mu)/denominator * weights[lib];
			info+=mu/denominator * weights[lib];
		}
		double step=dl/info;
		cur_beta+=step;
		if (fabs(step)<tolerance) {
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
void fit_one_group_mat (cmx *y, cmx *offsets, cmx *disp, cmx *weights, int maxit, double tol, double *beta, double *coef, int *conv) 
{
    int ntag = (y->nrow), nlib = (y->ncol);

    // Preparing for possible Poisson sums.
    int disp_is_zero=1, weight_is_one=1, tag_start=0;
    double sum_lib=0, zero=0, one=1;

    /* row vectors for y offsets disp weights */
    double *yptr = R_Calloc(nlib,double);
    double *optr = R_Calloc(nlib,double);
    double *wptr = R_Calloc(nlib,double);
    double *dptr = R_Calloc(nlib,double);
    
    if ((offsets->type) >= 2) {
        get_row(offsets,tag_start,optr);
        for (int lib=0; lib<nlib; ++lib) {
            sum_lib+=exp(optr[lib]);
        }
    }

    if ((disp->type) >= 2) {
        disp_is_zero=check_row_scalar(disp,tag_start,zero);
    }

    if ((weights->type) >= 2) {
        weight_is_one=check_row_scalar(weights,tag_start,one);
    }

    // Iterating through tags and fitting.
	for (int tag=0; tag<ntag; ++tag) {
        get_row4(y,offsets,disp,weights,tag,yptr,optr,dptr,wptr);

        // Checking for the Poisson special case with all-unity weights and all-zero dispersions.
        if ((disp->type <= 1)) {
            disp_is_zero=check_row_scalar(disp,tag,zero);
        }
        if ((weights->type) <= 1) {
            weight_is_one=check_row_scalar(weights,tag,one);
        }

        if (disp_is_zero && weight_is_one) {
            if ((offsets->type) <= 1) {
                // Only recalculate sum of library sizes if it has changed.
                sum_lib=0;
                for (int lib=0; lib<nlib; ++lib) { sum_lib += exp(optr[lib]); }
            }

            double sum_counts=0;
            for(int lib=0;lib<nlib;++lib){
                sum_counts += yptr[lib];
            } 

            if (sum_counts==0) {
                coef[tag]=R_NegInf;
            } else {
                coef[tag]=log(sum_counts/sum_lib);
            }
            conv[tag]=1;
        } 
        else 
        {
            // Otherwise going through NR iterations.
            double ocoef;
            int oconv;
            
            glm_one_group_vec(nlib, yptr, optr, dptr, wptr, maxit, tol, beta[tag], &ocoef, &oconv);
            
            coef[tag]=ocoef;
            conv[tag]=oconv;
        }
	}

    R_Free(yptr);
    R_Free(optr);
    R_Free(wptr);
    R_Free(dptr);

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

void get_one_way_fit(cmx *beta, cmx *offsets, int *group, double *mu)
{
    int ntag=(offsets->nrow), nlib=(offsets->ncol), nbeta=(beta->ncol);

    double *optr = R_Calloc(nlib, double);
    double *bptr = R_Calloc(nbeta, double);

    double *uptr;
    for(int tag=0;tag<ntag;++tag){
        get_row(offsets,tag,optr);
        get_row(beta,tag,bptr);

        /*
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            mu[ii]      = exp(optr[lib]+bptr[group[lib]]);
        }
        */
        uptr=mu+tag;
        for(int lib=0;lib<nlib;++lib,uptr+=ntag){
            (*uptr) = exp(optr[lib]+bptr[group[lib]]);
        }
    }

    R_Free(optr);
    R_Free(bptr);

    return;
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
void get_leven_start (cmx *y, cmx *offsets, cmx *disp, cmx *weights, cmx *design, int use_null, double *beta)
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

    double *xdpt    = R_Calloc(nlib*ncoef, double);
    double *tau     = R_Calloc(ncoef, double);
    double *effects = R_Calloc(nlib, double);

    // Setting up the workspace for dgeqrf and dormqr with optimal WORK
    double tmpwork;
    F77_CALL(dgeqrf)(&nlib, &ncoef, xdpt, &nlib, tau, &tmpwork, &lwork_geqrf, &info);
    lwork_geqrf=(int)(tmpwork+0.5);
    if (lwork_geqrf < 1) { lwork_geqrf = 1; }

    F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, xdpt, &nlib, tau, effects, &nlib, &tmpwork, &lwork_ormqr, &info FCONE FCONE);
    lwork_ormqr=(int)(tmpwork+0.5);
    if (lwork_ormqr < 1) { lwork_ormqr = 1; } 

    double *work_geqrf = R_Calloc(lwork_geqrf, double);
    double *work_ormqr = R_Calloc(lwork_ormqr, double);       

    // row vectors for y offsets disp weights 
    double *yptr = R_Calloc(nlib,double);
    double *optr = R_Calloc(nlib,double);
    double *wptr = R_Calloc(nlib,double);
    double *dptr = R_Calloc(nlib,double);

    double *bptr;
    if(use_null){
        // make a copy of design matrix
        for(int i=0;i<nlib*ncoef;++i) xdpt[i] = (design->dmat)[i];
        
        /* DGEQRF computes a QR factorization of a real M-by-N matrix A:
         * A = Q * R.
         *
         * https://www.netlib.org/lapack/lapack-3.1.1/html/dgeqrf.f.html
         */
        F77_CALL(dgeqrf)(&nlib, &ncoef, xdpt, &nlib, tau, work_geqrf, &lwork_geqrf, &info);
        if (info) { error("QR decomposition failed"); }  

        for(int tag=0;tag<ntag;++tag){
            get_row4(y,offsets,disp,weights,tag,yptr,optr,dptr,wptr);
            
            // Computing weighted average of the count:library size ratios.
            double sum_weight=0, sum_exprs=0;
            for (int lib=0; lib<nlib; ++lib) {
                double curN=exp(optr[lib]);
                double curweight=wptr[lib]*curN/(1 + dptr[lib] * curN);
                sum_exprs  += yptr[lib] * curweight / curN;
                sum_weight += curweight;
            }

            for (int lib=0; lib<nlib; ++lib) {
                effects[lib]=log(sum_exprs/sum_weight);
            }
            /*  DORMQR overwrites the general real M-by-N matrix C with
             *
             *                  SIDE = 'L'     SIDE = 'R'
             *  TRANS = 'N':      Q * C          C * Q
             *  TRANS = 'T':      Q**T * C       C * Q**T
             *
             *  where Q is a real orthogonal matrix defined as the product of k
             *  elementary reflectors
             *
             *        Q = H(1) H(2) . . . H(k)
             *
             *  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
             *  if SIDE = 'R'.
             *
             *  https://netlib.org/lapack/explore-3.1.1-html/dormqr.f.html
             */
            F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, xdpt, &nlib, tau, effects, &nlib, work_ormqr, &lwork_ormqr, &info FCONE FCONE);
            if (info) { error("Q**T multiplication failed"); }

            /*  DTRTRS solves a triangular system of the form
             *
             *     A * X = B  or  A**T * X = B,
             *
             *  where A is a triangular matrix of order N, and B is an N-by-NRHS
             *  matrix.  A check is made to verify that A is nonsingular.
             *
             * https://netlib.org/lapack/explore-3.1.1-html/dtrtrs.f.html
             */
            F77_CALL(dtrtrs)(&uplo, &trans_t, &diag, &ncoef, &unity, xdpt, &nlib, effects, &nlib, &info FCONE FCONE FCONE);
            if (info) { error("failed to solve the triangular system"); }

            bptr = beta+tag;
            for(int var=0;var<ncoef;++var,bptr+=ntag){
                (*bptr) = effects[var];
            }
        } 
    }
    else{
        // Finding the delta
        double delta = max_cmx(y);
        delta = fmin(delta, 1.0/6);

        // check whether weights is row repeated
        if((weights->type) >= 2){
            int tag_start = 0;
            get_row(weights,tag_start,wptr);
            
            for(int i=0;i<nlib*ncoef;++i) xdpt[i]=(design->dmat)[i]*sqrt(wptr[i % ncoef]);          
            F77_CALL(dgeqrf)(&nlib, &ncoef, xdpt, &nlib, tau, work_geqrf, &lwork_geqrf, &info);
            if (info) { error("QR decomposition failed"); }  
        }

        for(int tag=0;tag<ntag;++tag){
            get_row4(y,offsets,disp,weights,tag,yptr,optr,dptr,wptr);

            if((weights->type) <= 1 && tag >= 1){
                for(int i=0;i<nlib*ncoef;++i) xdpt[i]=(design->dmat)[i]*sqrt(wptr[i % ncoef]);          
                F77_CALL(dgeqrf)(&nlib, &ncoef, xdpt, &nlib, tau, work_geqrf, &lwork_geqrf, &info);
                if (info) { error("QR decomposition failed"); }  
            }
            
            // Computing normalized log-expression values.
            for (int lib=0; lib<nlib; ++lib) {
                yptr[lib]=log(fmax(delta, yptr[lib])) - optr[lib];
            }

            for (int lib=0; lib<nlib; ++lib) {
                effects[lib]=yptr[lib]*sqrt(wptr[lib]);
            }

            F77_CALL(dormqr)(&side, &trans_o, &nlib, &unity, &ncoef, xdpt, &nlib, tau, effects, &nlib, work_ormqr, &lwork_ormqr, &info FCONE FCONE);
            if (info) { error("Q**T multiplication failed"); }

            F77_CALL(dtrtrs)(&uplo, &trans_t, &diag, &ncoef, &unity, xdpt, &nlib, effects, &nlib, &info FCONE FCONE FCONE);
            if (info) { error("failed to solve the triangular system"); }

            bptr=beta+tag;
            for(int coef=0;coef<ncoef;++coef,bptr+=ntag){
                (*bptr) = effects[coef];
            }
        }
    }

    R_Free(xdpt);
    R_Free(tau);
    R_Free(effects);

    R_Free(work_geqrf);
    R_Free(work_ormqr);

    R_Free(yptr);
    R_Free(optr);
    R_Free(wptr);
    R_Free(dptr);

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
                double *mu, double *mbeta, double *dev, int *iter, int *failed) 
{
    int ntag = (y->nrow), nlib = (y->ncol), ncoef = (design->ncol);
    double *dm; dm = (design->dmat);

    // row vectors for y offsets disp weights 
    double *yptr = R_Calloc(nlib,double);
    double *optr = R_Calloc(nlib,double);
    double *wptr = R_Calloc(nlib,double);
    double *dptr = R_Calloc(nlib,double);

    // obt
    double *bptr  = R_Calloc(ncoef,double);
    // nbt
    double *nbeta = R_Calloc(ncoef,double);
    // dl
    double *dl    = R_Calloc(ncoef,double);
    // db
    double *dbeta = R_Calloc(ncoef,double);

    // working matrix XtWX and copy
    double *xtwx  = R_Calloc(ncoef*ncoef,double);
    double *xtwxc = R_Calloc(ncoef*ncoef,double);

    // omu
    double *uptr = R_Calloc(nlib,double);
    // nmu
    double *nmu  = R_Calloc(nlib,double);
    
    // working weights
    double *zwpt = R_Calloc(nlib,double);
    // derivative
    double *drvt = R_Calloc(nlib,double);

    int oiter, ofail;
    double odev;
    double *uupt, *bbpt;
    for (int tag=0; tag<ntag; ++tag) {
        get_row4(y,offsets,disp,weights,tag,yptr,optr,dptr,wptr);
        get_row(beta,tag,bptr);

        fit_leven_vec(nlib,yptr,optr,dptr,wptr,ncoef,dm,maxit,tol,zwpt,drvt,dl,dbeta,xtwx,xtwxc,nbeta,nmu,bptr,uptr,&odev,&oiter,&ofail);                

        /*
        if (ofail) { error("solution using Cholesky decomposition failed"); } 
        */

        /*      
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            mu[ii] = uptr[lib];
        }
        
        for(int coef=0;coef<ncoef;++coef){
            R_xlen_t ii = (R_xlen_t)(ntag)*coef+tag;
            mbeta[ii] = bptr[coef];
        }
        */  

        uupt=mu+tag;
        for(int lib=0;lib<nlib;++lib,uupt+=ntag){
            (*uupt) = uptr[lib];
        }
        bbpt=mbeta+tag;
        for(int coef=0;coef<ncoef;++coef,bbpt+=ntag){
            (*bbpt) = bptr[coef];
        }

		dev[tag]    = odev;
		iter[tag]   = oiter;
		failed[tag] = ofail;
    }
    R_Free(yptr);
    R_Free(optr);
    R_Free(wptr);
    R_Free(dptr);

    R_Free(bptr);
    R_Free(nbeta);
    R_Free(dl);
    R_Free(dbeta);

    R_Free(xtwx);
    R_Free(xtwxc);
    R_Free(uptr);
    R_Free(nmu);

    R_Free(zwpt);
    R_Free(drvt);

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

void fit_leven_vec(int nlib, double *y, double *offset, double *disp, double *w, int ncoef, double *dm, int maxit, double tol,
                   double *zwpt, double *drvt, double *dl, double *db, double *xtwx, double *xtwc, double *nbt, double *nmu, 
                   double *obt, double *omu, double *odev, int *oiter, int *ofail) 
{
    // const values
    const double low_value = 1e-10;
    const double one_millionth = 1e-6;
    const double supremely_low_value = 1e-13;
    const double ridiculously_low_value = 1e-100;

    const char uplo='U';
    const int nrhs=1;

	// We expect 'beta' to be supplied. We then check the maximum value of the counts.
    double ymax=0;
    for (int lib=0; lib<nlib; ++lib) { 
        ymax = (y[lib]>ymax)? y[lib] : ymax;
 	}	

    // If we start off with all entries at zero, there's really no point continuing. 
    if (ymax<low_value) {
        for(int coef=0;coef<ncoef; ++coef){
            obt[coef] = NA_REAL;
        }
        for(int lib=0;lib<nlib; ++lib){
            omu[lib] = 0; 
        }
        (*odev)  = 0;
        (*oiter) = 0;
        (*ofail) = 0;
        return;
    }
    
	// Otherwise, we compute 'mu' based on 'beta'. Returning if there are no coefficients!
    fit_leven_autofill(ncoef,obt,nlib,offset,omu,dm);

    double dev=0;
    for(int lib=0;lib<nlib;++lib){
        dev += w[lib]*compute_unit_nb_deviance(y[lib],omu[lib],disp[lib]);
    }
   
    /* it is possible to pass a null design?
    if(ncoef == 0) { 
        (*ofail)=0;
        return; 
    }
    */

    // Iterating using reweighted least squares; setting up assorted temporary objects.
    double max_info=-1, lambda=0;
    int info=0, iter=0, failed=0;
    while ((++iter) <= maxit) {

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
        for (int lib=0; lib<nlib; ++lib) {
            double cur_mu=omu[lib];
			double denom=(1+cur_mu*disp[lib]);
            zwpt[lib]=cur_mu/denom*w[lib];
            drvt[lib]=(y[lib]-cur_mu)/denom*w[lib];
        }

        compute_xtwx(nlib, ncoef, dm, zwpt, xtwx);

        double *dmc, *xtwxIt, *xtwcIt; dmc=dm, xtwxIt=xtwx;
        for (int coef=0; coef<ncoef; ++coef, dmc+=nlib, xtwxIt+=ncoef) {
            dl[coef]=0;
            for(int lib=0;lib<nlib;++lib) { dl[coef] += drvt[lib]*dmc[lib]; }
            if (xtwxIt[coef]>max_info) { max_info=xtwxIt[coef]; }
        }
        if (iter==1) {
            lambda=max_info*one_millionth;
            if (lambda < supremely_low_value) { lambda=supremely_low_value; } 
        }

        /* Levenberg/Marquardt damping reduces step size until the deviance increases or no 
         * step can be found that increases the deviance. In short, increases in the deviance
         * are enforced to avoid problems with convergence.
         */ 
        int lev=0, low_dev=0;
        while (++lev) {
			do {
             	/* We need to set up copies as the decomposition routine overwrites the originals, and 
 				 * we want the originals in case we don't like the latest step. For efficiency, we only 
	 			 * refer to the upper triangular for the XtWX copy (as it should be symmetrical). We also add 
	 			 * 'lambda' to the diagonals. This reduces the step size as the second derivative is increased.
        	     */
                xtwxIt=xtwx, xtwcIt=xtwc;
         		for (int coef1=0; coef1<ncoef; ++coef1, xtwxIt+=ncoef, xtwcIt+=ncoef) {
                    for(int coef2=0;coef2<=coef1;++coef2){
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
                 */
                F77_CALL(dpotrf)(&uplo, &ncoef, xtwc, &ncoef, &info FCONE);
                if (info!=0) {
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
                	if (lambda <= 0) { lambda=ridiculously_low_value; } // Just to make sure it actually increases.
                } else { 
                    break; 
                }
            } while (1);

            for(int coef=0;coef<ncoef;++coef) { db[coef]=dl[coef]; }

            /*  DPOTRS solves a system of linear equations A*X = B with a symmetric
             *  positive definite matrix A using the Cholesky factorization
             *  A = U**T*U or A = L*L**T computed by DPOTRF.
             *
             *  https://www.netlib.org/lapack/lapack-3.1.1/html/dpotrs.f.html
             */
            F77_CALL(dpotrs)(&uplo, &ncoef, &nrhs, xtwc, &ncoef, db, &ncoef, &info FCONE);
            if (info!=0) { error("solution using the Cholesky decomposition failed"); }

            // Updating beta and the means. 'dbeta' stores 'Y' from the solution of (X*VX)Y=dl, corresponding to a NR step.
            for (int coef=0; coef<ncoef; ++coef) { 
                nbt[coef]=obt[coef]+db[coef]; 
            }

            fit_leven_autofill(ncoef,nbt,nlib,offset,nmu,dm);

            /* Checking if the deviance has decreased or if it's too small to care about. Either case is good
             * and means that we'll be using the updated fitted values and coefficients. Otherwise, if we have
             * to repeat the inner loop, then we want to do so from the original values (as we'll be scaling
             * lambda up so we want to retake the step from where we were before). This is why we don't modify the values
             * in-place until we're sure we want to take the step.
             */
            double ndev=0;  
            for(int lib=0;lib<nlib;++lib){ 
                ndev += w[lib]*compute_unit_nb_deviance(y[lib],nmu[lib],disp[lib]);
            }

            if (ndev/ymax < supremely_low_value) { low_dev=1; }
            if (ndev <= dev || low_dev) {
                for(int coef=0;coef<ncoef;++coef) { 
                    obt[coef]=nbt[coef]; 
                }
                for(int lib=0;lib<nlib;++lib) { 
                    omu[lib]=nmu[lib]; 
                }
                dev=ndev;
                break; 
            }
            
            // Increasing lambda, to increase damping. Again, we have to make sure it's not zero.
            lambda*=2;
            if (lambda <= 0) { lambda=ridiculously_low_value; }

            // Excessive damping; steps get so small that it's pointless to continue.
            if (lambda/max_info > 1/supremely_low_value) { 
            	failed=1; 
            	break; 
            }
        } 

        /* Terminating if we failed, if divergence from the exact solution is acceptably low 
         * (cross-product of dbeta with the log-likelihood derivative) or if the actual deviance 
         * of the fit is acceptably low.
         */
        double divergence=0;
        for(int coef=0;coef<ncoef;++coef) {divergence += dl[coef]*db[coef]; } 

        if ( failed || low_dev || (divergence<tol)) { 
            (*odev)  = dev;
            (*oiter) = iter;
            (*ofail) = failed;
            break; 
        }

        /* If we quit the inner levenberg loop immediately and survived all the break conditions above, that means that deviance is decreasing
 		 * substantially. Thus, we need larger steps to get there faster. To do so, we decrease the damping factor. Note that this only applies 
 		 * if we didn't decrease the damping factor in the inner levenberg loop, as that would indicate that we need to slow down. 
         */
        if (lev==1) { lambda/=10; }
    }
	return;
}

// Computing updated mean = beta %*% design + offset.
void fit_leven_autofill(int ncoef, double *beta, int nlib, double *offset, double *mu, double *dm) 
{
    for(int lib=0; lib<nlib; ++lib) { mu[lib]=offset[lib]; }
    
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
     */
    F77_CALL(dgemv)(&trans, &nlib, &ncoef, &first_scaling, dm, &nlib, beta, &incx, &second_scaling, mu, &incy FCONE);
	
    for (int lib=0; lib<nlib; ++lib) {
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

void check_poi_bound(cmx *mu, cmx *disp, cmx *s2, int *out)
{
    int ntag=(mu->nrow), nlib=(mu->ncol);

    double *dptr = R_Calloc(nlib, double);
    double *sptr = R_Calloc(nlib, double);

    double *uptr;
    for(int tag=0;tag<ntag;++tag){
        get_row(disp,tag,dptr);
        get_row(s2,tag,sptr);
        
        /* a gene is below poisson bound if existing i,
         * such that s2*(1+u_i*d_i) < 1, that is
         * s2 < (1+u_i*d_i)^{-1}, which means
         * the quasi dispersion is too small
         */
        out[tag]=0;
        uptr=(mu->dmat)+tag;
        for(int lib=0;lib<nlib;++lib,uptr+=ntag){
            if(sptr[lib]*(1+dptr[lib]*(*uptr)) < 1){
                out[tag]=1;
                break;
            }    
        }
    }

    R_Free(dptr);
    R_Free(sptr);

    return;
}
