#include "edgeR.h"

/* function to compute the unit deviance and the sum
 *
 * inputs:
 * y       row vector or matrix of counts (double)
 * mu      row vector or matrix of means
 * disp    compressedMatrix of dispersion
 * weights compressedMatrix of weights
 * 
 * output:
 * dev     row vector or matrix of deviance
 * 
 * comment:
 * this is converted from R_compute_nbdev.cpp written by Aaron
 */

void compute_nbdev_sum(cmx *y, cmx *mu, cmx *disp, cmx *weights, double *dev)
{
    int ntag=(y->nrow), nlib=(y->ncol);
    
    /* row vectors for y mu disp weights */
    double *dptr = R_Calloc(nlib,double);
    double *wptr = R_Calloc(nlib,double);

    double *yptr, *uptr;
    for(int tag=0;tag<ntag;++tag){
        get_row(disp,tag,dptr);
        get_row(weights,tag,wptr);

        /*
        dev[tag]=0;
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            dev[tag]   += compute_unit_nb_deviance((y->dmat)[ii],(mu->dmat)[ii],dptr[lib])*wptr[lib];
        }
        */

        yptr=(y->dmat)+tag;
        uptr=(mu->dmat)+tag;
        dev[tag]=0;
        for(int lib=0;lib<nlib;++lib,yptr+=ntag,uptr+=ntag){
            dev[tag]+=compute_unit_nb_deviance((*yptr),(*uptr),dptr[lib])*wptr[lib];
        }
    }

    R_Free(wptr);
    R_Free(dptr);

    return;
}

void compute_nbdev_unit(cmx *y, cmx *mu, cmx *disp, double *dev)
{
    int ntag=(y->nrow), nlib=(y->ncol);
    
    /* row vectors for y mu disp weights */
    double *dptr = R_Calloc(nlib,double);

    double *yptr, *uptr, *vptr;
    for(int tag=0;tag<ntag;++tag){
        get_row(disp,tag,dptr);

        /*
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            dev[ii]     = compute_unit_nb_deviance((y->dmat)[ii],(mu->dmat)[ii],dptr[lib]);
        }
        */
        
        yptr=(y->dmat)+tag;
        uptr=(mu->dmat)+tag;
        vptr=dev+tag;
        for(int lib=0;lib<nlib;++lib,yptr+=ntag,uptr+=ntag,vptr+=ntag){
            (*vptr)=compute_unit_nb_deviance((*yptr),(*uptr),dptr[lib]);
        }
    }

    R_Free(dptr);

    return;
}

/* Function to calculate the deviance. Note the protection for very large mu*phi (where we
 * use a gamma instead) or very small mu*phi (where we use the Poisson instead). This
 * approximation protects against numerical instability introduced by subtracting
 * a very large log value in (log mu) with another very large logarithm (log mu+1/phi).
 * We need to consider the 'phi' as the approximation is only good when the product is
 * very big or very small.
 */

/* This C version is converted from nbdev.cpp written by Aaron */

double compute_unit_nb_deviance (double y, double mu, double phi) 
{
    const double one_tenthousandth=1e-4, mildly_low_value=1e-8, one_million=1e6;

    double out;
	  
    // We add a small value to protect against zero during division and logging.
    y+=mildly_low_value;
    mu+=mildly_low_value;

    /* Calculating the deviance using either the Poisson (small phi*mu), the Gamma (large) or NB (everything else).
     * Some additional work is put in to make the transitions between families smooth.
     */
    if (phi < one_tenthousandth) {
		const double resid = y - mu;
		out = 2 * ( y * log(y/mu) - resid - 0.5*resid*resid*phi*(1+phi*(2/3*resid-y)) );
    } else {
		const double product=mu*phi;
		if (product > one_million) {
            out = 2 * ( (y - mu)/mu - log(y/mu) ) * mu/(1+product);
        } else {
			const double invphi=1/phi;
            out = 2 * (y * log( y/mu ) + (y + invphi) * log( (mu + invphi)/(y + invphi) ) );
        }
	}
    out = fmax2(out,0);

    return out;
}