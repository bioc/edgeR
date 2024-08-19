#include "edgeR.h"

/* function to compute cpm or logCPM
 * 
 * inputs:
 * y        matrix of fitted values
 * libsizes compressedMatrix of library sizes
 * priors   compressedMatrix of priors
 * 
 * output:
 * cpm      matrix of cpm or logCPM 
 * 
 * comment:
 * this is converted from R_calculate_cpm.cpp written by Aaron
 */


void calc_cpm_log(cmx *y, cmx *libsizes, cmx *priors, double *cpm)
{
    int ntag=(y->nrow), nlib=(y->ncol);
    int log_in=0, log_out=1;    
    
    const double LNtwo = log(2), LNmillion = log(1e6);
    
    double *lptr = R_Calloc(nlib, double);
    double *pptr = R_Calloc(nlib, double);

    /*
    for(int tag=0;tag<ntag;++tag){
        get_row(libsizes,tag,lptr);
        get_row(priors,tag,pptr);

        compute_offsets(priors,libsizes,tag,log_in,log_out,pptr,lptr);

        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag) * lib + tag;
            cpm[ii] += pptr[lib];
            cpm[ii] = (cpm[ii]>0)? (log(cpm[ii])-lptr[lib]+LNmillion)/LNtwo : R_NaN; 
        }
    }
    */
    
    double *cptr;
    for(int tag=0;tag<ntag;++tag,++cpm){
        get_row(libsizes,tag,lptr);
        get_row(priors,tag,pptr);

        compute_offsets(priors,libsizes,tag,log_in,log_out,pptr,lptr);

        cptr=cpm;
        for(int lib=0;lib<nlib;++lib,cptr+=ntag){
            (*cptr)+= pptr[lib];
            (*cptr) = ((*cptr)>0)? (log(*cptr)-lptr[lib]+LNmillion)/LNtwo : R_NaN; 
        }
    }

    R_Free(lptr);
    R_Free(pptr);

    return;
}

void calc_cpm_raw(cmx *y, cmx *libsizes, double *cpm)
{
    const double one_million=1e6;

    int ntag=(y->nrow), nlib=(y->ncol);

    double *lptr = R_Calloc(nlib, double);

    /*
    for(int tag=0;tag<ntag;++tag){
        get_row(libsizes,tag,lptr);

        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag) * lib + tag;
            cpm[ii] = cpm[ii]*one_million/lptr[lib];
        }
    }
    */
    
    double *cptr;
    for(int tag=0;tag<ntag;++tag,++cpm){
        get_row(libsizes,tag,lptr);
        cptr=cpm;
        for(int lib=0;lib<nlib;++lib,cptr+=ntag){
            (*cptr) = (*cptr)*one_million/lptr[lib];
        }
    }

    R_Free(lptr);

    return;
}

/* function to compute average logCPM
 * 
 * inputs:
 * y       matrix of fitted values
 * offsets compressedMatrix of offsets 
 * priors  compressedMatrix of priors
 * disp    compressedMatrix of dispersion
 * weights compressedMatrix of weights
 * 
 * passing to glm_one_group_vec():
 * maxit   maximal iteration
 * tol     tolerance
 * 
 * output:
 * avecpm  matrix of average logCPM 
 * 
 * comment:
 * this is converted from R_ave_log_cpm.cpp written by Aaron
 */
void average_log_cpm(cmx *y, cmx *offsets, cmx *priors, cmx *disp, cmx *weights, int maxit, double tol, double *avecpm) 
{
    const double LNmillion=log(1e6), LNtwo=log(2);

    int ntag = (y->nrow), nlib = (y->ncol);
    int log_in=1, log_out=1;

    /* row vectors for y offsets, weights, dispersion, priors */
    double *yptr = R_Calloc(nlib,double);
    double *optr = R_Calloc(nlib,double);
    double *wptr = R_Calloc(nlib,double);
    double *dptr = R_Calloc(nlib,double);
    double *pptr = R_Calloc(nlib,double);

    /* whether both offsets and priors are row repeated */
    int repeat_row = ((offsets->type) >= 2) && ((priors->type) >=2);

    if(repeat_row)
    {
        int tag_start = 0;
        compute_offsets(priors,offsets,tag_start,log_in,log_out,pptr,optr);
    }

    // Returning average log-cpm
    double ocoef;
    int oconv;

    for (int tag=0; tag<ntag; ++tag) {
        get_row3(y,disp,weights,tag,yptr,dptr,wptr);

        if((tag >= 1) && (!repeat_row)){
            compute_offsets(priors,offsets,tag,log_in,log_out,pptr,optr);
        }
           
        // Adding the current set of priors.
        for (int lib=0; lib<nlib; ++lib) {
            yptr[lib] += pptr[lib];
        }

        // Fitting a one-way layout.
        glm_one_group_vec(nlib, yptr, optr, dptr, wptr, maxit, tol, NA_REAL, &ocoef, &oconv);
        avecpm[tag]=(ocoef + LNmillion)/LNtwo;
    }

    R_Free(yptr);
    R_Free(optr);
    R_Free(wptr);
    R_Free(dptr);
    R_Free(pptr);
    
    return;
}
