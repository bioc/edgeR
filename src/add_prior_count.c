#include "edgeR.h"

/* this function computes the adjusted offsets, adjusted counts by prior
 *
 * inputs:
 * y       count matrix
 * offsets compressedMatrix 
 * priors  compressedMatrix
 * 
 * outputs:
 * yy      adjusted count matrix
 * offset  vector or matrix 
 * 
 * comment:
 * This is converted from R_add_count.cpp written by Aaron
*/

void add_prior_count_vec(cmx *y, cmx *offsets, cmx *priors, double *yy, double *offset)
{
    int ntag=(y->nrow), nlib=(y->ncol), tag_start=0;
    int log_in=1, log_out=1;

    double *yptr = R_Calloc(nlib, double);      
    double *pptr = R_Calloc(nlib, double);

    /* compute adjusted prior and offset */
    compute_offsets(priors,offsets,tag_start,log_in,log_out,pptr,offset);

    /* add prior to the count matrix */
    
    /*
    for(int tag=0;tag<ntag;++tag){
        get_row(y,tag,yptr);
        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag)*lib+tag;
            yy[ii]      = yptr[lib]+pptr[lib];
        }
    }
    */
    
    // effificent way to loop matrix
    double *yypt; 
    for(int tag=0;tag<ntag;++tag,++yy){
        get_row(y,tag,yptr);
        yypt=yy;
        for(int lib=0;lib<nlib;++lib,yypt+=ntag){
           (*yypt) = yptr[lib]+pptr[lib];
        }
    }

    R_Free(yptr);
    R_Free(pptr);

    return;
}

void add_prior_count_mat(cmx *y, cmx *offsets, cmx *priors, double *yy, double *offset)
{
    int ntag = (y->nrow), nlib = (y->ncol);
    int log_in=1, log_out=1;

    double *yptr = R_Calloc(nlib, double);      
    double *pptr = R_Calloc(nlib, double);
    double *optr = R_Calloc(nlib, double);

    double *yypt, *oopt;
    for(int tag=0;tag<ntag;++tag,++yy,++offset){
        get_row(y,tag,yptr);
        compute_offsets(priors,offsets,tag,log_in,log_out,pptr,optr);

        yypt=yy, oopt=offset;
        for(int lib=0;lib<nlib;++lib,yypt+=ntag,oopt+=ntag){
            (*yypt) = yptr[lib]+pptr[lib];
            (*oopt) = optr[lib];
        }
    }

    R_Free(yptr);
    R_Free(pptr);
    R_Free(optr);

    return;
}

/* compute prior and library sizes
 * this is converted from add_prior::compute()
 */
void compute_offsets (cmx *priors, cmx *offsets, int tag, int log_in, int log_out, double *prior, double *offset)
{
    int nlib = (priors->ncol);

    double *optr = R_Calloc(nlib, double);
    double *pptr = R_Calloc(nlib, double);

    /* get row tag from lib sizes, if logged, recover it*/
    get_row(offsets, tag, optr);
    if(log_in){
        for(int lib=0;lib<nlib;++lib) offset[lib]=exp(optr[lib]);
    }else{
        for(int lib=0;lib<nlib;++lib) offset[lib]=optr[lib];
    }

    /* compute average lib size*/
    double ave_lib=0;
    for(int lib=0;lib<nlib;++lib) ave_lib += offset[lib];
    ave_lib=ave_lib/nlib;

    /* compute the adjusted prior count for each library */
    get_row(priors,tag, pptr);
    for(int lib=0;lib<nlib;++lib) prior[lib]=pptr[lib]*offset[lib]/ave_lib;

    /* add it twice back to the library sizes */
    for(int lib=0;lib<nlib;++lib) offset[lib]+=2*prior[lib];
    if(log_out){
        for(int lib=0;lib<nlib;++lib) offset[lib]=log(offset[lib]);
    }

    R_Free(optr);
    R_Free(pptr);
    
    return;
}