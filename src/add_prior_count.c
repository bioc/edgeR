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

void add_prior_count_vec(cmx *y, cmx *offsets, cmx *priors, double *yy, double *offset, int nthreads)
{
    int ntag=(y->nrow), nlib=(y->ncol), tag_start=0;
    int log_in=1, log_out=1;

    int nth = clamp_threads(nthreads);

    /* prior is constant across tags here: compute it once (serial), then read
     * it inside the parallel loop. so/sp are scratch for compute_offsets. */
    double *pptr = R_Calloc(nlib, double);
    double *so   = R_Calloc(nlib, double);
    double *sp   = R_Calloc(nlib, double);

    /* compute adjusted prior and offset */
    compute_offsets(priors,offsets,tag_start,log_in,log_out,pptr,offset,so,sp);
    R_Free(so);
    R_Free(sp);

    /* one y row buffer per thread */
    double **yptr_t = R_Calloc(nth, double*);
    for(int t=0;t<nth;++t)
    {
        yptr_t[t] = R_Calloc(nlib, double);
    }

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

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0;tag<ntag;++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *yptr = yptr_t[tid];
        get_row(y,tag,yptr);
        double *yypt=yy+tag;
        for(int lib=0;lib<nlib;++lib,yypt+=ntag)
        {
            (*yypt) = yptr[lib]+pptr[lib];
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(yptr_t[t]);
    }
    R_Free(yptr_t);
    R_Free(pptr);

    return;
}

/* per-thread scratch for add_prior_count_mat: y row, prior/offset outputs and
 * the two buffers compute_offsets needs */
typedef struct
{
    double *yptr, *pptr, *optr, *so, *sp;  /* y row; adjusted prior; adjusted offset; two compute_offsets scratch buffers */
} apc_ws;

/* Add prior counts to a count matrix per gene, parallelised over tags.
 * The prior and offset are recomputed for each tag, so the adjusted offset is
 * written as a full matrix (contrast add_prior_count_vec, which shares one
 * prior across all tags).
 *
 * inputs:
 *   y        count matrix (compressed matrix)
 *   offsets  compressed matrix of offsets
 *   priors   compressed matrix of prior counts
 *
 * outputs (caller-allocated):
 *   yy      adjusted count matrix (ntag x nlib)
 *   offset  adjusted offset matrix (ntag x nlib)
 *
 * No return value; writes only into the caller-owned output arrays.
 */
void add_prior_count_mat(cmx *y, cmx *offsets, cmx *priors, double *yy, double *offset, int nthreads)
{
    int ntag = (y->nrow), nlib = (y->ncol);
    int log_in=1, log_out=1;

    int nth = clamp_threads(nthreads);

    apc_ws *ws = R_Calloc(nth, apc_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].yptr = R_Calloc(nlib, double);
        ws[t].pptr = R_Calloc(nlib, double);
        ws[t].optr = R_Calloc(nlib, double);
        ws[t].so   = R_Calloc(nlib, double);
        ws[t].sp   = R_Calloc(nlib, double);
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
        apc_ws *w = &ws[tid];
        get_row(y,tag,w->yptr);
        compute_offsets(priors,offsets,tag,log_in,log_out,w->pptr,w->optr,w->so,w->sp);

        double *yypt=yy+tag, *oopt=offset+tag;
        for(int lib=0;lib<nlib;++lib,yypt+=ntag,oopt+=ntag)
        {
            (*yypt) = w->yptr[lib]+w->pptr[lib];
            (*oopt) = w->optr[lib];
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].pptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].so);
        R_Free(ws[t].sp);
    }
    R_Free(ws);

    return;
}

/* compute prior and library sizes
 * this is converted from add_prior::compute()
 *
 * optr and pptr are caller-owned nlib scratch buffers, so this helper does
 * no allocation of its own and is safe to call inside a parallel region
 */
void compute_offsets (cmx *priors, cmx *offsets, int tag, int log_in, int log_out, double *prior, double *offset, double *optr, double *pptr)
{
    int nlib = (priors->ncol);

    /* get row tag from lib sizes, if logged, recover it*/
    get_row(offsets, tag, optr);
    if(log_in)
    {
        for(int lib=0;lib<nlib;++lib)
        {
            offset[lib]=exp(optr[lib]);
        }
    }
    else
    {
        for(int lib=0;lib<nlib;++lib)
        {
            offset[lib]=optr[lib];
        }
    }

    /* compute average lib size*/
    double ave_lib=0;
    for(int lib=0;lib<nlib;++lib)
    {
        ave_lib += offset[lib];
    }
    ave_lib=ave_lib/nlib;

    /* compute the adjusted prior count for each library */
    get_row(priors,tag, pptr);
    for(int lib=0;lib<nlib;++lib)
    {
        prior[lib]=pptr[lib]*offset[lib]/ave_lib;
    }

    /* add it twice back to the library sizes */
    for(int lib=0;lib<nlib;++lib)
    {
        offset[lib]+=2*prior[lib];
    }
    if(log_out)
    {
        for(int lib=0;lib<nlib;++lib)
        {
            offset[lib]=log(offset[lib]);
        }
    }

    return;
}