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


/* per-thread scratch for calc_cpm_log: offset/prior outputs and the two
 * buffers compute_offsets needs */
typedef struct {
    double *lptr, *pptr, *so, *sp;
} clog_ws;

void calc_cpm_log(cmx *y, cmx *libsizes, cmx *priors, double *cpm, int nthreads)
{
    int ntag=(y->nrow), nlib=(y->ncol);
    int log_in=0, log_out=1;

    const double LNtwo = log(2), LNmillion = log(1e6);

    int nth = clamp_threads(nthreads);

    clog_ws *ws = R_Calloc(nth, clog_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].lptr = R_Calloc(nlib, double);
        ws[t].pptr = R_Calloc(nlib, double);
        ws[t].so   = R_Calloc(nlib, double);
        ws[t].sp   = R_Calloc(nlib, double);
    }

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

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for(int tag=0;tag<ntag;++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        clog_ws *w = &ws[tid];

        compute_offsets(priors,libsizes,tag,log_in,log_out,w->pptr,w->lptr,w->so,w->sp);

        double *cptr=cpm+tag;
        for(int lib=0;lib<nlib;++lib,cptr+=ntag)
        {
            (*cptr)+= w->pptr[lib];
            (*cptr) = ((*cptr)>0)? (log(*cptr)-w->lptr[lib]+LNmillion)/LNtwo : R_NaN;
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].lptr);
        R_Free(ws[t].pptr);
        R_Free(ws[t].so);
        R_Free(ws[t].sp);
    }
    R_Free(ws);

    return;
}

/* function to compute raw cpm
 *
 * inputs:
 * y        matrix of fitted values
 * libsizes compressedMatrix of library sizes
 *
 * output:
 * cpm      matrix of cpm, scaled in place by one million / library size
 *
 * comment:
 * no return value; results are written into the caller-owned cpm array
 */
void calc_cpm_raw(cmx *y, cmx *libsizes, double *cpm, int nthreads)
{
    const double one_million=1e6;

    int ntag=(y->nrow), nlib=(y->ncol);

    int nth = clamp_threads(nthreads);

    /* one libsize row buffer per thread */
    double **lptr_t = R_Calloc(nth, double*);
    for(int t=0;t<nth;++t)
    {
        lptr_t[t] = R_Calloc(nlib, double);
    }

    /*
    for(int tag=0;tag<ntag;++tag){
        get_row(libsizes,tag,lptr);

        for(int lib=0;lib<nlib;++lib){
            R_xlen_t ii = (R_xlen_t)(ntag) * lib + tag;
            cpm[ii] = cpm[ii]*one_million/lptr[lib];
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
        double *lptr = lptr_t[tid];
        get_row(libsizes,tag,lptr);
        double *cptr=cpm+tag;
        for(int lib=0;lib<nlib;++lib,cptr+=ntag)
        {
            (*cptr) = (*cptr)*one_million/lptr[lib];
        }
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(lptr_t[t]);
    }
    R_Free(lptr_t);

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
/* per-thread scratch for average_log_cpm */
typedef struct {
    double *yptr, *optr, *wptr, *dptr, *pptr;   /* nlib row buffers          */
    double *so, *sp;                            /* compute_offsets scratch    */
} avecpm_ws;

void average_log_cpm(cmx *y, cmx *offsets, cmx *priors, cmx *disp, cmx *weights, int maxit, double tol, double *avecpm, int nthreads)
{
    const double LNmillion=log(1e6), LNtwo=log(2);

    int ntag = (y->nrow), nlib = (y->ncol);
    int log_in=1, log_out=1;

    int nth = clamp_threads(nthreads);

    /* row vectors per thread (optr/pptr are R_Calloc'd to zero, matching the
     * serial pre-loop state used by the first tag when offsets/priors are not
     * row-repeated) */
    avecpm_ws *ws = R_Calloc(nth, avecpm_ws);
    for(int t=0;t<nth;++t)
    {
        ws[t].yptr = R_Calloc(nlib,double);
        ws[t].optr = R_Calloc(nlib,double);
        ws[t].wptr = R_Calloc(nlib,double);
        ws[t].dptr = R_Calloc(nlib,double);
        ws[t].pptr = R_Calloc(nlib,double);
        ws[t].so   = R_Calloc(nlib,double);
        ws[t].sp   = R_Calloc(nlib,double);
    }

    /* whether both offsets and priors are row repeated */
    int repeat_row = ((offsets->type) >= 2) && ((priors->type) >=2);

    if(repeat_row)
    {
        int tag_start = 0;
        /* compute the constant prior/offset once and replicate to every thread */
        compute_offsets(priors,offsets,tag_start,log_in,log_out,ws[0].pptr,ws[0].optr,ws[0].so,ws[0].sp);
        for(int t=1;t<nth;++t)
        {
            for(int lib=0;lib<nlib;++lib)
            {
                ws[t].pptr[lib]=ws[0].pptr[lib];
                ws[t].optr[lib]=ws[0].optr[lib];
            }
        }
    }

    // Returning average log-cpm
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(static)
    #endif
    for (int tag=0; tag<ntag; ++tag)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        avecpm_ws *w = &ws[tid];
        double ocoef;
        int oconv;
        get_row3(y,disp,weights,tag,w->yptr,w->dptr,w->wptr);

        if(!repeat_row)
        {
            compute_offsets(priors,offsets,tag,log_in,log_out,w->pptr,w->optr,w->so,w->sp);
        }

        // Adding the current set of priors.
        for (int lib=0; lib<nlib; ++lib)
        {
            w->yptr[lib] += w->pptr[lib];
        }

        // Fitting a one-way layout.
        glm_one_group_vec(nlib, w->yptr, w->optr, w->dptr, w->wptr, maxit, tol, NA_REAL, &ocoef, &oconv);
        avecpm[tag]=(ocoef + LNmillion)/LNtwo;
    }

    for(int t=0;t<nth;++t)
    {
        R_Free(ws[t].yptr);
        R_Free(ws[t].optr);
        R_Free(ws[t].wptr);
        R_Free(ws[t].dptr);
        R_Free(ws[t].pptr);
        R_Free(ws[t].so);
        R_Free(ws[t].sp);
    }
    R_Free(ws);

    return;
}
