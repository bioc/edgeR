#include "edgeR.h"

/* convert matrix to struct cmx */
cmx SEXPtocmx1 (SEXP obj)
{
    cmx res;
    if(TYPEOF(obj)==REALSXP){
        res.dmat=REAL(obj);
        res.isint=0;
    }
    else{
        res.imat=INTEGER(obj);
        res.isint=1;
    }

    SEXP ans;
    ans = getAttrib(obj, R_DimSymbol);
    int *dims = INTEGER(ans);
    res.nrow=dims[0], res.ncol=dims[1];
    res.type=0;

    return res;
}

/* convert CompressMatrix to struct cmx */
cmx SEXPtocmx2 (SEXP obj)
{
    cmx res;
    res.dmat=REAL(obj);
    res.isint=0;

    SEXP ans;
    char dim_cmx[] = "Dims";
    ans = getAttrib(obj, install(dim_cmx));
    int *dims = INTEGER(ans);
    res.nrow=dims[0], res.ncol=dims[1];

    char rep_row[] = "repeat.row", rep_col[] = "repeat.col";
    int is_rep_row, is_rep_col;
        
    ans = getAttrib(obj,install(rep_row));
    is_rep_row = asLogical(ans);  
    ans = getAttrib(obj,install(rep_col));
    is_rep_col = asLogical(ans);
    res.type=2*is_rep_row+is_rep_col;
    
    return res;
}

/* get a row vector from struct cmx */
void get_row (cmx *cmat, int k, double *y)
{
    int m=cmat->nrow, n=cmat->ncol;

    switch (cmat->type)
    {
    case 0:
        /*
        if(cmat->isint){
            for(int i=0;i<n;++i){
                R_xlen_t ii = (R_xlen_t)(m) * i+k;
                y[i]=(double)(cmat->imat)[ii];
            }
        }else{
            for(int i=0;i<n;++i){
                R_xlen_t ii = (R_xlen_t)(m) * i+k;
                y[i]=(cmat->dmat)[ii];
            }
        }        
        */

        // loop matrix using pointers
        if(cmat->isint){
            int *yptr; yptr = (cmat->imat)+k;
            for(int i=0;i<n;++i,yptr+=m){
                y[i]=(double)(*yptr);
            }
        }else{
            double *yptr; yptr = (cmat->dmat)+k;
            for(int i=0;i<n;++i,yptr+=m){
                y[i]=(*yptr);
            }
        }
        break;
    case 1:
        for(int i=0;i<n;++i) y[i]=(cmat->dmat)[k];
        break;
    case 2:
        for(int i=0;i<n;++i) y[i]=(cmat->dmat)[i];
        break;
    case 3:
        for(int i=0;i<n;++i) y[i]=(cmat->dmat)[0];
        break;
    }

    return;
}

/* get 3 or 4 rows from 3 or 4 cmx objects */
void get_row3 (cmx *cmat1, cmx *cmat2, cmx *cmat3, int k, double *y1, double *y2, double *y3)
{
    get_row(cmat1, k, y1);
    get_row(cmat2, k, y2);
    get_row(cmat3, k, y3);

    return;
}
void get_row4 (cmx *cmat1, cmx *cmat2, cmx *cmat3, cmx *cmat4, int k, double *y1, double *y2, double *y3, double *y4)
{
    get_row(cmat1, k, y1);
    get_row(cmat2, k, y2);
    get_row(cmat3, k, y3);
    get_row(cmat4, k, y4);

    return;
}

/* check whether a row of dmax from cmx is equal to 0 or 1 
 * work for dispersion and weights
 */
int check_row_scalar (cmx *cmat, int k, double y)
{
    int m=cmat->nrow, n=cmat->ncol, ans=1;

    switch (cmat->type)
    {
    case 0:;
        /*
        for(int i=0;i<n;++i){
            R_xlen_t ii = (R_xlen_t)(m) * i+k;
            if((cmat->dmat)[ii] != y)
            {
                ans=0;
                break;
            }
        }
        */
        double *yptr; yptr = (cmat->dmat)+k;
        for(int i=0;i<n;++i,yptr+=m){
            if((*yptr) != y)
            {
                ans=0;
                break;
            }
        }
        break;
    case 1:
        if((cmat->dmat)[k] != y) {
            ans=0;
        }
        break;
    case 2:
        for(int i=0;i<n;++i){
            if((cmat->dmat)[i] != y)
            {
                ans=0;
                break;
            }            
        }
        break;
    case 3:
        if((cmat->dmat)[0] != y){ 
            ans=0; 
        }
        break;
    }

    return ans; 
}

/* find the maximum of the cmx 
 * work for the count matrix y
 */

double max_cmx (cmx *cmat)
{
    int m=(cmat->nrow), n=(cmat->ncol);
    double ans=0;

    switch (cmat->type)
    {
    case 0:
        /*
        R_xlen_t nn = (R_xlen_t)(m) * n;
        if(cmat ->isint){
            int cmat_max = (cmat->imat)[0];
            for(R_xlen_t ii=0;ii<nn;++ii){
                if((cmat->imat)[ii]>cmat_max){
                    cmat_max=(cmat->imat)[ii];
                }
            }
            ans = (double)(cmat_max);
        }
        else
        {
            double cmat_max = (cmat->dmat)[0];
            for(R_xlen_t ii=0;ii<nn;++ii){
                if((cmat->dmat)[ii]>cmat_max){
                    cmat_max=(cmat->dmat)[ii];
                }
            }
            ans = cmat_max;
        }
        */
        if(cmat ->isint){
            int cmat_max = (cmat->imat)[0];
            int *cptr; cptr = (cmat->imat);
            for(int i=0;i<m;++i,++cptr){
                for(int j=0;j<n;++j,++cptr){
                    cmat_max=((*cptr)>cmat_max)? (*cptr):cmat_max;
                }
            }
            ans = (double)(cmat_max);
        }
        else
        {
            double cmat_max = (cmat->dmat)[0];
            double *cptr; cptr = (cmat->dmat);
            for(int i=0;i<m;++i,++cptr){
                for(int j=0;j<n;++j,++cptr){
                    cmat_max=((*cptr)>cmat_max)? (*cptr):cmat_max;
                }
            }
            ans = cmat_max;
        }
        break;
    case 1:
        ans = (cmat->dmat)[0];
        for(int i=0;i<m;++i){
            if((cmat->dmat)[i]>ans){
                ans=(cmat->dmat)[i];
            }
        }
        break;
    case 2:
        ans = (cmat->dmat)[0]; 
        for(int i=0;i<n;++i){
            if((cmat->dmat)[i] > ans)
            {
                ans=(cmat->dmat)[i];
            }            
        }
        break;
    case 3:
        ans = (cmat->dmat)[0]; 
        break;
    }

    return ans; 

}