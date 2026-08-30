#include "edgeR.h"

/* convert matrix to struct cmx
 *
 * input:
 *   obj  a numeric (REALSXP) or integer (INTSXP) R matrix
 *
 * output:
 *   a cmx whose dmat or imat aliases obj's data (no copy), with isint set
 *   accordingly, nrow/ncol taken from the dim attribute, and type 0
 *   (common, uncompressed matrix)
 *
 * return value: the populated cmx (returned by value)
 * side effects: none; obj is read only, not modified
 */
cmx SEXPtocmx1 (SEXP obj)
{
    cmx res;
    if(TYPEOF(obj)==REALSXP)
    {
        res.dmat=REAL(obj);
        res.isint=0;
    }
    else
    {
        res.imat=INTEGER(obj);
        res.isint=1;
    }

    SEXP ans;
    ans = getAttrib(obj, R_DimSymbol);
    int *dims = INTEGER(ans);
    res.nrow=dims[0];
    res.ncol=dims[1];
    res.type=0;
    res.pnrow=res.nrow;
    res.rowsel=NULL;

    return res;
}

/* convert CompressMatrix to struct cmx
 *
 * input:
 *   obj  a compressedMatrix stored as REALSXP, carrying the attributes
 *        "Dims", "repeat.row" and "repeat.col"
 *
 * output:
 *   a cmx whose dmat aliases obj's data (no copy), isint 0, nrow/ncol from
 *   the "Dims" attribute, and type = 2*repeat.row + repeat.col encoding the
 *   compression scheme (0 common, 1 by row, 2 by column, 3 by both)
 *
 * return value: the populated cmx (returned by value)
 * side effects: none; obj is read only, not modified
 */
cmx SEXPtocmx2 (SEXP obj)
{
    cmx res;
    res.dmat=REAL(obj);
    res.isint=0;

    SEXP ans;
    char dim_cmx[] = "Dims";
    ans = getAttrib(obj, install(dim_cmx));
    int *dims = INTEGER(ans);
    res.nrow=dims[0];
    res.ncol=dims[1];

    char rep_row[] = "repeat.row", rep_col[] = "repeat.col";
    int is_rep_row, is_rep_col;
        
    ans = getAttrib(obj,install(rep_row));
    is_rep_row = asLogical(ans);  
    ans = getAttrib(obj,install(rep_col));
    is_rep_col = asLogical(ans);
    res.type=2*is_rep_row+is_rep_col;
    res.pnrow=res.nrow;
    res.rowsel=NULL;

    return res;
}

/* build a cmx view over a caller-owned buffer
 *
 * inputs:
 *   dmat   pointer to column-major doubles (aliased, not copied), or NULL
 *   imat   pointer to column-major ints    (aliased, not copied), or NULL
 *   nrow   number of rows
 *   ncol   number of columns
 *   type   compression code: 0 common, 1 by row, 2 by column, 3 by both
 *   isint  1 if the data is integer (read imat), 0 if double (read dmat)
 *
 * output:
 *   a cmx carrying both pointers, the dims, type and isint verbatim
 *
 * return value: the populated cmx (returned by value)
 * side effects: none; the buffers are aliased, not owned
 */
cmx make_cmx (double *dmat, int *imat, int nrow, int ncol, int type, int isint)
{
    cmx res;
    res.dmat  = dmat;
    res.imat  = imat;
    res.nrow  = nrow;
    res.ncol  = ncol;
    res.type  = type;
    res.isint = isint;
    res.pnrow = nrow;
    res.rowsel = NULL;

    return res;
}

/* build a zero-copy row-subset VIEW of a full cmx (src.rowsel must be NULL).
 *
 * inputs:
 *   src    a full cmx (not itself a subset)
 *   rows   length-nrows physical row indices (0-based), caller-owned and outliving the view
 *   nrows  the subset size (becomes the view's logical nrow)
 *
 * output:
 *   a cmx sharing src's data/type/isint and physical stride (pnrow), whose logical
 *   nrow is nrows and whose row access is remapped through rows[] by the accessors.
 *
 * side effects: none; rows and src's data are aliased, not copied
 */
cmx subset_cmx (cmx src, const int *rows, int nrows)
{
    src.nrow   = nrows;
    src.rowsel = rows;   /* src.pnrow (physical stride) and data pointers are kept */
    return src;
}

/* turn a logical indicator into physical row indices for subset_cmx.
 *
 * inputs:
 *   keep  length-n array; a row is selected where keep[i] != 0
 *   n     length of keep
 *   idx   caller-owned output buffer of length >= n
 *
 * output: idx[0..count-1] holds the positions where keep is true
 * return value: the number of selected rows (count)
 */
int which_index (const int *keep, int n, int *idx)
{
    int c = 0;
    for(int i=0;i<n;++i)
    {
        if(keep[i]) idx[c++] = i;
    }
    return c;
}

/* get a row vector from struct cmx
 *
 * inputs:
 *   cmat  source compressed matrix
 *   k     zero-based row index to extract
 *   y     caller-owned output buffer of length cmat->ncol
 *
 * output:
 *   y[0..ncol-1] is filled with row k expanded to double according to the
 *   compression type (integer data is cast to double)
 *
 * return value: none
 * side effects: writes y; cmat is read only
 */
void get_row (cmx *cmat, int k, double *y)
{
    int m=cmat->pnrow, n=cmat->ncol;
    int kk = cmat->rowsel ? cmat->rowsel[k] : k;

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
        if(cmat->isint)
        {
            int *yptr;
            yptr = (cmat->imat)+kk;
            for(int i=0;i<n;++i,yptr+=m)
            {
                y[i]=(double)(*yptr);
            }
        }
        else
        {
            double *yptr;
            yptr = (cmat->dmat)+kk;
            for(int i=0;i<n;++i,yptr+=m)
            {
                y[i]=(*yptr);
            }
        }
        break;
    case 1:
        for(int i=0;i<n;++i)
        {
            y[i]=(cmat->dmat)[kk];
        }
        break;
    case 2:
        for(int i=0;i<n;++i)
        {
            y[i]=(cmat->dmat)[i];
        }
        break;
    case 3:
        for(int i=0;i<n;++i)
        {
            y[i]=(cmat->dmat)[0];
        }
        break;
    }

    return;
}

/* get 3 or 4 rows from 3 or 4 cmx objects
 *
 * wrapper over get_row: extracts row k from cmat1, cmat2 and cmat3 into the
 * caller-owned buffers y1, y2 and y3 (each of length ncol).
 *
 * return value: none
 * side effects: writes y1, y2, y3; the cmx inputs are read only
 */
void get_row3 (cmx *cmat1, cmx *cmat2, cmx *cmat3, int k, double *y1, double *y2, double *y3)
{
    get_row(cmat1, k, y1);
    get_row(cmat2, k, y2);
    get_row(cmat3, k, y3);

    return;
}

/* get the same row from 4 cmx objects
 *
 * wrapper over get_row: extracts row k from cmat1..cmat4 into the caller-owned
 * buffers y1..y4 (each of length ncol).
 *
 * return value: none
 * side effects: writes y1, y2, y3, y4; the cmx inputs are read only
 */
void get_row4 (cmx *cmat1, cmx *cmat2, cmx *cmat3, cmx *cmat4, int k, double *y1, double *y2, double *y3, double *y4)
{
    get_row(cmat1, k, y1);
    get_row(cmat2, k, y2);
    get_row(cmat3, k, y3);
    get_row(cmat4, k, y4);

    return;
}

/* check whether a row of dmax from cmx is equal to 0 or 1 
 * work for offsets, dispersion and weights
 *
 * inputs:
 *   cmat  source compressed matrix
 *   k     zero-based row index to test
 *   y     scalar to compare every entry of row k against
 *
 * return value: 1 if every entry of row k equals y, 0 otherwise
 * side effects: none; cmat is read only
 */
int check_row_scalar (cmx *cmat, int k, double y)
{
    int m=cmat->pnrow, n=cmat->ncol, ans=1;
    int kk = cmat->rowsel ? cmat->rowsel[k] : k;

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
        double *yptr;
        yptr = (cmat->dmat)+kk;
        for(int i=0;i<n;++i,yptr+=m)
        {
            if((*yptr) != y)
            {
                ans=0;
                break;
            }
        }
        break;
    case 1:
        if((cmat->dmat)[kk] != y)
        {
            ans=0;
        }
        break;
    case 2:
        for(int i=0;i<n;++i)
        {
            if((cmat->dmat)[i] != y)
            {
                ans=0;
                break;
            }
        }
        break;
    case 3:
        if((cmat->dmat)[0] != y)
        {
            ans=0;
        }
        break;
    }

    return ans; 
}

/* find the maximum of the cmx 
 * work for the count matrix y
 *
 * input:
 *   cmat  source compressed matrix
 *
 * return value: the largest entry of cmat as a double (integer data is cast
 *               to double)
 * side effects: none; cmat is read only
 */

double max_cmx (cmx *cmat)
{
    int m=(cmat->pnrow), n=(cmat->ncol), nr=(cmat->nrow);
    const int *sel=cmat->rowsel;
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
        if(cmat ->isint)
        {
            const int *d = cmat->imat;
            int cmat_max = sel ? d[sel[0]] : d[0];
            if(sel)
            {
                for(int k=0;k<nr;++k)
                {
                    int kk=sel[k];
                    for(int j=0;j<n;++j)
                    {
                        int v=d[kk + (R_xlen_t)j*m];
                        if(v>cmat_max) cmat_max=v;
                    }
                }
            }
            else
            {
                const int *cptr = d;
                for(int i=0;i<m;++i)
                {
                    for(int j=0;j<n;++j,++cptr)
                    {
                        if((*cptr)>cmat_max) cmat_max=(*cptr);
                    }
                }
            }
            ans = (double)(cmat_max);
        }
        else
        {
            const double *d = cmat->dmat;
            double cmat_max = sel ? d[sel[0]] : d[0];
            if(sel)
            {
                for(int k=0;k<nr;++k)
                {
                    int kk=sel[k];
                    for(int j=0;j<n;++j)
                    {
                        double v=d[kk + (R_xlen_t)j*m];
                        if(v>cmat_max) cmat_max=v;
                    }
                }
            }
            else
            {
                const double *cptr = d;
                for(int i=0;i<m;++i)
                {
                    for(int j=0;j<n;++j,++cptr)
                    {
                        if((*cptr)>cmat_max) cmat_max=(*cptr);
                    }
                }
            }
            ans = cmat_max;
        }
        break;
    case 1:
        if(sel)
        {
            ans = (cmat->dmat)[sel[0]];
            for(int k=0;k<nr;++k)
            {
                if((cmat->dmat)[sel[k]]>ans) ans=(cmat->dmat)[sel[k]];
            }
        }
        else
        {
            ans = (cmat->dmat)[0];
            for(int i=0;i<m;++i)
            {
                if((cmat->dmat)[i]>ans) ans=(cmat->dmat)[i];
            }
        }
        break;
    case 2:
        ans = (cmat->dmat)[0];
        for(int i=0;i<n;++i)
        {
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