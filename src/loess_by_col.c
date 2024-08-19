#include "edgeR.h"

/* this function fits loess by column of y  
 *
 * inputs:
 * xptr covariate vector     
 * y    response matrix
 * span width of smoothing window
 * 
 * outputs:
 * yy   smoothed response matrix
 * lv   leverage vector
 * 
 * comments:
 * 1. This is converted from R_loess_by_col.cpp written by Aaron
 * 2. The code is simplified without pointer
 * 3. sort x in C and simplify R code
*/

/* 'y' is a m by n matrix where each column corresponds to a different
 * dataset and each row corresponds to a point in 'x'. 'x' is a 'm'-long
 * vector containing sorted x-coordinates. 's' describes the span.
 */

void loess_by_column(double *xptr, cmx *y, int span, double *yy, double *lv) 
{
    int total=(y->nrow), ncols=(y->ncol);
    double low_value=1e-10;
    double *yptr, *yypt;

    /* sort x with rsort_with_index */
	int *ind = R_Calloc(total, int);
	for(int i=0;i<total;++i) ind[i]=i;
	rsort_with_index(xptr,ind,total);

    /* First we determine which of the x-axis values are closest together. This means
     * that we go through all points to determine which 'frame' brings gets the closest
     * set of points i.e. the maximum distance to any point is at a minimum. The frame
     * represents 50% of the points (i.e. the span) and the boundaries of the frame can
     * be used to compute the minimum distance to the current point.
     */
    int frame_end=span-1;
    for (int cur_p=0; cur_p<total; ++cur_p) {
        if (cur_p>frame_end) { 
            frame_end=cur_p; 
        }
        double back_dist  = xptr[cur_p]-xptr[frame_end-span+1];
        double front_dist = xptr[frame_end]-xptr[cur_p];
        double max_dist   = fmax2(back_dist, front_dist);

        while (frame_end < total-1 && cur_p+span-1>frame_end) {
            /* Every time we advance, we twiddle with the ends of the frame to see if we can't get 
                * a better fit. The frame will always advance in the forward direction. This is because the 
                * current frame is optimal with respect to the previous tag. If the previous maximal distance 
                * was at the back, shifting the frame backward will increase the back distance with respect to 
                * the current tag (and thus increase the maximal distance).
                *
                * If the previous maximal distance was at the front, shifting the frame backward may 
                * decrease the front distance with respect to the current tag. However, we note that 
                * because of optimality, having a previous maximal distance at the front must mean
                * that a back-shifted frame will result in an even larger previous maximal distance at 
                * the back (otherwise the optimal frame would be located further back to start with). In 
                * short, shifting the frame backwards will flip the maximal distance to that of the back
                * distance which is even larger than the non-shifted forward distance.
                *
                * Thus, the frame can only go forwards. Note that below, the frame is defined by 
                * the 'end' position which notes the end point of the current frame. The start
                * point is inherently defined by revolving around the minimum point.
                */
            back_dist  = xptr[cur_p]-xptr[frame_end-span+2];
            front_dist = xptr[frame_end+1]-xptr[cur_p];
            double next_max = fmax2(back_dist, front_dist);

            /* This bit provides some protection against near-equal values, by forcing the frame
                * forward provided that the difference between the lowest maximum distance and
                * the maximum distance at any other frame is less than a low_value. This ensures
                * that values following a stretch of identical x-coordinates are accessible
                * to the algorithm (rather than being blocked off by inequalities introduced by
                * double imprecision).
                */
            const double diff=(next_max-max_dist)/max_dist;
            if (diff > low_value) { 
                break; 
            } else if (diff < 0) {
                max_dist=next_max;                                
            }
            ++frame_end;
        }

        /* Now that we've located our optimal window, we can calculate the weighted average
            * across the points in the window (weighted according to distance from the current point).
            * and we can calculate the leverages. Unfortunately, we have to loop over the points in the 
            * window because each weight must be recomputed according to its new distance and new maximal
            * distance.
            */
        double total_weight=0;
        
        /*
        for (int i=0; i<ncols; ++i) {
            R_xlen_t ii = (R_xlen_t)(total)*i+ind[cur_p];
            yy[ii]=0; 
        }
        */

        yypt=yy+ind[cur_p];
        for (int i=0; i<ncols; ++i, yypt+=total) {
            (*yypt)=0; 
        }

        /* For non-zero maximum distances, we can compute the relative distance; otherwise, we set it to zero.
            * This means that all observations will have the same weight (user specifications aside). This makes
            * sense as they all lie on the same x-coordinate. Note that funny calculations might happen with the
            * leverage as there are multiple valid frames with the same minimum distance when many x-coordinates
            * are equal. 
            *
            * Note that we have to look for more than just the 'span' number of points. Consider the series
            * A,B,C,C where each is a value and A < B < C and C - B > B - A. The algorithm above will move the 
            * frame to [1,3] when calculating the maximum distance for B. This is the same as [0, 2] in terms
            * of distance, but only using the frame components to calculate the mean will miss out on element 0.
            * So, the computation should work from [0, 3]. There's no need to worry about the extra 'C' as it
            * will have weight zero.
            */
        for (int m=frame_end; m>=0; --m) { 
            double rel_dist = (max_dist > low_value) ? fabs(xptr[cur_p]-xptr[m])/max_dist : 0;
            double weight   = fcube(1.0-fcube(rel_dist));

            if (weight < 0) { 
                continue; 
            }

            total_weight+=weight;

            /*
            for (int i=0; i<ncols; ++i) {
                R_xlen_t ii = (R_xlen_t)(total)*i+ind[cur_p];
                R_xlen_t jj = (R_xlen_t)(total)*i+ind[m];
                yy[ii]     += weight*yptr[jj]; 
            }
            */            
            
            yptr=(y->dmat)+ind[m];
            yypt=yy+ind[cur_p];
            for (int i=0; i<ncols; ++i, yptr+=total, yypt+=total) {
                (*yypt) += weight*(*yptr); 
            }

            if (m==cur_p) { 
                lv[ind[cur_p]]=weight; 
            }
        }

        // Normalizing by the total weight.
        lv[ind[cur_p]]/=total_weight;

        /*
        for (int i=0; i<ncols; ++i) {
            R_xlen_t ii = (R_xlen_t)(total)*i+ind[cur_p];
            yy[ii]/=total_weight; 
        }
        */
        
        yypt=yy+ind[cur_p];
        for (int i=0; i<ncols; ++i, yypt+=total) {
            (*yypt)/=total_weight; 
        }
    }
	
    R_Free(ind);
    
    return;
}
