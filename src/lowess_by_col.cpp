#include "edgeR.h"

extern "C" {

/* 'y' is a m by n matrix where each column corresponds to a different
 * dataset and each row corresponds to a point in 'x'. 'x' is a 'm'-long
 * vector containing sorted x-coordinates. 's' describes the span.
 */

SEXP lowess_by_col(SEXP x, SEXP y, SEXP n_cols, SEXP s) {
    // Setting up input data structures.
    PROTECT(x=AS_NUMERIC(x));
    PROTECT(y=AS_NUMERIC(y));
    const long total=LENGTH(x);
    const long span=INTEGER_VALUE(s);
    if (span>total) {
        error("Number of smoothing points should less than the total number of points.");
    } else if (span<=0) {
        error("Number of smoothing points should be positive.");
    }
    double* x_ptr=NUMERIC_POINTER(x);
    const int ncols=INTEGER_VALUE(n_cols);
    if (LENGTH(y)!=ncols*total) {
        error("Supplied dimensions for matrix 'y' are not consistent.");
    }
    std::deque<double*> y_ptrs(ncols);
    for (int i=0; i<ncols; ++i) {
        if (i==0) { 
            y_ptrs[i]=NUMERIC_POINTER(y);
        } else {
            y_ptrs[i]=y_ptrs[i-1]+total;
        }
    }

    // Setting up output vectors.
    SEXP output;
    PROTECT(output=NEW_LIST(2));
    SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, total, ncols));
    SET_VECTOR_ELT(output, 1, NEW_NUMERIC(total));
    std::deque<double*> f_ptrs(ncols);
    for (int i=0; i<ncols; ++i) {
        if (i==0) {
            f_ptrs[i]=NUMERIC_POINTER(VECTOR_ELT(output, 0));
        } else {
            f_ptrs[i]=f_ptrs[i-1]+total;
        }
    }
    double* w_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 1));

    /* First we determine which of the x-axis values are closest together. This means
     * that we go through all points to determine which 'frame' brings gets the closest
     * set of points i.e. the maximum distance to any point is at a minimum.  
     */
    long frame_end=span-1;
    for (long cur_p=0; cur_p<total; ++cur_p) {
        if (cur_p>frame_end) { frame_end=cur_p; }
        const double& cur_point=x_ptr[cur_p];
        double back_dist=cur_point-x_ptr[frame_end-span+1], front_dist=x_ptr[frame_end]-cur_point,
            max_dist=(back_dist > front_dist ? back_dist : front_dist);

        // Checking out the ideal place for the frame to be.
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
            back_dist=cur_point-x_ptr[frame_end-span+2];
            front_dist=x_ptr[frame_end+1]-cur_point;
            const double& next_max=(back_dist > front_dist ? back_dist : front_dist);
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

		/* If there aren't any columns, there isn't much point doing all the hoopla; so
		 * we'll just return the maximum distance in the 'weights' column and break.
		 */
		if (ncols==0) { 
			w_ptr[cur_p]=max_dist;
			continue; 
		}

        /* Now that we've located our optimal window, we can calculate the weighted average
         * across the points in the window (weighted according to distance from the current point).
         * and we can calculate the leverages. Unfortunately, we have to loop over the points in the 
         * window because each weight must be recomputed according to its new distance and new maximal
         * distance.
         */
        double total_weight=0;
        double& out_leverage=(w_ptr[cur_p]=-1);
        for (int i=0; i<ncols; ++i) { f_ptrs[i][cur_p]=0; }

        /* For non-zero maximum distances, we can compute the relative distance; otherwise, we set it to zero.
         * This means that all observations will have the same weight (user specifications aside). This makes
         * sense as they all lie on the same x-coordinate. Note that funny calculations might happen with the
         * leverage as there are multiple valid frames with the same minimum distance when many x-coordinates
         * are equal.
         */
        for (long m=frame_end-span+1; m<=frame_end; ++m) {
            const double rel_dist=(max_dist > low_value ? std::abs(x_ptr[m]-cur_point)/max_dist : 0);
            const double weight=std::pow(1-std::pow(rel_dist, 3.0), 3.0);
            total_weight+=weight;            
            
            for (int i=0; i<ncols; ++i) { f_ptrs[i][cur_p]+=weight*y_ptrs[i][m]; }
            if (m==cur_p) { out_leverage=weight; }
        }

        // Normalizing by the total weight.
        out_leverage/=total_weight;
        for (int i=0; i<ncols; ++i) { f_ptrs[i][cur_p]/=total_weight; }
    }

    UNPROTECT(3);
    return output;
}


}
