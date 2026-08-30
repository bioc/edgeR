#include "edgeR.h"

/* This function is converted from interpolant.cpp / interpolant.h
 * and R_maximize_interpolant.cpp by Aaron
 *
 * The code is simplified by Lizhong
 */

/* Splines a la Forsythe Malcolm and Moler (from splines.c in the stats package).
 *
 * 'n' is the number of points, 'x' is the vector of x-coordinates, 'y' is the
 * vector of y-coordinates (unchanged and equal to the constant for the interpolating
 * cubic spline upon output), 'b' is the coefficient degree 1, 'c' is the coefficient
 * of degree '2', 'd' is the coefficient degree 3. We need to access the coefficients
 * directly, which makes evaluation from R undesirable.
 */


/************************************
 *
 * It fits the spline and grabs the coefficients of each segment.
 * It then calculates the maxima at the segments neighbouring
 * the maximally highest point. This avoids NR optimisation
 * as well as the need to call R's splinefun's from within C.
 *
 ***********************************/

/* find_max: return the abscissa that maximises the interpolating cubic spline of (x, y).
 *
 * Fits the natural cubic spline through the npts knots, then inspects the two segments
 * neighbouring the highest knot and solves each segment cubic for its interior
 * stationary point, keeping whichever yields the largest interpolated value. This
 * avoids Newton-Raphson optimisation and calling R's splinefun from C.
 *
 * inputs:
 * npts   number of interpolation points / knots
 * x      abscissas of the knots, strictly increasing (length npts)
 * y      ordinates of the knots (length npts, unchanged on output)
 *
 * outputs (caller-allocated scratch, overwritten):
 * b      degree-1 spline coefficients (length npts)
 * c      degree-2 spline coefficients (length npts)
 * d      degree-3 spline coefficients (length npts)
 *
 * return value:
 * the x-coordinate at which the fitted spline attains its maximum
 *
 * side effects: fills b, c, d via F77_CALL(spline); does not modify x or y.
 */
static double find_max (int npts, double *x, double *y,double *b, double *c, double *d)
{
    double maxed = -1;
    int maxed_at = -1;

    for (int i=0; i<npts; ++i)
    {
        // Getting a good initial guess for the MLE.
        if (maxed_at < 0 || y[i] > maxed)
        {
            maxed=y[i];
            maxed_at=i;
        }
    }
    double x_max=x[maxed_at];

    /* spline: Forsythe-Malcolm-Moler cubic spline coefficients (spline.f, the Fortran
     * ancestor of R stats' splines.c). Fills b, c, d with the degree-1/2/3 coefficients
     * of the natural interpolating cubic spline through (x, y).
     * Called locally so find_max can read the per-segment coefficients directly and
     * locate the interpolant maximum without R's splinefun.
     * Args -> Fortran: &npts=n (in); x,y=knots (in); b,c,d=coeffs (out).
     * spline.f is edgeR's own bundled Fortran (the FMM natural cubic spline, the same
     * algorithm as stats::splinefun(..., "fmm")); it is not called back through R here.
     * Per-argument mapping: &npts -> n, the number of knots (input, unchanged); x -> x(n)
     * and y -> y(n), the input abscissas (strictly increasing) and ordinates (both left
     * unchanged, y also being the degree-0 coefficient); b -> b(n), c -> c(n), d -> d(n)
     * are the output degree-1/2/3 coefficient arrays, which spline.f also uses in place as
     * the tridiagonal-solve workspace, so no separate workspace or derivative array passes.
     * Result: on the two segments flanking the highest knot these coefficients let find_max
     * solve each cubic's stationary point in closed form and keep the interpolant maximum,
     * i.e. the optimal dispersion location returned via max_interpolant.
     * Equivalent R operation: coefficients of stats::splinefun(x, y, "fmm").
     * Netlib references: Forsythe, Malcolm & Moler (1977); netlib fmm/spline.f.
     */
    F77_CALL(spline)(&npts, x, y, b, c, d);

    // First we have a look at the segment on the left and see if it contains the maximum.
    if (maxed_at>0)
    {
        double ld=d[maxed_at-1];
        double lc=c[maxed_at-1];
        double lb=b[maxed_at-1];

        double delta = fsquare(lc)-3*ld*lb;
        int solvable = (delta < 0)? 0 : 1;

        if (solvable)
        {
            /* Using the solution with the maximum (not minimum). If the curve is mostly increasing, the
             * maximal point is located at the smaller solution (i.e. sol1 for a>0). If the curve is mostly
             * decreasing, the maximal is located at the larger solution (i.e., sol1 for a<0).
             */
            double numerator = -lc-sqrt(delta);
            double chosen_sol=(ld!=0)? numerator/(3*ld) : -1.0;

            /* The spline coefficients are designed such that 'x' in 'y + b*x + c*x^2 + d*x^3' is
             * equal to 'x_t - x_l' where 'x_l' is the left limit of that spline segment and 'x_t'
             * is where you want to get an interpolated value. This is necessary in 'splinefun' to
             * ensure that you get 'y' (i.e. the original data point) when 'x=0'. For our purposes,
             * the actual MLE corresponds to 'x_t' and is equal to 'solution + x_0'.
             */
            if (chosen_sol > 0 && chosen_sol < x[maxed_at]-x[maxed_at-1])
            {
                double temp=((ld*chosen_sol+lc)*chosen_sol+lb)*chosen_sol+y[maxed_at-1];
                if (temp > maxed)
                {
                    maxed=temp;
                    x_max=chosen_sol+x[maxed_at-1];
                }
            }
        }
    }

    // Repeating for the segment on the right.
    if (maxed_at<npts-1)
    {
        double rd=d[maxed_at];
        double rc=c[maxed_at];
        double rb=b[maxed_at];

        double delta = fsquare(rc)-3*rd*rb;
        int solvable = (delta < 0)? 0 : 1;

        if (solvable)
        {
            double numerator = -rc-sqrt(delta);
            double chosen_sol=(rd!=0)? numerator/(3*rd) : -1.0;
            if (chosen_sol > 0 && chosen_sol < x[maxed_at+1]-x[maxed_at])
            {
                double temp=((rd*chosen_sol+rc)*chosen_sol+rb)*chosen_sol+y[maxed_at];
                if (temp>maxed)
                {
                    maxed=temp;
                    x_max=chosen_sol+x[maxed_at];
                }
            }
        }
    }

    return x_max;
}

/* wrapper function for find_max and apply it to each row
 *
 * For each row (tag) of the log-likelihood matrix ll, copies that row out and calls
 * find_max at the shared spline abscissas spts, recording the maximising abscissa.
 *
 * inputs:
 * spts   spline abscissas shared by every row (length ll->ncol)
 * ll     compressed matrix of interpolation ordinates, ll->nrow rows by ll->ncol cols
 *
 * outputs:
 * out    caller-allocated, the maximising abscissa for each row (length ll->nrow)
 *
 * side effects: allocates then frees the per-row spline scratch (lptr, b, c, d);
 * writes only into out.
 */
void max_interpolant(double *spts, cmx *ll, double *out)
{
    int npts=(ll->ncol), ntag=(ll->nrow);

    double *lptr = R_Calloc(npts, double);
    double *b    = R_Calloc(npts, double);
    double *c    = R_Calloc(npts, double);
    double *d    = R_Calloc(npts, double);

    for (int tag=0; tag<ntag; ++tag)
    {
        get_row(ll,tag,lptr);
        out[tag]=find_max(npts, spts, lptr, b, c, d);
    }

    R_Free(lptr);
    R_Free(b);
    R_Free(c);
    R_Free(d);

    return;
}
