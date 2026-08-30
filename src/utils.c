#include "edgeR.h"
#include <float.h>
#include <math.h>

/* Small general-purpose helpers shared across the package. These were previously
 * scattered in domain files (fsquare/fcube in ql_glm.c, clamp_threads inline in
 * edgeR.h, brent_fmin in coxreid.c); they are collected here because none is
 * specific to the module that used to host it. Prototypes are in edgeR.h. */

/* fsquare, fcube: small math helpers (used by ql_glm.c, interpolator.c, loess_by_col.c). */
double fsquare(double x)
{
    return x * x;
}

double fcube(double x)
{
    return x * x * x;
}

/* clamp a requested thread count to the [1, available] range */
int clamp_threads(int n)
{
#ifdef _OPENMP
    int hi = omp_get_max_threads();
    if (n < 1) n = 1;
    if (n > hi) n = hi;
    return n;
#else
    (void) n;
    return 1;
#endif
}

/* Brent's method for 1-D minimization on [ax, bx]. Port of the classic
 * Forsythe/Malcolm/Moler fmin used by R's optimize(). Returns the minimizer. */
double brent_fmin(double ax, double bx, double (*f)(double, void*), void *info, double tol)
{
    /* c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    const double eps = sqrt(DBL_EPSILON);

    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, tol1, tol3;

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.;
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

    for(;;)
    {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;

        if (fabs(x - xm) <= t2 - (b - a) * .5)
        {
            break;
        }
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1)   /* fit parabola */
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.)
            {
                p = -p;
            }
            else
            {
                q = -q;
            }
            r = e;
            e = d;
        }

        if (fabs(p) >= fabs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x))   /* golden-section step */
        {
            e = (x < xm) ? b - x : a - x;
            d = c * e;
        }
        else   /* parabolic interpolation step */
        {
            d = p / q;
            u = x + d;
            if (u - a < t2 || b - u < t2)
            {
                d = tol1;
                if (x >= xm)
                {
                    d = -d;
                }
            }
        }

        if (fabs(d) >= tol1)
        {
            u = x + d;
        }
        else if (d > 0.)
        {
            u = x + tol1;
        }
        else
        {
            u = x - tol1;
        }

        fu = (*f)(u, info);

        if (fu <= fx)
        {
            if (u < x)
            {
                b = x;
            }
            else
            {
                a = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if (u < x)
            {
                a = u;
            }
            else
            {
                b = u;
            }
            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    return x;
}
