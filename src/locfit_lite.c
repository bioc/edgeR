/* locfit_lite.c -- local regression smoothing of the columns of a matrix (locfitByCol)
 *
 * Each column of y is smoothed against the same covariate x, with the same
 * prior weights w, span alpha and degree. The method is locfit's default 1-D
 * local regression (Loader, 1999, Local Regression and Likelihood, Springer):
 *
 *   - At a point t the bandwidth h is the distance to the k-th nearest x,
 *     k = floor(n * alpha), and the kernel is tricube. The local fit is the
 *     Gaussian local likelihood: weighted least squares of degree 0 to 3,
 *     solved by a scaled eigen pseudo-inverse.
 *   - The fit is computed only at the vertices of an adaptive tree over
 *     [min x, max x]. A cell is split at its midpoint while it is wider than
 *     0.8 times its vertex bandwidths. A new vertex in a cell already narrower
 *     than that is a pseudo vertex: its value comes from the cubic through its
 *     parents rather than from a local fit.
 *   - Fitted values are interpolated within each leaf cell: linearly for
 *     degree 0, by cubic Hermite with the vertex slopes otherwise.
 *
 * The results agree with fitted(locfit(y ~ x, weights=w, alpha=alpha,
 * deg=degree)) to rounding error. They differ in one case: a vertex whose
 * window has no positive weight (ties in x can cause it) gets its value from
 * its parents, as a pseudo vertex does.
 *
 * The tree, the windows and the local solves depend on x, w, alpha and degree
 * only, and every fit is linear in y. So the setup builds, once, the tree and
 * for each fitted vertex the kernel rows giving its value and slope as
 * weighted sums of y. Each column is then a few dot products and an
 * interpolation, run on OpenMP threads.
 */
#include <math.h>

#include "edgeR.h"

#define CUT 0.8           /* split a cell wider than CUT x its bandwidths */
#define MAXDEG 3
#define MAXP (MAXDEG + 1)

/* status of smooth_setup, and the error message for each failure */
#define SETUP_OK 0
#define SETUP_NVERT 1
#define SETUP_EMPTY 2
#define SETUP_EIGEN 3
static const char *setup_error[] = {"", "too many vertices in the adaptive tree", "a local window has no positive weights", "eigen decomposition failed"};

typedef struct
{
    int n, deg, nv, nleaf;
    int nfill;        /* fitted vertices with an empty window, filled from their parents */
    int *ord;         /* x[ord[i]] is the i-th smallest x */
    double *xs, *ws;  /* x and w in that order */
    /* vertex v: position, bandwidth, parents (for v >= 2) */
    double *vx, *vh;
    int *pa, *pb;
    /* fitted vertex: kernel rows over sorted points lo .. lo + len - 1,
     * value row at ker + koff, slope row (deg > 0) after it; len 0 if the
     * vertex takes its value from its parents */
    int *lo, *len;
    size_t *koff;
    double *ker;
    int *leaf;        /* leaf j: the cell between vertices leaf[2j] and leaf[2j + 1] */
} smoother;

/* tricube kernel weight of a point at distance d >= 0 for bandwidth h; with
 * h = 0, weight 1 at distance 0 and 0 elsewhere */
static double tricube(double d, double h)
{
    double u;
    if (h == 0)
        return d == 0;
    u = d / h;
    if (u >= 1)
        return 0;
    u = 1 - u * u * u;
    return u * u * u;
}

/* bandwidth at t
 *
 * inputs:
 * xs  sorted x, length n
 * k   floor(n * alpha), the number of nearest neighbours
 * t   the fitting point
 *
 * outputs:
 * lo, len  the run xs[lo .. lo + len - 1] of points with positive kernel weight
 *
 * returns the distance to the k-th nearest point, or, when k >= n, the largest
 * distance times k / n */
static double bandwidth(const double *xs, int n, double k, double t, int *lo, int *len)
{
    int l, r, j, a, b;
    double h = 0;

    if (k >= n)
    {
        h = fmax(t - xs[0], xs[n - 1] - t) * (k / n);
        l = 0;
        r = n - 1;
    }
    else
    {
        /* a: the first point at or right of t */
        a = 0;
        b = n;
        while (a < b)
        {
            j = (a + b) / 2;
            if (xs[j] < t)
                a = j + 1;
            else
                b = j;
        }
        /* take the k nearest, merging outwards from t */
        l = a - 1;
        r = a;
        for (j = 0; j < k; j++)
        {
            if (l >= 0 && (r >= n || t - xs[l] <= xs[r] - t))
                h = t - xs[l--];
            else
                h = xs[r++] - t;
        }
        l++;
        r--;
        /* with h = 0, every point at t has weight 1, not only the k taken */
        while (l > 0 && tricube(t - xs[l - 1], h) > 0)
            l--;
        while (r < n - 1 && tricube(xs[r + 1] - t, h) > 0)
            r++;
    }
    /* points at distance h have weight 0 */
    while (l <= r && tricube(fabs(xs[l] - t), h) == 0)
        l++;
    while (r >= l && tricube(fabs(xs[r] - t), h) == 0)
        r--;
    *lo = l;
    *len = r - l + 1;
    return h;
}

/* equivalent kernel of the local fit at t
 *
 * The local fit is weighted least squares on 1, (x - t), ..., (x - t)^deg,
 * started at the weighted mean and moved by one Newton step with the
 * pseudo-inverse of X'WX. It is linear in y, so its value and slope at t are
 * weighted sums of the y of the window.
 *
 * inputs:
 * xs, ws  the m sorted x and prior weights of the window
 * t, h    fitting point and bandwidth
 * deg     degree, 0 to 3
 * W       workspace, m doubles (kernel x prior weights)
 *
 * outputs:
 * ker     m x (1 + (deg > 0)): the value row, then for deg > 0 the slope row
 *
 * returns 1, or 0 if the window has no weight, or -1 if dsyev fails */
static int local_kernel(const double *xs, const double *ws, int m, double t, double h, int deg, double *ker, double *W)
{
    int p = deg + 1, nc = (deg > 0) + 1, i, a, b, j, info, lwork = 64;
    double Z[MAXP * MAXP] = {0}, S[MAXP * MAXP], M[MAXP * MAXP], D[MAXP], ev[MAXP], work[64];
    double pw[MAXP], c[2], s0 = 0, tol, sum;

    for (i = 0; i < m; i++)
    {
        W[i] = tricube(fabs(xs[i] - t), h) * ws[i];
        s0 += W[i];
        pw[0] = 1;
        for (a = 1; a < p; a++)
            pw[a] = pw[a - 1] * (xs[i] - t);
        for (a = 0; a < p; a++)
            for (b = 0; b < p; b++)
                Z[a + p * b] += W[i] * pw[a] * pw[b];
    }
    if (s0 == 0)
        return 0;

    /* M = D Q diag(1/ev) Q' D, where D Z D = Q diag(ev) Q' has unit diagonal;
     * eigenvalues at or below 1e-8 of the largest count as zero */
    for (a = 0; a < p; a++)
        D[a] = (Z[a * (p + 1)] > 0) ? 1 / sqrt(Z[a * (p + 1)]) : 0;
    for (a = 0; a < p; a++)
        for (b = 0; b < p; b++)
            S[a + p * b] = D[a] * Z[a + p * b] * D[b];

    /* dsyev (LAPACK, symmetric eigenproblem): all eigenvalues and eigenvectors
     * of the p x p symmetric matrix S = D Z D.
     *
     * edgeR calls it to form the pseudo-inverse of the scaled local Gram
     * matrix X'WX: the eigenvalues are needed to drop the directions below the
     * 1e-8 relative cutoff, which keeps the fit defined when the window's
     * design is singular (e.g. tied x with deg > 0).
     *
     * Argument mapping:
     *   JOBZ = "V" (eigenvectors too), UPLO = "U" (upper triangle of S read),
     *   N = &p (order, at most 4), A = S (p x p, overwritten by the
     *   orthonormal eigenvectors, one per column), LDA = &p (leading
     *   dimension), W = ev (eigenvalues, ascending), WORK = work,
     *   LWORK = &lwork (64 >= 3p - 1, so no workspace query), INFO = &info
     *   (0 on success, > 0 if the QR iteration did not converge).
     *   FCONE FCONE pass the hidden Fortran lengths of JOBZ and UPLO.
     *
     * On return S holds Q and ev the eigenvalues, with D Z D = Q diag(ev) Q'.
     *
     * Equivalent R operation:
     *   e <- eigen(S, symmetric = TRUE)
     *   Q <- e$vectors
     *   ev <- e$values
     *   (R gives them in decreasing order).
     *
     * Netlib references: https://netlib.org/lapack/explore-html/ (dsyev). */
    F77_CALL(dsyev)("V", "U", &p, S, &p, ev, work, &lwork, &info FCONE FCONE);
    if (info != 0)
        return -1;
    tol = 1e-8 * ev[p - 1];
    for (a = 0; a < p; a++)
        for (b = 0; b < p; b++)
        {
            sum = 0;
            for (j = 0; j < p; j++)
                if (ev[j] > tol)
                    sum += S[a + p * j] * S[b + p * j] / ev[j];
            M[a + p * b] = D[a] * D[b] * sum;
        }

    /* coefficients = M X'W y + (e1 - M Z e1) * (weighted mean of y); the
     * second term is zero unless X'WX is singular */
    for (a = 0; a < nc; a++)
    {
        c[a] = (a == 0);
        for (b = 0; b < p; b++)
            c[a] -= M[a + p * b] * Z[b];
        c[a] /= s0;
    }
    for (i = 0; i < m; i++)
    {
        pw[0] = 1;
        for (a = 1; a < p; a++)
            pw[a] = pw[a - 1] * (xs[i] - t);
        for (a = 0; a < nc; a++)
        {
            sum = c[a];
            for (b = 0; b < p; b++)
                sum += M[a + p * b] * pw[b];
            ker[a * m + i] = W[i] * sum;
        }
    }
    return 1;
}

/* the vertex space locfit allows for span alpha */
static int max_vertices(double alpha)
{
    double a0 = (alpha > 1) ? 1 : 1 / alpha;
    return (int)((5 * a0 / CUT + 1) * 2);
}

/* the setup shared by all columns
 *
 * inputs:
 * x, w   covariate and prior weights, length n
 * alpha  span; deg  degree, 0 to 3
 *
 * outputs:
 * s      the sorted data, the tree, the kernels of its fitted vertices and the
 *        leaves; its arrays are R_Calloc'ed and freed by smooth_free, also when
 *        the setup fails
 *
 * returns SETUP_OK, or the SETUP_* code of the failure */
static int smooth_setup(smoother *s, const double *x, const double *w, int n, double alpha, int deg)
{
    int i, v, a, b, top, nvm = max_vertices(alpha), nc = (deg > 0) + 1, fit, status = SETUP_OK;
    double k = floor(n * alpha + 1e-12), le, hmin, score;
    size_t nker;

    s->n = n;
    s->deg = deg;
    s->nfill = 0;
    s->ker = NULL;
    s->ord = R_Calloc(n, int);
    s->xs = R_Calloc(n, double);
    s->ws = R_Calloc(n, double);
    for (i = 0; i < n; i++)
    {
        s->ord[i] = i;
        s->xs[i] = x[i];
    }
    rsort_with_index(s->xs, s->ord, n);
    for (i = 0; i < n; i++)
        s->ws[i] = w[s->ord[i]];

    s->vx = R_Calloc(nvm, double);
    s->vh = R_Calloc(nvm, double);
    s->pa = R_Calloc(nvm, int);
    s->pb = R_Calloc(nvm, int);
    s->lo = R_Calloc(nvm, int);
    s->len = R_Calloc(nvm, int);
    s->koff = R_Calloc(nvm, size_t);
    s->leaf = R_Calloc(2 * (size_t)nvm, int);
    int *stack = R_Calloc(2 * (size_t)nvm, int);
    char *pseudo = R_Calloc(nvm, char);
    double *W = R_Calloc(n, double);

    /* the tree: vertices at the ends, then split cells depth first, left
     * half first, so the leaves come out in order of x */
    for (v = 0; v < 2; v++)
    {
        s->vx[v] = s->xs[v == 0 ? 0 : n - 1];
        s->vh[v] = bandwidth(s->xs, n, k, s->vx[v], &s->lo[v], &s->len[v]);
        s->pa[v] = s->pb[v] = -1;
        pseudo[v] = 0;
    }
    s->nv = 2;
    s->nleaf = 0;
    a = 0;
    b = 1;
    top = 0;
    while (1)
    {
        /* split if the cell is wider than CUT x its smaller positive
         * bandwidth, or, with none, CUT / 2 x the range of x */
        le = s->vx[b] - s->vx[a];
        hmin = fmin(s->vh[a], s->vh[b]);
        if (hmin <= 0)
            hmin = fmax(s->vh[a], s->vh[b]);
        score = (hmin > 0) ? le / hmin : 2 * le / (s->vx[1] - s->vx[0]);
        if (score > CUT)
        {
            if (s->nv == nvm)
            {
                status = SETUP_NVERT;
                break;
            }
            v = s->nv++;
            s->vx[v] = (s->vx[a] + s->vx[b]) / 2;
            s->pa[v] = a;
            s->pb[v] = b;
            pseudo[v] = le < CUT * fmin(s->vh[a], s->vh[b]);
            if (pseudo[v])
            {
                s->vh[v] = (s->vh[a] + s->vh[b]) / 2;
                s->len[v] = 0;
            }
            else
                s->vh[v] = bandwidth(s->xs, n, k, s->vx[v], &s->lo[v], &s->len[v]);
            stack[top++] = v;
            stack[top++] = b;
            b = v;
            continue;
        }
        s->leaf[2 * s->nleaf] = a;
        s->leaf[2 * s->nleaf + 1] = b;
        s->nleaf++;
        if (top == 0)
            break;
        b = stack[--top];
        a = stack[--top];
    }

    /* kernels of the fitted vertices; an empty window is filled from the
     * parents, except at the two end vertices, which have none */
    if (status == SETUP_OK)
    {
        nker = 0;
        for (v = 0; v < s->nv; v++)
        {
            s->koff[v] = nker;
            nker += (size_t)nc * s->len[v];
        }
        s->ker = R_Calloc(nker, double);
        for (v = 0; v < s->nv && status == SETUP_OK; v++)
        {
            if (pseudo[v])
                continue;
            fit = 0;
            if (s->len[v] > 0)
                fit = local_kernel(s->xs + s->lo[v], s->ws + s->lo[v], s->len[v], s->vx[v], s->vh[v], deg, s->ker + s->koff[v], W);
            if (fit == -1)
                status = SETUP_EIGEN;
            else if (fit == 0 && v < 2)
                status = SETUP_EMPTY;
            else if (fit == 0)
            {
                s->len[v] = 0;
                s->nfill++;
            }
        }
    }

    R_Free(stack);
    R_Free(pseudo);
    R_Free(W);
    return status;
}

static void smooth_free(smoother *s)
{
    R_Free(s->ord);
    R_Free(s->xs);
    R_Free(s->ws);
    R_Free(s->vx);
    R_Free(s->vh);
    R_Free(s->pa);
    R_Free(s->pb);
    R_Free(s->lo);
    R_Free(s->len);
    R_Free(s->koff);
    R_Free(s->leaf);
    R_Free(s->ker);
}

/* one column
 *
 * inputs:
 * s      the setup
 * y      the column, length n
 * ys, f, df  workspace: n doubles, and nv doubles each for the vertex values
 *        and slopes
 *
 * outputs:
 * y      overwritten by its fitted values */
static void smooth_column(const smoother *s, double *y, double *ys, double *f, double *df)
{
    int i, j, v, a, b, n = s->n, deg = s->deg;
    double d, u, h1, sum;

    for (i = 0; i < n; i++)
        ys[i] = y[s->ord[i]];

    /* vertices in creation order, so parents come before children */
    for (v = 0; v < s->nv; v++)
    {
        int m = s->len[v];
        if (m > 0)
        {
            const double *kptr = s->ker + s->koff[v], *yptr = ys + s->lo[v];
            sum = 0;
            for (i = 0; i < m; i++)
                sum += kptr[i] * yptr[i];
            f[v] = sum;
            if (deg > 0)
            {
                sum = 0;
                for (i = 0; i < m; i++)
                    sum += kptr[m + i] * yptr[i];
                df[v] = sum;
            }
        }
        else
        {
            /* value and slope at the midpoint of the cubic through the parents */
            a = s->pa[v];
            b = s->pb[v];
            d = s->vx[b] - s->vx[a];
            f[v] = (f[a] + f[b]) / 2;
            if (deg > 0)
            {
                f[v] += d * (df[a] - df[b]) / 8;
                df[v] = 1.5 * (f[b] - f[a]) / d - (df[a] + df[b]) / 4;
            }
        }
    }

    /* the points in order of x, walking the leaves; a point on a cell
     * boundary belongs to the right-hand cell */
    j = 0;
    for (i = 0; i < n; i++)
    {
        while (j < s->nleaf - 1 && s->xs[i] >= s->vx[s->leaf[2 * j + 1]])
            j++;
        a = s->leaf[2 * j];
        b = s->leaf[2 * j + 1];
        d = s->vx[b] - s->vx[a];
        if (d == 0)
            sum = f[a];
        else
        {
            u = (s->xs[i] - s->vx[a]) / d;
            if (deg == 0)
                sum = (1 - u) * f[a] + u * f[b];
            else
            {
                h1 = u * u * (3 - 2 * u);
                sum = (1 - h1) * f[a] + h1 * f[b] + d * u * (1 - u) * ((1 - u) * df[a] - u * df[b]);
            }
        }
        y[s->ord[i]] = sum;
    }
}

/* edgeR locfitByCol after its argument checks (degree 0 to 3, ncol >= 1)
 *
 * inputs:
 * y         column-major ntags x ncol matrix
 * x, w      covariate and prior weights, length ntags
 * span      alpha; degree  0 to 3
 * nthreads  threads over columns
 *
 * outputs:
 * y         each column replaced by its fitted values
 *
 * Errors if the setup fails; warns if vertices were filled from their parents. */
void locfit_by_col_mat(double *y, int ntags, int ncol, const double *x, const double *w, double span, int degree, int nthreads)
{
    smoother s;
    int nth = clamp_threads(nthreads), status, nfill;

    if (nth > ncol)
        nth = ncol;

    status = smooth_setup(&s, x, w, ntags, span, degree);
    if (status != SETUP_OK)
    {
        smooth_free(&s);
        Rf_error("locfitByCol: %s", setup_error[status]);
    }

    /* per thread: sorted y, vertex values and slopes */
    size_t per = (size_t)ntags + 2 * (size_t)s.nv;
    double *buf = R_Calloc(nth * per, double);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nth) schedule(dynamic)
    #endif
    for (int j = 0; j < ncol; j++)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *bptr = buf + tid * per;
        smooth_column(&s, y + (size_t)j * ntags, bptr, bptr + ntags, bptr + ntags + s.nv);
    }

    nfill = s.nfill;
    R_Free(buf);
    smooth_free(&s);
    if (nfill > 0)
        Rf_warning("locfitByCol: %d local fits had no points with positive weight and were interpolated", nfill);
}
