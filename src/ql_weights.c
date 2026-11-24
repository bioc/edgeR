#include <stdint.h>
#include <string.h>
#include "edgeR.h"

/*
  Weights alpha = 2E/V and kappa = 2E^2/V of the unit deviance, from which the QL code takes the
  mean E = kappa/alpha and variance V = 2 kappa/alpha^2 of the unit deviance, for the Poisson,
  negative binomial and binomial distributions.
*/

/* low bound, the same for the Poisson, NB and binomial weights: zero weights for mu < 1e-10. The Poisson
 * and binomial tables start at the octave 2^-34 < 1e-10; the NB small-mu piece is fitted on [1e-32, 2^-14]
 * and used from 1e-10. */
static const double low_bound = 1e-10;

/* top of the NB octave cells, 2^16; the gamma closure (one log and one exp per value) is used from here on */
static const double nb_top = 65536.0;

/* number of NB mu-blocks: the small-mu piece and the octaves e = -14..15 (block b = e + 15) */
static const int nb_blocks = 31;

/* top of the small-mu piece of the NB cells, 2^-14 */
static const double nb_small = 6.103515625e-05;

/* Exponent and local variable of an octave piece. */
static inline double octave_x(double v, int *e)
{
    uint64_t bits;
    double m;
    memcpy(&bits, &v, sizeof(double));
    *e = (int) ((bits >> 52) & 0x7ff) - 1023;
    bits = (bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL;
    memcpy(&m, &bits, sizeof(double));
    return 2 * m - 3;
}

/* Chebyshev series c[0] T_0(x) + ... + c[5] T_5(x) by the Clenshaw recurrence, unrolled. */
static inline double cheb6(const double *c, double x)
{
    double x2 = 2 * x;
    double b5 = c[5];
    double b4 = c[4] + x2 * b5;
    double b3 = c[3] + x2 * b4 - b5;
    double b2 = c[2] + x2 * b3 - b4;
    double b1 = c[1] + x2 * b2 - b3;
    return c[0] + x * b1 - b2;
}

/* Chebyshev series c[0] + c[1] T_1(y) + c[2] T_2(y), with T_2(y) = 2y^2 - 1. */
static inline double cheb3(const double *c, double y)
{
    return c[0] + y * (c[1] + 2 * y * c[2]) - c[2];
}

/* Tensor series of one 6 x 3 cell: sum_i sum_j C[6i + j] T_j(x) T_i(y), rows in x first. */
static inline double cell6x3(const double *C, double x, double y)
{
    double r[3];
    r[0] = cheb6(C, x);
    r[1] = cheb6(C + 6, x);
    r[2] = cheb6(C + 12, x);
    return cheb3(r, y);
}

/* Both Poisson QL weights for large mu, from the large-mu forms
 * alpha = 1 - 1/(6mu) - 1/(2mu^2) and kappa = 1 - 1/(2.5mu^2).
 *
 * inputs:
 *   mu   Poisson mean, mu >= 2^14 (checked by the caller)
 *
 * outputs (caller-allocated):
 *   wpt  length-2 array: wpt[0] alpha weight, wpt[1] kappa weight
 */
static inline void pois_weight_large(double mu, double *wpt)
{
    double r = 1 / mu;
    wpt[0] = 1 - r / 6 - 0.5 * r * r;
    wpt[1] = 1 - 0.4 * r * r;
}

/* Dispersion piece of the NB cells.
 *
 * inputs:
 *   phi  dispersion; values below 0 are treated as 0 and values of 4 or more as 4 (the cap of the
 *        previous approximation). The first test is written so that a NaN phi takes the bottom piece
 *        (index 0, NaN y) rather than an out-of-range octave.
 *
 * outputs:
 *   *y   local variable in [-1, 1]
 *
 * returns:
 *   p = 0 for phi < 2^-7 (y = 256 phi - 1), else p = k + 8 for the octave [2^k, 2^(k+1)), k = -7..1.
 */
static inline int nb_phi_piece(double phi, double *y)
{
    int k;
    if (!(phi >= 0.0078125))
    {
        if (phi < 0)
        {
            phi = 0;
        }
        *y = 256 * phi - 1;
        return 0;
    }
    if (phi >= 4)
    {
        *y = 1;
        return 9;
    }
    *y = octave_x(phi, &k);
    return k + 8;
}

/* Gamma closure for mu >= 2^16: W = W_inf + (W_top - W_inf) (2^16/mu)^(1/phi), where W_inf is the
 * gamma limit of the weight and W_top its value at mu = 2^16. phi is the capped dispersion. */
static inline double nb_closure(double mu, double phi, double winf, double wtop)
{
    double f;
    if (phi > 0)
    {
        f = exp(log(nb_top / mu) / phi);
    }
    else
    {
        f = (mu == nb_top) ? 1.0 : 0.0;
    }
    return winf + (wtop - winf) * f;
}

/* Effective dispersion used by the cells and the closure: phi clamped to [0, 4]. */
static inline double nb_phi_eff(double phi)
{
    if (phi < 0)
    {
        return 0;
    }
    return fmin2(phi, 4.0);
}

/* Computes both negative binomial QL weights from the octave cells.
 *
 * inputs:
 *   u      raw fitted mean
 *   phi    dispersion (capped at 4, as before)
 *   prior  divisor applied to the mean (mu = u/prior)
 *
 * outputs (caller-allocated):
 *   wpt    length-2 array: wpt[0] alpha weight, wpt[1] kappa weight
 *
 * Zero weights for mu < 1e-10. For mu >= 2^16 the gamma closure;
 * for mu < 2^-14 the small-mu piece, where the kappa table holds kappa/mu.
 * No return value; writes only into the caller-owned wpt array.
 */
void compute_weight_negbin(double u, double phi, double prior, double *wpt)
{
    double mu = u / prior;
    double x;
    double y;
    double uu;
    const double *ca;
    const double *ck;
    int b;
    int e;
    int p;
    if (!(mu >= low_bound))
    {
        wpt[0] = 0;
        wpt[1] = 0;
        return;
    }
    p = nb_phi_piece(phi, &y);
    if (mu >= nb_top)
    {
        ca = ql_nb_alpha + ((nb_blocks - 1) * 10 + p) * 18;
        ck = ql_nb_kappa + ((nb_blocks - 1) * 10 + p) * 18;
        wpt[0] = nb_closure(mu, nb_phi_eff(phi), cheb3(ql_nb_winf_alpha + 3 * p, y), cell6x3(ca, 1.0, y));
        wpt[1] = nb_closure(mu, nb_phi_eff(phi), cheb3(ql_nb_winf_kappa + 3 * p, y), cell6x3(ck, 1.0, y));
        return;
    }
    if (mu < nb_small)
    {
        uu = 1 / (5 - log(mu));
        x = (2 * uu - ql_nb_ua - ql_nb_ub) / (ql_nb_ub - ql_nb_ua);
        b = 0;
    }
    else
    {
        x = octave_x(mu, &e);
        b = e + 15;
    }
    wpt[0] = cell6x3(ql_nb_alpha + (b * 10 + p) * 18, x, y);
    wpt[1] = cell6x3(ql_nb_kappa + (b * 10 + p) * 18, x, y);
    if (b == 0)
    {
        wpt[1] *= mu;
    }
    return;
}

/* Binomial QL weights. With mu = m min(pi, 1 - pi), the moments depend on (mu, m) only:
 *   2 <= m < 16:      exact sum over y = 0..m;
 *   m >= 1024, or m >= 16 and mu < 1/16:
 *                     Poisson weights at mu times (1 + h(mu)/m), h a 6-term series per mu-octave;
 *                     the Poisson closed forms alone for mu >= 2^14;
 *   otherwise (16 <= m < 1024, 1/16 <= mu < 512):
 *                     octave cells, 6 terms in x (mu-octave) times 6 in y, where y = 32/m - 1 (1/m on
 *                     [0, 1/16]) for mu < 8 and y = 4 pi - 1 (pi on [0, 1/2]) for mu >= 8.
 * The 1/m correction is accurate to about 3e-4 from m = 1024 (1e-4 when pi = mu/m < 1/8 or m >= 4096),
 * and for every m >= 16 when mu < 1/16 (1e-4 when mu < 2^-18).
 * The binomial is the NB with phi = -1/m, whose branch points phi = -1/y pile up at 1/m = 0; this is
 * why the series in 1/m (or pi) needs six terms and the smallest m are summed exactly. */

/* the exact sum is used for m below bin_exact_top, the Poisson correction from bin_large on and for
 * means below bin_small_mu */
static const int bin_exact_top = 16;
static const int bin_large = 1024;
static const double bin_small_mu = 0.0625;

/* Exact binomial QL weights for 2 <= mm < 16 and 0 < pp <= 1/2, mu = mm pp >= low_bound.
 *
 * inputs:
 *   pp   success probability, at most 1/2
 *   mm   number of trials, 2..15
 *
 * outputs (caller-allocated):
 *   wpt  length-2 array: wpt[0] alpha weight, wpt[1] kappa weight
 *
 * The probabilities follow p(0) = (1 - pp)^mm = exp(mm log1p(-pp)) and the recurrence
 * p(y+1) = p(y) (mm - y)/(y + 1) pp/(1 - pp). The unit deviance is written as
 * d(y) = 2{y (log y - log mu) + (mm - y)(log(mm - y) - log mm - log1p(-pp))} with log k from ql_logk,
 * which keeps d(0) = -2 mm log1p(-pp) accurate for tiny pp. The mean and variance are summed in two
 * passes. One log, one log1p and one exp per call; no allocation.
 */
static void binom_exact(double pp, int mm, double *wpt)
{
    double pr[16];
    double de[16];
    double lq = log1p(-pp);
    double lmu = ql_logk[mm] + log(pp);
    double r = pp / (1 - pp);
    double dd = 0;
    double vv = 0;
    double t;
    int k;
    pr[0] = exp(mm * lq);
    for (k = 0; k < mm; k++)
    {
        pr[k + 1] = pr[k] * (mm - k) / (k + 1) * r;
    }
    de[0] = -2 * mm * lq;
    de[mm] = 2 * mm * (ql_logk[mm] - lmu);
    for (k = 1; k < mm; k++)
    {
        de[k] = 2 * (k * (ql_logk[k] - lmu) + (mm - k) * (ql_logk[mm - k] - ql_logk[mm] - lq));
    }
    for (k = 0; k <= mm; k++)
    {
        dd += pr[k] * de[k];
    }
    for (k = 0; k <= mm; k++)
    {
        t = de[k] - dd;
        vv += pr[k] * t * t;
    }
    wpt[0] = 2 * dd / vv;
    wpt[1] = 2 * dd * dd / vv;
}

/* Tensor series of one binomial 6 x 6 cell: sum_k sum_j C[6k + j] T_j(x) T_k(y), rows in x first. */
static inline double cell6x6(const double *C, double x, double y)
{
    double r[6];
    r[0] = cheb6(C, x);
    r[1] = cheb6(C + 6, x);
    r[2] = cheb6(C + 12, x);
    r[3] = cheb6(C + 18, x);
    r[4] = cheb6(C + 24, x);
    r[5] = cheb6(C + 30, x);
    return cheb6(r, y);
}

/* Computes both binomial QL weights for a proportion y over mm trials.
 * Works with the smaller tail probability by symmetry (pp = min(y, 1-y)) and mu = pp * mm.
 * Returns zeros for mm < 2 (a single trial carries no information on the dispersion) and for
 * mu < 1e-10 (the test is written so that a NaN or negative mu also gives zeros); otherwise uses the
 * exact sum, the corrected Poisson weights or the octave cells described above.
 *
 * inputs:
 *   y    fitted proportion
 *   mm   number of binomial trials
 *
 * outputs (caller-allocated):
 *   wpt  length-2 array: wpt[0] alpha weight, wpt[1] kappa weight
 *
 * No return value; writes only into the caller-owned wpt array.
 */
void compute_weight_binomial (double y, int mm, double *wpt)
{
    double pp;
    double mu;
    double x;
    double yv;
    double rm;
    int e;
    pp = (y < 0.5) ? y : 1 - y;
    mu = pp * mm;
    if (!(mu >= low_bound) || mm < 2)
    {
        wpt[0] = 0;
        wpt[1] = 0;
        return;
    }
    if (mm < bin_exact_top)
    {
        binom_exact(pp, mm, wpt);
        return;
    }
    if (mm >= bin_large || mu < bin_small_mu)
    {
        if (mu >= 16384)
        {
            pois_weight_large(mu, wpt);
            return;
        }
        /* Poisson octave series and h on the same octave, so the octave is found once */
        x = octave_x(mu, &e);
        rm = 1.0 / mm;
        wpt[0] = cheb6(ql_pois_alpha + (e + 34) * 6, x) * (1 + cheb6(ql_binh_alpha + (e + 34) * 6, x) * rm);
        wpt[1] = cheb6(ql_pois_kappa + (e + 34) * 6, x) * (1 + cheb6(ql_binh_kappa + (e + 34) * 6, x) * rm);
        return;
    }
    /* 16 <= mm < 1024 and 1/16 <= mu < 512, so the octave e lies in -4..8 */
    x = octave_x(mu, &e);
    yv = (e <= 2) ? 32.0 / mm - 1 : 4 * pp - 1;
    wpt[0] = cell6x6(ql_bin_alpha + (e + 4) * 36, x, yv);
    wpt[1] = cell6x6(ql_bin_kappa + (e + 4) * 36, x, yv);
    return;
}