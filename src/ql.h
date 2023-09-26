#ifndef QL_H
#define QL_H

/* constants */
extern const double low_bound;

/* compute weight functions */
void compute_weight(double, double, double, double*);

/* compute unit deviance */
double compute_unit_nb_deviance (double, double, double);

/* qr decomposition*/
void QR_hat (double*, int, int, double*);

#endif
