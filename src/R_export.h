#ifndef REXPORT_H
#define REXPORT_H

#include "edgeR.h"

/* Defining all R-accessible functions. */

SEXP compute_nbdev(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP exact_test_by_deviance(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP fit_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_levenberg_start (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP loess_by_col(SEXP, SEXP, SEXP);

SEXP maximize_interpolant(SEXP, SEXP);

SEXP fit_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_one_way_fitted (SEXP, SEXP, SEXP);

SEXP simple_good_turing (SEXP, SEXP, SEXP);

SEXP add_prior_count (SEXP, SEXP, SEXP);

SEXP calculate_cpm_log (SEXP, SEXP, SEXP);

SEXP calculate_cpm_raw (SEXP, SEXP);

SEXP ave_log_cpm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP check_poisson_bound (SEXP, SEXP, SEXP);

SEXP compute_adj_vec (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_adj_mat (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP compute_ave_qd (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void processHairpinReads(int *, int *, char**, char**, int*,
		char**, char**, int*, int*, int*, int*, int*, int*,
		int*, int*, int*, int*, int*, int*, int *,
		int *, char**, int*);

#endif
