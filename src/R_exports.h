#ifndef REXPORT_H
#define REXPORT_H

#include "edgeR.h"

/* Defining all R-accessible functions. */

/* genewise negative-binomial residual deviances */
SEXP compute_nbdev(SEXP, SEXP, SEXP, SEXP, SEXP);

/* adjusted profile likelihood of the NB dispersion, genewise */
SEXP compute_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* fit NB GLMs and return the adjusted profile likelihood */
SEXP fit_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* Cox-Reid common-dispersion optimization (coxreid.c) */
SEXP coxreid_disp (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP coxreid_disp_top (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* exact negative-binomial test p-values by the deviance method */
SEXP exact_test_by_deviance(SEXP, SEXP, SEXP, SEXP, SEXP);

/* genewise NB GLM fit by the Levenberg-Marquardt algorithm */
SEXP fit_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* starting coefficients for the Levenberg-Marquardt GLM fit */
SEXP get_levenberg_start (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* column-wise lowess smoothing of a matrix */
SEXP loess_by_col(SEXP, SEXP, SEXP);

/* maximize a spline interpolant to locate the optimal dispersion */
SEXP maximize_interpolant(SEXP, SEXP);

/* fit a single-group NB GLM, genewise */
SEXP fit_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* fitted values for a one-way layout from group coefficients */
SEXP get_one_way_fitted (SEXP, SEXP, SEXP, SEXP);

/* Good-Turing frequency estimation */
SEXP simple_good_turing (SEXP, SEXP, SEXP);

/* add scaled prior counts and return the corresponding offsets */
SEXP add_prior_count (SEXP, SEXP, SEXP, SEXP);

/* genewise log2 counts-per-million */
SEXP calculate_cpm_log (SEXP, SEXP, SEXP, SEXP);

/* genewise raw counts-per-million */
SEXP calculate_cpm_raw (SEXP, SEXP, SEXP);

/* average log2 counts-per-million per gene */
SEXP ave_log_cpm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* flag genes whose QL variance falls below the Poisson bound */
SEXP check_poisson_bound (SEXP, SEXP, SEXP, SEXP);

/* NB QL bias-adjusted deviance, residual df and dispersion (vectors) */
SEXP compute_adj_vec (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* NB QL adjusted unit deviance, df and leverage (matrices) */
SEXP sample_weights (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* average NB QL dispersion across genes */
SEXP compute_ave_qd (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* one-group binomial fit, genewise */
SEXP bin_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* binomial GLM fit by IWLS for a general design */
SEXP bin_fit_iwls (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* unified binomial GLM fit (oneway shortcut or IWLS, chosen in C) */
SEXP bin_fit (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* single-call QL binomial fit (fit + QL adjustment in one call) */
SEXP bin_ql_fit (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* binomial QL adjusted unit deviance, df and leverage (matrices) */
SEXP sample_weights_bin (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* single-call binomial QL F-test with TREAT/UPSHOT p-values */
SEXP bin_ftest (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* single-call negative-binomial QL F-test with TREAT/UPSHOT p-values */
SEXP glm_ftest (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* extra fitters retained in this version's C backend */
/* one-way layout NB GLM fit, genewise */
SEXP fit_one_way (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* unified NB GLM fit (oneway shortcut or Levenberg, chosen in C) */
SEXP fit_glm (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* single-call predictive fold-change fit (prior-count augmentation + GLM fit) */
SEXP pred_fc (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* genewise differential-splicing GLM fit and F-tests */
SEXP fit_diff_splice (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* parse and tabulate hairpin-barcode reads (R_process_hairpin_reads.c) */
void processHairpinReads(int *, int *, char**, char**, int*,
		char**, char**, int*, int*, int*, int*, int*, int*,
		int*, int*, int*, int*, int*, int*, int *,
		int *, char**, int*);

#endif
