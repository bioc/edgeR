#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "R_exports.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(compute_nbdev, 5),
	CALLDEF(compute_apl, 7),
	CALLDEF(fit_apl, 10),
	CALLDEF(coxreid_disp, 9),
	CALLDEF(coxreid_disp_top, 10),
	CALLDEF(exact_test_by_deviance, 5),
	CALLDEF(loess_by_col, 3),
	CALLDEF(maximize_interpolant, 2),

    CALLDEF(fit_levenberg, 9),
	CALLDEF(get_levenberg_start, 7),
	CALLDEF(fit_one_group, 8),
	CALLDEF(fit_one_way, 10),
	CALLDEF(fit_glm, 9),
	CALLDEF(pred_fc, 9),
	CALLDEF(fit_diff_splice, 14),
	CALLDEF(get_one_way_fitted, 4),
	CALLDEF(simple_good_turing, 3),

    CALLDEF(add_prior_count, 4),
    CALLDEF(calculate_cpm_log, 4),
    CALLDEF(calculate_cpm_raw, 3),
    CALLDEF(ave_log_cpm, 8),

    CALLDEF(check_poisson_bound, 4),

    CALLDEF(compute_adj_vec, 7),
    CALLDEF(sample_weights, 8),
    CALLDEF(compute_ave_qd,  7),

    CALLDEF(bin_one_group,  6),
    CALLDEF(bin_fit_iwls,   7),
    CALLDEF(bin_fit,        9),
    CALLDEF(bin_ql_fit,     10),

    CALLDEF(sample_weights_bin, 6),

    CALLDEF(bin_ftest, 12),
    CALLDEF(glm_ftest, 13),

	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {"processHairpinReads", (DL_FUNC) &processHairpinReads, 22},
    {NULL, NULL, 0}
  };

void attribute_visible R_init_edgeR(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
