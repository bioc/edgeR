#include "glm.h"
#include "objects.h"

SEXP compute_apl(SEXP y, SEXP means, SEXP disps, SEXP weights, SEXP adjust, SEXP design) {
    BEGIN_RCPP
    any_numeric_matrix counts(y);
    const int num_tags=counts.get_nrow();
    const int num_libs=counts.get_ncol();

    // Setting up the means.
    Rcpp::NumericMatrix Means(means);
    if (Means.nrow()!=num_tags || Means.ncol()!=num_libs) {
        throw std::runtime_error("fitted value and count matrices must have the same dimensions");
    }

    // Setting up the dispersions and weights.
    compressed_matrix alld=check_CM_dims(disps, num_tags, num_libs, "dispersion", "count");
    compressed_matrix allw=check_CM_dims(weights, num_tags, num_libs, "weight", "count");

    // Determining whether we want to do the adjustment.
    bool do_adjust=check_logical_scalar(adjust, "adjustment specifier");

    // Setting up the CR adjustment object.
    int num_coefs=0;
    std::vector<double> X=check_design_matrix(design, num_libs, num_coefs);
    adj_coxreid acr(num_libs, num_coefs, X.data());

    // Generating output values.
    Rcpp::NumericVector output(num_tags);
    std::vector<double> working_weights(num_libs), current(num_libs);
    for (int tag=0; tag<num_tags; ++tag) {

        double& sum_loglik=output[tag];
        counts.fill_row(tag, current.data());
        auto curmeans=Means.row(tag);
        const double* dptr=alld.get_row(tag);
        const double* wptr=allw.get_row(tag);
    
        /* First computing the log-likelihood. */
        auto cmIt=curmeans.begin();
        for (int lib=0; lib<num_libs; ++lib, ++cmIt) { 
            const double& cury=current[lib];
            const double& curmu=(*cmIt);
            if (curmu==0) { 
                continue; // Should only be possible if count is zero, where the log-likelihood would then be 0.
            }

            double loglik=0;
            if (dptr[lib] > 0) {
                // same as loglik <- rowSums(weights*dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
                const double r=1/dptr[lib];
                const double logmur=std::log(curmu+r);
                loglik = cury*std::log(curmu) - cury*logmur + r*std::log(r) - r*logmur + lgamma(cury+r) - lgamma(cury+1) - lgamma(r); 
            } else {
                // same as loglik <- rowSums(weights*dpois(y,lambda=mu,log = TRUE))
                loglik = cury*std::log(curmu) - curmu - lgamma(cury+1);
            }
            sum_loglik += loglik * wptr[lib]; // with weights.
        }
        
        if (do_adjust) {
            /* Computing 'W', the matrix of negative binomial probabilities. 
             * The class computes 'XtWX' and performs an LDL decomposition 
             * to compute the Cox-Reid adjustment factor.
             */
            auto cmIt=curmeans.begin();
            for (int lib=0; lib<num_libs; ++lib, ++cmIt) {
                const double& curmu=*cmIt;
                working_weights[lib] = wptr[lib] * curmu/(1 + dptr[lib] * curmu);
            }

            double adj=0;
            if (num_coefs==1) {
                adj=std::accumulate(working_weights.begin(), working_weights.end(), 0.0);
                adj=std::log(std::abs(adj))/2;
            } else {
                std::pair<double, bool> x=acr.compute(working_weights.data());
                if (!x.second) { 
                    std::stringstream errout;
                    errout << "LDL factorization failed for tag " << tag+1;
                    throw std::runtime_error(errout.str());
                }
                adj=x.first;
            }
		    sum_loglik-=adj;
		} 
    }

    return output;
    END_RCPP
}

