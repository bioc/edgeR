#include "edgeR.h"

/* this function performs exact test by deviance 
 *
 * inputs:
 * s1    rowsum vector for group 1
 * s2    rowsum vector for group 2
 * ntag  number of rows
 * n1    number of samples in group 1
 * n2    number of samples in group 2
 * disp  dispersion vector
 * 
 * output:
 * pval  p values vector
 * 
 * comment:
 * This is converted from R_exact_test_by_deviance.cpp written by Aaron
*/

void exact_test_by_dev(int *s1, int *s2, int ntag, int n1, int n2, double *disp, double *pval) 
{	
	int nlib = n1+n2;
    for (int tag=0; tag<ntag; ++tag) {
        const int stotal=s1[tag]+s2[tag];

		// Computing current means and sizes for each library (probability is the same).
		double mu = stotal/nlib;
		double mu1=mu*n1, mu2=mu*n2, r1=n1/disp[tag], r2=n2/disp[tag];
        double p = r1/(r1+mu1);

		/* The aim is to sum conditional probabilities for all partitions of the total sum with deviances 
 		 * greater than that observed for the current partition. We start computing from the extremes
 		 * in both cases.
 		 */
	    double phi1=1/r1, phi2=1/r2;
	    double obsdev=compute_unit_nb_deviance(s1[tag], mu1, phi1)+compute_unit_nb_deviance(s2[tag], mu2, phi2);
		
        pval[tag]=0;

		// Going from the left.	
		int j=0;
		while (j <= stotal) {
			if (obsdev <= compute_unit_nb_deviance(j, mu1, phi1)+compute_unit_nb_deviance(stotal-j, mu2, phi2)) { 
				pval[tag]+=dnbinom(j, r1, p, 0) * dnbinom(stotal-j, r2, p, 0);
			} else { break; }
			++j;
		}

		// Going from the right, or what's left of it.
		for (int k=0; k<=stotal-j; ++k) {
			if (obsdev <= compute_unit_nb_deviance(k, mu2, phi2)+compute_unit_nb_deviance(stotal-k, mu1, phi1)) { 
				pval[tag]+=dnbinom(k, r2, p, 0) * dnbinom(stotal-k, r1, p, 0);
			} else { break; }
		}

    	double totalr=r1+r2;
	    pval[tag] /= dnbinom(stotal, totalr, totalr/(totalr+mu1+mu2), 0);
    }
    
    return;
}