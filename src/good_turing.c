#include "edgeR.h"

/* comments from Aaron:
 * Implements the simple version of the Good-Turing frequency estimator in C++.
 * This is based on the C code written by Geoffrey Sampson from Sussex University.
 * It takes in a vector of observed frequencies and another vector of the same
 * length of frequencies (of observed frequencies). The first vector must be 
 * sorted in ascending order. It also takes a numeric scalar which describes
 * the confidence factor.
 */

/* this function performs simple good turing  
 *
 * inputs:
 * Obs   observation vector     
 * Freq  frequence vector
 * nrows length of observation
 * conf  confidence factor
 * 
 * outputs:
 * pzero zero proportion
 * out   estimated proportion
 * 
 * comments:
 * This is converted from R_simple_good_turing.cpp written by Aaron
 */

void good_turing (int *Obs, int *Freq, int nrows, double conf, double *pzero, double *out) 
{
	// Computing constant values.
	double bigN=0;
	double XYs=0, meanX=0, meanY=0, Xsquares=0;

    double *log_obs = R_Calloc(nrows, double);
    for(int i=0;i<nrows;++i) log_obs[i] = log(Obs[i]);

	int last=nrows-1;

    for (int i=0; i<nrows; ++i) { 
		bigN += Obs[i] * Freq[i];

		// Computing log data.
		int x = (i==0)? 0 : Obs[i-1];
		double logO = log_obs[i];

        int xx = (i==last)? 2*(Obs[i]-x) : Obs[i+1]-x;
		double logZ = log(2*Freq[i]) - log(xx);

		meanX    += logO;
		meanY    += logZ;
		XYs      += logO*logZ;
		Xsquares += logO*logO;
	}

	meanX    /= nrows;
	meanY    /= nrows;
	XYs      -= meanX*meanY*nrows;
	Xsquares -= meanX*meanX*nrows;

	double slope = XYs/Xsquares;

	// Computing other various bits and pieces.
	pzero[0] = (nrows==0 || Obs[0]!=1) ? 0 : Freq[0]/bigN;

    // Collecting results.
    double bigNprime=0;
	int indiffValsSeen=0;

    for (int i=0; i<nrows; ++i) {
		int next_obs=Obs[i]+1;
		double y = next_obs*exp(slope*(log(next_obs)-log_obs[i])); // don't need intercept, cancels out.
		if (i==last || Obs[i+1]!=next_obs) { 
            indiffValsSeen=1; 
        }

		if (1-indiffValsSeen) {
            double x = (double)(next_obs)*Freq[i+1]/Freq[i];
            if (fabs(x - y) <= conf * x * sqrt(1.0/Freq[i+1] + 1.0/Freq[i])) { // Simplified expression.
                indiffValsSeen=1;
            } else { 
                out[i]=x;
			}
		} 
		if (indiffValsSeen) { 
			out[i]=y;
        }
        bigNprime += out[i]*Freq[i];
	}

	// Running through them to compute the remaining bit.
	double factor=(1.0-pzero[0])/bigNprime;       
    for (int i=0;i<nrows;++i) { 
        out[i]*=factor; 
    }

    R_Free(log_obs);
	
    return; 
}
