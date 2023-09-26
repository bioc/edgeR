#include <Rmath.h>

const double one_tenthousandth=1e-4, mildly_low_value=1e-8, one_million=1e6;

/* Function to calculate the deviance. Note the protection for very large mu*phi (where we
 * use a gamma instead) or very small mu*phi (where we use the Poisson instead). This
 * approximation protects against numerical instability introduced by subtracting
 * a very large log value in (log mu) with another very large logarithm (log mu+1/phi).
 * We need to consider the 'phi' as the approximation is only good when the product is
 * very big or very small.
 */

double compute_unit_nb_deviance (double y, double mu, double phi) {
	// We add a small value to protect against zero during division and logging.
    y+=mildly_low_value;
    mu+=mildly_low_value;

    /* Calculating the deviance using either the Poisson (small phi*mu), the Gamma (large) or NB (everything else).
     * Some additional work is put in to make the transitions between families smooth.
     */
    if (phi < one_tenthousandth) {
		const double resid = y - mu;
		return 2 * ( y * log(y/mu) - resid - 0.5*resid*resid*phi*(1+phi*(2/3*resid-y)) );
    } else {
		const double product=mu*phi;
		if (product > one_million) {
            return 2 * ( (y - mu)/mu - log(y/mu) ) * mu/(1+product);
        } else {
			const double invphi=1/phi;
            return 2 * (y * log( y/mu ) + (y + invphi) * log( (mu + invphi)/(y + invphi) ) );
        }
	}
}
