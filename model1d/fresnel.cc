// This program is public domain.

/* fzw */
/** file
 *  calculating the reflectivity with fresnel_amplitude.
 */

#include <complex>
#include "reflcalc.h"

// Note: I've verified that at least for some values, the reflectivity
// in an extremely thick layer with absorption matches the reflectivity
// returned by fresnel_amplitude.  Therefore, the substrate can have
// an absorption term.  Obviously, the vacuum cannot since any absorption
// through a semi-infinite medium will have zero transmission.  However,
// since this is being used for both positive and negative Q values we
// will accept an absorption term for both, and you will have to apply
// an absorption term separately based on the details of your substrate.

extern "C" void 
fresnel_reflectivity(const double vrho, const double srho, 
		     const double vmu, const double smu, 
		     int points, const double Q[], double R[],
		     const double lambda)
{
  const double cutoff = 1e-10;
  const std::complex<double> Sp(16*M_PI*(srho-vrho),8*M_PI*smu/lambda);
  const std::complex<double> Sm(16*M_PI*(vrho-srho),8*M_PI*vmu/lambda);

  for (int i=0; i < points; i++) {
    if (fabs(Q[i]) <= cutoff) {
      R[i] = 1.;
    } else if (Q[i] > 0.) {
      const std::complex<double> f = sqrt(Q[i]*Q[i] - Sp);
      const std::complex<double> amp = (Q[i] - f) / (Q[i] + f);
      R[i] = real(amp*conj(amp));
    } else {
      const std::complex<double> f = sqrt(Q[i]*Q[i] - Sm);
      // Note: The following is inverted relative to q>0 case
      // because we are not using |Q|.  It is equivalent to
      // (-Q - f) / (-Q + f).
      const std::complex<double> amp = (Q[i] + f) / (Q[i] - f);
      R[i] = real(amp*conj(amp));
    }
  }
}


extern "C" void 
fresnel_amplitude(const double vrho, const double srho, 
		  const double vmu, const double smu, 
		  int points, const double Q[], std::complex<double> R[],
		  const double lambda)
{
  const double cutoff = 1e-10;
  const std::complex<double> Sp(16*M_PI*(srho-vrho),8*M_PI*smu/lambda);
  const std::complex<double> Sm(16*M_PI*(vrho-srho),8*M_PI*vmu/lambda);

  for (int i=0; i < points; i++) {
    if (fabs(Q[i]) <= cutoff) {
      R[i] = std::complex<double>(-1.,0.);
    } else if (Q[i] > 0.) {
      const std::complex<double> f = sqrt(Q[i]*Q[i] - Sp);
      R[i] = (Q[i] - f) / (Q[i] + f);
    } else {
      const std::complex<double> f = sqrt(Q[i]*Q[i] - Sm);
      // Note: The following is inverted relative to q>0 case
      // because we are not using |Q|.  It is equivalent to
      // (-Q - f) / (-Q + f).
      R[i] = (Q[i] + f) / (Q[i] - f);
    }
  }
}

// $Id$
