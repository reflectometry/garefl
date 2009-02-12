// This program is public domain.

#include <complex>
#include "reflconfig.h"
#include "reflcalc.h"

#define Fr4xa   F77_FUNC(r4xa,R4XA)

extern "C" void
Fr4xa(const int &n, const double d[], const double rho[],
      const double mu[], const double &lambda, const double P[],
      const std::complex<double> expth[], const double &aguide,
      const double &q, std::complex<double> &A, std::complex<double> &B,
      std::complex<double> &C, std::complex<double> &D);

extern "C" void
magnetic_amplitude(const int layers, const double d[],
		   const double rho[], const double mu[], const double L,
		   const double alignment,
		   const double P[], const std::complex<double> expth[],
		   const double Aguide, const int points, const double Q[],
		   std::complex<double> Ra[], std::complex<double> Rb[],
		   std::complex<double> Rc[], std::complex<double> Rd[])
{
  for (int i=0; i < points; i++) {
    Fr4xa(layers,d,rho,mu,L,P,expth,Aguide,
          adjust_alignment(Q[i],alignment,L),
          Ra[i],Rb[i],Rc[i],Rd[i]);
  }
}

extern "C" void
magnetic_reflectivity(const int layers, const double d[],
		      const double rho[], const double mu[], const double L,
		      const double alignment,
		      const double P[], const std::complex<double> expth[],
		      const double Aguide, const int points, const double Q[],
		      double Ra[], double Rb[], double Rc[], double Rd[])
{
  for (int i=0; i < points; i++) {
    std::complex<double> A,B,C,D;
    Fr4xa(layers,d,rho,mu,L,P,expth,Aguide,
          adjust_alignment(Q[i],alignment,L),A,B,C,D);
    Ra[i] = real(A * conj(A));
    Rb[i] = real(B * conj(B));
    Rc[i] = real(C * conj(C));
    Rd[i] = real(D * conj(D));
  }
}

// $Id$
