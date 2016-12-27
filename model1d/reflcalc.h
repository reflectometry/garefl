/* This program is public domain. */

#ifndef _REFLCALC_H
#define _REFLCALC_H

/*


Reflectivity calculations on rectangular slabs.

   Parameters:
   - M:      Number of slabs
   - d[M]:   depth of each slab; first and last ignored (semi-infinite)
   - rho[M]: scattering length density in number density units
   - mu[M]:  absorption; incident layer ignored
   - lambda: wavelength (required if absorption != 0)
   - P[M]:   magnetic scattering length density
   - exptheta[2*M]: cos(theta),sin(theta) pairs for magnetic scattering angle
   - Aguide  polarization angle relative to the guide (usually -90)
   - N:      Number of Q points
   - Q[N]:   Q points; negative for back reflectivity, which is implemented
             by reversing the order of the layers
   - R[N]:   returned reflectivity = |r|^2
   - r[N]:   returned complex reflectivity amplitude
   - real_r[N]: returned real portion of complex reflectivity amplitude

   Functions:
   reflectivity(M,d,rho,mu,lambda,N,Q,R)
   reflectivity_amplitude(M,d,rho,mu,lambda,N,Q,r)
   reflectivity_real(M,d,rho,mu,lambda,N,Q,real_r)
   reflectivity_imag(M,d,rho,mu,lambda,N,Q,real_r)
   magnetic_reflectivity(M,d,rho,mu,lambda,P,exptheta,Aguide,N,Q,R)
   magnetic_amplitude(M,d,rho,mu,lambda,P,exptheta,Aguide,N,Q,r)


Fresnel reflectivity from a single interface.

   Parameters:
   - Vrho:   incident medium rho
   - Srho:   substrate rho
   - Vmu:    incident medium absorption
   - Smu:    substrate absorption
   - lambda: wavelength (required if absorption != 0)
   - N:      number of points
   - f[N]:   Fresnel complex reflectivity amplitude
   - F[N]:   Fresnel reflectivity = |f|^2

   Functions:
   fresnel_reflectivity(Vrho, Srho, Vmu, Smu, N, Q, F, lambda)
   fresnel_amplitude(Vrho, Srho, Vmu, Smu, N, Q, f, lambda)

Resolution convolution

   Parameters:
   - M          number of points in model
   - Qin[M]     Q values of points in model
   - Rin[M]     R values of points in model
   - N          number of points in data
   - Q[N]       Q values of points in data
   - dQ[N]      dQ uncertainty in Q values
   - R[N]       returned estimate of reflectivity

   Functions:
   resolution(M, Qin, Rin, N, Q, dQ, R)

Resolution estimate to determine the width of the gaussian dQ to associate
with each point Q

   Parameters:
   - s1         slit 1 opening
   - s2         slit 2 opening
   - d          distance between slits
   - A3         angle at which s1,s2 were measured
   - dT         angular divergence for fixed slits
   - dToT       angular divergence for opening slits
   - L          wavelength
   - dLoL       wavelength divergence
   - N          number of points in data
   - Q[N]       Q values of points in data
   - dQ[N]      returned estimate of uncertainty


   Functions:
   dT = resolution_dT(s1,s2,d)
   dToT = resolution_dToT(s1,s2,d,A3)
   resolution_fixed(L,dLoL,dT,N,Q,dQ)
   resolution_varying(L,dLoL,dToT,N,Q,dQ)
   constant_resolution(res, N, dQ)

Resolution padding to determine the number of additional steps beyond
Q of a given step size in order to reach 0.001 on a gaussian of width dQ

   Parameters:
   - w          number of additional steps needed to compute reflectivity
   - step       step size

   Functions:
   w = resolution_padding(dQ, step)
*/

#ifdef __cplusplus
#include <complex>
#include <cmath>

extern "C" {
  typedef std::complex<Real> refl_complex;
#else
#include <math.h>
  typedef Real refl_complex;
#endif

#ifdef NEED_ERF
double erf(double);
#endif

#ifndef M_PI
# define M_PI 3.141592653589793
#endif

Real
adjust_alignment(Real Q,Real lambda,Real alignment);

void
reflectivity(const int layers, const Real d[],
	     const Real rho[], const Real mu[], const Real lambda,
             const Real alignment,
             const int points, const Real Q[], Real R[]);
void
reflectivity_amplitude(const int layers, const Real d[],
		       const Real rho[], const Real mu[], const Real L,
		       const Real alignment,
		       const int points, const Real Q[], refl_complex R[]);
void
reflectivity_real(const int layers, const Real d[],
		  const Real rho[], const Real mu[], const Real lambda,
                  const Real alignment,
		  const int points, const Real Q[], Real R[]);

void
reflectivity_imag(const int layers, const Real d[],
		  const Real rho[], const Real mu[], const Real lambda,
                  const Real alignment,
		  const int points, const Real Q[], Real R[]);

void
reflrough(const int layers, const Real d[], const Real sigma[],
	  const Real rho[], const Real mu[], const Real lambda,
          const Real alignment,
	  const int points, const Real Q[], Real R[]);
void
reflrough_amplitude(const int layers, const Real d[], const Real sigma[],
		    const Real rho[], const Real mu[], const Real lambda,
                    const Real alignment,
		    const int points, const Real Q[], refl_complex R[]);

void
magnetic_reflectivity(const int layers, const Real d[],
		      const Real rho[], const Real mu[], const Real L,
              const Real alignment,
		      const Real P[], const refl_complex expth[],
		      const Real Aguide, const int points, const Real Q[],
		      Real Ra[], Real Rb[], Real Rc[], Real Rd[]);
void
magnetic_amplitude(const int layers, const Real d[],
		   const Real rho[], const Real mu[], const Real L,
           const Real alignment,
		   const Real P[], const refl_complex expth[],
		   const Real Aguide, const int points, const Real Q[],
		   refl_complex Ra[], refl_complex Rb[],
		   refl_complex Rc[], refl_complex Rd[]);

void
fresnel_reflectivity(const Real vrho, const Real srho,
		     const Real vmu, const Real smu,
		     int points, const Real Q[], Real R[],
		     const Real lambda);
void
fresnel_amplitude(const Real vrho, const Real srho,
		  const Real vmu, const Real smu,
		  int points, const Real Q[], refl_complex R[],
		  const Real lambda);

#define T2Q(L,T) (4.*M_PI/(L)*sin((T)*M_PI/180.))
#define Q2T(L,Q) 180./M_PI*asin(Q*L/(4.*M_PI))
void
resolution(int Nin, const Real Qin[], const Real Rin[],
	   int N, const Real Q[], const Real dQ[], Real R[]);
Real
resolution_dT(Real s1,Real s2,Real d);
Real
resolution_dToT(Real s1,Real s2,Real d,Real A3);

void
constant_resolution(Real res, int n, Real dQ[]);

void
resolution_fixed(Real L, Real dLoL, Real dT,
		 int n, const Real Q[], Real dQ[]);

void
resolution_varying(Real L, Real dLoL, Real dToT,
		   int n, const Real Q[], Real dQ[]);

void
resolution_dQoQ(Real dQoQ, int n, const Real Q[], Real dQ[]);

int
resolution_padding(Real dQ, Real step);


#ifdef __cplusplus
}
#endif

#endif /* _REFLCALC_H */

/* $Id$ */
