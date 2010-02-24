/** \file
 * Prototypes and definitions for Levenberg - Marquardt fitting routines.
 */

#ifndef _NLLS_H

#ifndef Real
# ifdef USE_SINGLE
#  define Real float
# else
#  define Real double
# endif
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef void (*nlls_fn)(const Real *p, Real *y,
			const int n, const int k, void *data);

int nlls_worksize(int n, int k);

void nlls(nlls_fn f, int n, int k, 
	  const Real x[], const Real y[], Real p[], Real covar[]);

void box_nlls(nlls_fn f, int n, int k, 
	      const Real x[], const Real y[],
	      const Real bounds[], Real p[], Real covar[]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NLLS_H */
