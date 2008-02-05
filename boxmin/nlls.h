/** \file
 * Prototypes and definitions for Levenberg - Marquardt fitting routines.
 */

#ifndef _NLLS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef void (*nlls_fn)(const double *p, double *y,
			const int n, const int k, void *data);

int nlls_worksize(int n, int k);

void nlls(nlls_fn f, int n, int k, 
	  const double x[], const double y[], double p[], double covar[]);

void box_nlls(nlls_fn f, int n, int k, 
	      const double x[], const double y[], 
	      const double bounds[], double p[], double covar[]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _NLLS_H */
