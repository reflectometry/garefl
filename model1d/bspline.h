#ifndef _BSPLINE_H
#define _BSPLINE_H

#ifdef __cplusplus
extern "C" {
#endif

/** \file
 * Fill a grid of n equally spaced values starting at t0 and stepping by
 * dt from a cubic bspline.
 *
 * Values beyond the range of the knots are either clamped to 0 or to
 * the value of the last control point, depending on how bspline was
 * compiled.
 *
 * The knot vector is of length N_knots.
 * The control vector is of length N_knots-4.
 * The t,v return vectors are of length n.
 */
void
bspline3(const int N_knots,
	 const double knot[], const double control[],
	 const double t_vals[],
	 const int n, double t[], double v[]);
void
bspline3_refine(int n, const double control[],
		double control_output[]);

void
bspline3_eval_derivs(const int N_knots, const double knot[],
	    const double control[], const double currT,
	    double result[]);

void
bspline3_eval_all_derivs(const int N_knots,
			 const double knot[], const double control[],
			 double t, double derivs[]);

int
bspline3_interpolate(const int N_data, double data[],
		     double result[], double work[]);


#ifdef __cplusplus
}
#endif

#endif /* _BSPLINE_H */

/* $Id$ */
